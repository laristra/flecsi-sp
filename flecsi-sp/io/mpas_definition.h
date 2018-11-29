/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

/// \file

// user includes
#include <flecsi/topology/mesh_definition.h>
#include <flecsi/utils/logging.h>

// thirdparty includes
extern "C" {
#include <hdf5.h>
}
//#include "H5Cpp.h"

// system includes
#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// flecsi-sp includes
#include "flecsi-sp/io/io_utils.h"

namespace flecsi_sp {
namespace io {

/* Notes
 *
 * The use of HDF5 here is fairly primitive. Any API failure will cause a
 * program abort, but there are likely ways to recover gracefully in many cases.
 *
 * It might be worth trying to manage what resources are open and need to be
 * closed with good C++ idioms, but I don't know if that really makes sense in a
 * multi-node context.  I'll have to answer that once I've got a better handle
 * on the parallel use case.
 * */

/*!
 * \brief Simple wrappers for some HDF5 functions that abort on error
 */
namespace h5 {

/*!
 * \brief Open an HDF5 file with the given flags and access properties
 *
 * This will fail hard if there are problems opening the file (e.g. not
 * readable/writeable, does not exist, is not an HDF5 file)
 *
 * \param[in] name file name to open
 *
 * \param[in] flags read/write/create flags to open the file with
 *
 * \param[in] fapl_id handle to the file access property list to use
 *
 * \return a handle to the open file
 */
hid_t open_file(const std::string &name, unsigned flags, hid_t fapl_id) {
  // Note: Third param (fapl_id) is the id for file access props For
  //   parallel access, the fapl_id will hold the communicator
  hid_t file = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file < 0) {
    // Try to give some meaningful feedback if open fails
    htri_t res = H5Fis_hdf5(name.c_str());

    if (res == 0) {  // res == 0 ==> exists but not in HDF5 format
      clog_fatal(name << " is not an HDF5 file\n");
    }  // if
    else {
      // TODO: At least check if the file exists
      clog_fatal("Couldn't open " << name << " ... not sure why");
    }  // else
  }    // if
}  // open_file

/*!
 * \brief wrapper for H5Fclose
 *
 * Note: API-level failure causes program abort.  What effect this may have on
 * the file being closed is unclear.
 *
 * \param[in] file Handle to an open file
 */
void close_file(hid_t file_id) {
  if (H5Fclose(file_id) < 0) {
    clog_fatal("Closing file with handle " << file_id << " failed!");
  }  // if
}  // close_file

/*!
 * \brief Open the specified dataset, failing hard if the open is
 * unsuccessful
 *
 * \param[in] loc the handle to an open location identifier
 *
 * \param[in] name name of the dataset to access
 *
 * \return a handle to the dataset
 */
hid_t open_dataset(hid_t loc, const std::string &name) {
  hid_t handle = H5Dopen(loc, name.c_str(), H5P_DEFAULT);
  clog_assert(handle >= 0, "Failed to open dataset: " << name);
  return handle;
}  // open_dataset

/*!
 * \brief wrapper for H5Dclose
 *
 * Note: API-level failure causes program abort.  What effect this may have on
 * the underlying file is unclear.
 *
 * \param[in] dset_id Handle to an open dataset
 */
void close_dataset(hid_t dset_id) {
  if (H5Dclose(dset_id) < 0) {
    clog_fatal("Closing dataset with handle " << dset_id << " failed!");
  }  // if
}  // close_dataset

/*!
 * \brief Get a handle to a copy of a dataset's dataspace, failing hard if the
 *        open is unsuccessful
 *
 * Note that callers are responsible for closing the copied dataset using
 * close_dataspace (or H5Sclose)
 *
 * \param[in] dset_id handle to the open dataspace
 *
 * \return a handle to the copied dataspace
 */
hid_t get_space(hid_t dset_id) {
  hid_t handle = H5Dget_space(dset_id);

  // TODO: use H5Iget_name() to give a more meaningful error message
  clog_assert(handle >= 0, "Failed to make a copy of dataspace");

  return handle;
}  // get_space

/*!
 * \brief wrapper for H5Sclose
 *
 * Note: API-level failure causes program abort.  What effect this may have on
 * the underlying file is unclear.
 *
 * \param[in] dspace_id Handle to an open dataspace
 */
void close_dataspace(hid_t dspace_id) {
  if (H5Sclose(dspace_id) < 0) {
    clog_fatal("Closing dataspace with handle " << dspace_id << " failed!");
  }  // if
}  // close_dataspace

/*!
 * \brief fails hard if the given dataset does not hold the expected class
 *
 * \param[in] dset handle to the dataset
 *
 * \param[in] expect_class expected class of data
 */
void assert_dataset_typeclass(hid_t dset, H5T_class_t expect_class) {
  hid_t dtype_id = H5Dget_type(dset);
  H5T_class_t real_class = H5Tget_class(dtype_id);
  clog_assert(real_class == expect_class, "Unexpected dataset type class "
                                              << real_class
                                              << " != " << expect_class);
  H5Tclose(dtype_id);
}

/*!
 * \brief fails hard if the given dataset does not hold the expected type
 *
 * Example usage: assert_dataset_type(dset_id, H5T_NATIVE_INT)
 *
 * \param[in] dset handle to the dataset
 *
 * \param[in] expect_type handle to hid_t of expected type of data
 */
void assert_dataset_type(hid_t dset, hid_t expect_type) {
  hid_t dtype_id = H5Dget_type(dset);

  // TODO: get type names for better message
  clog_assert(H5Tequal(dtype_id, expect_type),
              "Dataset does not hold expected type!");

  H5Tclose(dtype_id);
}  // assert_dataset_type

/*! \brief Get the rank of a dataset, failing hard if the HDF5 API call is
 *         unsuccessful
 *
 *  \param[in] dset_id handle to the dataset
 *
 *  \return the rank of the dataset (number of dimensions)
 */
int get_rank(hid_t dset_id) {
  hid_t space = get_space(dset_id);
  int rank = H5Sget_simple_extent_ndims(space);

  // TODO: use H5Iget_name() to give a more meaningful error message
  clog_assert(rank >= 0, "Failed to retrieve dataset rank");

  H5Sclose(space);
  return rank;
}  // get_rank

/*!
 * \brief Wrapper for H5Sget_simple_extent_dims to retrieve dimension
 * information of a dataset
 *
 * The caller must allocate appropriate space to hold the dimension extents
 *
 * \param[in] dset_id dataset handle
 *
 * \param[out] dims pointer to array
 *
 * \return the number of dimesions (rank) of the dataset
 */
int get_simple_dims(hid_t dset_id, hsize_t *dims) {
  hid_t space = get_space(dset_id);
  int rank = H5Sget_simple_extent_dims(space, dims, NULL);
  H5Sclose(space);
  return rank;
}  // get_simple_dims

/*!
 * \brief Read the coordinates for a set of cells
 *
 * \param[in] file_id handle to open HDF5 file
 *
 * \param[in] name name of the dataset to read from
 *
 * \param[out] entities vector to add coordinates to
 *
 * \return number of entities read
 */
template <typename T>
size_t read_coord_dset(hid_t file_id, const std::string &name,
                       std::vector<T> &entities) {
  hid_t dset = open_dataset(file_id, name);

  // Probably only need the latter of these two checks, but redundancy doesn't
  // hurt for now
  assert_dataset_typeclass(dset, H5T_FLOAT);
  assert_dataset_type(dset, H5T_NATIVE_DOUBLE);

  // Expecting one-dimensional data
  clog_assert(get_rank(dset) == 1,
              "Expected single dimension for coordinate data");

  // since the rank is known, we can use the stack here.  Ordinarily you'd need
  // to allocate `dims` dynamically
  hsize_t dims[1];
  get_simple_dims(dset, dims);
  size_t num_entities = dims[0];

  // Since the rank of these datasets is 1, a simple array is OK.  For
  // higher-dimensional data, a more dynamic allocation is probably necessary
  T coords[num_entities];
  for (size_t i = 0; i < num_entities; i++) coords[i] = 0;

  // Read native doubles, all spaces to all spaces with default props
  herr_t status =
      H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);

  for (size_t i = 0; i < num_entities; i++) entities.push_back(coords[i]);

  // Done with the dataset, close it.
  close_dataset(dset);

  return num_entities;
}  // read_coord_set

}  // namespace h5

namespace detail {

// reading entity coordinates
template <typename T>
size_t read_coordinates(const std::string &file_name,
                        const std::string &x_dataset_name,
                        const std::string &y_dataset_name,
                        std::vector<T> &entity) {
  // handle to file to open - equivalent to root Group
  hid_t file = h5::open_file(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // read the coordinate datasets individually, adding the data to the entity
  // vector
  size_t xsize = h5::read_coord_dset(file, x_dataset_name, entity);
  size_t ysize = h5::read_coord_dset(file, y_dataset_name, entity);

  clog_assert(xsize == ysize,
              "Expected x and y coordinate datasets to have the same size");

  // TODO: This could fail - do something about it?
  H5Fclose(file);

  return xsize;
}  // read_coordinates

// reading connectivity information
void read_connectivity(const std::string &file_name,
                       const std::string &dataset_name,
                       std::vector<std::vector<size_t>> &connectivity) {
  hid_t file = h5::open_file(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dset = h5::open_dataset(file, dataset_name);

  h5::assert_dataset_typeclass(dset, H5T_INTEGER);
  h5::assert_dataset_type(dset, H5T_NATIVE_LONG);

  // Expecting two-dimensional data
  clog_assert(h5::get_rank(dset) == 2,
              "Expected two dimensions for coordinate data");

  // Get the dimension size of each dimension in the dataspace and
  // display them.
  hsize_t dims[2];
  h5::get_simple_dims(dset, dims);

  size_t dim1_size = dims[0];
  size_t dim2_size = dims[1];

  size_t conn[dims[0]][dims[1]];

  for (size_t i = 0; i < dim1_size; i++)
    for (size_t j = 0; j < dim2_size; j++) conn[i][j] = 0;

  H5Dread(dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, conn);

  // Adding connectivity information from the HDF5 array to FleCSI
  // connectivity
  for (size_t i = 0; i < dim1_size; i++) {
    std::vector<size_t> tmp;
    for (size_t j = 0; j < dim2_size; j++) {
      // there is no connectivity if  conn[i][j]==0
      if (conn[i][j] != 0)
        // in HDF5 file entitiy ordering starts from 1 where in flecsi we
        // start orfering from 0
        tmp.push_back(conn[i][j] - 1);
    }
    connectivity.push_back(tmp);
  }
}  // read_connectivity

void dump_connectivity(std::vector<std::vector<size_t>> &connectivity) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "dump connectivity";
    for (size_t i = 0; i < connectivity.size(); i++) {
      std::cout << "conn[" << i << "] = " << std::endl;
      auto &tmp = connectivity[i];
      for (size_t j = 0; j < tmp.size(); j++) std::cout << tmp[j] << "   ";
      std::cout << std::endl;
    }  // fo

  }  // if
}  // dump_connectivity

}  // namespace detail

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        MPAS HDF5 file format.
////////////////////////////////////////////////////////////////////////////////
template <int D, typename T>
class mpas_base__ {
 public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the size type
  using size_t = std::size_t;
  //! \brief the counter type
  using counter_t = flecsi::utils::counter_t;
  //! \brief the floating point type
  using real_t = T;
  //! \brief the type used for indexing arrays
  using index_t = std::size_t;

  //! \brief an alias for the vector class
  template <typename U>
  using vector = typename std::vector<U>;

  //! \brief an alias for the matrix class
  template <typename U>
  using sparse_matrix = vector<vector<U>>;

  //! \brief the data type for an index vector
  using index_vector_t = vector<index_t>;

  //! \brief the data type for connectivity
  using connectivity_t = sparse_matrix<index_t>;

  //! the number of dimensions
  static constexpr size_t num_dims = D;

  enum class block_t { tri, quad, polygon, tet, hex, polyhedron, unknown };
};
////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        MPAS HDF5 file format.
////////////////////////////////////////////////////////////////////////////////
template <int D, typename T>
class mpas_definition__ {};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        MPAS HDF5 file format.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
class mpas_definition__<2, T> : public flecsi::topology::mesh_definition__<2> {
 public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = mpas_base__<2, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition__<2>;

  //! the number of dimensions
  using mesh_definition_t::dimension;

  //! the floating point type
  using real_t = typename base_t::real_t;
  //! the index type
  using index_t = typename base_t::index_t;

  //! the vector type
  template <typename U>
  using vector = typename base_t::template vector<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;

  using point_t = mesh_definition_t::point_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Default constructor
  mpas_definition__() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  mpas_definition__(const std::string &filename) { read(filename); }

  /// Copy constructor (disabled)
  mpas_definition__(const mpas_definition__ &) = delete;

  /// Assignment operator (disabled)
  mpas_definition__ &operator=(const mpas_definition__ &) = delete;

  /// Destructor
  ~mpas_definition__() = default;

  //============================================================================
  //! \brief Implementation of mpas mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //============================================================================
  void read(const std::string &name) {
    clog(info) << "Reading mesh from: " << name << std::endl;

    //--------------------------------------------------------------------------
    // read coordinates

    {  // cells
      const std::string x_dataset_name("xCell");
      const std::string y_dataset_name("yCell");
      num_cells_ = detail::read_coordinates(name, x_dataset_name,
                                            y_dataset_name, cells_);
    }  // scope

    {  // vertices
      const std::string x_dataset_name("xVertex");
      const std::string y_dataset_name("yVertex");
      num_vertices_ = detail::read_coordinates(name, x_dataset_name,
                                               y_dataset_name, vertices_);
    }  // scope

    /*
      std::cout<<"num_cells = "<<num_cells_<< " , num_vert = "<<
      num_vertices_<<std::endl;

      for (size_t i = 0; i< 2*num_cells_; i++)
      std::cout<< "cells_["<<i<<"] = %g"<<cells_[i]<<std::endl;

      for (size_t i = 0; i< 2*num_vertices_; i++)
      std::cout<< "vertices_["<<i<<"] = "<<vertices_[i]<<std::endl;
    */

    // read connectivity information
    const std::string vertOnCell_NAME("verticesOnCell");
    detail::read_connectivity(name, vertOnCell_NAME, entities_[2][0]);

    const std::string edgesOnCell_NAME("edgesOnCell");
    detail::read_connectivity(name, edgesOnCell_NAME, entities_[2][1]);

    const std::string vertOnEdge_NAME("verticesOnEdge");
    detail::read_connectivity(name, vertOnEdge_NAME, entities_[1][0]);

    const std::string cellsOnVertex_NAME("cellsOnVertex");
    detail::read_connectivity(name, cellsOnVertex_NAME, entities_[0][2]);

    const std::string cellsOnEdge_NAME("cellsOnEdge");
    detail::read_connectivity(name, cellsOnEdge_NAME, entities_[1][2]);

    const std::string edgesOnVertex_NAME("edgesOnVertex");
    detail::read_connectivity(name, edgesOnVertex_NAME, entities_[0][1]);

    const std::string cellsOnCell_NAME("cellsOnCell");
    detail::read_connectivity(name, cellsOnCell_NAME, entities_[2][2]);

    const std::string edgesOnEdge_NAME("edgesOnEdge");
    detail::read_connectivity(name, edgesOnEdge_NAME, entities_[1][1]);
  }

  //============================================================================
  //! \brief Implementation of mpas mesh write for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //============================================================================
  template <typename U = int>
  void write(
      const std::string &name,
      const std::initializer_list<std::pair<const char *, std::vector<U>>>
          &element_sets = {},
      const std::initializer_list<std::pair<const char *, std::vector<U>>>
          &node_sets = {}) const {
    clog(info) << "Writing mesh to: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file
  }

  //============================================================================
  // Required Overrides
  //============================================================================

  /// Return the number of entities of a particular dimension
  /// \param [in] dim  The entity dimension to query.
  size_t num_entities(size_t dim) const override {
    switch (dim) {
      case 0:
        return vertices_.size() / dimension();
      case 1:
      case 2:
        return entities_.at(dim).at(0).size();
      default:
        clog_fatal("Dimension out of range: 0 < " << dim << " </ "
                                                  << dimension());
        return 0;
    }
  }

  /// Return the set of vertices of a particular entity.
  /// \param [in] dimension  The entity dimension to query.
  /// \param [in] entity_id  The id of the entity in question.
  const auto &entities(size_t from_dim, size_t to_dim) const {
    return entities_.at(from_dim).at(to_dim);
  }  // vertices

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in question.
  std::vector<size_t> entities(size_t from_dim, size_t to_dim,
                               size_t from_id) const override {
    return entities_.at(from_dim).at(to_dim).at(from_id);
  }  // vertices

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  template <typename POINT_TYPE>
  auto vertex(size_t vertex_id) const {
    auto num_vertices = vertices_.size() / dimension();
    POINT_TYPE p;
    for (int i = 0; i < dimension(); ++i)
      p[i] = vertices_[i * num_vertices + vertex_id];
    return p;
  }  // vertex

 private:
  //============================================================================
  // Private data
  //============================================================================

  size_t num_vertices_ = 0;
  size_t num_cells_ = 0;
  //  size_t num_edges_ = 0;

  //! \brief storage for element verts
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for cells coordinates
  vector<real_t> cells_;

  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;
};

}  // namespace io
}  // namespace flecsi_sp
