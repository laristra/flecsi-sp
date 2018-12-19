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
#include <hdf5.h>
#include "H5Cpp.h"

// system includes
#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>


//flecsi-sp includes
#include "flecsi-sp/io/io_utils.h"

namespace flecsi_sp {
namespace io {

namespace detail{
//reading entity coordinates
template <typename T>
size_t read_coordinates (
  const std::string & file_name,
  const H5std_string &x_dataset_name,
  const H5std_string &y_dataset_name,
  std::vector<T> &entity)
{
      using namespace H5;

      H5File file(file_name, H5F_ACC_RDONLY);

      DataSet x_dataset = file.openDataSet(x_dataset_name);

      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class = x_dataset.getTypeClass();

      assert(type_class == H5T_FLOAT);
      // Get dataspace of the dataset.
      DataSpace x_dataspace = x_dataset.getSpace();

      // Get the number of dimensions in the dataspace.
      int x_dimension = x_dataspace.getSimpleExtentNdims();

      assert (x_dimension==1);

      // Get the dimension size of each dimension in the dataspace and
      // display them.
      hsize_t dims_out[2];
      int ndims = x_dataspace.getSimpleExtentDims( dims_out, NULL);

      size_t num_entities=dims_out[0];

      //std::cout << "num_cells =  " << num_cells_ << std::endl;

      T xCoord[num_entities];
      for (size_t i=0; i<num_entities; i++)
        xCoord[i]=0;

      x_dataset.read( xCoord, PredType::NATIVE_DOUBLE, x_dataspace );

      T yCoord[num_entities];
      for (size_t i=0; i<num_entities; i++)
        yCoord[i]=0;

      DataSet y_dataset = file.openDataSet(y_dataset_name);

      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class2 = y_dataset.getTypeClass();

      assert(type_class2 == H5T_FLOAT);
      // Get dataspace of the dataset.
      DataSpace y_dataspace = y_dataset.getSpace();

      y_dataset.read( yCoord, PredType::NATIVE_DOUBLE, y_dataspace );

      for (size_t i=0; i<num_entities; i++)
        entity.push_back(xCoord[i]);
      for (size_t i=0; i<num_entities; i++)
        entity.push_back(yCoord[i]);

      return num_entities;

}

//reading connectivity information
void read_connectivity(
  const std::string & file_name,
	const H5std_string &dataset_name,
  std::vector<std::vector<size_t>> &connectivity
)
{
      using namespace H5;

      H5File file(file_name, H5F_ACC_RDONLY);
      
      //open dataset
      DataSet conn_dataset = file.openDataSet(dataset_name);

      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class = conn_dataset.getTypeClass();
      //connectivity type should be Integer
      assert(type_class == H5T_INTEGER);
      // Get dataspace of the dataset.
      DataSpace conn_dataspace = conn_dataset.getSpace();

      // Get the number of dimensions in the dataspace.
      int dimension = conn_dataspace.getSimpleExtentNdims();

      assert (dimension==2);

      // Get the dimension size of each dimension in the dataspace and
      // display them.
      hsize_t dims_out[2];
      int ndims = conn_dataspace.getSimpleExtentDims( dims_out, NULL);

      size_t num_entities1=dims_out[0];
      size_t num_entities2=dims_out[1];

      size_t conn[dims_out[0]][dims_out[1]];
      for (size_t i=0; i<num_entities1; i++)
       for (size_t j=0; j<num_entities2; j++)
         conn[i][j]=0;

     //conn_dataset.read( conn, H5T_INTEGER, conn_dataspace );
      conn_dataset.read( conn, PredType::NATIVE_LONG, conn_dataspace );


      // Adding connectivity information from the HDF5 array to FleCSI
      // connectivity
      for (size_t i=0; i<num_entities1; i++){
       std::vector<size_t> tmp;
       for (size_t j=0; j<num_entities2; j++){
         //there is no connectivity if  conn[i][j]==0 
         if (conn[i][j]!=0)
           //in HDF5 file entitiy ordering starts from 1 where in flecsi we
           //start orfering from 0
           tmp.push_back(conn[i][j]-1);
       }
       connectivity.push_back(tmp);
      }
}//read_connectivity

void dump_connectivity(
std::vector<std::vector<size_t>> &connectivity
)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank ==0){
    std::cout<<"dump connectivity";
    for (size_t i=0; i<connectivity.size(); i++)
    {
     std::cout<<"conn["<< i <<"] = " <<std::endl;
     auto &tmp = connectivity[i];
     for (size_t j=0; j< tmp.size(); j++)
			std::cout <<tmp[j]<<"   ";
     std::cout<<std::endl;
    }//fo
    
  }//if 
}//dump_connectivity



}//end details

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class hdf5_base__ {

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
  template<typename U>
  using vector = typename std::vector<U>;

  //! \brief an alias for the matrix class
  template<typename U>
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
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class hdf5_definition__ {};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class hdf5_definition__<2, T> : public flecsi::topology::mesh_definition_u<2>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = hdf5_base__<2, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition_u<2>;

  //! the number of dimensions
  using mesh_definition_t::dimension;

  //! the floating point type
  using real_t = typename base_t::real_t;
  //! the index type
  using index_t = typename base_t::index_t;

  //! the vector type
  template<typename U>
  using vector = typename base_t::template vector<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;

  using point_t = mesh_definition_t::point_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Default constructor
  hdf5_definition__() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  hdf5_definition__(const std::string & filename) {
    read(filename);
  }

  /// Copy constructor (disabled)
  hdf5_definition__(const hdf5_definition__ &) = delete;

  /// Assignment operator (disabled)
  hdf5_definition__ & operator=(const hdf5_definition__ &) = delete;

  /// Destructor
  ~hdf5_definition__() = default;

  //============================================================================
  //! \brief Implementation of hdf5 mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void read(const std::string & name) {

    clog(info) << "Reading mesh from: " << name << std::endl;

    using namespace H5;
  
    //--------------------------------------------------------------------------
    try
	  {
      //read coordinates

      {//cells
        const H5std_string x_dataset_name( "xCell" );
        const H5std_string y_dataset_name( "yCell" );
        num_cells_ =
          detail::read_coordinates ( name , x_dataset_name,
          y_dataset_name, cells_);
      }//scope

      {//vertices
        const H5std_string x_dataset_name( "xVertex" );
        const H5std_string y_dataset_name( "yVertex" );
        num_vertices_ =
          detail::read_coordinates ( name , x_dataset_name,
          y_dataset_name, vertices_);
       }//scope

/* 
std::cout<<"num_cells = "<<num_cells_<< " , num_vert = "<<
num_vertices_<<std::endl;

for (size_t i = 0; i< 2*num_cells_; i++)
std::cout<< "cells_["<<i<<"] = %g"<<cells_[i]<<std::endl;

for (size_t i = 0; i< 2*num_vertices_; i++)
std::cout<< "vertices_["<<i<<"] = "<<vertices_[i]<<std::endl;
*/

      //read connectivity information
      const H5std_string vertOnCell_NAME( "verticesOnCell" );
      detail::read_connectivity(name, vertOnCell_NAME, entities_[2][0]);

      const H5std_string edgesOnCell_NAME( "edgesOnCell" );
      detail::read_connectivity(name, edgesOnCell_NAME, entities_[2][1]);

      const H5std_string vertOnEdge_NAME( "verticesOnEdge" );
      detail::read_connectivity(name, vertOnEdge_NAME, entities_[1][0]);

      const H5std_string cellsOnVertex_NAME( "cellsOnVertex" );
      detail::read_connectivity(name, cellsOnVertex_NAME, entities_[0][2]);

      const H5std_string cellsOnEdge_NAME( "cellsOnEdge" );
      detail::read_connectivity(name, cellsOnEdge_NAME, entities_[1][2]);

      const H5std_string edgesOnVertex_NAME( "edgesOnVertex" );
      detail::read_connectivity(name, edgesOnVertex_NAME, entities_[0][1]);

      const H5std_string cellsOnCell_NAME( "cellsOnCell" );
      detail::read_connectivity(name, cellsOnCell_NAME, entities_[2][2]);

      const H5std_string edgesOnEdge_NAME( "edgesOnEdge" );
      detail::read_connectivity(name, edgesOnEdge_NAME, entities_[1][1]);

		}//end try
    
    //catch failure caused by the H5File operations
    catch(FileIException error)
    {
      //error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printError();
   }
   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printError();
   }
   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printError();
   }



  }

  //============================================================================
  //! \brief Implementation of hdf5 mesh write for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  template<typename U = int>
  void write(
      const std::string & name,
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
          element_sets = {},
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
          node_sets = {}) const {
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
        clog_fatal(
            "Dimension out of range: 0 < " << dim << " </ " << dimension());
        return 0;
    } 
  }

  /// Return the set of vertices of a particular entity.
  /// \param [in] dimension  The entity dimension to query.
  /// \param [in] entity_id  The id of the entity in question.
  const auto & entities(size_t from_dim, size_t to_dim) const {
    return entities_.at(from_dim).at(to_dim);
  } // vertices

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in question.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override {
    return entities_.at(from_dim).at(to_dim).at(from_id);
  } // vertices

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  template<typename POINT_TYPE>
  auto vertex(size_t vertex_id) const {
    auto num_vertices = vertices_.size() / dimension();
    POINT_TYPE p;
    for (int i = 0; i < dimension(); ++i)
      p[i] = vertices_[i * num_vertices + vertex_id];
    return p;
  } // vertex

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


} // namespace io
} // namespace flecsi
