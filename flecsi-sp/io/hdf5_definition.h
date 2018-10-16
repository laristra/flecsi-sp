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
class hdf5_definition__<2, T> : public flecsi::topology::mesh_definition__<2>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = hdf5_base__<2, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition__<2>;

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
    // Open file

    // open the hdf5 file and read cells
    
    try{
      // Open an existing file and dataset.
			H5File file(name, H5F_ACC_RDONLY);
      const H5std_string DATASET_NAME( "xCell" );

			DataSet xCells_dataset = file.openDataSet(DATASET_NAME); 
      
      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class = xCells_dataset.getTypeClass();

      assert(type_class == H5T_FLOAT);
      // Get dataspace of the dataset.
      DataSpace xCells_dataspace = xCells_dataset.getSpace();

      // Get the number of dimensions in the dataspace.
      int dimension = xCells_dataspace.getSimpleExtentNdims();

      assert (dimension==1);

      // Get the dimension size of each dimension in the dataspace and
      // display them.
      hsize_t dims_out[2];
      int ndims = xCells_dataspace.getSimpleExtentDims( dims_out, NULL);

      num_cells_=dims_out[0];

      std::cout << "num_cells =  " << num_cells_ << std::endl;
    
      real_t xCells[num_cells_];
      for (size_t i=0; i<num_cells_; i++)
        xCells[i]=0;

      xCells_dataset.read( xCells, PredType::NATIVE_DOUBLE, xCells_dataspace );
/*
//      dataset.read( xCells, PredType::NATIVE_DOUBLE,dataspace );
std::cout<<"xCels = "<<xCells[0]<<" , "<<xCells[1]<<" , "<<
xCells[2]<<" , "<< xCells[3]<<" , "<< xCells[4]<<" , "<<
xCells[5]<<" , "<< xCells[6]<<" , "<< xCells[7]<<" , "<<std::endl;
*/
      real_t yCells[num_cells_];
      for (size_t i=0; i<num_cells_; i++)
        yCells[i]=0;

      const H5std_string DATASET_NAME2( "yCell" );
   
      DataSet yCells_dataset = file.openDataSet(DATASET_NAME2);

      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class2 = yCells_dataset.getTypeClass();

      assert(type_class2 == H5T_FLOAT);
      // Get dataspace of the dataset.
      DataSpace yCells_dataspace = yCells_dataset.getSpace(); 

      yCells_dataset.read( yCells, PredType::NATIVE_DOUBLE, yCells_dataspace );

/*
std::cout<<"yCels = "<<yCells[0]<<" , "<<yCells[1]<<" , "<<
yCells[2]<<" , "<< yCells[3]<<" , "<< yCells[4]<<" , "<<
yCells[5]<<" , "<< yCells[6]<<" , "<< yCells[7]<<" , "<<std::endl;
*/
      for (size_t i=0; i<num_cells_; i++)
        cells_.push_back(xCells[i]);
       for (size_t i=0; i<num_cells_; i++)
        cells_.push_back(yCells[i]);

    }
    //catch failure caused by the H5File operations
    catch(FileIException error)
    {
    	error.printErrorStack();
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

   //read vertex information:

   try{
      // Open an existing file and dataset.
      H5File file(name, H5F_ACC_RDONLY);
      const H5std_string DATASET_NAME( "xVertex" );

      DataSet xVertex_dataset = file.openDataSet(DATASET_NAME);
   
      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class = xVertex_dataset.getTypeClass();

      assert(type_class == H5T_FLOAT);
      // Get dataspace of the dataset.
      DataSpace xVertex_dataspace = xVertex_dataset.getSpace();

      // Get the number of dimensions in the dataspace.
      int dimension = xVertex_dataspace.getSimpleExtentNdims();

      assert (dimension==1);

      // Get the dimension size of each dimension in the dataspace and
      // display them.
      hsize_t dims_out[2];
      int ndims = xVertex_dataspace.getSimpleExtentDims( dims_out, NULL);

      num_vertices_=dims_out[0];

      std::cout << "num_vertices =  " << num_vertices_ << std::endl;
   
      real_t xVertex[num_vertices_];
      for (size_t i=0; i<num_vertices_; i++)
        xVertex[i]=0;

      xVertex_dataset.read( xVertex, PredType::NATIVE_DOUBLE,
				xVertex_dataspace );

//      dataset.read( xCells, PredType::NATIVE_DOUBLE,dataspace );
std::cout<<"xVertex = "<<xVertex[0]<<" , "<<xVertex[1]<<" , "<<
xVertex[2]<<" , "<< xVertex[3]<<" , "<< xVertex[4]<<" , "<<
xVertex[5]<<" , "<< xVertex[6]<<" , "<< xVertex[7]<<" , "<<std::endl;

      real_t yVertex[num_vertices_];
      for (size_t i=0; i<num_vertices_; i++)
        yVertex[i]=0;

      const H5std_string DATASET_NAME2( "yVertex" );

      DataSet yVertex_dataset = file.openDataSet(DATASET_NAME2);

      //Get the class of the datatype that is used by the dataset.
      H5T_class_t type_class2 = yVertex_dataset.getTypeClass();

      assert(type_class2 == H5T_FLOAT);
      // Get dataspace of the dataset.
      DataSpace yVertex_dataspace = yVertex_dataset.getSpace();

      yVertex_dataset.read( yVertex, PredType::NATIVE_DOUBLE,
				yVertex_dataspace );
    /*
std::cout<<"yVertex = "<<yVertex[0]<<" , "<<yVertex[1]<<" , "<<
yVertex[2]<<" , "<< yVertex[3]<<" , "<< yVertex[4]<<" , "<<
yVertex[5]<<" , "<< yVertex[6]<<" , "<< yVertex[7]<<" , "<<std::endl;
*/
      for (size_t i=0; i<num_vertices_; i++)
        vertices_.push_back(xVertex[i]);
       for (size_t i=0; i<num_vertices_; i++)
        vertices_.push_back(yVertex[i]);
    }//end try

    //catch failure caused by the H5File operations
    catch(FileIException error)
    {
      error.printErrorStack();
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

//  size_t dimension()
//  {
//    return 2;
//  }

private:
  //============================================================================
  // Private data
  //============================================================================

  size_t num_vertices_;
  size_t num_cells_;
  size_t num_edges_;

  //! \brief storage for element verts
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for cells coordinates
  vector<real_t> cells_;
  
  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;
  
};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class hdf5_definition__<3, T> : public flecsi::topology::mesh_definition__<3>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = hdf5_base__<3, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition__<3>;

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

    //--------------------------------------------------------------------------
    // Open file and read the initialization paeameters

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


    //--------------------------------------------------------------------------
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
      default:
        return entities_.at(dim).at(0).size();
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

  //! \brief storage for element verts
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;
};

} // namespace io
} // namespace flecsi
