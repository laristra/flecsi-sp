/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes

// library includes
#include <portage/support/portage.h>
//#include <ristra/math/array.h>

// system includes
#include <utility>
#include <cstring>
#include <string>
#include <map>
#include <cstdarg>

namespace flecsi_sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
/// \brief Provides access to data stored in Flecsi_State
////////////////////////////////////////////////////////////////////////////////
template< typename M >
class flecsi_state_t {

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief The mesh type
  using mesh_t = M;
  //! \brief the size type
  using size_t = typename mesh_t::size_t;
  //! \brief the real type
  using real_t = typename mesh_t::real_t;
  
  using vector_t = typename mesh_t::vector_t;


  //! \brief The entity kind type
  using entity_kind_t = Portage::Entity_kind;
  //! \brief The entity type 
  using entity_type_t = Portage::Entity_type;
  //! \brief The field type
  using field_type_t = Portage::Field_type;

public:


  //============================================================================
  // Member Variables
  //============================================================================

  std::map < std::string, double*> var_map;
  std::map < std::string, std::string> type_map;
  int number_materials=0;

  //============================================================================
  // Constructors
  //============================================================================

  //!  \brief Default constructor.
  //!  \param[in] mesh The minimum coordinates of the domain.
  //!  \param[in] mesh The minimum coordinates of the domain.
  explicit flecsi_state_t(mesh_t & mesh) : mesh_(&mesh)
    {}

  //! Default constructor deleted
  flecsi_state_t() = default;

  //! Default copy constructor
  flecsi_state_t(const flecsi_state_t &) = default;

  //! Default assignment operator
  flecsi_state_t & operator=(const flecsi_state_t &) = default;

  //============================================================================
  // Public Members
  //============================================================================

  //! \brief Add a field that needs to be remapped to the variable map
  void add_field( std::string var_name, double* data, std::string type_name) {

    var_map.insert( std::pair <std::string, double*> (var_name, data));
    type_map.insert( std::pair < std::string, std::string> (var_name, type_name));

  }

  //! \brief Number of materials in problem
  int num_materials() const {
    return number_materials;
  }

  void num_materials(int num_mats){
    number_materials = num_mats;
  }

  //! \brief Name of material
  std::string material_name(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return "UNKNOWN";
  }

  //! \brief Get number of cells containing a particular material
  //! \param matid    Index of material (0, num_materials()-1)
  //! \return         Number of cells containing material 'matid'
  int mat_get_num_cells(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return 0;
  }

  //! \brief Get cell indices containing a particular material
  //! \param matid    Index of material (0, num_materials()-1)
  //! \param matcells Cells containing material 'matid'
  void mat_get_cells(int matid, std::vector<int> *matcells) const {
    assert(matid >= 0 && matid < num_materials());
    matcells->clear();
  }

  //! \brief Get number of materials contained in a cell
  //! \param cellid  Index of cell in mesh
  //! \return        Number of materials in cell
  int cell_get_num_mats(int cellid) const {
    return 0;
  }

  //! \brief Get the IDs of materials in a cell
  //! \param cellid    Index of cell in mesh
  //! \param cellmats  Indices of materials in cell
  void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
    cellmats->clear();
  }

  //! \brief Get the local index of mesh cell in material cell list
  //! \param meshcell    Mesh cell ID
  //! \param matid       Material ID
  //! \return             Local cell index in material cell list
  int cell_index_in_material(int meshcell, int matid) const {
    return -1;
  }
  
  //! \brief Type of field (MESH_FIELD or MULTIMATERIAL_FIELD)
  //! \param[in] onwhat   Entity_kind that field is defined on
  //! \param[in] varname  Name of field
  //! \return             Field_type
  field_type_t field_type(entity_kind_t on_what, std::string const& var_name)
      const {
    return field_type_t::MESH_FIELD;  // MULTI-MATERIAL FIELDS NOT ACCESSED YET
  }
  

  //! \brief Get the entity type on which the given field is defined
  //! \param[in] var_name The string name of the data field
  //! \return The Entity_kind enum for the entity type on which the field is defined
  //!
  //! \todo  THIS ASSUMES ONLY DOUBLE VECTORS - WE HAVE TO ACCOUNT FOR OTHER TYPES
  //!        OR WE HAVE TO GENERALIZE THE FIND FUNCTION!!!
  //! \todo  THIS ALSO DOES NOT CHECK FOR OTHER ENTITY TYPES LIKE EDGE, FACE,
  //!        SIDE, WEDGE AND CORNER
  entity_kind_t get_entity(std::string const& var_name) const 
  {
    
    auto type = type_map.at(var_name);

    if (type == "CELL"){
      return entity_kind_t::CELL;
    } else if (type == "NODE"){
      return entity_kind_t::NODE;
    } else {
      return entity_kind_t::UNKNOWN_KIND;
    }

  }


  //! \brief Get pointer to scalar data
  //! \param[in] on_what The entity type on which to get the data
  //! \param[in] var_name The string name of the data field
  //! \param[in,out] data A pointer to an array of data
  template <class T>
  void mesh_get_data(entity_kind_t on_what, std::string const& var_name, 
                     T ** data) const {

    *data = (var_map.at(var_name));
    return;
    // Ignore on_what here - the state manager knows where it lives
    // based on its name
    auto mesh_handle = flecsi_get_client_handle(mesh_t, meshes, mesh0);

    // first check cells
    auto cell_field_list = flecsi_get_handle(
					     mesh_handle, hydro, remap_data, real_t, dense, 0);
    if ( cell_field_list.index_space == burton_mesh_t::index_spaces_t::cells) {
      
      *data = (var_map.at(var_name));
      return;

    }

    // now check nodes
    auto nodal_field_list = flecsi_get_handle(
					      mesh_handle, hydro, remap_data, real_t, dense, 0);

    if ( nodal_field_list.index_space == burton_mesh_t::index_spaces_t::vertices) {
      return;
    }
    // if we got here, there is something wrong
    throw_runtime_error( "Could not find variable to ReMAP!" );
  }

  void check_map(std::string const & var_name) const {
    auto search = var_map.find(var_name);
    if (search != var_map.end()) {
	    std::cout << "Found " << search->first << " " << search->second << '\n';
    } else {
	    std::cout << "Not found\n";
    }
  }


  //! \brief Get pointer to read-only scalar cell data for a particular material
  //! \param[in] var_name The string name of the data field
  //! \param[in] matid   Index (not unique identifier) of the material
  //! \param[out] data   vector containing the values corresponding to cells in the material
  void mat_get_celldata( std::string const& var_name, int matid,
    double const **data) const 
  {}

  // TEMPORARY: UNTIL THIS GETS TEMPLATED ON TYPE OF DATA
  void mat_get_celldata( std::string const& var_name, int matid, 
      Wonton::Point<2> const **data) const
  {}
  void mat_get_celldata(std::string const& var_name, int matid,
      Wonton::Point<3> const **data) const
  {}


  //! \brief Get pointer to read-write scalar data for a particular material
  //! \param[in] on_what The entity type on which to get the data
  //! \param[in] var_name The string name of the data field
  //! \param[in] matid   Index (not unique identifier) of the material
  //! \param[out] data   vector containing the values corresponding to cells in the material
  void mat_get_celldata(std::string const& var_name, int matid, double **data)
  {}

  //! \brief Get a pointer to data from the state manager with a given
  //! variable @c name and on @c on_what mesh entities.
  //! \param[in] on_what The Entity_kind (e.g. CELL) on which the data lives.
  //! \param[in] name The name of the variable.
  //! \param data A @c pointer to the const data array.  If the requested
  //! data is not found in the state manager, a @c nullptr is returned.
  void mesh_add_data(entity_kind_t on_what, std::string const& name,
      double const **data) const
  {}

  //! \brief Add a scalar multi-valued data field on cells and initialize its
  //! material data to a single value
  //! \param[in] var_name The name of the data field
  //! \param[in] value Initialize with this value

  //! The 2D array will be read and values copied according to which materials
  //! are contained in which cells. If a material+cell combination is not active
  //! providing a value for the array will have no effect. 
  void mat_add_celldata(std::string const& var_name, double value) {
  }


  //! \brief Add a scalar multi-valued data field on cells and initialize its
  //! material data according to a 2D array
  //! \param[in] var_name The name of the data field
  //! \param[in] layout  Whether 2D array is laid out with first index being
  //! the cell (CELL_CENRIC) or material (MATERIAL CENTRIC)
  //! \param[in] value Initialize with this value
  //!
  //! The 2D array will be read and values copied according to which
  //! materials are contained in which cells. If a material+cell
  //! combination is not active providing a value for the array will
  //! have no effect.

  void mat_add_celldata( std::string const& var_name, 
      double const * const *values = nullptr, 
      Portage::Data_layout layout = Portage::Data_layout::MATERIAL_CENTRIC )
  {}


  //! \brief Add a scalar multi-valued data field on cells and add
  //! data to one of its materials
  //! \param[in] var_name The name of the data field
  //! \param[in] matid  Index of material in the problem
  //! \param[in] layout Data layout - 
  //! \param[in] values Initialize with this array of values
  //!
  //! Subsequent calls to this function with the same name will find the added
  //! field and just add the data.

  void mat_add_celldata(std::string const& var_name, int matid,
      double const * values)
  {}


  //! \brief Add a scalar multi-valued data field on cells and initialize one
  //! of its material data to a uniform value
  //! \param[in] var_name The name of the data field
  //! \param[in] matid Index of material in the problem
  //! \param[in] value Initialize with this value
  //! Subsequent calls to this function with the same name will find the added
  //! field and just add the data.
  template <typename T>
  void mat_add_celldata(std::string const& var_name, int matid, T const * const value)
  {}

  //! \brief Add cells to material (or add material to cells)
  //! \param[in] matid  Material ID
  //! \param[in] newcells Vector of new cells in material
  void mat_add_cells(int matid, std::vector<int> const& newcells)
  {}


  //! \brief Remove cells from material (or remove material from cells)
  //! \param[in] matid  Material ID
  //! \param[in] matcells Vector of to be removed cells
  void mat_rem_cells(int matid, std::vector<int> const& delcells)
  {}


  //! \brief Add a material to state
  //! \param[in] matname  Name of material
  //! \param[in] matcells Cells containing the material
  void add_material(std::string const& matname, std::vector<int> const& matcells)
  {}

  int get_data_size(entity_kind_t on_what, std::string const& var_name) const
  {
    return mesh_->num_cells();
    //    return (var_map_copy.at(var_name)).size();
    //    throw_runtime_error( "get_data_size not implemented yet!" );
  }



  //!
  //! @brief Get the data type of the given field
  //! @param[in] var_name The string name of the data field
  //! @return A reference to the type_info struct for the field's data type
  const std::type_info& get_data_type(std::string const& var_name) const {
    return typeid(double);  // thats the only type we can represent
  }
  
 private:
  
  //! \brief the flecsi mesh pointer
  mesh_t * mesh_ = nullptr;
  
 };  // Flecsi_State_Wrapper
 
} // namespace burton
} // namespace flecsi_sp
