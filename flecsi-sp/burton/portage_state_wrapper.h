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
class portage_state_wrapper_t {

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
  std::map < std::string, std::vector<int> > mat_map;
  std::map < std::string, entity_kind_t> entity_map;
  std::map < std::string, field_type_t> type_map;
  std::vector<std::vector<int>> mat_cells_;
  std::vector<std::vector<int>> cell_mats_;
  std::vector<int> mat_offsets_;
  std::vector<int> mat_ids_;
  std::vector<int> offsets_;
  std::vector<int> index_vals_;
  std::map < int, std::vector<double>> volfrac_map;
  int number_materials=0;

  //============================================================================
  // Constructors
  //============================================================================

  //!  \brief Default constructor.
  //!  \param[in] mesh The minimum coordinates of the domain.
  //!  \param[in] mesh The minimum coordinates of the domain.
  explicit portage_state_wrapper_t(mesh_t & mesh) : mesh_(&mesh)
    {}

  //! Default constructor deleted
  portage_state_wrapper_t() = default;

  //! Default copy constructor
  portage_state_wrapper_t(const portage_state_wrapper_t &) = default;

  //! Default assignment operator
  portage_state_wrapper_t & operator=(const portage_state_wrapper_t &) = default;

  //============================================================================
  // Public Members
  //============================================================================

  /* //! \brief Add a material variable of entity type (cell)  */
  /* //! field that needs to be remapped to the variable map */
  /* void add_mat_cell_field( std::string var_name, int matid, double* data) { */

  /*   //    var_name = var_name + std::to_string(matid); */
  /*   var_map.insert( std::pair <std::string, double*> (var_name, data)); */
  /*   type_map.insert( std::pair < std::string, std::string> (var_name, "CELL")); */

  /* } */

  //! \brief Add a variable of entity type (cell) 
  //! field that needs to be remapped to the variable map
  void add_cell_field( std::string var_name, double* data, 
                       entity_kind_t entity=entity_kind_t::CELL,
                       field_type_t type=field_type_t::MESH_FIELD) {
    var_map.insert( std::pair <std::string, double*> (var_name, data));
    entity_map.insert( std::pair < std::string, entity_kind_t> 
                       (var_name, entity));
    type_map.insert( std::pair < std::string, field_type_t> (var_name, type));
  }

  void add_field_type( std::string var_name,
                       entity_kind_t entity=entity_kind_t::CELL,
                       field_type_t type=field_type_t::MESH_FIELD) {
    entity_map.insert( std::pair < std::string, entity_kind_t> (var_name, entity));
    type_map.insert( std::pair < std::string, field_type_t> (var_name, type));
  }

  //! \brief Set the offset for the beginning of a material in a state
  void set_mat_offsets(std::vector<int> mat_offsets){
    mat_offsets_ = mat_offsets;
  }

  //! \brief Number of materials in problem
  int num_materials() const {
    return number_materials;
  }

  //! \brief Set the number of materials in problem
  void num_materials(int num_mats){
    number_materials = num_mats;
  }

  //! \brief Name of material
  std::string material_name(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return std::to_string(matid);
  }

  //! \brief Get number of cells containing a particular material
  //! \param matid    Index of material (0, num_materials()-1)
  //! \return         Number of cells containing material 'matid'
  int mat_get_num_cells(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return mat_cells_[matid].size();
  }

  //! \brief Get cell indices containing a particular material
  //! \param matid    Index of material (0, num_materials()-1)
  //! \param matcells Cells containing material 'matid'
  void mat_get_cells(int matid, std::vector<int> *matcells) const{
    assert(matid >= 0 && matid < num_materials());
    matcells->clear();
    for (auto i: mat_cells_[matid]){
      matcells->push_back(i);
    }
  }

  //! \brief Set cell indices for the entire domain
  //! \param matcells  a vector of vectors for the cells with a given material
  void mat_set_cells(std::vector<std::vector<int>> mat_cells) {
    mat_cells_ = mat_cells;
  }

  //! \brief Set materials within a cell for the entire domain
  //! \param matcells  a vector of vectors for the materials in a given cell
  void cell_set_mats(std::vector<std::vector<int>> cell_mats) {
    cell_mats_ = cell_mats;
  }

  //! \brief Get number of materials contained in a cell
  //! \param cellid  Index of cell in mesh
  //! \return        Number of materials in cell
  int cell_get_num_mats(int cellid) const {
    return cell_mats_[cellid].size();
  }

  //! \brief Get the IDs of materials in a cell
  //! \param cellid    Index of cell in mesh
  //! \param cellmats  Indices of materials in cell
  void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
    cellmats->clear();
    for (auto element: cell_mats_[cellid])
      cellmats->push_back(element);
  }

  //! \brief Get the local index of mesh cell in material cell list
  //! \param meshcell    Mesh cell ID
  //! \param matid       Material ID
  //! \return             Local cell index in material cell list
  int cell_index_in_material(int meshcell, int matid) const {
    auto offset_val_0 = offsets_[meshcell];
    auto offset_val_1 = offsets_[meshcell+1];
    auto mat_iter = std::find(mat_ids_.begin() + offset_val_0,
                              mat_ids_.begin() + offset_val_1,
                              matid);

    if (mat_iter == mat_ids_.end()){
      std::cerr << "MAT_ITER not found " << meshcell << " " 
                << matid << std::endl;
    }

    int mat_iter_offset = std::distance(mat_ids_.begin() + offset_val_0,
                                        mat_iter);
    int index = index_vals_[offset_val_0 + mat_iter_offset];
    return index;
  }

  //! \brief Support function for cell_index_in_material
  //! \param offsets  the offset value for a given cell into the state data
  void add_offsets(std::vector<int> offsets){
    offsets_ = offsets;
  }

  //! \brief Support function for cell_index_in_material
  //! \param mat_ids  the material corresponding to a specific offset value
  void add_mat_ids(std::vector<int> mat_ids){
    mat_ids_ = mat_ids;
  }

  //! \brief Support function for cell_index_in_material
  //! \para index_vals the cell id for a material w.r.t a specific offset value
  void add_index_vals(std::vector<int> index_vals){
    index_vals_ = index_vals;
  }
  
  //! \brief Type of field (MESH_FIELD or MULTIMATERIAL_FIELD)
  //! \param[in] onwhat   Entity_kind that field is defined on
  //! \param[in] varname  Name of field
  //! \return             Field_type
  field_type_t field_type(entity_kind_t on_what, std::string const& var_name)
      const {
  
    auto mat_name = var_name + std::to_string(this->number_materials-1);
    if ( type_map.find(var_name) != type_map.end() ){
      return type_map.at(var_name);
    } else if ( type_map.find(mat_name) != type_map.end() ){
      return type_map.at(mat_name);
    } else {
      std::cerr << " Could not find state variable field type " << var_name
                << " or " << mat_name << std::endl;
      //ASSUME MULTIMATERIAL_FIELD IF NOT FOUND (NO UNKNOWN_TYPE EXISTS)
      return field_type_t::MULTIMATERIAL_FIELD;
    }
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
    auto mat_name = var_name + std::to_string(this->number_materials-1);   
    if ( entity_map.find(var_name) != entity_map.end()){
      return entity_map.at(var_name);
    } else if (entity_map.find(mat_name) != entity_map.end()) {
      return entity_map.at(mat_name);
    } else {
      std::cerr << " Could not find state variable entity kind " << var_name << " "
                << "or " << mat_name <<  std::endl;
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
  }

  void check_map(std::string const & var_name) const {
    auto search = var_map.find(var_name);
    if (search != var_map.end()) {
	    std::cout << "Found " << search->first << " " << search->second << '\n';
    } else {
	    std::cout << var_name << " Not found\n" ;
    }
  }

  void check_map_state(std::string const & var_name, int size) const {
    auto search = var_map.find(var_name);
    if (search != var_map.end()) {
	    std::cout << "Found " << search->first << " " << search->second << '\n';
	    for ( int i=0; i < size; ++i ) {
	      std::cout << (search->second)[i] << " " << &(search->second)[i] << std::endl;
	    }
    } else {
	    std::cout << var_name << " Not found\n";
    }
  }

  void check_mat_map(std::string const & var_name) const {
    auto search = mat_map.find(var_name);
    if (search != mat_map.end()) {
      std::cout << "Found " << search->first << " " << (search->second)[0] << '\n';
    } else {
	    std::cout << var_name << " Not found\n";
    }
  }


  //! \brief Get pointer to read-only scalar cell data for a particular material
  //! \param[in] var_name The string name of the data field
  //! \param[in] matid   Index (not unique identifier) of the material
  //! \param[out] data   vector containing the values corresponding to cells in the material
  void mat_get_celldata( std::string const& var_name, int matid,
    double const **data) const 
  {
    auto name = var_name + std::to_string(matid);   
    if ( var_map.find(name) != var_map.end() ){
      *data = var_map.at(name);
    } else if (var_map.find(var_name) != var_map.end()) {
      *data = var_map.at(var_name) + mat_offsets_[matid];
    } else {
      std::cerr << " Could not find state variable " << name << std::endl;
      *data = nullptr;
    }
  }

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
  {
    auto name = var_name + std::to_string(matid);   
    if ( var_map.find(name) != var_map.end() ){
      *data = var_map.at(name);
    } else if (var_map.find(var_name) != var_map.end()) {
      *data = var_map.at(var_name);
    } else {
      std::cerr << " Could not find state variable " << name << std::endl;
      *data = nullptr;
    }
  }

  //! \brief Get the volume fractions for a particular material
  //! \param[in] matid  Material ID
  //! \param[in] data   vector in which the volume fractions are placed
  void mat_get_volfracs(int matid, double **data)
  {
    if ( volfrac_map.find(matid) != volfrac_map.end() ){
      *data = (volfrac_map.at(matid)).data();
    } else {
      std::cerr << " Could not find state variable " << matid << std::endl;
      *data = nullptr;
    }
  }

  //! \brief Get a pointer to data from the state manager with a given
  //! variable @c name and on @c on_what mesh entities.
  //! \param[in] on_what The Entity_kind (e.g. CELL) on which the data lives.
  //! \param[in] name The name of the variable.
  //! \param data A @c pointer to the const data array.  If the requested
  //! data is not found in the state manager, a @c nullptr is returned.
  void mesh_add_data(entity_kind_t on_what, std::string const& name,
      double const **data) const
  {
    THROW_RUNTIME_ERROR("mesh_add_data not implemented");
  }

  //! \brief Add a scalar multi-valued data field on cells and initialize its
  //! material data to a single value
  //! \param[in] var_name The name of the data field
  //! \param[in] value Initialize with this value

  //! The 2D array will be read and values copied according to which materials
  //! are contained in which cells. If a material+cell combination is not active
  //! providing a value for the array will have no effect. 
  /* void mat_add_celldata(std::string const& var_name, double value) */
  /* { */
  /*   if ( var_map.find(var_name) == var_map.end() ){ */
  /*     this->add_cell_field( var_name, &value); */
  /*   } else { */
  /*     var_map[var_name] = &value; */
  /*     //      std::cerr << " Attempted to add duplicate variable vector. Ignoring" << std::endl;  */
  /*   } */
  /* } */

  /* void mat_add_celldata(std::string const& var_name, double& value) */
  /* { */
  /*   if ( var_map.find(var_name) == var_map.end() ){ */
  /*     this->add_cell_field( var_name, value); */
  /*   } else { */
  /*     var_map[var_name] = value; */
  /*     // std::cerr << " Attempted to add duplicate variable vector. Ignoring" << std::endl;  */
  /*   } */
  /* } */

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
  {
    auto name = var_name + std::to_string(matid);
    auto vals = const_cast<double*>(values);	
    auto size_max = mat_cells_[matid].size();
    std::vector<double> volfracs_;
    if ( var_map.find(name) == var_map.end() ){
      this->add_cell_field(name, vals);
    } else {
      var_map[name] = vals;
      for (size_t i=0; i < size_max; ++i){
        volfracs_.push_back(*(vals+i));
      }
      volfrac_map.insert( std::pair <int, std::vector<double>> (matid, volfracs_));
    }
  }
  
  /* //! \brief Add a scalar multi-valued data field on cells and initialize one */
  /* //! of its material data to a uniform value */
  /* //! \param[in] var_name The name of the data field */
  /* //! \param[in] matid Index of material in the problem */
  /* //! \param[in] value Initialize with this value */
  /* //! Subsequent calls to this function with the same name will find the added */
  /* //! field and just add the data. */
  template <typename T>
  void mat_add_celldata(std::string const& var_name, int matid, T const * const value)
  {}

  //! \brief Add cells to material (or add material to cells)
  //! \param[in] matid  Material ID
  //! \param[in] newcells Vector of new cells in material
  void mat_add_cells(int matid, std::vector<int> /* const */& newcells)
  {
    mat_cells_[matid].clear();
    for (auto element: newcells){
      mat_cells_[matid].push_back(element);
    }
  }


  //! \brief Remove cells from material (or remove material from cells)
  //! \param[in] matid  Material ID
  //! \param[in] matcells Vector of to be removed cells
  void mat_rem_cells(int matid, std::vector<int> const& delcells)
  {
    std::cerr << "mat_rem_cells --- NOT IMPLEMENTED" << std::endl;
  }


  //! \brief Add a material to state
  //! \param[in] matname  Name of material
  //! \param[in] matcells Cells containing the material
  void add_material(std::string const& matname, std::vector<int> const& matcells)
  {
    if (mat_map.find(matname) == mat_map.end()){
      mat_map.insert( std::pair<std::string, std::vector<int>> (matname, matcells));
    } else {
      std::cerr << "Attempted to add duplicate material. Ignoring" << std::endl;
    }
  }

  //! Get the size of a particular quantity
  //! WARNING: this may not be the correct value for multi-material quantities
  int get_data_size(entity_kind_t on_what, std::string const& var_name) const
  {
    return mesh_->num_cells();
  }

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
