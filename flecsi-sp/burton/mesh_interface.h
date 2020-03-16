/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tasks related to solving full hydro solutions.
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <ristra/initialization/input.h>

#include <flecsi/execution/context.h>
#include <flecsi/execution/execution.h>

#include <flecsi-sp/utils/types.h>
#include <flecsi-sp/burton/burton_mesh.h>

namespace flecsi_sp {
namespace burton {


// mesh and some underlying data types
using mesh_t = flecsi_sp::burton::burton_mesh_t;
using real_t = mesh_t::real_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief Update mesh geometry
//!
//! \param [in] mesh the mesh object
////////////////////////////////////////////////////////////////////////////////
void update_geometry( 
  flecsi_sp::utils::client_handle_r<mesh_t> mesh,
  int alpha_,
  int axis_
) {
  
  mesh.update_geometry(alpha_, axis_);

}

////////////////////////////////////////////////////////////////////////////////
//! \brief Check if the mesh is correct
//!
//! \param [in] mesh the mesh object
////////////////////////////////////////////////////////////////////////////////
void validate_mesh( 
  flecsi_sp::utils::client_handle_r<mesh_t> mesh
) {
  
  mesh.is_valid();

}

////////////////////////////////////////////////////////////////////////////////
// TASK REGISTRATION
////////////////////////////////////////////////////////////////////////////////

flecsi_register_task(update_geometry, flecsi_sp::burton, loc, index|flecsi::leaf);
flecsi_register_task(validate_mesh, flecsi_sp::burton, loc, index|flecsi::leaf);


///////////////////////////////////////////////////////////////////////////////
// Clent Registration happens here because the specialization initialization
// needs to know which mesh to access
///////////////////////////////////////////////////////////////////////////////
flecsi_register_data_client(mesh_t, meshes, mesh0);


////////////////////////////////////////////////////////////////////////////////
// MESH INTERFACE
////////////////////////////////////////////////////////////////////////////////

class mesh_interface_t {

  using mesh_handle_t = decltype(flecsi_get_client_handle(mesh_t, meshes, mesh0));

  mesh_handle_t handle_;

public:

  mesh_interface_t()
    : handle_( flecsi_get_client_handle(mesh_t, meshes, mesh0) )
  {}

  void setup(size_t time_cnt, real_t soln_time, ristra::initialization::input_t & input)
  {
   
    auto mesh_input = input["mesh"];
    int alpha_= 0;
    int axis_ = 1;
    if (!mesh_input.empty()){
      auto geom = mesh_input["geom"].as<std::string>();
      auto alpha_axis_vec = initialize_geometry(geom);
      alpha_ = alpha_axis_vec[0];
      axis_ = alpha_axis_vec[1];
    }

    auto f = flecsi_execute_task( 
      update_geometry, 
      flecsi_sp::burton,
      index, 
      handle_,
      alpha_,
      axis_
    );
    f.wait();  // DONT GO FORWARD UNTIL DONE!
  
    flecsi_execute_task( 
      validate_mesh, 
      flecsi_sp::burton,
      index, 
      handle_
  );
  }

  auto get_handle_ptr() { return &handle_; }

  //==============================================================================
  //! \brief set the geometry variables
  //! cartesian -> cartesian geometry, no assumptions of symmetry
  //! cylindrical_x -> cylindrical geometry, x-axis of symmetry
  //! cylindrical_y -> cylindrical geometry, y-axis of symmetry
  //! \param [in] geom  The chosen geometry for the problem (string)
  //==============================================================================
  static std::vector<int> initialize_geometry(std::string geom){

    // Create an enumerated list of the geometries
    //  explicitly assign enumerated values
    enum geom_enum{cartesian=0,
                   cylindrical_x = 1,
                   cylindrical_y = 2};
    // Create a map of the possible geometry strings
    //  same values as with the enumeration
    std::map<std::string, int> geom_choices_ = {
      {"cartesian", 0},
      {"cylindrical_x", 1},
      {"cylindrical_y", 2}};

    // Check if the geometry is available
    //  assign value for runtime error if geom not available
    int geom_val_ = -1;
    if (geom_choices_.find(geom) != geom_choices_.end())
      geom_val_ = geom_choices_.at(geom);
    
    // Create container and fill the geometry values
    // alpha -> zeroth value, corresponds to switch from cartesian to cylindrical
    // axis -> first value, choice of symmetry axis, irrelevant for cartesian
    std::vector<int> alpha_axis_vec;
    switch(geom_val_)
      {
      case cartesian:
        alpha_axis_vec.push_back(0); //alpha
        alpha_axis_vec.push_back(0); //axis
        break;
      case cylindrical_x:
        alpha_axis_vec.push_back(1); //alpha
        alpha_axis_vec.push_back(0); //axis
        break;
      case cylindrical_y:
        alpha_axis_vec.push_back(1); //alpha
        alpha_axis_vec.push_back(1); //axis
        break;
      default:
        THROW_RUNTIME_ERROR(
          "Unknown geometry type provided: \""<< geom << "\""
          );
        break;
      }
    return alpha_axis_vec;
  }

};

} // namespace burton
} // namespace flecsi-sp
