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
  flecsi_sp::utils::client_handle_r<mesh_t> mesh
) {
  
  mesh.update_geometry();

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
    auto f = flecsi_execute_task( 
      update_geometry, 
      flecsi_sp::burton,
      index, 
      handle_
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

};

} // namespace burton
} // namespace flecsi-sp
