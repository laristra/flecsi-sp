/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef driver_h
#define driver_h

#include <cinchlog.h>

#include "config.h"
#include "state.h"

constexpr size_t num_mats = 3;

void driver(int argc, char ** argv) {

  clog_init("all");

  // Mesh object
	minimal_mesh_t mesh;

  // Initialize the mesh
	init_mesh(mesh);

  // Register material state
  register_data(mesh, hydro, mat_cell_state, mat_state_t, sparse, 2, cells,
    num_mats);

  // Register cell state
  register_data(mesh, hydro, cell_state, cell_state_t, dense, 2, cells);

  for(auto z: mesh.cells()) {
    
  } // for

} // driver

#endif // driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
