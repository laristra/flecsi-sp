/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef mmstate_h
#define mmstate_h

#include <set>

#include "config.h"

///
// The driver task is the top-level user task for a FleCSI application. The
// driver task is executed after runtime and specailization initialization
// have been executed. Semantically, the driver task is sequential. However,
// when built against a parallel backend to FleCSI, it will execute in a
// distributed-memory mode. The fact of parallel execution is not exposed to
// the user.
///
void driver(int argc, char ** argv) {

  // Mesh object
	minimal_mesh_t mesh;

  // Initialize the mesh
	init_mesh(mesh);

  ///
  //
  ///
	register_data(mesh, solver, pressure, double, dense, 1, cells);

  // iterate over the interior cells
  {
  auto cell = get_accessor(mesh, solver, pressure, double, dense, 0);

  for(auto c: mesh.cells(boundary)) {
    cell[c] = 1.0;
  } // for

  } // scope

  // iterate over all cells
  {
  auto cell = get_accessor(mesh, solver, pressure, double, dense, 0);

  for(auto c: mesh.cells()) {
    std::cout << c->id<0>() << " " << cell[c] << std::endl;
  } // for

  } // scope

} // driver

#endif // mmstate_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
