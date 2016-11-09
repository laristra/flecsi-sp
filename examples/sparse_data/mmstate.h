/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef mmstate_h
#define mmstate_h

#include "config.h"

static constexpr size_t num_mats = 10;

void driver(int argc, char ** argv) {

	minimal_mesh_t mesh;

  // Initialize the mesh
	init_mesh(mesh);

  // Register sparse material variable 'mats'
	register_data(mesh, solver, mats, double, sparse, 1, vertices, num_mats);

  // Initialize material data
  {
  auto ap = get_mutator(mesh, solver, mats, double, sparse, 0, 3);

  for(auto v: mesh.vertices()) {
    double value = 1.0;
    std::cout << "vertex: " << v->id<0>() << std::endl;
    for(size_t i(0); i<3; ++i) {
      ap(v, i) = value + 1.0;
    } // for
  } // for

  } // scope

  auto ma = get_accessor(mesh, solver, mats, double, sparse, 0);

  // Algorithm 1: Iterate through all indices with non-zero entries
  for(auto m: ma.indices()) {
    std::cout << "index: " << m << std::endl;
    std::cout << "value: " << ma(m) << std::endl;
  } // for

} // driver

#endif // mmstate_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
