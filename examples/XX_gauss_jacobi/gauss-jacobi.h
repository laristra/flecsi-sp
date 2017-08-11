/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef driver_h
#define driver_h

///
//
///

#include <iostream>

// Files from flecsi
#include <flecsi/data/data.h>

// Files from specializations
#include <flecsi-sp/minimal/minimal_mesh.h>

#include "../common/init_mesh.h"

using namespace flecsi;
using namespace flecsi::data;

void driver(int argc, char ** argv) {
	minimal_mesh_t m;

	init_mesh(m);

	register_data(m, solver, unknowns, double, dense, 2, vertices);
	register_data(m, solver, p, double, sparse, 2, vertices, 3);

  {
  auto ap = get_mutator(m, solver, p, double, sparse, 0, 3);

  for(auto v: m.vertices()) {
    for(size_t i(0); i<3; ++i) {
      ap(v, i) = 1.0;
    } // for
  } // for

  } // scope

  auto u = get_accessor(m, solver, unknowns, double, dense, 0);

  for(auto v: m.vertices()) {
    std::cout << v->coordinates() << std::endl;
    u[v] = 0.0;
  } // for

} // driver

#endif // driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
