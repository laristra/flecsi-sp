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

// Files from flecsi-sp
#include <flecsi-sp/minimal/minimal_mesh.h>

using namespace flecsi;
using namespace flecsi::data;
using namespace flecsi::sp;
using namespace flecsi::sp::minimal;
using vertex_t = minimal_mesh_t::vertex_t;

static constexpr size_t N = 8;

///
// Initialize a basic mesh.
///
void init_mesh(minimal_mesh_t & m) {
  std::vector<vertex_t *> vs;

    for(size_t j(0); j<N+1; ++j) {
      for(size_t i(0); i<N+1; ++i) {
        vs.push_back(m.make_vertex({double(i), double(j)}));
      } // for
    } // for

    size_t width = N+1;

    for(size_t j(0); j<N; ++j) {
      for(size_t i(0); i<N; ++i) {
        m.make_cell({
          vs[ i    + ( j    * width)],
          vs[(i+1) + ( j    * width)],
          vs[(i+1) + ((j+1) * width)],
          vs[ i    + ((j+1) * width)]
        });
      } // for
    } // for

    m.init();

} // init_mesh

void driver(int argc, char ** argv) {
	minimal_mesh_t m;

	// Initialize the mesh
	//
	// For real programs, mesh initialization will be handled through
	// an I/O object.
	init_mesh(m);

	register_data(m, sovler, unknowns, double, dense, 2, vertices);

  auto u = get_accessor(m, solver, unknowns, double, dense, 0);

  for(auto v: m.vertices()) {
    u[v] = 0.0;
  } // for

} // driver

#endif // driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
