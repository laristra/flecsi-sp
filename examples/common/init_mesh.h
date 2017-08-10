/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef init_mesh_h
#define init_mesh_h

#include <flecsi-sp/minimal/mesh.h>

//using namespace flecsi;
//using namespace flecsi::data;
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
        bool is_domain_boundary = i==0 || j==0 || i==(N-1) || j==(N-1);

        m.make_cell({
          vs[ i    + ( j    * width)],
          vs[(i+1) + ( j    * width)],
          vs[(i+1) + ((j+1) * width)],
          vs[ i    + ((j+1) * width)]
          },
          is_domain_boundary ? cell_type_t::domain_boundary :
          cell_type_t::unknown);
      } // for
    } // for

    m.init();

} // init_mesh

#endif // init_mesh_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
