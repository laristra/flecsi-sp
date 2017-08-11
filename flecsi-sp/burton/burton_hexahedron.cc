/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsi-sp/burton/burton_mesh_topology.h"
#include "flecsi-sp/burton/burton_hexahedron.h"


namespace flecsi {
namespace sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
// 3D hexahedron
////////////////////////////////////////////////////////////////////////////////

// the centroid
burton_hexahedron_t::point_t burton_hexahedron_t::centroid() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geometry::shapes::hexahedron::centroid( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates(),
    vs[4]->coordinates(), vs[5]->coordinates(),
    vs[6]->coordinates(), vs[7]->coordinates() );
}

// the midpoint
burton_hexahedron_t::point_t burton_hexahedron_t::midpoint() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geometry::shapes::hexahedron::midpoint( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates(),
    vs[4]->coordinates(), vs[5]->coordinates(),
    vs[6]->coordinates(), vs[7]->coordinates() );
}


// the area of the cell
burton_hexahedron_t::real_t burton_hexahedron_t::volume() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geometry::shapes::hexahedron::volume( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates(),
    vs[4]->coordinates(), vs[5]->coordinates(),
    vs[6]->coordinates(), vs[7]->coordinates() );
}

} // namespace
} // namespace
} // namespace
