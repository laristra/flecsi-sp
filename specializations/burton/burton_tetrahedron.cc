/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsi-sp/burton/burton_mesh_topology.h"
#include "flecsi-sp/burton/burton_tetrahedron.h"


namespace flecsi {
namespace sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
// 3D tetrahedron
////////////////////////////////////////////////////////////////////////////////

// the centroid
burton_tetrahedron_t::point_t burton_tetrahedron_t::centroid() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geometry::shapes::tetrahedron::centroid( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

// the midpoint
burton_tetrahedron_t::point_t burton_tetrahedron_t::midpoint() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geometry::shapes::tetrahedron::midpoint( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}


// the area of the cell
burton_tetrahedron_t::real_t burton_tetrahedron_t::volume() const
{
  auto msh = static_cast<const burton_3d_mesh_topology_t *>(mesh()); 
  auto vs = msh->template entities<vertex_t::dimension, vertex_t::domain>(this);
  return geometry::shapes::tetrahedron::volume( 
    vs[0]->coordinates(), vs[1]->coordinates(), 
    vs[2]->coordinates(), vs[3]->coordinates() );
}

} // namespace
} // namespace
} // namespace
