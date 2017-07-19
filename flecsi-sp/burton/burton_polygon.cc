/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

// user includes
#include "flecsi-sp/burton/burton_mesh_topology.h"
#include "flecsi-sp/burton/burton_polygon.h"


namespace flecsi {
namespace sp {
namespace burton {

// some type aliases
using burton_2d_polygon_t = burton_polygon_t<2>;
using burton_3d_polygon_t = burton_polygon_t<3>;

////////////////////////////////////////////////////////////////////////////////
// 2D polygon
////////////////////////////////////////////////////////////////////////////////


// the centroid
burton_2d_polygon_t::point_t burton_2d_polygon_t::centroid() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::centroid( coords );
}

// the midpoint
burton_2d_polygon_t::point_t burton_2d_polygon_t::midpoint() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::midpoint( coords );
}


// the area of the cell
burton_2d_polygon_t::real_t burton_2d_polygon_t::area() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::area( coords );
}


////////////////////////////////////////////////////////////////////////////////
// 3D polygon
////////////////////////////////////////////////////////////////////////////////


// the centroid
burton_3d_polygon_t::point_t burton_3d_polygon_t::centroid() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::centroid( coords );
}

// the midpoint
burton_3d_polygon_t::point_t burton_3d_polygon_t::midpoint() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::midpoint( coords );
}

// the normal
burton_3d_polygon_t::vector_t burton_3d_polygon_t::normal() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::normal( coords );
}

// the area of the cell
burton_3d_polygon_t::real_t burton_3d_polygon_t::area() const
{
  auto coords = coordinates();
  return geometry::shapes::polygon<num_dimensions>::area( coords );
}


} // namespace
} // namespace
} // namespace
