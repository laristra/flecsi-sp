/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef felcsisp_minimal_mesh_h
#define felcsisp_minimal_mesh_h

#include "flecsi-sp/geometry/point.h"
#include "flecsi-sp/minimal/minimal_types.h"

///
// \file minimal_mesh.h
// \authors bergen
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace sp {

///
// Use this namespace to expose enumerations and types.
///
namespace minimal {

  enum minimal_index_spaces_t : size_t {
    vertices,
    edges,
    faces,
    cells
  }; // enum minimal_index_spaces_t

  enum minimal_cell_index_spaces_t : size_t {
    interior,
    boundary
  }; // enum minimal_cell_index_spaces_t

} // namespace minimal

///
// \class minimal_mesh_t minimal_mesh.h
// \brief minimal_mesh_t provides...
///
class minimal_mesh_t
  : public topology::mesh_topology_t<minimal_types_t>
{
public:

  using base_t = topology::mesh_topology_t<minimal_types_t>;

  static constexpr size_t dimension = minimal_config_t::num_dimensions;
  
  using vertex_t = minimal_types_t::vertex_t;
  using cell_t = minimal_types_t::cell_t;

  using point_t = minimal_config_t::point_t;

  /// Default constructor
  minimal_mesh_t()
    : base_t() {}

  /// Copy constructor (disabled)
  minimal_mesh_t(const minimal_mesh_t &) = delete;

  /// Assignment operator (disabled)
  minimal_mesh_t & operator = (const minimal_mesh_t &) = delete;

  /// Destructor
   ~minimal_mesh_t() {}

  ///
  // Initialize the mesh.
  ///
  void
  init()
  {
    // Initialize domain 0 of the mesh topology.
    base_t::init<0>();

    interior_cells_ =
      base_t::entities<dimension, 0>().filter(is_interior);
    boundary_cells_ =
      base_t::entities<dimension, 0>().filter(is_domain_boundary);
  } // init

  ///
  // Add a vertex to the mesh topology.
  ///
  vertex_t *
  make_vertex(
    const point_t & pos
  )
  {
    auto v = base_t::make<vertex_t>(*this, pos);
    base_t::add_entity<0, 0>(v);
    return v;
  } // make_vertex

  ///
  // Add a cell to the mesh topology
  ///
  cell_t *
  make_cell(
    const std::initializer_list<vertex_t *> & vertices,
    cell_type_t type
  )
  {
    auto c = base_t::make<cell_t>(*this, type);
    base_t::add_entity<dimension, 0>(c);
    base_t::init_entity<0, dimension, 0>(c, vertices);
    return c;
  } // make_cell

  ///
  //
  ///
  size_t
  indices(
    size_t index_space_id
  ) override
  {
    switch(index_space_id) {
      case minimal::vertices:
        return base_t::num_entities(0);
      case minimal::cells:
        return base_t::num_entities(dimension);
      default:
        assert(false && "unknown index space");
    } // switch
  } // indices

  ///
  //
  ///
  auto
  vertices()
  {
    return base_t::entities<0, 0>();
  } // vertices

  ///
  //
  ///
  template<
    typename E
  >
  auto
  vertices(
    E * e
  )
  {
    return base_t::entities<0, 0>(e);
  } // vertices

  ///
  //
  ///
  auto
  cells()
  {
    return base_t::entities<dimension, 0>();
  } // cells

  ///
  //
  ///
  template<
    typename E
  >
  auto
  cells(E * e)
  {
    return base_t::entities<dimension, 0>(e);
  } // cells

  auto
  cells(
    size_t is
  )
  {
    switch(is) {
      case minimal::interior:
        return interior_cells_;
      case minimal::boundary:
        return boundary_cells_;
    } // switch
  } // cells

private:

  ///
  // Predicate function to create index space for accessing
  // domain boundary cells.
  ///
  static
  bool
  is_domain_boundary(
    cell_t * c
  )
  {
    return c->type() == cell_type_t::domain_boundary;
  } // is_domain_boundary

  ///
  // Predicate function to create index space for accessing
  // interior cells.
  ///
  static
  bool
  is_interior(
    cell_t * c
  )
  {
    return !is_domain_boundary(c);
  } // is_interior

  topology::index_space<
    topology::domain_entity<0, cell_t>, false, true, false> interior_cells_;
  topology::index_space<
    topology::domain_entity<0, cell_t>, false, true, false> boundary_cells_;

}; // class minimal_mesh_t

} // namespace sp
} // namespace flecsi

#endif // felcsisp_minimal_mesh_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
