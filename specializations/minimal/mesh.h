/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef felcsisp_minimal_mesh_h
#define felcsisp_minimal_mesh_h

#include "flecsi-sp/geometry/point.h"
#include "flecsi-sp/minimal/types.h"

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

  enum minimal_edge_index_spaces_t : size_t {
    interior_edge,
    boundary_reflec,
    boundary_src,
    boundary_vacuum
  };//enum minimal_edge_index_space_t
  
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
  using edge_t = minimal_types_t::edge_t;
#if FLECSI_MESH_DIMENSION == 3
  using face_t = minimal_types_t::face_t;
#endif  
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

    // use a predicate function to create interior edge/faces.
    interior_edges_ =
      base_t::entities<dimension-1,0>().filter(is_interior_edge);
    // similar filter for reflec bc
    reflec_bc_edges_ =
      base_t::entities<dimension-1,0>().filter(is_reflec_bc);
    // for src bc.
    src_bc_edges_ =
      base_t::entities<dimension-1,0>().filter(is_src_bc);
    // for vacuum
    vacuum_bc_edges_=
      base_t::entities<dimension-1,0>().filter(is_vacuum_bc);
    
    // Use a predicate function to create the interior cells
    // index space
    interior_cells_ =
      base_t::entities<dimension, 0>().filter(is_interior);

    // Use a predicate function to create the domain boundary cells
    // index space
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
  // add edge to the mesh topology.
  ///
  edge_t *
  make_edge(
    const std::initializer_list<vertex_t *> &vertices,
    edge_type_t type
  )
  {
    auto e = base_t::make<edge_t>(*this,type);
    base_t::add_entity<dimension-1,0>(e);
    base_t::init_entity<0,dimension-1,0>(e, vertices);
    return e;
  }
  
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
<<<<<<< HEAD
  ) const override
=======
  )
  const
  override
>>>>>>> upstream/master
  {
    switch(index_space_id) {
      case minimal::vertices:
        return base_t::num_entities(vertex_t::dimension);
      case minimal::edges:
        return base_t::num_entities(edge_t::dimension);
      case minimal::cells:
        return base_t::num_entities(cell_t::dimension);
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

  auto
  edges()
  {
    return base_t::entities<dimension-1,0>();
  }
  
  template<
    typename E
  >
  auto
  edges(
    E * e
  )
  {
    return base_t::entities<dimension-1,0>(e);
  }

  auto
  edges(
    size_t is
  )
  {
    switch(is) {
      case minimal::interior_edge:
        return interior_edges_;
      case minimal::boundary_reflec:
        return reflec_bc_edges_;
      case minimal::boundary_src:
        return src_bc_edges_;
      case minimal::boundary_vacuum:
        return vacuum_bc_edges_;
      default:
        assert(false && "unknown index space");
    } // switch
  } // edges

  
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
      default:
        assert(false && "unknown index space");
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

  ///
  //predicate functions to create interior/boundary faces/edges.
  static
  bool
  is_reflec_bc(
    edge_t *e
  )
  {
    return e->type() == edge_type_t::reflec_bc;
  }//is_reflec_bc

  static
  bool
  is_src_bc(
    edge_t *e
  )
  {
    return e->type() == edge_type_t::src_bc;
  }

  static
  bool
  is_vacuum_bc(
    edge_t *e
  )
  {
    return e->type() == edge_type_t::vacuum_bc;
  }
  
  static
  bool
  is_interior_edge(
    edge_t *e
  )
  {
    return e->type() == edge_type_t::interior_edge;
  }
  
  
  topology::index_space<
    topology::domain_entity<0, cell_t>, false, true, false> interior_cells_;
  topology::index_space<
    topology::domain_entity<0, cell_t>, false, true, false> boundary_cells_;

  
  topology::index_space<
    topology::domain_entity<0, edge_t>, false, true, false> interior_edges_;
  topology::index_space<
    topology::domain_entity<0, edge_t>, false, true, false> reflec_bc_edges_;
  topology::index_space<
    topology::domain_entity<0, edge_t>, false, true, false> src_bc_edges_;
  topology::index_space<
    topology::domain_entity<0, edge_t>, false, true, false> vacuum_bc_edges_;
}; // class minimal_mesh_t

} // namespace sp
} // namespace flecsi

#endif // felcsisp_minimal_mesh_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
