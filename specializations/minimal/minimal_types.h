/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_sp_minimal_types_h
#define flecsi_sp_minimal_types_h

#include <flecsi/topology/mesh_topology.h>

#include "flecsi-sp/geometry/point.h"
#include "flecsi-sp/minimal/minimal_entity_types.h"

///
// \file minimal_types.h
// \authors bergen bird
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace sp {

///
// \class minimal_types_t minimal_types.h
// \brief minimal_types_t provides...
///
struct minimal_types_t
{
  //--------------------------------------------------------------------------//
  // Define local traits to satisfy mesh_topology requirements.
  //--------------------------------------------------------------------------//

  /// The dimension of the mesh
  static constexpr size_t num_dimensions =
    minimal_config_t::num_dimensions;

  /// The number of domains
  static constexpr size_t num_domains =
    minimal_config_t::num_domains;

  //--------------------------------------------------------------------------//
  // Define basic types.
  //--------------------------------------------------------------------------//

  /// Mesh vertex type
  using vertex_t = minimal_vertex_t;

  /// Mesh edge type
  using edge_t = minimal_edge_t;

#if FLECSI_MESH_DIMENSION == 3
  /// Mesh face type
  using face_t = minimal_face_t;
#endif

  /// Mesh cell type
  using cell_t = minimal_cell_t;

  /// Convenience type
  template<size_t D>
  using domain_ = topology::domain_<D>;

  ///
  // Definitions of burton mesh entities and their domain.
  // clang-format off
  ///
  using entity_types =
      std::tuple<
        std::pair<domain_<0>, vertex_t>,
        std::pair<domain_<0>, edge_t>,
#if FLECSI_MESH_DIMENSION == 3
        std::pair<domain_<0>, face_t>,
#endif
        std::pair<domain_<0>, cell_t>
      >;

  ///
  // Connectivities are adjacencies of entities within a single domain.
  ///
  using connectivities =
    std::tuple<
      std::tuple<domain_<0>, vertex_t, cell_t>,
      std::tuple<domain_<0>, cell_t, vertex_t>
    >;
#if 0
  using connectivities =
    std::tuple<
      std::tuple<domain_<0>, vertex_t, edge_t>,
      std::tuple<domain_<0>, vertex_t, cell_t>,
      std::tuple<domain_<0>, edge_t, vertex_t>,
      std::tuple<domain_<0>, edge_t, cell_t>,
#if FLECSI_MESH_DIMENSION == 3
      std::tuple<domain_<0>, vertex_t, face_t>,
      std::tuple<domain_<0>, cell_t, face_t>,
      std::tuple<domain_<0>, face_t, cell_t>,
#endif
      std::tuple<domain_<0>, cell_t, vertex_t>,
      std::tuple<domain_<0>, cell_t, edge_t>
    >;
#endif

  ///
  // Bindings are adjacencies of entities across two domains.
  ///
  using bindings = std::tuple<>;

  //-------------------------------------------------------------------------//
  //
  //-------------------------------------------------------------------------//

  ///
  // \tparam M The topological domain.
  // \tparam D The topological dimension for which to create an entity.
  ///
  template<
    size_t M,
    size_t D
  >
  static topology::mesh_entity_base_t<num_domains> *
  create_entity(
    topology::mesh_topology_base_t * mesh,
    size_t num_vertices
  )
  {
    switch(M) {

      case 0:
      {
        switch(D) {
          case 1:
            return mesh->make<edge_t>(*mesh);

          default:
            assert(false && "invalid topological dimension");
        } // switch

        break;
      } // case

      default:
        assert(false && "invalid domain");

    } // switch
  } // create_entity

}; // class minimal_types_t

} // namespace sp
} // namespace flecsi

#endif // flecsi_sp_minimal_types_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
