/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsi/data/data.h>
#include <flecsi/execution/task.h>
#include <flecsi-sp/burton/burton_mesh_topology.h>
#include <flecsi-sp/burton/burton_types.h>
#include <ristra/assertions/errors.h>
#include <ristra/compatibility/type_traits.h>
#include <ristra/utils/type_traits.h>

// system includes
#include <set>
#include <string>
#include <sstream>


namespace flecsi_sp {
namespace burton {

//! This namespace is used to expose enumerations and types.
namespace attributes {
  
////////////////////////////////////////////////////////////////////////////////
/// \brief Attributes for flecsi.
////////////////////////////////////////////////////////////////////////////////
enum data_attributes_t : size_t {
  persistent
};

} // namespace attributes

////////////////////////////////////////////////////////////////////////////////
/// \brief A specialization of the flecsi low-level mesh topology, state and
///   execution models.
/// \tparam N The number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N, bool Extra_Elements = false >
class burton_mesh  : public burton_mesh_topology_t<N, Extra_Elements>

{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the base type
  using base_t = burton_mesh_topology_t<N, Extra_Elements>;

  //! \brief the mesh types
  using types_t = burton_types_t<N, Extra_Elements>;

  //! \brief the mesh traits
  using config_t = typename types_t::config_t;

  //! \brief the number of dimensions
  static constexpr auto num_dimensions = config_t::num_dimensions;
  //! \brief the number of domains
  static constexpr auto num_domains = config_t::num_domains;

  //! a compile string type
  using const_string_t = typename config_t::const_string_t;
  
  //! a bitfield type
  using bitfield_t = typename config_t::bitfield_t;

  //! Integer data type.
  using integer_t = typename config_t::integer_t;

  //! Floating point data type.
  using real_t = typename config_t::real_t;

  //! The size type.
  using size_t = typename config_t::size_t;

  //! The type used for loop indexing
  using counter_t = typename config_t::counter_t;

  //! Point data type.
  using point_t = typename config_t::point_t;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //! Vertex type.
  using vertex_t = typename types_t::vertex_t;

  //! Edge type.
  using edge_t = typename types_t::edge_t;

  //! Cell type.
  using face_t = typename types_t::face_t;

  //! Cell type.
  using cell_t = typename types_t::cell_t;

  //! Wedge type.
  using wedge_t = typename types_t::wedge_t;

  //! Corner type.
  using corner_t = typename types_t::corner_t;

  //! Side type.
  using side_t = typename types_t::side_t;


  //! the index spaces type
  using index_spaces_t = typename types_t::index_spaces_t;
  //! special subspace for 
  using index_subspaces_t = typename types_t::index_subspaces_t;

  //! \brief The locations of different bits that we set as flags
  using bits = typename config_t::bits;

  //! \brief the type of id for marking boundaries
  using tag_t = typename config_t::tag_t;
  using tag_list_t = typename config_t::tag_list_t;

  //! Shape data type.
  using shape_t = typename config_t::shape_t;

  //! the ownership ( exclusive, shared, ghost ) types
  using partition_t = flecsi::partition_t;

  //! other special subsets that we have defined
  enum class subset_t {
    overlapping
  };

  //============================================================================
  // Constructors
  //============================================================================

  //! Default constructor
  burton_mesh() = default;

  //! \brief Assignment operator (default)
  burton_mesh & operator=(const burton_mesh &) = default;

  //! \brief Copy constructor

  burton_mesh(const burton_mesh &src) = default;

  //! \brief allow move construction
  burton_mesh( burton_mesh && ) = default;

  //! \brief Copy constructor for data client handle
  burton_mesh(const burton_mesh& m, bool dummy) : base_t(m, dummy)
  {}

  //! Destructor
  virtual ~burton_mesh() {};

  //============================================================================
  // Accessors
  //============================================================================


  //============================================================================
  // Vertex Interface
  //============================================================================

  //! \brief Return number of vertices in the burton mesh.
  //! \return The number of vertices in the burton mesh.
  auto num_vertices() const
  {
    return 
      base_t::template num_entities<vertex_t::dimension, vertex_t::domain>();
  }

  auto num_vertices( partition_t subset ) const
  {
    return 
      base_t::template num_entities<vertex_t::dimension, vertex_t::domain>(
        subset
      );
  }

  auto num_vertices( subset_t ) const
  {
    // only return overlapping set right now 
    return base_t::template num_subentities<index_subspaces_t::overlapping_vertices>();
  }

  //! \brief Return all vertices in the burton mesh.
  //! \return Return all vertices in the burton mesh as a sequence for use,
  //!   e.g., in range based for loops.
  decltype(auto) vertices() const 
  { 
    return base_t::template entities<vertex_t::dimension, vertex_t::domain>(); 
  }

  decltype(auto) vertices( partition_t subset ) const 
  { 
    return base_t::template entities<vertex_t::dimension, vertex_t::domain>(
      subset
    ); 
  }

  decltype(auto) vertices( subset_t )
  { 
    // only returns overlapping set right now
    return base_t::template subentities<index_subspaces_t::overlapping_vertices>();
  }

  //! \brief Return connected \e EOUT entities for \e EIN entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam EIN Entity type to get EOUTs for.
  //! \tparam EOUT Entity type to return.
  //!
  //! \param[in] e Entity to get EOUTs for.
  //!
  //! \return Corners for entity \e e in domain \e M.
  template<
    size_t M,
    class EIN,
    class EOUT,
    bool Enabled = ( !ristra::compatibility::is_same_v<EIN, EOUT> ),
    typename std::enable_if_t< Enabled >* = nullptr
  >
  FLECSI_INLINE_TARGET decltype(auto) query_entities(const flecsi::topology::domain_entity_u<M, EIN> & e) const
  {
    // normal case:  EIN and EOUT are different, query topology
    return
      base_t::template entities<EOUT::dimension, M, EOUT::domain>(
        e.entity()
      );
  }

  template<
    size_t M,
    class EIN,
    class EOUT,
    bool Enabled = ( ristra::compatibility::is_same_v<EIN, EOUT> ),
    typename std::enable_if_t< Enabled >** = nullptr
  >
  FLECSI_INLINE_TARGET decltype(auto) query_entities(const flecsi::topology::domain_entity_u<M, EIN> & e) const
  {
    // degenerate case:  EIN and EOUT are the same, return trivial set
    using etype = flecsi::topology::domain_entity_u<M, EIN>;
    // TODO:  use index_space instead of array?
    // TODO:  figure out how to get rid of the const_cast
    return std::array<etype, 1>{const_cast<etype &>(e)};
  }

  //! \brief Return vertices associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertices for.
  //!
  //! \param[in] e instance of entity to return vertices for.
  //!
  //! \return Return vertices associated with entity instance \e e as a
  //!    sequence.
  template <class E>
  decltype(auto) vertices(E * e) const
  {
    return base_t::template entities<vertex_t::dimension, vertex_t::domain>(e);
  }

  //! \brief Return vertices for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get vertices for.
  //!
  //! \param[in] e Entity to get vertices for.
  //!
  //! \return Vertices for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) vertices(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, vertex_t>(e);
  }


  //! \brief Return vertices associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertices for.
  //!
  //! \param[in] e instance of entity to return vertices for.
  //!
  //! \return Return vertices associated with entity instance \e e as a
  //!    sequence.
  template <
    typename P,
    typename = typename std::enable_if_t< ristra::utils::is_callable_v<P> >
  >
  decltype(auto) vertices( P && p ) const
  {
    
    auto vs = vertices();
    return vs.filter( std::forward<P>(p) );
  }

  //! \brief Return ids for all vertices in the burton mesh.
  //!
  //! \return Ids for all vertices in the burton mesh.
  decltype(auto) vertex_ids() const
  {
    return base_t::template entity_ids<vertex_t::dimension, vertex_t::domain>();
  }

  //! \brief Return vertex ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return vertex ids for.
  //!
  //! \param[in] e instance of entity to return vertex ids for.
  //!
  //! \return Return vertex ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) vertex_ids(E * e) const
  {
    return 
      base_t::template entity_ids<vertex_t::dimension, vertex_t::domain>(e);
  }

  //! \brief Return vertices for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get vertices for.
  //!
  //! \param[in] e Entity to get vertices for.
  //!
  //! \return Vertices for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) 
  vertex_ids(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return 
      base_t::template entity_ids<vertex_t::dimension, M, vertex_t::domain>(
        e.entity()
      );
  }

  //! \brief Return boundary  vertices in the burton mesh.
  //! \return Return all boundary vertices in the burton mesh as a sequence for 
  //!         use, e.g., in range based for loops.
  decltype(auto) boundary_vertices() const
  { 
    auto verts = vertices();
    auto bnd_verts = verts.filter(
      [](auto v){ return v->is_boundary(); }
    );
    return bnd_verts;
  }


  //============================================================================
  // Edge Interface
  //============================================================================

  //! \brief Return the number of burton mesh edges.
  //! \return The number of burton mesh edges.
  auto num_edges() const
  {
    return base_t::template num_entities<edge_t::dimension, edge_t::domain>();
  }

  auto num_edges(partition_t subset) const
  {
    return base_t::template num_entities<edge_t::dimension, edge_t::domain>(
      subset
    );
  }

  auto num_edges( subset_t ) const
  {
    // only returns overlapping right now
    return base_t::template num_subentities<index_subspaces_t::overlapping_edges>();
  }

  //! \brief Return all edges in the burton mesh.
  //! \return Return all edges in the burton mesh as a sequence for use, e.g., 
  //!         in range based for loops.
  decltype(auto) edges() const
  { 
    return base_t::template entities<edge_t::dimension, edge_t::domain>(); 
  }

  decltype(auto) edges( partition_t subset ) const 
  { 
    return base_t::template entities<edge_t::dimension, edge_t::domain>(subset); 
  }
  
  decltype(auto) edges( subset_t ) const 
  { 
    // only returns overlapping right now
    return base_t::template subentities<index_subspaces_t::overlapping_edges>();
  }


  //! \brief Return edges for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get edges for.
  //!
  //! \param[in] e Entity to get edges for.
  //!
  //! \return Edges for entity \e e in domain \e M.
  template <size_t M, class E>
  decltype(auto) edges(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, edge_t>(e);
  }

  //! \brief Return ids for all edges in the burton mesh.
  //!
  //! \return Ids for all edges in the burton mesh.
  decltype(auto) edge_ids() const
  {
    return base_t::template entity_ids<edge_t::dimension, edge_t::domain>();
  }

  //! \brief Return edge ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return edge ids for.
  //!
  //! \param[in] e instance of entity to return edge ids for.
  //!
  //! \return Return edge ids associated with entity instance \e e as a 
  //!         sequence.
  template <class E>
  decltype(auto) edge_ids(E * e) const
  {
    return base_t::template entity_ids<edge_t::dimension, edge_t::domain>(e);
  }

  //! \brief Return boundary edges in the burton mesh.
  //! \return Return all boundary vertices in the burton mesh as a sequence for 
  //!         use, e.g., in range based for loops.
  decltype(auto) boundary_edges() const
  { 
    auto es = edges();
    auto bnd_es = es.filter(
      [](auto e){ return e->is_boundary(); }
    );
    return bnd_es;
  }

  //============================================================================
  // Face Interface
  //============================================================================

  //! \brief Return the number of faces in the burton mesh.
  //! \return The number of faces in the burton mesh.
  auto num_faces() const
  {
    return base_t::template num_entities<face_t::dimension, face_t::domain>();
  } // num_faces

  auto num_faces(partition_t subset) const
  {
    return base_t::template num_entities<face_t::dimension, face_t::domain>(
      subset
    );
  } // num_faces

  auto num_faces( subset_t ) const
  {
    // only returns overlapping right now
    return base_t::template num_subentities<index_subspaces_t::overlapping_faces>();
  }


  //! \brief Return all faces in the burton mesh.
  //!
  //! \return Return all faces in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  FLECSI_INLINE_TARGET decltype(auto) faces() const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>();
  }

  FLECSI_INLINE_TARGET decltype(auto) faces(partition_t subset) const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>(
      subset
    );
  }
  
  FLECSI_INLINE_TARGET decltype(auto) faces( subset_t ) const 
  { 
    // only returns overlapping right now
    return base_t::template subentities<index_subspaces_t::overlapping_faces>();
  }


  //! \brief Return all faces in the burton mesh.
  //! \return Return all faces in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  FLECSI_INLINE_TARGET decltype(auto) faces() // FIXME const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>();
  }

  //! \brief Return faces associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return faces for.
  //!
  //! \param[in] e instance of entity to return faces for.
  //!
  //! \return Return faces associated with entity instance \e e as a sequence.
  template <class E>
  FLECSI_INLINE_TARGET decltype(auto) faces(E * e) const
  {
    return base_t::template entities<face_t::dimension, face_t::domain>(e);
  }

  //! \brief Return faces for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get faces for.
  //!
  //! \param[in] e Entity to get faces for.
  //!
  //! \return Faces for entity \e e in domain \e M.
  template <size_t M, class E>
  FLECSI_INLINE_TARGET decltype(auto) faces(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, face_t>(e);
  }

  //! \brief Return ids for all faces in the burton mesh.
  //! \return Ids for all faces in the burton mesh.
  decltype(auto) face_ids() const
  {
    return base_t::template entity_ids<face_t::dimension, face_t::domain>();
  }

  //! \brief Return face ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return face ids for.
  //!
  //! \param[in] e instance of entity to return face ids for.
  //!
  //! \return Return face ids associated with entity instance \e e as a 
  //!         sequence.
  template <class E>
  decltype(auto) face_ids(E * e) const
  {
    return base_t::template entity_ids<face_t::dimension, face_t::domain>(e);
  }

  //============================================================================
  // Cell Interface
  //============================================================================

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  FLECSI_INLINE_TARGET size_t num_cells() const
  {
    return base_t::template num_entities<cell_t::dimension, cell_t::domain>();
  }

  FLECSI_INLINE_TARGET size_t num_cells(partition_t subset) const
  {
    return base_t::template num_entities<cell_t::dimension, cell_t::domain>(
      subset
    );
  }

  //! \brief Return all cells in the burton mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  FLECSI_INLINE_TARGET decltype(auto) cells() const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>();
  }

  FLECSI_INLINE_TARGET decltype(auto) cells(partition_t subset) const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>(
      subset
    );
  }

  //! \brief Return all cells in the burton mesh.
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  FLECSI_INLINE_TARGET decltype(auto) cells() // FIXME const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>();
  }

  //! \brief Return cells associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return cells for.
  //!
  //! \param[in] e instance of entity to return cells for.
  //!
  //! \return Return cells associated with entity instance \e e as a sequence.
  template <class E>
  FLECSI_INLINE_TARGET decltype(auto) cells(E * e) const
  {
    return base_t::template entities<cell_t::dimension, cell_t::domain>(e);
  }

  //! \brief Return cells for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get cells for.
  //!
  //! \param[in] e Entity to get cells for.
  //!
  //! \return Cells for entity \e e in domain \e M.
  template <size_t M, class E>
  FLECSI_INLINE_TARGET decltype(auto) cells(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, cell_t>(e);
  }

  //! \brief Return ids for all cells in the burton mesh.
  //! \return Ids for all cells in the burton mesh.
  decltype(auto) cell_ids() const
  {
    return base_t::template entity_ids<cell_t::dimension, cell_t::domain>();
  }

  //! \brief Return cell ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return cell ids for.
  //!
  //! \param[in] e instance of entity to return cell ids for.
  //!
  //! \return Return cell ids associated with entity instance \e e as a 
  //!         sequence.
  template <class E>
  decltype(auto) cell_ids(E * e) const
  {
    return base_t::template entity_ids<cell_t::dimension, cell_t::domain>(e);
  }

  //! \brief Return the number of cells in the burton mesh.
  //! \return The number of cells in the burton mesh.
  decltype(auto) cell_types() const
  {
    auto cs = cells();
    using cell_type_t = decltype( cs[0]->type() );
    std::set< cell_type_t > cell_types;
    for ( auto c : cs ) cell_types.insert( c->type() );
    return cell_types;
  }

  //============================================================================
  // Wedge Interface
  //============================================================================

  //! \brief Return number of wedges in the burton mesh.
  //! \return The number of wedges in the burton mesh.
  size_t num_wedges() const
  {
    return base_t::template num_entities<wedge_t::dimension, wedge_t::domain>();
  }

  size_t num_wedges(partition_t subset) const
  {
    return base_t::template num_entities<wedge_t::dimension, wedge_t::domain>(subset);
  }

  //! \brief Return all wedges in the burton mesh.
  //!
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) wedges() const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>();
  }

  decltype(auto) wedges(partition_t subset) const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>(subset);
  }

  //! \brief Return all wedges in the burton mesh.
  //! \return Return all wedges in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) wedges() // FIXME const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>();
  }

  //! \brief Return wedges associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return wedges for.
  //!
  //! \param[in] e instance of entity to return wedges for.
  //!
  //! \return Return wedges associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) wedges(E * e) const
  {
    return base_t::template entities<wedge_t::dimension, wedge_t::domain>(e);
  }

  //! \brief Return wedges for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get wedges for.
  //!
  //! \param[in] e Entity to get wedges for.
  //!
  //! \return Wedges for entity \e e in domain \e M.
  template<size_t M, class E>
  decltype(auto) wedges(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, wedge_t>(e);
  }

  //! \brief Return ids for all wedges in the burton mesh.
  //! \return Ids for all wedges in the burton mesh.
  decltype(auto) wedge_ids() const
  {
    return base_t::template entity_ids<wedge_t::dimension, wedge_t::domain>();
  }

  //! \brief Return wedge ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return wedge ids for.
  //!
  //! \param[in] e instance of entity to return wedge ids for.
  //!
  //! \return Return wedge ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) wedge_ids(E * e) const
  {
    return base_t::template entity_ids<wedge_t::dimension, wedge_t::domain>(e);
  }

  //============================================================================
  // Corner Interface
  //============================================================================

  //! \brief Return number of corners in the burton mesh.
  //! \return The number of corners in the burton mesh.
  size_t num_corners() const
  {
    return 
      base_t::template num_entities<corner_t::dimension, corner_t::domain>();
  }

  size_t num_corners(partition_t subset) const
  {
    return 
      base_t::template num_entities<corner_t::dimension, corner_t::domain>(subset);
  }

  //! \brief Return all corners in the burton mesh.
  //!
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) corners() const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>();
  }

  decltype(auto) corners(partition_t subset) const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>(subset);
  }

  //! \brief Return all corners in the burton mesh.
  //! \return Return all corners in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) corners() // FIXME const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>();
  }


  //! \brief Return corners associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return corners for.
  //!
  //! \param[in] e instance of entity to return corners for.
  //!
  //! \return Return corners associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) corners(E * e) const
  {
    return base_t::template entities<corner_t::dimension, corner_t::domain>(e);
  }

  //! \brief Return corners for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get corners for.
  //!
  //! \param[in] e Entity to get corners for.
  //!
  //! \return Corners for entity \e e in domain \e M.
  template<size_t M, class E>
  decltype(auto) corners(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, corner_t>(e);
  }

  //! \brief Return ids for all corners in the burton mesh.
  //! \return Ids for all corners in the burton mesh.
  decltype(auto) corner_ids() const
  {
    return base_t::template entity_ids<corner_t::dimension, corner_t::domain>();
  }

  //! \brief Return corner ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return corner ids for.
  //!
  //! \param[in] e instance of entity to return corner ids for.
  //!
  //! \return Return corner ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) corner_ids(E * e) const
  {
    return 
      base_t::template entity_ids<corner_t::dimension, corner_t::domain>(e);
  }


  //============================================================================
  // Side Interface
  //============================================================================

  //! \brief Return number of sides in the burton mesh.
  //! \return The number of sides in the burton mesh.
  size_t num_sides() const
  {
    return 
      base_t::template num_entities<side_t::dimension, side_t::domain>();
  }

  size_t num_sides(partition_t subset) const
  {
    return 
      base_t::template num_entities<side_t::dimension, side_t::domain>(subset);
  }

  //! \brief Return all sides in the burton mesh.
  //!
  //! \return Return all sides in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) sides() const
  {
    return base_t::template entities<side_t::dimension, side_t::domain>();
  }

  decltype(auto) sides(partition_t subset) const
  {
    return base_t::template entities<side_t::dimension, side_t::domain>(subset);
  }

  //! \brief Return all sides in the burton mesh.
  //! \return Return all sides in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) sides() // FIXME const
  {
    return base_t::template entities<side_t::dimension, side_t::domain>();
  }


  //! \brief Return sides associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return sides for.
  //!
  //! \param[in] e instance of entity to return sides for.
  //!
  //! \return Return sides associated with entity instance \e e as a sequence.
  template <class E>
  decltype(auto) sides(E * e) const
  {
    return base_t::template entities<side_t::dimension, side_t::domain>(e);
  }

  //! \brief Return sides for entity \e e in domain \e M.
  //!
  //! \tparam M Domain.
  //! \tparam E Entity type to get sides for.
  //!
  //! \param[in] e Entity to get sides for.
  //!
  //! \return Sides for entity \e e in domain \e M.
  template<size_t M, class E>
  decltype(auto) sides(const flecsi::topology::domain_entity_u<M, E> & e) const
  {
    return query_entities<M, E, side_t>(e);
  }

  //! \brief Return ids for all sides in the burton mesh.
  //! \return Ids for all sides in the burton mesh.
  decltype(auto) side_ids() const
  {
    return base_t::template entity_ids<side_t::dimension, side_t::domain>();
  }

  //! \brief Return side ids associated with entity instance of type \e E.
  //!
  //! \tparam E entity type of instance to return side ids for.
  //!
  //! \param[in] e instance of entity to return side ids for.
  //!
  //! \return Return side ids associated with entity instance \e e as a
  //!   sequence.
  template <class E>
  decltype(auto) side_ids(E * e) const
  {
    return 
      base_t::template entity_ids<side_t::dimension, side_t::domain>(e);
  }


  //============================================================================
  // Region Interface
  //============================================================================

  //! \brief Return the number of regions in the burton mesh.
  //! \return The number of regions in the burton mesh.
  auto num_regions() const
  {
    return 1;
  }

  //! \brief Return the number of regions in the burton mesh.
  //! \param [in]  n  The number of regions in the burton mesh.
  template< typename T >
  void set_regions(T * region_ids)
  {
    for ( auto c : cells() )
      c->region() = region_ids[c.id()];
  }

void set_regions(std::vector<int> &region_ids)
{
    for(auto c: cells())
    {
        int id = region_ids[c.id()];
        c->set_region(id);
    }
    
}

  //! \brief Return all cells in the regions mesh.
  //!
  //! \return Return all cells in the burton mesh as a sequence for use, e.g.,
  //!   in range based for loops.
  decltype(auto) regions()
  {
    // get the number of regions
    auto n = num_regions();
    // get the cells
    auto cs = cells();

    // now filter out the cells of each region
    auto region_cells = cs.bin( [](const auto & c){ return c->region(); } );

    return region_cells;
  }

  //============================================================================
  // Element Creation
  //============================================================================


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell(
    V && verts,
    std::enable_if_t< 
      ristra::compatibility::is_same_v<
        typename std::decay_t<V>::value_type, vertex_t* > &&
      std::remove_pointer_t<typename
        std::decay_t<V>::value_type>::num_dimensions == 1
    >** = nullptr ) 
  {
    return create_1d_element_from_verts_<cell_t>( std::forward<V>(verts) );
  }

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell( 
    std::initializer_list<V*> verts,
    std::enable_if_t< 
      ristra::compatibility::is_same_v<V, vertex_t> && V::num_dimensions == 1 
    >** = nullptr ) 
  {
    return create_1d_element_from_verts_<cell_t>( verts );
  }


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell(
    V && verts,
    std::enable_if_t< 
      ristra::compatibility::is_same_v<
        typename std::decay_t<V>::value_type, vertex_t* > &&
      std::remove_pointer_t<typename
        std::decay_t<V>::value_type>::num_dimensions == 2
    >* = nullptr ) 
  {
    return create_2d_element_from_verts_<cell_t>( std::forward<V>(verts) );
  }

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto  create_cell( 
    std::initializer_list<V*> verts,
    std::enable_if_t< 
      ristra::compatibility::is_same_v<V, vertex_t> && V::num_dimensions == 2 
    >* = nullptr ) 
  {
    return create_2d_element_from_verts_<cell_t>( verts );
  }


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_cell(
    V && verts,
    typename std::enable_if_t< 
      ristra::compatibility::is_same_v< typename std::decay_t<V>::value_type, vertex_t* > &&
      std::remove_pointer_t<typename std::decay_t<V>::value_type>::num_dimensions == 3
    >* = nullptr ) 
  {
    return create_3d_element_from_verts_( std::forward<V>(verts) );
  }

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto  create_cell( 
    std::initializer_list<V*> verts,
    typename std::enable_if_t< 
      ristra::compatibility::is_same_v<V, vertex_t> && V::num_dimensions == 3
    >* = nullptr ) 
  {
    return create_3d_element_from_verts_( verts );
  }

  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  template< typename F >
  auto create_cell(
    F && faces,
    typename std::enable_if_t< 
      ! ristra::compatibility::is_same_v< face_t, vertex_t > &&
      ristra::compatibility::is_same_v< typename std::decay_t<F>::value_type, face_t* >  &&
      std::remove_pointer_t<typename std::decay_t<F>::value_type>::num_dimensions == 3
    >* = nullptr ) 
  {
    return create_3d_element_from_faces_( std::forward<F>(faces) );
  }
  
  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  auto  create_cell( std::initializer_list<face_t *> faces ) {
    return create_3d_element_from_faces_( faces );
  }



  //! \brief Create a face in the burton mesh.
  //! \param[in] verts The vertices defining the face.
  //! \return Pointer to cell created with \e faces.
  template< typename V >
  auto create_face(V && verts)
  {
    return create_2d_element_from_verts_<face_t>( std::forward<V>(verts) );
  } // create_cell

  
  //! \brief Create a face in the burton mesh.
  //! \param[in] verts The vertices defining the face.
  //! \return Pointer to cell created with \e faces.
  auto  create_face( std::initializer_list<vertex_t *> verts ) {
    return create_2d_element_from_verts_<face_t>( verts );
  }


  //! \brief Create a vertex in the burton mesh.
  //!
  //! \param[in] pos The position (coordinates) for the vertex.
  //!
  //! \return Pointer to a vertex created at \e pos.
  template< typename... ARGS >
  vertex_t * create_vertex( ARGS &&... args )
  {
    auto v = 
      base_t::template make<vertex_t>( std::forward<ARGS>(args)... );
    return v;
  }

  //============================================================================
  // Mesh Creation
  //============================================================================

  //! \brief Dump the burton mesh to standard out.
  std::ostream & dump( std::ostream & stream )
  {
    return base_t::dump( stream );
  }

  //! \brief Initialize burton mesh state for the number of \e vertices.
  //!
  //! \param[in] vertices The number of \e vertices to initialize the burton 
  //!                     mesh with.
  void init_parameters(size_t num_nodes)
  { }

  //!---------------------------------------------------------------------------
  //! \brief Initialize the burton mesh.
  //!---------------------------------------------------------------------------
  void init()
  {

    base_t::template init<0>();
    
    base_t::template init_bindings<1>();

    // now set the boundary flags.
    for ( auto f : faces() ) {
	    auto cs = cells(f);
      f->set_owner(cs[0]);
      f->set_boundary( (cs.size() == 1) );
      // if there is only one cell, it is a boundary
      if ( f->is_boundary() ) {
        if constexpr ( num_dimensions >= 2 ) {
	        // cell flags
	        for ( auto c : cs )
	          c->set_touching_boundary( true );
          // point flags are only for 2d and 3d
          auto ps = vertices(f);
          for ( auto p : ps )
            p->set_boundary( true );
        }
        // edge flags are only for 3d
        if constexpr ( num_dimensions == 3 ) {
          auto es = edges(f);
          for ( auto e : es ) 
            e->set_boundary( true );
        }
      } // is_boundary
    } // for
    
#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    // set boundary wedges
    for ( auto w : wedges() ) 
      w->set_boundary( edges(w).front()->is_boundary() );
#endif

    // identify the cell regions
    //for ( auto c : cells() ) c->region() = 0;

    // update the geometry
    // update_geometry();
  }


  //!---------------------------------------------------------------------------
  //! \brief Check the burton mesh.
  //!---------------------------------------------------------------------------
  bool is_valid( bool raise_on_error = true )
  {
    // some includes
    using ristra::math::dot_product;

    // a lambda function for raising errors or returning 
    // false
    auto raise_or_return = [=]( std::ostream & msg )
      {
        if ( raise_on_error )
          THROW_RUNTIME_ERROR( msg.rdbuf() );
        else 
          std::cerr << msg.rdbuf() << std::endl;
        return false;
      };
    
    // we use a stringstream to construct messages and pass them 
    // simultanesously
    std::stringstream ss;

    //--------------------------------------------------------------------------
    // make sure face normal points out from first cell
    
    auto & context = flecsi::execution::context_t::instance();
    auto rank = context.color();
    auto & vertex_map = context.index_map( index_spaces_t::vertices );
    auto & face_map = context.index_map( index_spaces_t::faces );
    auto & cell_map = context.index_map( index_spaces_t::cells );

    for( auto f : faces() ) {
      auto n = f->normal();
      auto fx = f->midpoint();
      auto c = cells(f).front();
      auto cx = c->midpoint();
      auto flipped = (f->is_flipped(c));
      auto delta = flipped ? cx - fx : fx - cx;
      auto dot = dot_product( n, delta );
      // std::cout << "Checking face with mid " << face_map[ f.id() ] << std::endl;
      // std::cout << "With cells : ";
      // for ( auto cl : cells(f) ) std::cout << cell_map[ cl.id() ] << ", ";
      // //for ( auto cl : cells(f) ) std::cout << cl.id() << ", ";
      // std::cout << std::endl;
      // std::cout << "And vertices : ";
      // for ( auto vt : vertices(f) ) std::cout << vertex_map[ vt.id() ] << ", ";
      // std::cout << std::endl;
      // std::cout << "Face has midpoint " << fx << std::endl;
      // std::cout << "Cell has midpoint " << c->midpoint() << std::endl;
      // std::cout << "Cell has centroid " << c->centroid() << std::endl;
      // std::cout << dot << ", " << std::boolalpha << flipped << std::endl;
      // std::cout << std::endl;
      if ( dot < 0 ) {
        ss << "Face " << face_map.at(f.id()) << " of rank " << rank << 
          " has opposite normal " << std::endl;
      }
    } 

    if ( ss.tellp() != 0 ) return raise_or_return( ss );

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS

    //--------------------------------------------------------------------------
    // check all the corners and wedges
    auto cnrs = corners();
    auto num_corners = cnrs.size();
    auto & wedge_map = context.index_map( index_spaces_t::wedges );

    //#omp parallel for reduction( || : bad_corner )
    for( counter_t cnid=0; cnid<num_corners; ++cnid ) {
  
      auto cn = cnrs[cnid];
      auto cs = cells(cn);
      auto fs = faces(cn);
      auto es = edges(cn);
      auto vs = vertices(cn);

      if ( cs.size() != 1 ) {
        ss << "Corner " << cn.id() << " has " << cs.size() << "/=1 cells" 
           << std::endl;
      }

      if ( fs.size() != num_dimensions ) {
        ss << "Corner " << cn.id() << " has " << fs.size() << "/=" 
           << num_dimensions << " faces" << std::endl;
      }

      if ( es.size() != num_dimensions ) {
        ss << "Corner " << cn.id() << " has " << es.size() << "/=" 
           << num_dimensions << " edges" << std::endl;
      }

      if ( vs.size() != 1 ) {
        ss << "Corner " << cn.id() << " has " << vs.size() << "/=1 vertices"
           << std::endl;
      }

      if ( num_dimensions == 1 ) continue;

      auto cl = cs.front();
      auto vt = vs.front();
      
      auto ws = wedges(cn);

      if ( ws.size() % 2 != 0 ) {
        if ( vs.size() != 1 ) {
          ss << "Corner " << cn.id() << " has " << ws.size() << "%2/=0 wedges"
             << std::endl;
        }
      }

      for ( auto wg : ws  ) {
        auto cls = cells( wg );
        auto fs = faces( wg );
        auto es = edges( wg );
        auto vs = vertices( wg );
        auto cns = corners( wg );
        if ( cls.size() != 1 ) {
          ss << "Wedge " << wg.id() << " has " << cls.size() 
             << "/=1 cells" << std::endl;
        }
        if ( fs.size() != 1 ) {
          ss << "Wedge " << wg.id() << " has " << fs.size() 
             << "/=1 faces" << std::endl;
        }
        if ( es.size() != 1 ) {
          ss << "Wedge " << wg.id() << " has " << es.size() 
             << "/=1 edges" << std::endl;
        }
        if ( vs.size() != 1 ) {
          ss << "Wedge " << wg.id() << " has " << vs.size() 
             << "/=1 vertices" << std::endl;
        }
        if ( cns.size() != 1 ) {
          ss << "Wedge " << wg.id() << " has " << cns.size() 
             << "/=1 corners" << std::endl;
        }
        auto vert = vs.front();
        auto cell = cls.front();
        auto corn = cns.front();
        if ( vert != vt ) {
          ss << "Wedge " << wg.id() << " has incorrect vertex " 
             << vert.id() << "!=" << vt.id() << std::endl;
        }
        if ( cell != cl ) {
          ss << "Wedge " << wg.id() << " has incorrect cell " 
             << cell.id() << "!=" << cl.id() << std::endl;
        }
        if ( corn != cn ) {
          ss << "Wedge " << wg.id() << " has incorrect corner " 
             << corn.id() << "!=" << cn.id() << std::endl;
        }
        auto cell_verts = vertices(cl);
        auto vert_exists = false;
        for ( auto v : vertices(cl) )
          if ( v == vt ) {
            vert_exists = true;
            break;
          }
        if ( !vert_exists ) {
          ss << "Wedge " << wg.id() << " has vertex " 
             << vt.id() << " that is not in cell" << std::endl;
        }
        auto fc = fs.front();            
        auto fx = fc->midpoint();
        auto cx = cl->centroid();
        auto delta = fx - cx;
        auto flipped = (fc->is_flipped(cl));
        real_t dot;
        auto n = wg->facet_normal();
        dot = dot_product( n, delta );
        if ( dot < 0 ) {
          ss << "Wedge " << wg.id() << " has opposite normal dot=" << dot
             << std::endl;
        }
      } // wedges
      
    } // corners


    if ( ss.tellp() != 0 ) return raise_or_return( ss );

#endif // FLECSI_SP_BURTON_MESH_EXTRAS

    return true;

  }

  //!---------------------------------------------------------------------------
  //! \brief Compute the goemetry.
  //!---------------------------------------------------------------------------
  void update_geometry()
  {
    // get the mesh info
    auto cs = cells();
    auto fs = faces();
    auto es = edges();  
    auto num_cells = cs.size();
    auto num_faces = fs.size();
    auto num_edges = es.size();

    #pragma omp parallel
    {

      //--------------------------------------------------------------------------
      // compute cell parameters (in 1D, this must come first;
      // edge/face update requires cell to be up-to-date)

      if ( num_dimensions == 1 ) {
        #pragma omp for
        for ( counter_t i=0; i<num_cells; i++ )
          cs[i]->update( this );
      }

      //--------------------------------------------------------------------------
      // compute edge parameters

      #pragma omp for nowait
      for ( counter_t i=0; i<num_edges; i++ )
        es[i]->update( this );

      //--------------------------------------------------------------------------
      // compute face parameters

      if ( num_dimensions == 3 ) {
        #pragma omp for nowait
        for ( counter_t i=0; i<num_faces; i++ )
          fs[i]->update( this );
      }

      //--------------------------------------------------------------------------
      // compute cell parameters (in 2D and 3D, this must follow
      // edges and faces)

      if ( num_dimensions > 1 ) {
        #pragma omp for
        for ( counter_t i=0; i<num_cells; i++ )
          cs[i]->update( this );
      }

      //--------------------------------------------------------------------------
      // compute wedge parameters

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS

      if ( num_dimensions == 1 ) {
        // loop through wedges directly
        auto ws = wedges();
        auto num_wedges = ws.size();
        #pragma omp for
        for ( counter_t i=0; i<num_wedges; ++i ) {
          ws[i]->update( this, true );
        }
      }  // num_dimensions == 1
      else {
        // loop through corners, then wedges per corner
        auto cnrs = corners();
        auto num_corners = cnrs.size();
        #pragma omp for
        for ( counter_t i=0; i<num_corners; ++i ) {
          auto cn = cnrs[i];
          auto ws = wedges(cn);
          // first compute the normals
          for ( auto wit = ws.begin(); wit != ws.end(); ++wit )
          {
            // get the first wedge normal
            (*wit)->update( this, false );
            // move to next wedge
            ++wit;
            assert( wit != ws.end() );
            // get the second wedge normal
            (*wit)->update( this, true );
          }
        }
      }  // else num_dimensions

#endif // FLECSI_SP_BURTON_MESH_EXTRAS

    } // end omp parallel

  }


  //============================================================================
  //! \brief Install a boundary and tag the relatex entities.
  //============================================================================
  template< typename P >
  void install_boundary( P && p, tag_t key ) 
  {
    std::vector<face_t*>   bnd_faces;
    std::vector<edge_t*>   bnd_edges;
    std::vector<vertex_t*> bnd_verts;

    // add the face tags and collect the attached edge and vertices
    for ( auto f : faces() )
      if ( p( f ) ) {
        // tag the face
        f->tag( key );
        bnd_faces.emplace_back( f );
        // tag the vertices in 2d or 3d
        if ( num_dimensions >= 2 ) {
          auto vs = vertices( f );
          bnd_verts.reserve( bnd_verts.size() + vs.size() );
          for ( auto v : vs ) bnd_verts.emplace_back( v );
        }
        // tag edges in 3d
        if ( num_dimensions == 3 ) {
          auto es = edges( f );
          bnd_edges.reserve( bnd_edges.size() + es.size() );
          for ( auto e : es ) bnd_edges.emplace_back( e );
        } // dims
      }

    // need to remove duplicates from edge and vertex lists
    std::sort( bnd_edges.begin(), bnd_edges.end() );
    std::sort( bnd_verts.begin(), bnd_verts.end() );

    bnd_edges.erase( 
      std::unique( bnd_edges.begin(), bnd_edges.end() ), bnd_edges.end() );
    bnd_verts.erase( 
      std::unique( bnd_verts.begin(), bnd_verts.end() ), bnd_verts.end() );

    // add the edge tags
    for ( auto e : bnd_edges ) e->tag( key );
    // add the vertex tags
    for ( auto v : bnd_verts ) v->tag( key );

  }

#if 0
  //============================================================================
  //! \brief Get the set of tagged vertices associated with a specific id
  //! \praram [in] id  The tag to lookup.
  //! \return The set of tagged vertices.
  //============================================================================
  const auto & tagged_vertices( tag_t id ) const noexcept
  {
    return vert_sets_[ id ];
  }
#endif

  auto cell_vertex_neighbors(
    const flecsi::topology::domain_entity_u<0, cell_t> & c0) const
  { 
 
    std::vector<cell_t*> cs;

    auto vs = vertices(c0);
    for ( auto v : vs )
      for ( auto c : cells(v) )
        if ( c != c0 ) 
          cs.emplace_back(c);
    
    std::sort( cs.begin(), cs.end() );
    auto last = std::unique( cs.begin(), cs.end() );
    cs.erase( last, cs.end() );
    return cs;
  }


  //============================================================================
  // Operators
  //============================================================================

  //! Print some statistics.
  template< std::size_t M >
  friend std::ostream& operator<< (std::ostream& stream, const burton_mesh<M>& mesh);



  //============================================================================
  // Private Members
  //============================================================================

 private:


  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename E, typename V >
  auto create_1d_element_from_verts_( V && verts  )
  {
    
    E * e = nullptr;

    // should be "segment", but that doesn't exist yet
    auto cell_type = shape_t::none;

    if ( verts.size() != 2 )
      THROW_RUNTIME_ERROR( "must have exactly 2 vertices" );

    e = base_t::template make<E>(cell_type);
    base_t::template init_entity<E::domain, E::dimension, vertex_t::dimension>(
      e, std::forward<V>(verts)
    );
    return e;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename E, typename V >
  auto create_2d_element_from_verts_( V && verts  )
  {
    
    E * e = nullptr;

    auto cell_type = shape_t::polygon;

    switch ( verts.size() ) {
    case (1,2):
      THROW_RUNTIME_ERROR( "can't have <3 vertices" );
    case (3):
      cell_type = shape_t::triangle;
      break;
    case (4):
      cell_type = shape_t::quadrilateral;
      break;
    }

    e = base_t::template make<E>(cell_type);
    base_t::template init_entity<E::domain, E::dimension, vertex_t::dimension>(
      e, std::forward<V>(verts)
    );
    return e;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] verts The vertices defining the cell.
  //! \return Pointer to cell created with \e verts.
  template< typename V >
  auto create_3d_element_from_verts_( V && verts )
  {
    
    cell_t * c = nullptr;
    
    auto cell_type = shape_t::none;

    switch ( verts.size() ) {
    case (1,2,3):
      THROW_RUNTIME_ERROR( "can't have <4 vertices" );
    case (4):
      cell_type = shape_t::tetrahedron;
      break;
    case (8):
      cell_type = shape_t::hexahedron;
      break;
    default:
      THROW_RUNTIME_ERROR( "can't build polyhedron from vertices alone" );
    }

    c = base_t::template make< cell_t >(cell_type);
    base_t::template init_cell<cell_t::domain>( c, std::forward<V>(verts) );

    return c;
  } // create_cell

  //! \brief Create a cell in the burton mesh.
  //! \param[in] faces The faces defining the cell.
  //! \return Pointer to cell created with \e faces.
  template< typename F >
  auto create_3d_element_from_faces_( F && faces )
  {
    
    cell_t * c = nullptr;

    switch ( faces.size() ) {
    case (1,2,3):
      THROW_RUNTIME_ERROR( "can't have <4 vertices" );
    }

    c = base_t::template make< cell_t >( shape_t::polyhedron );
    base_t::template init_entity<
      cell_t::domain, cell_t::dimension, face_t::dimension
    >( c, std::forward<F>(faces) );
    return c;
  } // create_cell


}; // class burton_mesh_t


////////////////////////////////////////////////////////////////////////////////
// External Class Definitions
////////////////////////////////////////////////////////////////////////////////

//==============================================================================
// Friends
//==============================================================================

//!  \brief Print some statistics.
//! 
//!  \param[in] stream The stream to print to.
//! 
//!  \param[in] mesh   The mesh object to print stats for.
//! 
//!  \return the stream operator.
template< std::size_t M >
inline
std::ostream& operator<< (std::ostream& stream, const burton_mesh<M>& mesh)
{
  using std::endl;
  stream << "Burton mesh:" << endl;
  stream << " + Num Points = " << mesh.num_vertices() << endl;
  stream << " + Num Edges  = " << mesh.num_edges() << endl;
  stream << " + Num Faces  = " << mesh.num_faces() << endl;
  stream << " + Num Cells  = " << mesh.num_cells() << endl;
  return stream;
}

//==============================================================================
// Final mesh type
//==============================================================================
#ifndef FLECSI_SP_BURTON_MESH_EXTRAS

#  ifdef FLECSI_SP_BURTON_MESH_DIMENSION
using burton_mesh_t = burton_mesh<FLECSI_SP_BURTON_MESH_DIMENSION>;
#  else
using burton_mesh_t = burton_mesh<2>;
#  endif

#else

#  ifdef FLECSI_SP_BURTON_MESH_DIMENSION
using burton_mesh_t = burton_mesh<FLECSI_SP_BURTON_MESH_DIMENSION,true>;
#  else
using burton_mesh_t = burton_mesh<2,true>;
#  endif

#endif

} // namespace burton
} // namespace flecsi-sp
