/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Provides the vertex type for burton_mesh_t. 
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsi/topology/mesh_types.h>
#include <flecsi-sp/burton/burton_config.h>
#include <ristra/math/general.h>
#include <ristra/geometry/shapes/geometric_shapes.h>
#include <ristra/assertions/errors.h>


namespace flecsi_sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//!
//! \tparam N The dimension of the vertex.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_vertex_t : public 
  flecsi::topology::mesh_entity__<0, burton_config_t<N>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<N>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = typename config_t::mesh_storage_t;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base__< mesh_storage_t >;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type.
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  template< typename...ARGS >
  burton_vertex_t(ARGS &&... args) 
    : coordinates_{ std::forward<ARGS>(args)... }
  {}

  // dissallow copying
  burton_vertex_t( burton_vertex_t & ) = delete;
  burton_vertex_t & operator=( burton_vertex_t & ) = delete;

  // dissallow moving
  burton_vertex_t( burton_vertex_t && ) = delete;
  burton_vertex_t & operator=( burton_vertex_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  const point_t & coordinates() const noexcept
  { return coordinates_; }

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  //! \remark this is the non const version
  point_t & coordinates() noexcept
  { return coordinates_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! return true if this is on a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the coordinates of the vertex
  point_t coordinates_ = 0;

  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

}; // class burton_vertex_t


////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_vertex_t type provides an interface for managing
//!   geometry and state associated with mesh vertices.
//!   This is the special 1D version.  Since vertices in 1D also serve
//!   as edges and faces, we give this class some trivial edge- and
//!   face-like methods.
//!
//! \tparam N The dimension of the vertex.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_vertex_t<1> : public
  flecsi::topology::mesh_entity__<0, burton_config_t<1>::num_domains>
{
public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<1>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = typename config_t::mesh_storage_t;
  //! the flecsi mesh topology type
  using mesh_topology_base_t =
    flecsi::topology::mesh_topology_base__< mesh_storage_t >;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type.
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  template< typename...ARGS >
  burton_vertex_t(ARGS &&... args)
    : coordinates_{ std::forward<ARGS>(args)... }
  {}

  // dissallow copying
  burton_vertex_t( burton_vertex_t & ) = delete;
  burton_vertex_t & operator=( burton_vertex_t & ) = delete;

  // dissallow moving
  burton_vertex_t( burton_vertex_t && ) = delete;
  burton_vertex_t & operator=( burton_vertex_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  const point_t & coordinates() const noexcept
  { return coordinates_; }

  //! \brief Get the coordinates at a vertex from the state handle.
  //! \return coordinates of vertex.
  //! \remark this is the non const version
  point_t & coordinates() noexcept
  { return coordinates_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! return true if this is on a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //============================================================================
  // Accessors / Modifiers for 1D edge/face
  //============================================================================

  //! the edge "midpoint" (in 1D, the vertex itself)
  const auto & midpoint() const
  { return coordinates_; }

  //! the edge "centroid" (in 1D, the vertex itself)
  const auto & centroid() const
  { return coordinates_; }

  //! the edge "length"
  auto length() const
  { return 1.; }

  //! the face "area"
  auto area() const
  { return 1.; }

  //! the face "normal"
  const auto & normal() const
  {
    return normal_;
  }

  //! \brief update the mesh geometry
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    // TODO:  This code recomputes the orientation on every update -
    //        can we do it just once at initialization?
    using ristra::math::sgn;
    auto cs = mesh->template entities<1, domain>(this);
    const auto & v = this->coordinates();
    const auto & c = cs[0]->midpoint();
    normal_ = { sgn(v[0] - c[0]) };
  }

  //============================================================================
  // Private Data
  //============================================================================

private:

  //! the coordinates of the vertex
  point_t coordinates_ = 0;

  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

  //! the "normal" vector
  vector_t normal_ = 0;

}; // class burton_vertex_t



} // namespace burton
} // namespace flecsi_sp

