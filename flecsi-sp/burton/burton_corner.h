/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This file defines the corner entity for the burton mesh.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsi-sp/burton/burton_vertex.h>
#include <flecsi-sp/burton/burton_element.h>


namespace flecsi_sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_corner_t type provides an interface for managing and
//!   geometry and state associated with mesh corners.
//!
//! \tparam N The domain of the corner.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_corner_t
  : public flecsi::topology::mesh_entity__<0, burton_config_t<N>::num_domains>
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
  static constexpr auto domain = 1;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = typename config_t::mesh_storage_t;
  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base__< mesh_storage_t >;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_corner_t() = default;

  // dissallow copying
  burton_corner_t( burton_corner_t & ) = delete;
  burton_corner_t & operator=( burton_corner_t & ) = delete;

  // dissallow moving
  burton_corner_t( burton_corner_t && ) = delete;
  burton_corner_t & operator=( burton_corner_t && ) = delete;

};

////////////////////////////////////////////////////////////////////////////////
//! \brief The burton_corner_t type provides an interface for managing and
//!   geometry and state associated with mesh corners.
//!   This is the special 1D version.  Since corners in 1D also serve
//!   as wedges, we give this class some trivial wedge-like methods.
//!
//! \tparam N The domain of the corner.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_corner_t<1>
  : public flecsi::topology::mesh_entity__<0, burton_config_t<1>::num_domains>
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
  static constexpr auto domain = 1;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = typename config_t::mesh_storage_t;
  //! the flecsi mesh topology type
  using mesh_topology_base_t =
    flecsi::topology::mesh_topology_base__< mesh_storage_t >;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_corner_t() = default;

  // dissallow copying
  burton_corner_t( burton_corner_t & ) = delete;
  burton_corner_t & operator=( burton_corner_t & ) = delete;

  // dissallow moving
  burton_corner_t( burton_corner_t && ) = delete;
  burton_corner_t & operator=( burton_corner_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers for 1D wedge
  //============================================================================

  //! \brief update the wedge geometry (in 1D, a no-op)
  //! \param [in] is_right  This wedge is the right orientation when true.
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, bool is_right ) {}

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  const auto & facet_normal() const
  { return facet_normal_; }

  //! \brief Is this wedge on the boundary.
  //! \return true if on boundary.
  auto is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  //============================================================================
  // Private Data
  //============================================================================

private:

  //! the "normal" of the outer facet
  const vector_t facet_normal_{1.};
  //! the entity flags
  bitfield_t flags_;

};


} // namespace burton
} // namespace flecsi_sp
