/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines the corner and wedge entities for the burton mesh.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsi-sp/burton/burton_vertex.h>
#include <flecsi-sp/burton/burton_element.h>
#include <ristra/math/general.h>
#include <ristra/geometry/shapes/triangle.h>


namespace flecsi_sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
//! \brief This type provides the base for all the dual entities (extras).
//!
//! \tparam NUM_DIMS  The total number of dimensions.
//! \tparam DIM       The actual dimension of the entity.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t NUM_DIMS, std::size_t DIM, typename Enable = void >
class burton_extras_t { };


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh corners and wedges.
//! \tparam N The total number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_extras_types_t {
  using corner_t = burton_extras_t<N,0>;
  using wedge_t = burton_extras_t<N,1>;
  using side_t = burton_extras_t<N,2>;
};

template<>
struct burton_extras_types_t<1> {
  using corner_t = burton_extras_t<1,0>;
  using wedge_t = burton_extras_t<1,0>;
  using side_t = burton_extras_t<1,0>;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh corners.
//! \tparam N The total number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_corner_t = typename burton_extras_types_t<N>::corner_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh wedges.
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_wedge_t = typename burton_extras_types_t<N>::wedge_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh sides.
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_side_t = typename burton_extras_types_t<N>::side_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!   two- or three-dimensional mesh corners.
//!
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////

template< std::size_t N >
class burton_extras_t<N,0>
  : public flecsi::topology::mesh_entity_u<0, burton_config_t<N>::num_domains>
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
    flecsi::topology::mesh_topology_base_u< mesh_storage_t >;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_extras_t() = default;

  // dissallow copying
  burton_extras_t( burton_extras_t & ) = delete;
  burton_extras_t & operator=( burton_extras_t & ) = delete;

  // dissallow moving
  burton_extras_t( burton_extras_t && ) = delete;
  burton_extras_t & operator=( burton_extras_t && ) = delete;

};

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!   one-dimensional mesh corners/wedges.
//!   Since corners in 1D also serve as wedges, we give this class some
//!   trivial wedge-like methods.
//!
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template<>
class burton_extras_t<1,0>
  : public flecsi::topology::mesh_entity_u<0, burton_config_t<1>::num_domains>
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
    flecsi::topology::mesh_topology_base_u< mesh_storage_t >;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //! Coordinate point type.
  using point_t = typename config_t::point_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_cell_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_extras_t() = default;

  // dissallow copying
  burton_extras_t( burton_extras_t & ) = delete;
  burton_extras_t & operator=( burton_extras_t & ) = delete;

  // dissallow moving
  burton_extras_t( burton_extras_t && ) = delete;
  burton_extras_t & operator=( burton_extras_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers for 1D wedge
  //============================================================================

  //! \brief update the wedge geometry (in 1D, a no-op)
  //! \param [in] is_right  Unused in 1D.
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, bool is_right );

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  const auto & facet_normal() const
  { return facet_normal_; }

  //! \brief the facet area
  auto facet_area() const
  { return 1.; }

  //! \brief Get the cell facet centroid
  const auto & facet_centroid() const
  { return facet_centroid_; }

  //! \brief Get the cell facet midpoint
  const auto & facet_midpoint() const
  { return facet_centroid_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

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

  //! "centroid" of the outer facet
  point_t facet_centroid_ = 0;
  //! the "normal" of the outer facet
  vector_t facet_normal_ = 0;
  //! the entity flags
  bitfield_t flags_;

};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!   two- or three-dimensional mesh wedges.
//!
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_extras_t<N,1,typename std::enable_if< N != 1 >::type>
  : public flecsi::topology::mesh_entity_u<1, burton_config_t<N>::num_domains>
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
    flecsi::topology::mesh_topology_base_u< mesh_storage_t >;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Physics vector type.
  using vector_t = typename config_t::vector_t;

  //! Coordinate point type.
  using point_t = typename config_t::point_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_cell_t<num_dimensions>;

  //============================================================================
  // Constructors
  //============================================================================

  // default constructor
  burton_extras_t() = default;

  // dissallow copying
  burton_extras_t( burton_extras_t & ) = delete;
  burton_extras_t & operator=( burton_extras_t & ) = delete;

  // dissallow moving
  burton_extras_t( burton_extras_t && ) = delete;
  burton_extras_t & operator=( burton_extras_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! \brief update the wedge geometry
  //! \param [in] is_right  This wedge is the right orientation when true.
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, bool is_right );

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  const auto & facet_normal() const
  { return facet_normal_; }

  //! \brief the facet area
  auto facet_area() const
  { return facet_area_; }

  //! \brief Get the cell facet normal for the wedge.
  //! \return Cell facet normal vector.
  const auto & internal_facet_normal() const
  { return internal_facet_normal_; }

  //! \brief the facet area
  auto internal_facet_area() const
  { return internal_facet_area_; }

  //! \brief Get the cell facet centroid
  const auto & facet_centroid() const
  { return facet_centroid_; }
  
  //! \brief Get the cell facet midpoint
  const auto & facet_midpoint() const
  { return facet_centroid_; }
  
  //! \brief Get the cell facet centroid
  const auto & internal_facet_centroid() const
  { return internal_facet_centroid_; }

  //! \brief Get the cell facet centroid
  const auto & internal_facet_midpoint() const
  { return internal_facet_centroid_; }


  //! \brief Get the corner-midpoint
  const auto & corner_midpoint() const
  { return corner_midpoint_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! \brief Is this wedge on the boundary.
  //! \return true if on boundary.
  auto is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }
  
  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  auto volume() const
  {
      return volume_;
  }
  auto is_right() const { return is_right_; }

  //============================================================================
  // Private Data
  //============================================================================

 private:

  //! centroid of the outer facetn
  point_t facet_centroid_ = 0;
  //! the normal of the outer facet
  vector_t facet_normal_ = 0;
  //! the area of the outer facet
  real_t facet_area_ = 0;
  //! the entity flags
  bitfield_t flags_;

  //! volume of wedge.
  real_t volume_=0;

  //! the normal of edge that connects mp and e
  vector_t internal_facet_normal_ = 0;
  //! the area of edge that connects mp and e.
  real_t internal_facet_area_ = 0;
  //! centroid of the internal facet
  point_t internal_facet_centroid_ = 0;

  //! this is midpoint between cell mid-point and vertex.
  point_t corner_midpoint_=0;
  bool is_right_ = false;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!   two- or three-dimensional mesh sides.
//!
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
class burton_extras_t<N,2,typename std::enable_if< N != 1 >::type>
  : public flecsi::topology::mesh_entity_u<2, burton_config_t<N>::num_domains>
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
    flecsi::topology::mesh_topology_base_u< mesh_storage_t >;

  //============================================================================
  // Constructors
  //============================================================================

  //! default constructor
  burton_extras_t() = default;

  // dissallow copying
  burton_extras_t( burton_extras_t & ) = delete;
  burton_extras_t & operator=( burton_extras_t & ) = delete;

  // dissallow moving
  burton_extras_t( burton_extras_t && ) = delete;
  burton_extras_t & operator=( burton_extras_t && ) = delete;

};

////////////////////////////////////////////////////////////////////////////////
// facet normal for 1d wedge
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TOPOLOGY >
void burton_extras_t<1,0>::update( const MESH_TOPOLOGY * mesh, bool is_right )
{
  // TODO:  This code recomputes the orientation on every update -
  //        can we do it just once at initialization?
  using ristra::math::sgn;
  auto vs = mesh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto cs = mesh->template entities<cell_t::dimension, domain, cell_t::domain>(this);
  assert( vs.size() == 1 );
  const auto & v = vs.front()->coordinates();
  const auto & c = cs.front()->midpoint();
  facet_normal_ = { sgn(v[0] - c[0]) };
  facet_centroid_ = v;
}

////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 2d wedge
////////////////////////////////////////////////////////////////////////////////
template<>
template< typename MESH_TOPOLOGY >
void burton_extras_t<2,1>::update( const MESH_TOPOLOGY * mesh, bool is_right )
{
  using ristra::math::abs;
  auto cs = mesh->template entities<cell_t::dimension, domain, cell_t::domain>(this);
  auto vs = mesh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = mesh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  assert( cs.size() == 1 );
  assert( vs.size() == 1 );
  assert( es.size() == 1 );

  //const auto & mp = cs.front()->midpoint();
  const auto & mp = cs.front()->centroid();

  const auto & e = es.front()->midpoint();
  const auto & v = vs.front()->coordinates();

  if ( is_right )
  {
    facet_normal_ = { e[1] - v[1], v[0] - e[0] };
    internal_facet_normal_ = {  mp[1]-e[1], e[0] - mp[0] };
  }
  else
  {
      
    facet_normal_ = { v[1] - e[1], e[0] - v[0] };
    internal_facet_normal_ = { e[1] - mp[1], mp[0] - e[0] };
  }

  facet_area_ = abs(facet_normal_);
  facet_normal_ /= facet_area_;
  facet_centroid_ = 0.5 * ( e + v );

  volume_ = std::abs((mp[0]*(e[1]-v[1])+e[0]*(v[1]-mp[1])+v[0]*(mp[1]-e[1]))/2);

  internal_facet_area_ = abs(internal_facet_normal_);
  internal_facet_normal_ /= internal_facet_area_;
  internal_facet_centroid_ = 0.5 * ( e + mp );

  corner_midpoint_ = 0.5*( mp + v );
}

////////////////////////////////////////////////////////////////////////////////
// left and right facet normals for 3d wedge
////////////////////////////////////////////////////////////////////////////////
template<>
template< typename MESH_TOPOLOGY >
void burton_extras_t<3,1>::update(const  MESH_TOPOLOGY* mesh, bool is_right)
{
  using ristra::math::abs;
  auto vs = mesh->template entities<vertex_t::dimension, domain, vertex_t::domain>(this);
  auto es = mesh->template entities<edge_t::dimension, domain, edge_t::domain>(this);
  auto fs = mesh->template entities<face_t::dimension, domain, face_t::domain>(this);
  assert( vs.size() == 1 );
  assert( es.size() == 1 );
  assert( fs.size() == 1 );
  auto e = es.front()->midpoint();
  auto v = vs.front()->coordinates();
  auto f = fs.front()->midpoint();

  auto cs = mesh->template entities<cell_t::dimension, domain, cell_t::domain>(this);
  const auto & mp = cs.front()->centroid();
  
  if ( !is_right )
  {
    
    facet_normal_ =
      ristra::geometry::shapes::triangle<num_dimensions>::normal( v, f, e );

    internal_facet_normal_ =
      ristra::geometry::shapes::triangle<num_dimensions>::normal( f, mp, e );
  }
  else
  {
    facet_normal_ =
      ristra::geometry::shapes::triangle<num_dimensions>::normal( v, e, f );

    internal_facet_normal_ =
      ristra::geometry::shapes::triangle<num_dimensions>::normal( f, e, mp );
  }
  
  facet_area_ = abs(facet_normal_);
  facet_normal_ /= facet_area_;

  internal_facet_area_ = abs(internal_facet_normal_);
  internal_facet_normal_ /= internal_facet_area_;

  facet_centroid_ =
    ristra::geometry::shapes::triangle<num_dimensions>::centroid( v, f, e );

  internal_facet_centroid_ =
    ristra::geometry::shapes::triangle<num_dimensions>::centroid( f, e, mp );

  volume_ = 
    ristra::geometry::shapes::tetrahedron::volume( v, f, e, mp );

  set_boundary( fs.front()->is_boundary() );
  is_right_ = is_right;

  // NOTE internal facet values are not defined....
}


} // namespace burton
} // namespace flecsi_sp
