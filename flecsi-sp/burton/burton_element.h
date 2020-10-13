/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines a base element type.
////////////////////////////////////////////////////////////////////////////////

#pragma once


// user includes
#include <flecsi/topology/mesh_types.h>
#include <flecsi-sp/burton/burton_config.h>
#include <ristra/assertions/errors.h>
#include <ristra/geometry/shapes/geometric_shapes.h>
#include <ristra/geometry/shapes/hexahedron.h>
#include <ristra/geometry/shapes/polygon.h>
#include <ristra/geometry/shapes/polyhedron.h>
#include <ristra/geometry/shapes/quadrilateral.h>
#include <ristra/geometry/shapes/tetrahedron.h>
#include <ristra/geometry/shapes/triangle.h>
#include <ristra/math/general.h>
#include <cmath>
#include <iostream>

namespace flecsi_sp {
namespace burton {

  
namespace detail {

////////////////////////////////////////////////////////////////////////////////
/// \brief compute the minimum length
////////////////////////////////////////////////////////////////////////////////
template< typename VS >
inline auto min_length( VS && vs ) {

  std::decay_t< decltype(vs[0]->coordinates()[0]) > len;
  bool first = true;

  // check each vertex combination
  for ( auto vi : vs ) {
    const auto & pi = vi->coordinates();
    for ( auto vj : vs ) {
      if ( vi == vj ) continue;
      const auto & pj = vj->coordinates();
      auto delta = pi - pj;
      if ( first ) {
        len = abs(delta);
        first = false;
      }
      else {
        len = std::min( abs(delta), len );
      }
    }
  }

  return len;
}

////////////////////////////////////////////////////////////////////////////////
// the list of actual coordinates
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TOPOLOGY, typename E >
auto coordinates( const MESH_TOPOLOGY * mesh, const E * ent )
{
  auto vs = mesh->template entities<0, 0>(ent);
  typename E::point_list_t coords;
  coords.reserve( vs.size() );
  for ( auto v : vs ) coords.emplace_back( v->coordinates() );
  return coords;
}
} // namespace detail



////////////////////////////////////////////////////////////////////////////////
//! \brief This type provides the base for all the primal entities except 
//!        for vertices.
//!
//! \tparam NUM_DIMS  The total number of dimensions.
//! \tparam DIM       The actual dimension of the entity.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t NUM_DIMS, std::size_t DIM  >
struct burton_element_t { };


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh edges, faces, and cells.
//! \tparam N The total number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_element_types_t {
  using edge_t = burton_element_t<N,1>;
  using face_t = burton_element_t<N,N-1>;
  using cell_t = burton_element_t<N,N>;
};

template<>
struct burton_element_types_t<1> {
  using edge_t = burton_vertex_t<1>;
  using face_t = burton_vertex_t<1>;
  using cell_t = burton_element_t<1,1>;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh edges.
//! \tparam N The total number of dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_edge_t = typename burton_element_types_t<N>::edge_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh faces.
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_face_t = typename burton_element_types_t<N>::face_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with
//!        multi-dimensional mesh cells.
//! \tparam N The total number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
using burton_cell_t = typename burton_element_types_t<N>::cell_t;


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        one-dimensional mesh cells.
//! \remark this is a 1d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<1,1>
  : public flecsi::topology::mesh_entity_u<1, burton_config_t<1>::num_domains>
{
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

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the shape type
  using shape_t = config_t::shape_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {};

  // Destructor
  ~burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid ( = midpoint )
  const auto & centroid() const 
  { return midpoint_; };

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; };

  //! the area of the element
  auto area() const
  { return area_; };

  //! the area of the element
  auto volume() const
  { return area_; };

  //! the minimum length in the element ( = area )
  auto min_length() const
  { return area_; }

  //! set the region id
  auto & region()
  { return region_; }

  //! get the region id
  auto region() const
  { return region_; }

  //! return the radius for cylindrical coordinates,
  //! just returns 1.0 for 1D systems
  auto get_radius() const
  { return 1.0;}

  //! the element type
  auto shape() const
  { return shape_; };


  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t>
  create_entities(
    const id_t & cell,
    size_t dim,
    const connectivity_t& conn,
    id_t * entities 
  ) {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  static 
  std::vector<size_t>
  create_bound_entities(
    size_t from_domain, 
    size_t to_domain, 
    size_t dim, 
    const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities
  ) {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
  {
    auto vs = mesh->template entities<0, domain>(this);

    const auto & a = vs[0]->coordinates();
    const auto & b = vs[1]->coordinates();
    midpoint_ = 0.5 * ( a + b );
    area_ = std::abs( b[0] - a[0] );
  }

private:

  //============================================================================
  // Private Data
  //============================================================================
  
  //! the region tag
  size_t region_ = 0;

  //! the geometry informtation
  real_t area_ = 0;
  point_t midpoint_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        two-dimensional mesh edges.
//! \remark this is a 2d edge
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<2,1> :
    public flecsi::topology::mesh_entity_u<1, burton_config_t<2>::num_domains>
{

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  // the constructor
  burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the list of vertices
  auto vertices() const;

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; }

  //! \brief in 2d, this doubles as a face, so the centroid is the same as
  //!    the midpoint
  //! \remark this is only enabled in 2d
  const auto & centroid() const
  { return midpoint_; }

  //! the edge length
  auto length() const
  { return length_; }

  //! in 2d, this doubles as a face, so the area is the same as the length
  //! \remark this is only enabled in 2d
  auto area() const
  { return length_; }

  //! true area of 2d element,differs from "area" when running rz.
  auto true_area() const
  {return true_area_;}

  //! the edge normal
  const auto & normal() const
  { return normal_; }

  //! is this a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //! check if flipped
  template<typename T>
  bool is_flipped(const T & c)
  {
    return owner_id_ != c->global_id().global();
  }
  
  template<typename T>
  void set_owner(const T & c) {
    owner_id_ = c->global_id().global();
  }

  //! \brief update the mesh geometry
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, int alpha=0 , int axis=1)
  {
    using ristra::math::sqr;
    using ristra::math::normal;
    auto vs = mesh->template entities<0, domain>(this);
    const auto & a = vs[0]->coordinates();
    const auto & b = vs[1]->coordinates();
    midpoint_[0] = 0.5*(a[0] + b[0]);
    midpoint_[1] = 0.5*(a[1] + b[1]);
    length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) );
    normal_ = normal( b, a );
    normal_ /= length_;

    if( alpha==0)
      true_area_ = length_;
    else
      true_area_ = 2*M_PI*length_*midpoint_[axis];
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:

  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

  //! the geometric information
  real_t length_ = 0;
  point_t midpoint_ = 0;
  vector_t normal_ = 0;

  real_t true_area_ =0;
  //! owner
  size_t owner_id_ = std::numeric_limits<size_t>::max();
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh edges.
//! \remark this is a 3d edge
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,1> :
    public flecsi::topology::mesh_entity_u<1, burton_config_t<3>::num_domains>
{

  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::array< point_t, 2 >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! a bitfield type
  using bitfield_t = typename config_t::bitfield_t;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type
  using tag_list_t = typename config_t::tag_list_t;

  //============================================================================
  // Constructors
  //============================================================================

  // the constructor
  burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; }

  //! the edge length
  auto length() const
  { return length_; }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! is this a boundary
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

  //! \brief update the mesh geometry
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, int alpha=0 , int axis=1)
  {
    using ristra::math::sqr;
    auto vs = mesh->template entities<0, domain>(this);
    const auto & a = vs[0]->coordinates();
    const auto & b = vs[1]->coordinates();
    midpoint_[0] = 0.5*(a[0] + b[0]);
    midpoint_[1] = 0.5*(a[1] + b[1]);
    midpoint_[2] = 0.5*(a[2] + b[2]);
    length_ = std::sqrt( sqr(a[0]-b[0]) + sqr(a[1]-b[1]) + sqr(a[2]-b[2]) );
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the entity tags
  tag_list_t tags_;

  //! the entity flags
  bitfield_t flags_;

  //! the geometric information
  real_t length_ = 0;
  point_t midpoint_ = 0;

}; // struct burton_element_t<3,1>


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        two-dimensional mesh cells.
//! \remark this is a 2d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<2,2>
  : public flecsi::topology::mesh_entity_u<2, burton_config_t<2>::num_domains>
{
  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<2>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t; 

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the shape type
  using shape_t = config_t::shape_t;

  //! a bitfield type
  using bitfield_t = typename config_t::bitfield_t;
  
  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {};

  // Destructor
  ~burton_element_t() = default;

  // dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  // dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;

  
  //============================================================================
  // Accessors / Modifiers
  //============================================================================
  //! the centroid
  FLECSI_INLINE_TARGET
  const auto & centroid() const 
  { return centroid_; };

  //! the edge midpoint
  FLECSI_INLINE_TARGET
  const auto & midpoint() const
  { return midpoint_; };

  //! the area of the element
  FLECSI_INLINE_TARGET
  auto area() const
  { return area_; };

  //! the area of the element
  FLECSI_INLINE_TARGET
  auto volume() const
  { return volume_; };
  
  FLECSI_INLINE_TARGET
  auto min_length() const
  { return min_length_; }

  //! set the region id
  auto & region()
  { return region_; }

  //! get the region id
  auto region() const
  { return region_; }

  //! return the radius for cylindrical coordinates
  auto get_radius() const
  { return volume_ / area_;}

  //! is this a cell that has a face on the boundary
  bool is_touching_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! \brief set whether this element is touching the boundary
  //! \param [in] is_boundary  True if touching the boundary.
  void set_touching_boundary( bool is_touching_boundary )
  { flags_.set(config_t::bits::boundary, is_touching_boundary); }


  //! the element type
  auto shape() const
  { return shape_; };


  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t>
  create_entities(
    const id_t & cell,
    size_t dim,
    const connectivity_t& conn,
    id_t * entities 
  ) {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  static 
  std::vector<size_t>
  create_bound_entities(
    size_t from_domain, 
    size_t to_domain, 
    size_t dim, 
    const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities
  ) {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //! \brief default geometry -> cartesian
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh, int alpha=0 , int axis=1)
  {
    // get general entity connectivity
    auto vs = mesh->template entities<0, domain>(this);
    auto es = mesh->template entities<1, domain>(this);
    switch (shape_) {

      // the element is a triangle
      case shape_t::triangle: {
        const auto & a = vs[0]->coordinates();
        const auto & b = vs[1]->coordinates();
        const auto & c = vs[2]->coordinates();
        centroid_ =
          ristra::geometry::shapes::triangle<num_dimensions>::centroid(
              a, b, c );
        midpoint_ = 
          ristra::geometry::shapes::triangle<num_dimensions>::midpoint(
              a, b, c );
        area_ = ristra::geometry::shapes::triangle<num_dimensions>::area(
            a, b, c );
        // check the edges first
        min_length_ = abs( a - b );
        min_length_ = std::min( abs( b - c ), min_length_ );
        min_length_ = std::min( abs( c - a ), min_length_ );
        break;
      }

      // the element is a quadrilateral
      case shape_t::quadrilateral: {
        const auto & a = vs[0]->coordinates();
        const auto & b = vs[1]->coordinates();
        const auto & c = vs[2]->coordinates();
        const auto & d = vs[3]->coordinates();
        centroid_ = 
          ristra::geometry::shapes::quadrilateral<num_dimensions>::centroid(
              a, b, c, d );
        midpoint_ = 
          ristra::geometry::shapes::quadrilateral<num_dimensions>::midpoint(
              a, b, c, d );
        area_ = 
          ristra::geometry::shapes::quadrilateral<num_dimensions>::area(
              a, b, c, d );
	
        // check the edges first
        min_length_ = abs( a - b );
        min_length_ = std::min( abs( b - c ), min_length_ );
        min_length_ = std::min( abs( c - d ), min_length_ );
        min_length_ = std::min( abs( d - a ), min_length_ );
        // now check the diagonal
        min_length_ = std::min( abs( a - c ), min_length_ );
        min_length_ = std::min( abs( b - d ), min_length_ );
        break;
      }

      // the element is a polygon
      case shape_t::polygon: {
        auto coords = detail::coordinates(mesh, this);
        centroid_ =
          ristra::geometry::shapes::polygon<num_dimensions>::centroid( coords );
        midpoint_ =
          ristra::geometry::shapes::polygon<num_dimensions>::midpoint( coords );
        area_ =
          ristra::geometry::shapes::polygon<num_dimensions>::area( coords );

        // now check min edge length
        auto vs = mesh->template entities<0, domain>(this);
        min_length_ = detail::min_length( vs );
        break;
      }

      // there should be no unmatched case
      default:
        THROW_RUNTIME_ERROR( "Unknown cell type" );

    } // switch
    //s:w
    volume_ = area_*2*M_PI*centroid_[axis];
    if (alpha == 0){
      volume_ = area_;
    }
  }
  
  //----------------------------------------------------------------------------
  // is  a point in a  cell
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  bool is_inside(const MESH_TOPOLOGY * mesh, const real_t * x)
  {

    auto coords = detail::coordinates(mesh, this);
    auto NumPoints  = coords.size();
    coords.emplace_back( coords.front() );

    for (unsigned i=0; i<NumPoints; ++i) {
      const auto & p0 = coords[i];
      const auto & p1 = coords[i+1];
      auto dx = p1[0] - p0[0];
      auto dy = p1[1] - p0[1];

      auto f = dx * ( x[1] - p0[1] ) - dy * ( x[0] - p0[0] );
      if (f<0) return false;
    }

    return true;
  }

private:

  //============================================================================
  // Private Data
  //============================================================================

  //! the entity flags
  bitfield_t flags_;

  //! the region tag
  size_t region_ = 0;

  //! the geometry informtation
  real_t area_ = 0;
  real_t volume_ = 0;
  point_t midpoint_ = 0;
  point_t centroid_ = 0;
  real_t min_length_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;
};


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh faces.
//! \remark this is a 3d face
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,2>
  : public flecsi::topology::mesh_entity_u<2, burton_config_t<3>::num_domains>
{
  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  //! The point list type.
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! The boundary id type
  using tag_t = typename config_t::tag_t;
  //! The boundary id list type.
  using tag_list_t = typename config_t::tag_list_t;

  //! the shape type
  using shape_t = config_t::shape_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {};

  //! Destructor
  ~burton_element_t() = default;

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! is this a boundary
  bool is_boundary() const
  { return flags_.test( config_t::bits::boundary ); }
  
  //! \brief set whether this element is on the boundary
  //! \param [in] is_boundary  True if on the boundary.
  void set_boundary( bool is_boundary )
  { flags_.set(config_t::bits::boundary, is_boundary); }

  //! return the bitfield flags
  const auto & flags() const { return flags_; }
  auto & flags() { return flags_; }

  //! get all entity tags
  const tag_list_t & tags() const
  { return tags_; }

  //! tag entity
  void tag(const tag_t & tag)
  { tags_.push_back(tag); }

  //! the centroid
  const auto & centroid() const
  { return centroid_; }

  //! the midpoint used in tesselating the element
  const auto & midpoint() const 
  { return midpoint_; }

  //! the normal
  const auto & normal() const 
  { return normal_; }

  //! the area of the element
  auto area() const
  { return area_; }

  auto true_area() const
  {return area_;}

  //! the minimum length in the element
  auto min_length() const
  { return min_length_; }

  //! the element type
  auto shape() const
  { return shape_; };

  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities
  )  {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  std::vector<id_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell_id,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities )
  { 
    THROW_IMPLEMENTED_ERROR(
      "burton_3d_face_t::create_bound_entities has not been implemented"
    );
  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update(  const MESH_TOPOLOGY * mesh, int alpha=0 , int axis=1)
  {
    using ristra::math::abs;
  
    auto vs = mesh->template entities<0, domain>(this);
  
    switch (shape_) {
  
      // the element is a triangle
      case shape_t::triangle: {
        const auto & a = vs[0]->coordinates();
        const auto & b = vs[1]->coordinates();
        const auto & c = vs[2]->coordinates();
        centroid_ =
          ristra::geometry::shapes::triangle<num_dimensions>::centroid(
              a, b, c );
        midpoint_ =
          ristra::geometry::shapes::triangle<num_dimensions>::midpoint(
              a, b, c );
        normal_ =
          ristra::geometry::shapes::triangle<num_dimensions>::normal( a, b, c );
        // check the edges first
        min_length_ = abs( a - b );
        min_length_ = std::min( abs( b - c ), min_length_ );
        min_length_ = std::min( abs( c - a ), min_length_ );
        break;
      }
  
      // the element is a quadrilateral
      case shape_t::quadrilateral: {
        const auto & a = vs[0]->coordinates();
        const auto & b = vs[1]->coordinates();
        const auto & c = vs[2]->coordinates();
        const auto & d = vs[3]->coordinates();
        centroid_ = 
          ristra::geometry::shapes::quadrilateral<num_dimensions>::centroid(
              a, b, c, d );
        midpoint_ = 
          ristra::geometry::shapes::quadrilateral<num_dimensions>::midpoint(
              a, b, c, d );
        normal_ = 
          ristra::geometry::shapes::quadrilateral<num_dimensions>::normal(
              a, b, c, d );
        // check the edges first
        min_length_ = abs( a - b );
        min_length_ = std::min( abs( b - c ), min_length_ );
        min_length_ = std::min( abs( c - d ), min_length_ );
        min_length_ = std::min( abs( d - a ), min_length_ );
        // now check the diagonal
        min_length_ = std::min( abs( a - c ), min_length_ );
        min_length_ = std::min( abs( b - d ), min_length_ );
        break;
      }
  
      // the element is a polygon
      case shape_t::polygon: {
        auto coords = coordinates( mesh );
        centroid_ =
          ristra::geometry::shapes::polygon<num_dimensions>::centroid( coords );
        midpoint_ =
          ristra::geometry::shapes::polygon<num_dimensions>::midpoint( coords );
        normal_ =
          ristra::geometry::shapes::polygon<num_dimensions>::normal( coords );
        min_length_ = detail::min_length( vs );
        break;
      }
  
      // there should be no unmatched case
      default:
        THROW_RUNTIME_ERROR( "Unknown cell type" );
  
    } // switch
    
    area_ = abs( normal_ );
    normal_ /= area_;
  }

  //----------------------------------------------------------------------------
  //! the list of actual coordinates
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  auto coordinates(
    const MESH_TOPOLOGY * mesh,
    bool reverse = false
  ) {
    auto vs = mesh->template entities<0, domain>(this);
    point_list_t coords;
    if ( reverse ) {
      coords.resize( vs.size() );
      size_t cnt = vs.size()-1;
      for ( auto v : vs ) 
        coords[cnt--] = v->coordinates();
    }
    else {
      coords.reserve( vs.size() );
      for ( auto v : vs ) 
        coords.emplace_back( v->coordinates() );     
    }
    return coords;
  }

  //! check if flipped
  template<typename T>
  bool is_flipped(const T & c)
  {
    return owner_id_ != c->global_id().global();
  }
  
  template<typename T>
  void set_owner(const T & c) {
    owner_id_ = c->global_id().global();
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:
  
  //! the entity tags
  tag_list_t tags_;
  
  //! the entity flags
  bitfield_t flags_;

  //! the geometric information
  real_t area_ = 0;
  real_t min_length_ = 0;
  vector_t normal_ = 0;
  point_t centroid_ = 0;
  point_t midpoint_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;

  //! owner
  size_t owner_id_ = std::numeric_limits<size_t>::max();

}; // struct burton_element_t<3,2>


////////////////////////////////////////////////////////////////////////////////
//! \brief An interface for managing geometry and state associated with 
//!        three-dimensional mesh cells.
//! \remark this is a 3d cell
////////////////////////////////////////////////////////////////////////////////
template<>
struct burton_element_t<3,3>
  : public flecsi::topology::mesh_entity_u<3, burton_config_t<3>::num_domains>
{
  //============================================================================
  // Typedefs
  //============================================================================

  //! the mesh traits
  using config_t = burton_config_t<3>;

  //! Number of domains in the burton mesh.
  static constexpr auto num_domains = config_t::num_domains;

  //! Number of domains in the burton mesh.
  static constexpr auto num_dimensions = config_t::num_dimensions;

  //! The domain of the entity
  static constexpr auto domain = 0;

  //! Type of floating point.
  using real_t = typename config_t::real_t;

  //! Type of floating point.
  using size_t = typename config_t::size_t;

  //! Type containing coordinates of the vertex.
  using point_t = typename config_t::point_t;
  using point_list_t = std::vector< point_t >;

  //! Type vector type.
  using vector_t = typename config_t::vector_t;

  //! the flecsi id type
  using id_t = typename config_t::id_t;

  //! The flecsi domain connectivity type.
  using connectivity_t = typename config_t::connectivity_t;

  //! the base vertex type
  using vertex_t = burton_vertex_t<num_dimensions>;

  //! the base edge type
  using edge_t = burton_edge_t<num_dimensions>;

  //! the base edge type
  using face_t = burton_face_t<num_dimensions>;

  //! the base cell type
  using cell_t = burton_element_t<3,3>;

  //! the shape type
  using shape_t = config_t::shape_t;

  //! The bitfield.
  using bitfield_t = typename config_t::bitfield_t;


  //============================================================================
  // Constructors
  //============================================================================

  //! Constructor
  burton_element_t(shape_t shape) : shape_(shape)
  {}

  //! Destructor
  ~burton_element_t() = default;

  //! dissallow copying
  burton_element_t( burton_element_t & ) = delete;
  burton_element_t & operator=( burton_element_t & ) = delete;

  //! dissallow moving
  burton_element_t( burton_element_t && ) = delete;
  burton_element_t & operator=( burton_element_t && ) = delete;


  //============================================================================
  // Accessors / Modifiers
  //============================================================================

  //! the centroid
  const auto & centroid() const 
  { return centroid_; };

  //! the midpoint used in tesselating the element
  const auto & midpoint() const 
  { return midpoint_; };

  //! the area of the element
  auto volume() const
  { return volume_; };

  //! the minimum length in the element
  auto min_length() const
  { return min_length_; }

  //! set the region id
  auto & region()
  { return region_; }

  //! get the region id
  auto region() const
  { return region_; }

  //! the element type
  auto type() const
  { return shape_; };

  //! return the radius for cylindrical coordinates,
  //! just returns 1.0 for 3D systems
  auto get_radius() const
  { return 1.0;}

  //! is this a cell that has a face on the boundary
  bool is_touching_boundary() const
  { return flags_.test( config_t::bits::boundary ); }

  //! \brief set whether this element is touching the boundary
  //! \param [in] is_boundary  True if touching the boundary.
  void set_touching_boundary( bool is_touching_boundary )
  { flags_.set(config_t::bits::boundary, is_touching_boundary); }

  //----------------------------------------------------------------------------
  //! \brief create_entities is a function that creates entities
  //!   of topological dimension dim, using vertices v, and puts the vertices
  //!   in e. See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] cell  The id of the entity in question.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] conn The currently build connectivity of the mesh.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of vertex collections making up the
  //!   entity and b) the number of vertices per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t> create_entities(
    const id_t & cell, size_t dim,
    const connectivity_t& conn,
    id_t * entities
  )  {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief create_bound_entities binds mesh entities across domains.
  //!   See, e.g., burton_quadrilateral_element_t for an implementation of
  //!   this pure virtual function.
  //!
  //! \param[in] from_domain  The domain index of the root entity.
  //! \param[in] to_domain  The domain index of the sub entities.
  //! \param[in] dim The topological dimension of the entity to create.
  //! \param[in] cell_id  The id of the entity in question.
  //! \param[in] primal_conn The currently built connectivity of the primal mesh.
  //! \param[in] domain_conn The currently built connectivity of the mesh for 
  //!                        this entities domain.
  //! \param[out] entities Vector to fill with ids of the lower level entities 
  //!                      that form this entity.
  //!
  //! \return A pair with a) the number of entity collections making up the
  //!   binding and b) the number of entities per collection.
  //----------------------------------------------------------------------------
  std::vector<size_t> create_bound_entities(
    size_t from_domain, size_t to_domain, size_t dim, const id_t & cell,
    const connectivity_t& primal_conn,
    const connectivity_t& domain_conn,
    id_t * entities
  ) {
    THROW_RUNTIME_ERROR("Should not get here");
    return {};
  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update(  const MESH_TOPOLOGY * mesh, int alpha=0 , int axis=1)
  {
    auto vs = mesh->template entities<0, domain>(this);
    auto fs = mesh->template entities<2, domain>(this);
  
    switch (shape_) {
  
      // the element is a tet
      case shape_t::tetrahedron: {
        const auto & v0 = vs[0]->coordinates();
        const auto & v1 = vs[1]->coordinates();
        const auto & v2 = vs[2]->coordinates();
        const auto & v3 = vs[3]->coordinates();
        centroid_ = 
          ristra::geometry::shapes::tetrahedron::centroid( v0, v1, v2, v3 );
        midpoint_ = 
          ristra::geometry::shapes::tetrahedron::midpoint( v0, v1, v2, v3 );
        volume_ = 
          ristra::geometry::shapes::tetrahedron::volume( v0, v1, v2, v3 );
        min_length_ = detail::min_length( vs );
        break;
      }
  
      // the element is a hex
      case shape_t::hexahedron: {
        ristra::geometry::shapes::polyhedron<point_t> poly;     
        const auto & v0 = vs[0]->coordinates();
        const auto & v1 = vs[1]->coordinates();
        const auto & v2 = vs[2]->coordinates();
        const auto & v3 = vs[3]->coordinates();
        const auto & v4 = vs[4]->coordinates();
        const auto & v5 = vs[5]->coordinates();
        const auto & v6 = vs[6]->coordinates();
        const auto & v7 = vs[7]->coordinates();
        poly.insert( {v0, v3, v2, v1} );
        poly.insert( {v1, v2, v6, v5} );
        poly.insert( {v2, v3, v7, v6} );
        poly.insert( {v3, v0, v4, v7} );
        poly.insert( {v0, v1, v5, v4} );
        poly.insert( {v4, v5, v6, v7} );
        centroid_ = poly.centroid();
        midpoint_ = poly.midpoint();
        volume_ = poly.volume();
        min_length_ = detail::min_length( vs );
        break;
      }
  
      // the element is a polyhedron
      case shape_t::polyhedron: {
        ristra::geometry::shapes::polyhedron<point_t> poly;     
        for ( auto f : fs ) {
          auto reverse = !f->is_flipped(this);
          auto coords = f->coordinates( mesh, reverse );
          poly.insert( coords );
        }
        centroid_ = poly.centroid();
        midpoint_ = poly.midpoint();
        volume_ = poly.volume();
        min_length_ = detail::min_length( vs );
        break;
      }
  
      // there should be no unmatched case
      default:
        THROW_RUNTIME_ERROR( "Unknown cell type" );
  
    } // switch
    
  }

  //----------------------------------------------------------------------------
  // is  a point in a  cell
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  bool is_inside(const MESH_TOPOLOGY * mesh, const real_t * x)
  {
    auto vs = mesh->template entities<0, domain>(this);
    auto fs = mesh->template entities<2, domain>(this);

    for ( auto f : fs ) {
      auto cs = mesh->template entities<3, domain>(f);
      auto reverse = !f->is_flipped(this);
      auto coords = f->coordinates( mesh, reverse );
    
      auto NumPoints = coords.size();
      coords.emplace_back(coords.front());

      const auto & fx = f->midpoint();
    
      for (unsigned i=0; i<NumPoints; ++i) {
        const auto & p0 = coords[i];
        const auto & p1 = coords[i+1];
        auto n = ristra::geometry::shapes::triangle<num_dimensions>::normal( p0, p1, fx );
        auto f = 
          n[0] * ( x[0] - p0[0] ) +
          n[1] * ( x[1] - p0[1] ) +
          n[2] * ( x[2] - p0[2] );
        if (f<0) return false;
      }
    
    }

    return true;
  }

  //============================================================================
  // Private Data
  //============================================================================
  
private:

  //! the entity flags
  bitfield_t flags_;  

  //! the region tag
  size_t region_ = 0;

  //! the geometry informtation
  real_t area_ = 0;
  real_t volume_ = 0;
  point_t centroid_ = 0;
  point_t midpoint_ = 0;
  real_t min_length_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;

}; // struct burton_element_t<3,3>


} // namespace burton
} // namespace flecsi_sp
