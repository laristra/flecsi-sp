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
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, /* dimension */ 0 );
    auto num_cell_verts = v.size();
    assert( num_cell_verts == 2 );

    entities[ 0 ] = v[ 0 ];
    entities[ 1 ] = v[ 1 ];

    return std::vector<size_t>(2, 1);
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

    // general connecitivity that is needed in all cases
    auto cell_verts = primal_conn.get_entity_vec( cell, /* dimension */ 0 );
    auto num_vertices = cell_verts.size();
    assert( num_vertices == 2 );

    // main counter
    size_t i = 0;


    switch (dim) {
      //------------------------------------------------------------------------
      // Corners ( = wedges in 1D )
    case 0: {
    
      entities[i++] = cell_verts[0];
      entities[i++] = cell_verts[1];

      return std::vector<size_t>( 2, 1 );
    }
      //------------------------------------------------------------------------
      // Failure
    default:
      THROW_RUNTIME_ERROR("Unknown bound entity type");
    } // switch

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
  FLECSI_INLINE_TARGET auto area() const
  { return length_; }

  //! the edge normal
  FLECSI_INLINE_TARGET const auto & normal() const
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
  void update( const MESH_TOPOLOGY * mesh )
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
  void update( const MESH_TOPOLOGY * mesh )
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
  const auto & centroid() const 
  { return centroid_; };

  //! the edge midpoint
  const auto & midpoint() const
  { return midpoint_; };

  //! the area of the element
  auto area() const
  { return area_; };

  //! the area of the element
  FLECSI_INLINE_TARGET auto volume() const
  { return area_; };

  //! the minimum length in the element
  auto min_length() const
  { return min_length_; }

  //! set the region id
  auto & region()
  { return region_; }

  //! get the region id
  auto region() const
  { return region_; }


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
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, /* dimension */ 0 );
    auto num_cell_verts = v.size();

    size_t ind=0;
    for ( int i=0; i<static_cast<int>(num_cell_verts)-1; i++ ) {
      auto vp = v[i];
      auto vn = v[i+1];
      entities[ ind++ ] = vp;
      entities[ ind++ ] = vn;
    }
    entities[ ind++ ] = v[ num_cell_verts-1 ];
    entities[ ind++ ] = v[ 0 ];

    return std::vector<size_t>(num_cell_verts, 2);
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

    // general connecitivity that is needed in all cases
    auto cell_verts = primal_conn.get_entity_vec( cell, /* dimension */ 0 );
    auto cell_edges = primal_conn.get_entity_vec( cell, /* dimension */ 1 );
    auto num_vertices = cell_verts.size();
    auto num_edges = cell_edges.size();
    
    // main counter
    size_t i = 0;

    // a lambda function to locate an edge connected to two points
    auto _find_edge = [&]( const auto & pa, const auto & pb  ) 
    {
      // locate the edge with the two vertices
      auto edge = std::find_if( 
        cell_edges.begin(), cell_edges.end(), 
        [&]( const auto & e ) 
        { 
          auto verts = primal_conn.get_entity_vec( e, /* dim */ 0 );
          assert( verts.size() == 2 && "should be two vertices per edge" );
          return ( (verts[0] == pa && verts[1] == pb) || 
                   (verts[0] == pb && verts[1] == pa) );
        } 
      );
      // make sure we found an edge
      assert( edge != cell_edges.end() && "should have found an edge");
      // return the edge
      return edge;
    };


    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      // The right edge is always first
    case 0: {
    
      // loop over each edge (pair of vertices)
      // there is a corner for each vertex
      auto p1 = std::prev( cell_verts.end() );
      auto p0 = std::prev( p1 );
      auto edge0 = _find_edge( *p0, *p1 );
      for ( auto p2=cell_verts.begin(); p2!=cell_verts.end(); ++p2 ) {
        // get the next edge
        auto edge1 = _find_edge( *p1, *p2 );
        // the first point, this is the right (even) one
        entities[i++] = *p1;
        entities[i++] = *edge1;
        entities[i++] = *edge0;
        // update the iterators
        p0 = p1;
        p1 = p2;
        edge0 = edge1;
      } // foreach vertex

      // rotate the results forward because we were looping over the vertices
      // in a slightly shifter manner
      std::rotate( entities, entities+3, entities+num_vertices*3 );

      return std::vector<size_t>( num_vertices, 3 );
    }
      //------------------------------------------------------------------------
      // wedges
      // The right wedge is always first
    case 1: {

      // get connectivity specific to wedges
      auto cell_corners = domain_conn.get_entity_vec( cell, /* dimension */ 0 ).vec();
      
      // sort cell corners for later intersection
      std::sort( cell_corners.begin(),  cell_corners.end() );

      // temporary storage for found corners
      std::vector<id_t> corners;
      corners.reserve(1);

      // loop over each edge (pair of vertices)
      // there is a wedge for each vertex -> edge combination
      auto p1 = std::prev( cell_verts.end() );
      auto p0 = std::prev( p1 );
      auto edge0 = _find_edge( *p0, *p1 );
      for ( auto p2=cell_verts.begin(); p2!=cell_verts.end(); ++p2 ) {
        // get the next edge
        auto edge1 = _find_edge( *p1, *p2 );
        // get the corner of this point, but also associated with this cell
        auto point_corners =
          domain_conn.get_entity_vec(  *p1, /* dim */ 0 ).vec();
        // sort the lists for intersectinos
        std::sort( point_corners.begin(), point_corners.end() );
        // get the intersections of the sets
        corners.clear();
        std::set_intersection( point_corners.begin(), point_corners.end(),
                               cell_corners.begin(), cell_corners.end(),
                               std::back_inserter(corners));
        assert( corners.size() == 1 );
        const auto & c1 = corners.front();
        // the first point, this is the right (even) one
        entities[i++] = *p1;
        entities[i++] = *edge1;
        entities[i++] = c1;
        // for the next point, this one is the left (odd) one
        entities[i++] = *p1;
        entities[i++] = *edge0;
        entities[i++] = c1;
        // update the iterators
        p0 = p1;
        p1 = p2;
        edge0 = edge1;
      } // foreach vertex
      
      // rotate the results forward because we were looping over the vertices
      // in a slightly shifter manner
      std::rotate( entities, entities+6, entities + 2*num_vertices*3 );


      return std::vector<size_t>(2*num_vertices, 3);

    }
      //------------------------------------------------------------------------
      // sides
    case 2: {

      // get connectivity specific to sides
      // list of wedges associated to this cell.
      auto cell_wedges = domain_conn.get_entity_vec( cell, /* dimension */ 1 ).vec();
      
      // sort cell wedges for later intersection
      std::sort( cell_wedges.begin(),  cell_wedges.end() );

      // temporary storage for found wedges
      std::vector<id_t> wedges;
      wedges.reserve(2);

      // loop over each edge (pair of vertices)
      for(auto e1 = cell_edges.begin(); e1 != cell_edges.end();++e1)
      {
          auto edge_wedges = domain_conn.get_entity_vec(*e1,/*dimension*/1).vec();
          
          // sort cell wedges for later intersection
          std::sort( edge_wedges.begin(),  edge_wedges.end() );

          wedges.clear();
          std::set_intersection( edge_wedges.begin(), edge_wedges.end(),
                                 cell_wedges.begin(), cell_wedges.end(),
                                 std::back_inserter(wedges));

          // it should have two wedges per edges..
          assert(wedges.size()==2);
          auto edge_verts = primal_conn.get_entity_vec( *e1, /* dimension */ 0 );

          for( auto v:edge_verts) 
              entities[i++]=v;
          entities[i++]=*e1;
          for(auto we: wedges)
          {
              entities[i++]= we;
          }//we
      }//e1

      return std::vector<size_t>(num_edges, 5);
    }
        
      //------------------------------------------------------------------------
      // Failure
    default:
      THROW_RUNTIME_ERROR("Unknown bound entity type");
    } // switch

  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
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
  FLECSI_INLINE_TARGET const auto & normal() const 
  { return normal_; }

  //! the area of the element
  FLECSI_INLINE_TARGET auto area() const
  { return area_; }

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
    assert( dim == 1 );

    auto v = conn.get_entity_vec( cell, /* dimension */ 0 );
    auto num_cell_verts = v.size();

    size_t ind=0;
    for ( int i=0; i<static_cast<int>(num_cell_verts)-1; i++ ) {
      auto vp = v[i];
      auto vn = v[i+1];
      entities[ ind++ ] = vp;
      entities[ ind++ ] = vn;
    }
    entities[ ind++ ] = v[ num_cell_verts-1 ];
    entities[ ind++ ] = v[ 0 ];

    return std::vector<size_t>(num_cell_verts, 2);
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
  void update( const MESH_TOPOLOGY * mesh )
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
  FLECSI_INLINE_TARGET auto volume() const
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
    // you should only be coming in here to build edges
    assert( dim == 1 ); 
  
    // get the cell entities
    auto cell_faces = conn.get_entity_vec( cell, /* dimension */ 2 );
    // make sure the faces exist
    assert( cell_faces.size() > 0 && "no cell faces yet" );
    
    // the list of edge pairs
    std::vector< std::pair< id_t, id_t > > cell_edges;
  
    // get the edges of each face
    for ( const auto & face : cell_faces ) { 
      // get all vertices in this face
      auto face_verts = conn.get_entity_vec( face, /* dimension */ 0 );
      assert( face_verts.size() > 2 && "not enough vertices for a valid face" );
      // reserve space 
      cell_edges.reserve( cell_edges.size() + face_verts.size() );
      // add each edge pair to the list
      auto vp = std::prev( face_verts.end() );
      assert( vp != face_verts.begin() && "no vertices in this face" ); 
      for ( auto vn = face_verts.begin(); vn != face_verts.end(); vn++ ) {
        assert( *vn != *vp && "edge has two equal vertices" );
        if ( *vn < *vp ) 
          cell_edges.emplace_back( std::make_pair( *vn, *vp ) );
        else
          cell_edges.emplace_back( std::make_pair( *vp, *vn ) );
        vp = vn;
      }          
    }
  
    // now sort the list of edges
    std::sort( 
      cell_edges.begin(), cell_edges.end(), 
      [](const auto & a, const auto & b) 
      { 
        id_t a1 = a.first;
        id_t b1 = b.first;
        if ( a1 == b1 )
          return ( a.second < b.second );
        else
          return (a1 < b1);
      }
    );
    // remove uniques
    auto end = std::unique( cell_edges.begin(), cell_edges.end() );
    auto num_edges = std::distance( cell_edges.begin(), end );
  
    // copy the unique results to the output array   
    size_t i = 0;
  
    std::for_each( 
      cell_edges.begin(), end, 
      [&](const auto & edge) 
      {
        entities[i++] = edge.first;
        entities[i++] = edge.second;
      }
    );
  
    return std::vector<size_t>( num_edges, 2 );
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
  
    size_t i = 0;
  
    switch (dim) {
      //------------------------------------------------------------------------
      // Corners
      //
      // Take your right hand, its origin is the vertex of the corner.  Curl 
      // your hand from the first edge to the second edge, with the third edge
      // aligned with your thumb.  You hand also curls from the first to the 
      // first to second face, with the third face on the bottom.
      //
    case 0: {
  
      std::vector<size_t> entity_count;
  
      // get the cell entities
      auto cell_verts = primal_conn.get_entity_vec( cell, /* dim */ 0 );
      auto cell_edges = primal_conn.get_entity_vec( cell, /* dim */ 1 ).vec();
      auto cell_faces = primal_conn.get_entity_vec( cell, /* dim */ 2 ).vec();
  
      // sort the edges and faces for intersections later      
      std::sort( cell_edges.begin(), cell_edges.end() );
      std::sort( cell_faces.begin(), cell_faces.end() );
  
      // temparary lists
      std::vector<id_t> edges, faces;
  
      for ( const auto & vert : cell_verts ) { 
        // clear temporary lits
        edges.clear();
        faces.clear();
        // get all entities attached to this vertex
        auto vert_edges = primal_conn.get_entity_vec( vert, /* dim */ 1 ).vec(); 
        auto vert_faces = primal_conn.get_entity_vec( vert, /* dim */ 2 ).vec(); 
        // sort the lists for intersectinos
        std::sort( vert_edges.begin(), vert_edges.end() );
        std::sort( vert_faces.begin(), vert_faces.end() );
        // get the intersections of the sets
        std::set_intersection( vert_edges.begin(), vert_edges.end(),
                               cell_edges.begin(), cell_edges.end(),
                               std::back_inserter(edges));
        std::set_intersection( vert_faces.begin(), vert_faces.end(),
                               cell_faces.begin(), cell_faces.end(),
                               std::back_inserter(faces));
        // add the entities to the list
        entities[i++] = vert;
        for ( auto e : edges ) entities[i++] = e;
        for ( auto f : faces ) entities[i++] = f;
        // add the final number of entities to the list
        entity_count.emplace_back( 1 + edges.size() + faces.size() );
        
      }  // foreach vertex
  
      return entity_count;
    }
      //------------------------------------------------------------------------
      // Wedges
      //
      // There is an even/odd ordering so we 
      // know which way to compute normals.  So edges are defined in pairs
      //
    case 1: {
  
      // get the higher dimensional entities
      auto cell_faces   = primal_conn.get_entity_vec( cell, /* dim */ 2 );
      auto cell_corners = domain_conn.get_entity_vec( cell, /* dim */ 0 ).vec();

      // sort cell corners for later intersections
      std::sort( cell_corners.begin(),  cell_corners.end() );
      
      // a counter for the number of wedges
      size_t num_wedges = 0;
  
      // loop over faces
      for ( const auto & face : cell_faces ) { 
  
        // get the vertices of the face
        auto face_verts = primal_conn.get_entity_vec( face, /* dim */ 0 );
        auto ccw_face_verts =
          std::vector<id_t>( face_verts.begin(), face_verts.end() );
        
        // get the edges of the face
        auto face_edges = primal_conn.get_entity_vec( face, /* dim */ 1 );
  
        // get the cells of this face
        auto face_cells = primal_conn.get_entity_vec( face, /* dim */ 3 );
        // reverse the list of vertices if this is backwards
        if ( face_cells[0] != cell ) 
          std::reverse( ccw_face_verts.begin(), ccw_face_verts.end() );
  
        // a lambda function to locate an edge connected to two points
        auto _find_edge = [&]( const auto & pa, const auto & pb  ) 
        {
          // locate the edge with the two vertices
          auto edge = std::find_if( 
            face_edges.begin(), face_edges.end(), 
            [&]( const auto & e ) 
            { 
              auto verts = primal_conn.get_entity_vec( e, /* dim */ 0 );
              assert( verts.size() == 2 && "should be two vertices per edge" );
              return ( (verts[0] == pa && verts[1] == pb) || 
                       (verts[0] == pb && verts[1] == pa) );
            } 
          );
          // make sure we found an edge
          assert( edge != face_edges.end() );
          // return the edge
          return edge;
        };
          
        // loop over each edge (pair of vertices of the face)
        // there is a wedge for each vertex -> edge -> face combination
        auto p1 = std::prev( ccw_face_verts.end() );
        auto p0 = std::prev( p1 );
        auto edge0 = _find_edge( *p0, *p1 );
        for ( auto p2=ccw_face_verts.begin(); p2!=ccw_face_verts.end(); p2++ ) {
          // get the next edge
          auto edge1 = _find_edge( *p1, *p2 );
          // get the corner of this point, but also associated with this cell
          auto point_corners =
            domain_conn.get_entity_vec(  *p1, /* dim */ 0 ).vec();
          // sort the lists for intersectinos
          std::sort( point_corners.begin(), point_corners.end() );
          // get the intersections of the sets
          std::vector<id_t> corners;
          corners.reserve(1);
          std::set_intersection( point_corners.begin(), point_corners.end(),
                                 cell_corners.begin(), cell_corners.end(),
                                 std::back_inserter(corners));
          assert( corners.size() == 1 );
          const auto & c1 = corners.front();
          // the first point, this is the right (even) one
          entities[i++] = *p1;
          entities[i++] = *edge0;
          entities[i++] = face;
          entities[i++] = c1;
          // for the next point, this one is the left (odd) one
          entities[i++] = *p1;
          entities[i++] = *edge1;
          entities[i++] = face;
          entities[i++] = c1;
          // update the iterators
          p0 = p1;
          p1 = p2;
          edge0 = edge1;
          // increment the wedge counter
          num_wedges += 2;
        }
      } // foreach vertex
  
      return std::vector<size_t>( num_wedges, 4 );
    }
      //------------------------------------------------------------------------
      // sides
    case 2: {

      /////////////////////////////////////////////////////
      // loop over each edge (pair of vertices of the face)
      auto num_sides(0);
      
      // get the cell entities
      auto cell_verts = primal_conn.get_entity_vec( cell, /* dim */ 0 );
      auto cell_edges = primal_conn.get_entity_vec( cell, /* dim */ 1 ).vec();
      auto cell_faces = primal_conn.get_entity_vec( cell, /* dim */ 2 ).vec();

      // get connectivity specific to sides
      // list of wedges associated to this cell.
      auto cell_wedges = domain_conn.get_entity_vec( cell, /* dimension */ 1 ).vec();
      
      // sort cell wedges for later intersection
      // somehow, sort is important.
      // without sort, the intersection fails..
      std::sort( cell_wedges.begin(),  cell_wedges.end() );

      // temporary storage for found wedges
      std::vector<id_t> wedges;
      wedges.reserve(2);

      std::vector<id_t> faces;
      faces.reserve(1);

      // loop over each face.
      for( auto f = cell_faces.begin(); f != cell_faces.end();++f)
      {
        auto face_edges = primal_conn.get_entity_vec( *f, /* dim */ 1 ).vec();

        // list of wedges that attached to this face.
        auto face_wedges = domain_conn.get_entity_vec(*f, 1).vec();
        std::sort(face_wedges.begin(), face_wedges.end());
        
        std::vector<id_t> cell_face_wedges;
        // find intersection, wedges that belongs to this cell and face.
        cell_face_wedges.clear();
        std::set_intersection( face_wedges.begin(), face_wedges.end(),
                               cell_wedges.begin(), cell_wedges.end(),
                               std::back_inserter(cell_face_wedges));

        std::sort(cell_face_wedges.begin(), cell_face_wedges.end());
        assert(cell_face_wedges.size() == face_edges.size()*2);
        
        
        // loop over edges that belongs to this face.
        for(auto e = face_edges.begin(); e != face_edges.end();++e)
        {
          // wedges attached to edge, e.
          auto edge_wedges = domain_conn.get_entity_vec(*e,/*dimension*/1).vec();
          std::sort(edge_wedges.begin(),edge_wedges.end());
          
          
          // find intersection, that belongs this cell/face/edge.
          std::vector<id_t> cell_face_edge_wedges;
          // find intersection, wedges that belongs to this cell and face.
          cell_face_edge_wedges.clear();
          std::set_intersection( cell_face_wedges.begin(), cell_face_wedges.end(),
                                 edge_wedges.begin(), edge_wedges.end(),
                                 std::back_inserter(cell_face_edge_wedges));
          assert(cell_face_edge_wedges.size() == 2);

          // now attach the entities.
          // sides consists of:
          // 2 vertices,
          // 1 edge,
          // 1 face,
          // 2 wedges,
          // 
          auto edge_verts  = primal_conn.get_entity_vec( *e, /* dimension */ 0 ).vec();

          // 2 vertex entries 
          for( auto v=edge_verts.begin(); v != edge_verts.end();++v) 
              entities[i++]=*v;
          
          //1 edge entries
          entities[i++]=*e;

          // 1 face entries
          entities[i++]=*f;

          // 2 wedges.
          for(auto we: cell_face_edge_wedges)
            entities[i++]= we;

          num_sides++;
        }//e
      }//f
      return std::vector<size_t>(num_sides, 6);
    }
        
      //------------------------------------------------------------------------
      // failure
    default:
      THROW_RUNTIME_ERROR("Unknown bound entity type");
      return {};
  
    } // switch
  }

  //----------------------------------------------------------------------------
  //! \brief update the mesh geometry
  //----------------------------------------------------------------------------
  template< typename MESH_TOPOLOGY >
  void update( const MESH_TOPOLOGY * mesh )
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
          auto cs = mesh->template entities<3, domain>(f);
          auto reverse = !f->is_flipped(cs[0]);
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

  //============================================================================
  // Private Data
  //============================================================================
  
private:

  //! the entity flags
  bitfield_t flags_;  

  //! the region tag
  size_t region_ = 0;

  //! the geometry informtation
  real_t volume_ = 0;
  point_t centroid_ = 0;
  point_t midpoint_ = 0;
  real_t min_length_ = 0;

  // the shape parameter
  shape_t shape_ = shape_t::none;

}; // struct burton_element_t<3,3>


} // namespace burton
} // namespace flecsi_sp
