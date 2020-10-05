/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// library includes
#include <portage/support/portage.h>
#include <wonton/support/Point.h>
#include <wonton/mesh/AuxMeshTopology.h>

// system includes
#include <map>
#include <memory>
#include <vector>
#include <iostream>

namespace flecsi_sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
/// \brief A utility function to convert a type to a portage point
////////////////////////////////////////////////////////////////////////////////
// @{
template<
  typename T, 
  template<typename,std::size_t> class A
>
Wonton::Point<3> make_point( const A<T,3> & a )
{
  return { a[0], a[1], a[2] };
}

template<
  typename T, 
  template<typename,std::size_t> class A
>
Wonton::Point<2> make_point( const A<T,2> & a )
{
  return { a[0], a[1] };
}

template<
  typename T, 
  template<typename,std::size_t> class A
>
Wonton::Point<1> make_point( const A<T,1> & a )
{
  return { a[0] };
}

template<typename T>
Wonton::Point<1> make_1d_point( T && a )
{
  return { std::forward<T>(a)[0] };
}

template<typename T>
Wonton::Point<2> make_2d_point( T && a )
{
  return { std::forward<T>(a)[0], std::forward<T>(a)[1] };
}

template<typename T>
Wonton::Point<3> make_3d_point( T && a )
{
  return {
    std::forward<T>(a)[0], std::forward<T>(a)[1], std::forward<T>(a)[2]
  };
}

// @}


////////////////////////////////////////////////////////////////////////////////
///  \brief Implements a mesh wrapper for Portage mesh queries.
////////////////////////////////////////////////////////////////////////////////
template< typename M >
class portage_mesh_wrapper_t : public Wonton::AuxMeshTopology<portage_mesh_wrapper_t<M>> 
{

public:

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief the auxiliary class
  using portage_mesh_aux_t = Wonton::AuxMeshTopology<portage_mesh_wrapper_t<M>>;

  //! \brief The mesh type
  using mesh_t = M;
  //! \brief the size type
  using size_t = typename mesh_t::size_t;
  //! \brief the real type
  using real_t = typename mesh_t::real_t;
  //! \brief the vector type
  using vector_t = typename mesh_t::vector_t;
  //! the gometric shape 
  using shape_t = typename mesh_t::shape_t;

  using counter_t = typename mesh_t::counter_t;

  using subset_t = typename mesh_t::subset_t;

  //! \brief The entity kind type
  using entity_kind_t = Wonton::Entity_kind;
  //! \brief The entity type 
  using entity_type_t = Wonton::Entity_type;
  //! \brief the portage point type
  using point_t = Wonton::Point< mesh_t::num_dimensions >;
  //! \brief the portage element type
  using element_type_t = Portage::Element_type;

  //! \brief a map type for mapping mesh shapes to portage shapes
  using shape_map_t = std::map< shape_t, element_type_t >;

  //! \brief the 2d portage point type
  using point_1d_t = Wonton::Point< 1 >;
  //! \brief the 2d portage point type
  using point_2d_t = Wonton::Point< 2 >;
  //! \brief the 3d portage point type
  using point_3d_t = Wonton::Point< 3 >;

  //============================================================================
  // Public Static Data
  //============================================================================

  //! \brief the map between 
  static const shape_map_t shapes_to_portage;

  //============================================================================
  // Constructors
  //============================================================================

  //!  \brief Constructor for creating a serial, 3D Cartesian mesh.
  //!  \param[in] mesh The minimum coordinates of the domain.
  explicit portage_mesh_wrapper_t(mesh_t & mesh) :
    cells_(mesh.cells()), faces_(mesh.faces()),
    vertices_(mesh.vertices()), mesh_(&mesh)
  {
    // base class (AuxMeshTopology) method that has to be called here
    // and not in the constructor of the base class because it needs
    // access to methods in this class which in turn need access to
    // its member variables. But these member vars don't get
    // initialized until the base class is constructed
    portage_mesh_aux_t::build_aux_entities(); 
    
    const auto & context = flecsi::execution::context_t::instance();
    vert_local_to_global_id_map_ = 
      &context.index_map( mesh_t::index_spaces_t::vertices );
    face_local_to_global_id_map_ = 
      &context.index_map( mesh_t::index_spaces_t::faces );
    cell_local_to_global_id_map_ =
      &context.index_map( mesh_t::index_spaces_t::cells );

  }

  //! Default constructor deleted
  portage_mesh_wrapper_t() = default;

  //! Default copy constructor
  portage_mesh_wrapper_t(const portage_mesh_wrapper_t &) = default;

  //! Default assignment operator
  portage_mesh_wrapper_t & operator=(const portage_mesh_wrapper_t &) = default;

  void set_new_coordinates(
    const vector_t * node_coordinates,
    const real_t * cell_volumes,
    const vector_t * cell_centroids = nullptr )
  {
    node_coordinates_ = node_coordinates;
    cell_volumes_ = cell_volumes;
    cell_centroids_ = cell_centroids;
  }

  //============================================================================
  // Public Members Required By Portage
  //============================================================================

  //! Dimension of space or mesh points
  constexpr auto space_dimension() const 
  {
    return mesh_t::num_dimensions;
  }

  //! Cell area/volume
  auto cell_volume(size_t cell_id) const
  {
    if ( cell_volumes_ )
      return cell_volumes_[cell_id];
    else
      return cells_[cell_id]->volume();
  }

  //! Number of owned cells in the mesh
  size_t num_owned_cells() const 
  {
    return mesh_->num_cells(flecsi::owned);
  }

  //! Number of owned faces in the mesh
  size_t num_owned_faces() const 
  {
    return mesh_->num_faces(subset_t::overlapping);
  }

  //! Number of owned edges in the mesh
  size_t num_owned_edges() const 
  {
    return mesh_->num_edges(subset_t::overlapping);
  }

  //! Number of owned nodes in the mesh
  size_t num_owned_nodes() const 
  {
    return mesh_->num_vertices(subset_t::overlapping);
  }

  //! Number of ghost cells in the mesh
  size_t num_ghost_cells() const 
  {
    return mesh_->num_cells() - mesh_->num_cells(flecsi::owned);
  }

  //! Number of ghost faces in the mesh
  size_t num_ghost_faces() const 
  {
    return mesh_->num_faces() - mesh_->num_faces(subset_t::overlapping);
  }

  //! Number of ghost edges in the mesh
  size_t num_ghost_edges() const 
  {
    return mesh_->num_edges() - mesh_->num_edges(subset_t::overlapping);
  }

  //! Number of ghost nodes in the mesh
  size_t num_ghost_nodes() const 
  {
    return mesh_->num_vertices() - mesh_->num_vertices(subset_t::overlapping);
  }

  //! Get list of nodes for a cell
  //! \param [in] cell_id  The id of the cell
  //! \param [in,out] nodes  The list of nodes to populate
  template< typename T >
  void cell_get_nodes(size_t cell_id, std::vector<T> * nodes) const 
  {
    auto c = cells_[cell_id];
    nodes->clear();
    for (auto v : mesh_->vertices(c))
      nodes->emplace_back( v.id() );
  }



  //! Get cell faces and the directions in which they are used
  //! \param [in] cell_id  The id of the cell in question
  //! \param [in,out] faces  The list of cell faces to populate
  //! \param [in,out] face_dirs  The list of face directions to populate
  template< typename T >
  void cell_get_faces_and_dirs(
    size_t cell_id, 
    std::vector<T> *faces,
    std::vector<T> *face_dirs
  ) const 
  {
    faces->clear();
    face_dirs->clear();
    auto c = cells_[cell_id];
    for (auto f : mesh_->faces(c)) {
      faces->emplace_back( f.id() );
      face_dirs->emplace_back( f->is_flipped(c) ? -1 : 1 );
    }
  }

  //! Get nodes of a face
  //! \param [in] face_id  The id of the face in question
  //! \param [in,out] nodes  The list of face nodes to populate
  template< typename T >
  void face_get_nodes(
    size_t face_id, 
    std::vector<T> *nodes
  ) const 
  {
    auto f = faces_[face_id];
    nodes->clear();
    for (auto v : mesh_->vertices(f))
      nodes->emplace_back( v.id() );
    //std::reverse( nodes->begin(), nodes->end() );
  }

  //! \brief Get connected cells of given node
  //!
  //! Get connected cells of given node 
  //!
  //! \param [in] node_id  The node index
  //! \param [in] type  The type of indexes to include (ghost, shared, 
  //!   all, etc...) 
  //! \param [in,out] adj_cells  The list of cell neighbors to populate
  //!
  //!  NOTE: For now we are not distinguishing between parallel types
  template< typename T >
  void node_get_cells(
    size_t node_id,
    entity_type_t type,
    std::vector<T> *adj_cells
  ) const 
  {
    adj_cells->clear();
    auto this_node = vertices_[node_id];
    for (auto cell : mesh_->cells(this_node))
      adj_cells->emplace_back(cell.id());
  }


  //!  \brief Check if entity is on the domain boundary
  //!  \param [in] type  The type of indexes to include (ghost, shared, 
  //!   all, etc...) 
  //!  \param[in] id  The ID of the entity.
  //!  \param[out] true/false  depending on if the entity is on the boundary
  //!
  bool on_exterior_boundary(
  			    entity_kind_t entity,
  			    size_t const id) const
  {
    switch(entity) {
    case entity_kind_t::NODE :
      {
  	auto v = vertices_[id];
  	if (v->is_boundary() == true) {
  	  return true;
  	} else {
  	  return false;
  	}
      }
    case entity_kind_t::CELL :
      {
  	auto c = cells_[id];
  	if (c->is_touching_boundary() == true) {
  	  return true;
  	} else {
  	  return false;
  	}
      }
    default :
      {
  	std::cerr << "Unknown Entity Type" << std::endl;
  	return 0;
      }
    }
  }
 

  //!  \brief Get the coords of a node
  //!  \param[in] node_id The ID of the node.
  //!  \param[in,out] pp The Wonton::Point object containing the coordinate
  //!    information.
  //!
  //!  \remark FleCSI Burton specialization doesn't currently fully support 1D.
  void node_get_coordinates(
    size_t node_id, 
    point_1d_t * pp
  ) const 
  {
    auto v = vertices_[node_id];
    if ( node_coordinates_ )
      *pp = make_1d_point( node_coordinates_[v] );
    else
      *pp = make_1d_point( v->coordinates() );
  }

  void node_get_coordinates(
    size_t node_id, 
    point_2d_t * pp
  ) const 
  {
    auto v = vertices_[node_id];
    if ( node_coordinates_ )
      *pp = make_2d_point( node_coordinates_[v] );
    else
      *pp = make_2d_point( v->coordinates() );
  }

  void node_get_coordinates(
    size_t node_id,
    point_3d_t * pp
  ) const
  {
    auto v = vertices_[node_id];
    if ( node_coordinates_ )
      *pp = make_3d_point( node_coordinates_[v] );
    else
      *pp = make_3d_point( v->coordinates() );
  }


  //! \brief Get the coodinates of the nodes of a cell.
  //! \param[in] cell_id The ID of the cell.
  //! \param[in,out] point_list The vector of Wonton::Point objects containing
  //!   the coordinates of a node.  The length of the vector is equal to the 
  //!   number of nodes in the cell with ID @c cellid.
  //! \remark FleCSI Burton specialization doesn't currently fully support 1D.
  void cell_get_coordinates(
    size_t const cell_id,
    std::vector<point_t> *point_list
  )  const 
  {
    // Get this cell object
    auto cell = cells_[cell_id];

    // Loop over vertices of this cell to get their coordinates
    point_list->clear();
    auto verts = mesh_->vertices(cell);

    // new node coordinates
    if ( node_coordinates_ ) {
      for (auto v : verts)
        point_list->emplace_back( make_point(node_coordinates_[v]) );
    }
    // original node coordinates
    else {
      for (auto v : verts)
        point_list->emplace_back( make_point(v->coordinates()) );
    }

  }

  //! \brief Centroid of a cell.
  //! \param[in]  cell_id The ID of the cell.
  //! \param[in,out] centroid The vector of coordinates of the cell @c cellid's
  //!   centroid.  The length of the vector is equal to the dimension of the 
  //!    mesh.
  void cell_centroid(
    size_t cell_id,
    point_t * centroid
  ) const
  {
    auto this_cell = cells_[cell_id];
    if ( cell_centroids_ )
      *centroid = make_point(cell_centroids_[this_cell]);
    else if (node_coordinates_) 
      THROW_RUNTIME_ERROR(
          "Your mesh wrapper has overriden coordinates, but no "
          " volumes have been provided!" );
    else
      *centroid = make_point(this_cell->centroid());
  }


  //============================================================================
  //! \brief Facet a cells face
  //============================================================================
  void cell_get_facetization(
      int const cellid,
      std::vector<std::vector<int>> *facetpoints,
      std::vector<point_3d_t> *points) const
  {
    auto c = cells_[cellid];
    for (auto f : mesh_->faces(c)) {

      const auto & vs = mesh_->vertices(f);
      auto nv = vs.size();
      auto add_one = (nv == 3) ? 0 : 1;
      points->resize( nv + add_one );

      // go backwards cause flipped
      if ( f->is_flipped(c) ) {
        for ( size_t i=nv, j=0; i --> 0; ++j)
          node_get_coordinates(vs[i]->id(), &points->operator[](j));
      }
      // iterate through points forward
      else {
        for ( size_t i=0; i<nv; ++i ) 
          node_get_coordinates(vs[i]->id(), &points->operator[](i));
      }

      // if there is only 3 vertices, it is a planar triangle
      if ( nv == 3 ) {
        facetpoints->resize(1);
        facetpoints->front().reserve(3);
        for ( size_t i=0; i<nv; ++i ) facetpoints->front().emplace_back(i);
      }

      // if there are more than 3 vertices, it might be non-planar,
      // so subdivide the face about the face midpoint
      else {

        // figureout the last point
        for ( size_t i=0; i<nv; ++i )
          points->operator[](nv) += points->operator[](i);
        points->operator[](nv) /= nv;

        // now create facets
        facetpoints->resize(nv);
        for ( int i=0; i<nv-1; ++i )
          facetpoints->operator[](i) = {i, i+1, int(nv)};
        facetpoints->operator[](nv-1) = {int(nv-1), 0, int(nv)};

      } // num verts
    
    } // face
  }

  //============================================================================
  //! \brief decompose a cell into tets
  //============================================================================
  void decompose_cell_into_tets(
    int cellid,
    std::vector<std::array<point_3d_t, 4>> *tcoords,
    const bool planar_hex) const
  {

    // The tet's vertices are ordered in the following way:
    //    3
    //  / | \
    // /  |  \
    // 2--|---1
    //  \ | /
    //    0
    
    point_3d_t xc;
    cell_centroid(cellid, &xc);
      
    auto c = cells_[cellid];
    for (auto f : mesh_->faces(c)) {

      const auto & vs = mesh_->vertices(f);
      auto nv = vs.size();

      // if there is only 3 vertices, it is a planar triangle
      if ( nv == 3 ) {
        // bump size by one
        tcoords->resize( tcoords->size() + 1 );
        auto & new_tet = tcoords->back();
        // set tet vertices, for flipped faces
        if ( f->is_flipped(c) ) {
          for ( size_t i=nv, j=0; i --> 0; ++j)
            node_get_coordinates(vs[i]->id(), &new_tet[j]);
        }
        // not flipped
        else {
          for ( int i=0; i<3; ++i )
            node_get_coordinates(vs[i]->id(), &new_tet[i]);
        }
        // now set last point
        new_tet[3] = xc;
      }

      // if there are more than 3 vertices, it might be non-planar,
      // so subdivide the face about the face midpoint
      else {
        // each edge has a new tet
        tcoords->reserve( tcoords->size() + nv );
        // need to compute face midpoint
        point_3d_t xm, xtmp1, xtmp2;
        for ( auto v : vs ) {
          node_get_coordinates(v->id(), &xtmp1);
          xm += xtmp1;
        }
        xm /= nv;
        // go backwards cause flipped
        if ( f->is_flipped(c) ) {
          node_get_coordinates(vs[nv-1]->id(), &xtmp1);
          for ( size_t i=nv; i --> 1; ) {
            node_get_coordinates(vs[i-1]->id(), &xtmp2);
            std::array<point_3d_t, 4> new_tet = {xtmp2, xtmp1, xm, xc};
            tcoords->emplace_back( std::move(new_tet) );
            xtmp1 = xtmp2;
          }
          node_get_coordinates(vs[nv-1]->id(), &xtmp2);
          std::array<point_3d_t, 4> new_tet = {xtmp2, xtmp1, xm, xc};
          tcoords->emplace_back( std::move(new_tet) );
        }
        // iterate through points forward
        else {
          node_get_coordinates(vs[0]->id(), &xtmp1);
          for ( size_t i=0; i<nv-1; ++i ) {
            node_get_coordinates(vs[i+1]->id(), &xtmp2);
            std::array<point_3d_t, 4> new_tet = {xtmp2, xtmp1, xm, xc};
            tcoords->emplace_back( std::move(new_tet) );
            xtmp1 = xtmp2;
          }
          node_get_coordinates(vs[0]->id(), &xtmp2);
          std::array<point_3d_t, 4> new_tet = {xtmp2, xtmp1, xm, xc};
          tcoords->emplace_back( std::move(new_tet) );
        }

      } // num vertices
    
    } // faces
  }

  //============================================================================
  // Public Members Required By Portage
  // These are only needed in 3D
  //============================================================================

  //! \brief Get the type of the cell - PARALLEL_OWNED or PARALLEL_GHOST
  //!
  //! Assumes a 1-1 correspondence between integer values of the
  //! enum types to avoid switch statements
  //! 
  //! \param [in] cell_id  The index of the cell in question.
  //!
  //! NOTE: Currently FleCSI is not exposing this info
  auto cell_get_type(size_t cell_id) const 
  {
    if ( cell_id < num_owned_cells() )
      return entity_type_t::PARALLEL_OWNED;
    else
      return entity_type_t::PARALLEL_GHOST;
  }

  //! Get the element type of a cell - TRI, QUAD, POLYGON, TET, HEX,
  //! PRISM OR POLYHEDRON

  //! \brief Get the element type of a cell
  //!
  //! Can be one of: TRI, QUAD, POLYGON, TET, HEX,
  //!   PRISM OR POLYHEDRON.
  //! 
  //! \param [in] cell_id  The index of the cell in question.
  auto cell_get_element_type(size_t cell_id) const 
  {
    // search for the type in the map
    /*auto tp = cells_[cell_id]->type();
    auto it = shapes_to_portage.find( tp );
    // if it was found, return the mapped type
    if ( it != shapes_to_portage.end() )
      return it->second;
    // otherwise we dont know what it is
    return element_type_t::UNKNOWN_TOPOLOGY;*/
 
    auto this_cell = cells_[cell_id];
    auto conn = mesh_->vertices(this_cell);
    if (conn.size() == 8)
       return element_type_t::HEX;
    else
      return element_type_t::UNKNOWN_TOPOLOGY;
  }

  //! \brief Get the type of the node - PARALLEL_OWNED or PARALLEL_GHOST
  //!
  //! Assumes a 1-1 correspondence between integer values of the
  //! enum types to avoid switch statements
  //! 
  //! \param [in] cell_id  The index of the node in question.
  //!
  //! NOTE: Currently FleCSI is not exposing this info
  auto node_get_type(size_t node_id) const 
  {
    return entity_type_t::PARALLEL_OWNED;
  }

  //! Get global id */
  int get_global_id(size_t id, entity_kind_t const kind) const
  {
     if (kind == entity_kind_t::NODE)
     {
       return vert_local_to_global_id_map_->at(id);
     }    
     else if (kind == entity_kind_t::FACE)
     { 
       return face_local_to_global_id_map_->at(id);
     }
     else if (kind == entity_kind_t::CELL)
     { 
       return cell_local_to_global_id_map_->at(id);
     }
     else {
       THROW_RUNTIME_ERROR("unknown kind");
       return 0;
     }
  }

  //============================================================================
  // Private Members
  //============================================================================

private:

  //! \brief a pointer to the mesh
  const mesh_t * mesh_ = nullptr;
  //! \brief the list of cells
  decltype( mesh_->cells() ) cells_;
  //! \brief the list of faces
  decltype( mesh_->faces() ) faces_;
  //! \brief the list of veritices
  decltype( mesh_->vertices() ) vertices_;
  
  const vector_t * node_coordinates_ = nullptr;
  const real_t * cell_volumes_ = nullptr;
  const vector_t * cell_centroids_ = nullptr;

  const std::map<size_t, size_t> * vert_local_to_global_id_map_ = nullptr;
  const std::map<size_t, size_t> * face_local_to_global_id_map_ = nullptr; 
  const std::map<size_t, size_t> * cell_local_to_global_id_map_ = nullptr; 
};


} // namespace burton
} // namespace flecsi_sp

