/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief This is the main initialization driver
///////////////////////////////////////////////////////////////////////////////
#pragma once

// user incldues
#include <cinch/logging/cinchlog.h>
#include <ristra/utils/algorithm.h>
#include <ristra/utils/string_utils.h>
#include <flecsi/coloring/dcrs_utils.h>
#include <flecsi/coloring/mpi_communicator.h>
#include <flecsi/coloring/parmetis_colorer.h>
#include <flecsi/execution/execution.h>
#include <flecsi/topology/closure_utils.h>
#include <flecsi-sp/burton/burton_mesh.h>
#include <flecsi-sp/io/exodus_definition.h>
#include <flecsi-sp/utils/char_array.h>
#include <flecsi-sp/utils/types.h>

#ifndef FLECSI_SP_ENABLE_EXODUS
#  error Exodus is needed to build burton specialization.
#endif

// system includes
#include <iostream>

namespace flecsi_sp {
namespace burton {

using dom_dim_t = std::pair<size_t, size_t>;

//! \brief holds some extra mesh info.  
//! This is used to pass information between tlt and spmd initializations.
//! The alternative would be to duplicate some computation.
struct extra_mesh_info_t {
  std::map<dom_dim_t, std::map<dom_dim_t, flecsi::coloring::crs_t>> connectivity;
  std::map<dom_dim_t, std::map<size_t, size_t>> global_to_local;
  std::map<dom_dim_t, std::vector<size_t>> local_to_global;

  std::map<dom_dim_t, std::map<dom_dim_t, flecsi::coloring::crs_t>> ghost_connectivity;
  std::map< dom_dim_t, std::vector<size_t> > ghost_ids;

  std::vector< size_t > cell_distribution;
};
  
using mesh_definition_t =
  flecsi_sp::io::mesh_definition<burton_mesh_t::num_dimensions, burton_mesh_t::real_t>;


namespace globals {

  std::unique_ptr<mesh_definition_t> mesh_def;

  std::unique_ptr<extra_mesh_info_t> extra_mesh_info;

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 1D and 2D version
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TYPE >
void create_cells(
  flecsi_sp::io::mesh_definition<2, burton_mesh_t::real_t> & mesh_def,
  extra_mesh_info_t & extra_mesh_info,
  MESH_TYPE & mesh )
{

  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------

  // some constant expressions
  constexpr auto num_dims = burton_mesh_t::num_dimensions;

  // Alias the index spaces type
  using index_spaces = typename burton_mesh_t::index_spaces_t;

  // alias some other types
  using point_t = typename burton_mesh_t::point_t;
  using vertex_t = typename burton_mesh_t::vertex_t;
  using edge_t = typename burton_mesh_t::edge_t;
  using cell_t = typename burton_mesh_t::cell_t;

  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  //----------------------------------------------------------------------------
  // create the vertices
  //----------------------------------------------------------------------------
  
  const auto & vertex_lid_to_mid = context.index_map( index_spaces::vertices );
  const auto & vertex_mid_to_lid = context.reverse_index_map(
    index_spaces::vertices
  );

  std::vector< vertex_t * > vertices;
  vertices.reserve( vertex_lid_to_mid.size() );
  
  const auto & vert_global2local = mesh_def.global_to_local(vertex_t::dimension);
  const auto & vert_local2global = mesh_def.local_to_global(vertex_t::dimension);

  point_t temp_point{0};
  for(auto & vm: vertex_lid_to_mid) {

    // search this ranks mesh definition for the matching offset
    auto it = vert_global2local.find( vm.second );
    // the vertex exists on this rank so create it
    if ( it != vert_global2local.end() ) {
      // get the point
      mesh_def.vertex(it->second, temp_point.data());
    }
    else {
      // doesnt matter but useful for debugging
      temp_point = 0;
    }
    // now create it, if the vertex did not exist, it doesnt matter what
    // data we give it
    auto v = mesh.create_vertex( temp_point );
    v->global_id().set_global(vm.second);
    vertices.emplace_back(v);
  } // for vertices
  
  //----------------------------------------------------------------------------
  // create the edges
  //----------------------------------------------------------------------------
  
  const auto & edge_coloring = context.coloring( index_spaces::edges );
  
  const auto & edge_lid_to_mid = context.index_map( index_spaces::edges );
  const auto & edge_mid_to_lid = context.reverse_index_map(
    index_spaces::edges
  );

  //const auto & edges = context.coloring( index_spaces::edges );
  std::vector< edge_t * > edges;
  edges.reserve( edge_lid_to_mid.size() );
  
  
  const auto & edge_global2local = mesh_def.global_to_local(edge_t::dimension);
  const auto & edge_local2global = mesh_def.local_to_global(edge_t::dimension);

  const auto & edge_to_vertices =
      mesh_def.entities_crs( edge_t::dimension, vertex_t::dimension );
  
  const auto & ghost_edges = extra_mesh_info.ghost_ids.at({0,edge_t::dimension});
  const auto & ghost_edge2vertex = extra_mesh_info.ghost_connectivity.
    at({0,edge_t::dimension}).at({0,vertex_t::dimension});
  
  // create a list of vertex pointers
  std::vector< vertex_t * > elem_vs;
  std::vector< edge_t * > elem_es;
  
  // assume numbering goes, exlusive, then shared, then ghost
  auto num_owned_edges =
    edge_coloring.exclusive.size() + edge_coloring.shared.size();


  // get the side vertices
  const auto & side_vertices = mesh_def.side_vertices();
  const auto & side_ids = mesh_def.side_ids();
  std::vector<size_t> sides;

  // create the edges
  for(auto & em: edge_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;

    // clear the list
    elem_vs.clear();
    sides.clear();

    // search this ranks mesh definition for the matching offset
    auto it = edge_global2local.find( mid );
    bool is_mine = lid < num_owned_edges;
    
    // the vertex exists on this rank so create it
    if ( is_mine && it != edge_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = edge_to_vertices.at(id);
      // create a list of vertex pointers
      elem_vs.resize( vs.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        elem_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return vertices[ flecsi_id ];
        }
      );

      // convert to global for search
      std::vector<size_t> sorted_vs( vs.begin(), vs.end() );
      for ( auto & v : sorted_vs ) v = vert_local2global[v];
      std::sort(sorted_vs.begin(), sorted_vs.end());

      // find matching sides
      for ( size_t is=0; is<side_vertices.size(); ++is ) {
        const auto & test_vs = side_vertices.at(is);
        bool match{false};
        if ( test_vs[0] < test_vs[1] ) {
          match = test_vs[0] == sorted_vs[0] && test_vs[1] == sorted_vs[1];
        }
        else {
          match = test_vs[1] == sorted_vs[0] && test_vs[0] == sorted_vs[1];
        }
        if (match) {
          sides.emplace_back( side_ids[is] );
        }
      }

    }
    
    // otherwise it is a ghost edge
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_edges.begin(), ghost_edges.end(), mid );
      assert( it != ghost_edges.end() && "Could not find ghost edge id" );
      auto i = std::distance( ghost_edges.begin(), it );
      // reserve space
      const auto & vs = ghost_edge2vertex.at(i);
      elem_vs.reserve( vs.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        elem_vs.emplace_back( vertices[id] );
      }
    }

    // make the edge
    auto new_edge = mesh.template make<edge_t>();
    new_edge->global_id().set_global(mid);
    edges.emplace_back(new_edge);
    // add the connectivity
    mesh.template init_entity<0, edge_t::dimension, vertex_t::dimension>(
        new_edge, elem_vs);
    // add any tags
    std::sort( sides.begin(), sides.end() );
    auto last = std::unique( sides.begin(), sides.end() );
    for ( auto it=sides.begin(); it != last; ++it ) {
      new_edge->tag(*it);
    }
  }


  //----------------------------------------------------------------------------
  // create the cells
  //----------------------------------------------------------------------------

  // get the list of vertices
  const auto & cell_to_vertices =
      mesh_def.entities_crs( cell_t::dimension, vertex_t::dimension );
  const auto & cell_to_edges =
      mesh_def.entities_crs( cell_t::dimension, edge_t::dimension );
  
  const auto & ghost_cells = extra_mesh_info.ghost_ids.at({0,cell_t::dimension});
  const auto & ghost_cell2vertex = extra_mesh_info.ghost_connectivity.
    at({0,cell_t::dimension}).at({0,vertex_t::dimension});
  const auto & ghost_cell2edge = extra_mesh_info.ghost_connectivity.
    at({0,cell_t::dimension}).at({0,edge_t::dimension});
  
  const auto & cell_lid_to_mid = context.index_map( index_spaces::cells );

  const auto & cell_region_ids = mesh_def.region_ids();

  // create the cells
  for(auto & cm: cell_lid_to_mid) {

    auto lid = cm.first;
    auto mid = cm.second;

    // clear the list
    elem_vs.clear();
    elem_es.clear();

    // search this ranks mesh definition for the matching offset
    auto r = flecsi::coloring::rank_owner( extra_mesh_info.cell_distribution, mid );

    // set default region id
    size_t region_id{0};

    // if i am the owner
    if ( r == rank ) {
      // figure out the original local id (note flecsi will have renumberd this)
      auto id = mid - extra_mesh_info.cell_distribution[rank];
      // get the list of vertices
      const auto & vs = cell_to_vertices.at(id);
      const auto & es = cell_to_edges.at(id);
      // create a list of vertex pointers
      elem_vs.resize( vs.size() );
      elem_es.resize( es.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        elem_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        elem_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return edges[ flecsi_id ];
        }
      );
      // get the region id
      region_id = cell_region_ids[id] - 1;
    }
    // otherwise, it is a ghost cell
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_cells.begin(), ghost_cells.end(), cm.second );
      assert( it != ghost_cells.end() && "Could not find ghost cell id" );
      auto i = std::distance( ghost_cells.begin(), it );
      // reserve space
      const auto & vs = ghost_cell2vertex.at(i);
      const auto & es = ghost_cell2edge.at(i);
      elem_vs.reserve( vs.size() );
      elem_es.reserve( es.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        elem_vs.emplace_back( vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        elem_es.emplace_back( edges[id] );
      }
    }

    // create the cell (only need the number of vertices rn)
    // if this is a ghost cell, it is junk for now
    auto c = mesh.create_cell( elem_vs );
    c->global_id().set_global(mid);
    c->region() = region_id;
    // add the connectivity
    mesh.template init_entity<cell_t::domain, edge_t::domain, 
      cell_t::dimension, edge_t::dimension>(c, elem_es);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 3D version
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TYPE >
void create_cells(
  flecsi_sp::io::mesh_definition<3, burton_mesh_t::real_t> & mesh_def,
  extra_mesh_info_t & extra_mesh_info,
  MESH_TYPE & mesh )
{

  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------

  // some constant expressions
  constexpr auto num_dims = burton_mesh_t::num_dimensions;

  // Alias the index spaces type
  using index_spaces = burton_mesh_t::index_spaces_t;

  // alias some other types
  using point_t = burton_mesh_t::point_t;
  using vertex_t = burton_mesh_t::vertex_t;
  using edge_t = burton_mesh_t::edge_t;
  using face_t = burton_mesh_t::face_t;
  using cell_t = burton_mesh_t::cell_t;

  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  //----------------------------------------------------------------------------
  // create the vertices
  //----------------------------------------------------------------------------
  
  const auto & vertex_lid_to_mid = context.index_map( index_spaces::vertices );
  const auto & vertex_mid_to_lid = context.reverse_index_map(
    index_spaces::vertices
  );

  std::vector< vertex_t * > vertices;
  vertices.reserve( vertex_lid_to_mid.size() );
  
  const auto & vert_global2local = mesh_def.global_to_local(vertex_t::dimension);
  const auto & vert_local2global = mesh_def.local_to_global(vertex_t::dimension);

  point_t temp_point{0};
  for(auto & vm: vertex_lid_to_mid) {

    // search this ranks mesh definition for the matching offset
    auto it = vert_global2local.find( vm.second );
    // the vertex exists on this rank so create it
    if ( it != vert_global2local.end() ) {
      // get the point
      mesh_def.vertex(it->second, temp_point.data());
    }
    else {
      // doesnt matter but useful for debugging
      temp_point = 0;
    }
    // now create it, if the vertex did not exist, it doesnt matter what
    // data we give it
    auto v = mesh.create_vertex( temp_point );
    v->global_id().set_global(vm.second);
    vertices.emplace_back(v);
  } // for vertices
  
  //----------------------------------------------------------------------------
  // create the edges
  //----------------------------------------------------------------------------
  
  const auto & edge_lid_to_mid = context.index_map( index_spaces::edges );
  
  const auto & edge_global2local = mesh_def.global_to_local(edge_t::dimension);
  const auto & edge_to_vertices =
      mesh_def.entities_crs( edge_t::dimension, vertex_t::dimension );
  
  const auto & ghost_edges = extra_mesh_info.ghost_ids.at({0,edge_t::dimension});
  const auto & ghost_edge2vertex = extra_mesh_info.ghost_connectivity.
    at({0,edge_t::dimension}).at({0,vertex_t::dimension});
  
  // create a list of vertex pointers
  std::vector< vertex_t * > elem_vs;

  // create the edges
  for(auto & em: edge_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;

    // clear the list
    elem_vs.clear();

    // search this ranks mesh definition for the matching offset
    auto it = edge_global2local.find( mid );
    
    // the vertex exists on this rank so create it
    if ( it != edge_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = edge_to_vertices.at(id);
      // create a list of vertex pointers
      elem_vs.resize( vs.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        elem_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return vertices[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost edge
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_edges.begin(), ghost_edges.end(), mid );
      assert( it != ghost_edges.end() && "Could not find ghost edge id" );
      auto i = std::distance( ghost_edges.begin(), it );
      // reserve space
      const auto & vs = ghost_edge2vertex.at(i);
      elem_vs.reserve( vs.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        elem_vs.emplace_back( vertices[id] );
      }
    }

    // make the edge
    auto new_edge = mesh.template make<edge_t>();
    new_edge->global_id().set_global(mid);
    // add the connectivity
    mesh.template init_entity<0, edge_t::dimension, vertex_t::dimension>(
        new_edge, elem_vs);
  }

  //----------------------------------------------------------------------------
  // create the faces
  //----------------------------------------------------------------------------

  const auto & face_coloring = context.coloring( index_spaces::faces );
  
  const auto & face_lid_to_mid = context.index_map( index_spaces::faces );
  
  const auto & face_global2local = mesh_def.global_to_local(face_t::dimension);
  const auto & face_to_vertices =
      mesh_def.entities_crs( face_t::dimension, vertex_t::dimension );
  
  const auto & ghost_faces = extra_mesh_info.ghost_ids.at({0,face_t::dimension});
  const auto & ghost_face2vertex = extra_mesh_info.ghost_connectivity.
    at({0,face_t::dimension}).at({0,vertex_t::dimension});
  
  const auto & face_owners = mesh_def.face_owners();
  
  std::vector< face_t * > faces;
  faces.reserve( face_lid_to_mid.size() );

  // assume numbering goes, exlusive, then shared, then ghost
  auto num_owned_faces =
    face_coloring.exclusive.size() + face_coloring.shared.size();
  
  // get the side vertices
  const auto & side_vertices = mesh_def.side_vertices();
  const auto & side_ids = mesh_def.side_ids();
  std::vector<size_t> sides;

  auto num_sides = side_vertices.size();
  std::vector< std::vector<size_t> > sorted_side_vertices( num_sides );
  for ( size_t i=0; i<num_sides; ++i ) {
    const auto & vs = side_vertices.at(i);
    sorted_side_vertices[i].assign( vs.begin(), vs.end() );
    std::sort( sorted_side_vertices[i].begin(), sorted_side_vertices[i].end() );
  }
  
  // create the faces
  for(auto & em: face_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;

    // clear the list
    elem_vs.clear();
    sides.clear();

    // search this ranks mesh definition for the matching offset
    auto it = face_global2local.find( mid );
    bool is_mine = lid < num_owned_faces;
    
    // the vertex exists on this rank so create it
    if ( is_mine && it != face_global2local.end() ) {

      // get the list of vertices
      auto id = it->second;
      const auto & vs = face_to_vertices.at(id);
      // create a list of vertex pointers
      elem_vs.resize( vs.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        elem_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return vertices[ flecsi_id ];
        }
      );
      // check if this cell owns this face, if it doesnt, rotate the
      // vertices
      //if ( face_owners[id] != cell_local2global[cell_id] )
      //  std::reverse( elem_vs.begin(), elem_vs.end() );
      
      // convert to global for search
      std::vector<size_t> sorted_vs( vs.begin(), vs.end() );
      for ( auto & v : sorted_vs ) v = vert_local2global[v];
      std::sort(sorted_vs.begin(), sorted_vs.end());
      
      // find matching sides
      for ( size_t is=0; is<sorted_side_vertices.size(); ++is ) {
        const auto & test_vs = sorted_side_vertices.at(is);
        if ( sorted_vs.size() == test_vs.size() &&
             std::equal( sorted_vs.begin(), sorted_vs.end(), test_vs.begin() ) )
          sides.emplace_back( side_ids[is] );
      }
    }
    
    // otherwise it is a ghost face
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_faces.begin(), ghost_faces.end(), mid );
      assert( it != ghost_faces.end() && "Could not find ghost face id" );
      auto i = std::distance( ghost_faces.begin(), it );
      // reserve space
      const auto & vs = ghost_face2vertex.at(i);
      elem_vs.reserve( vs.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        elem_vs.emplace_back( vertices[id] );
      }
    }

    // make the face
    auto f = mesh.create_face(elem_vs);
    f->global_id().set_global(mid);
    faces.emplace_back( f );
    // add any tags
    std::sort( sides.begin(), sides.end() );
    auto last = std::unique( sides.begin(), sides.end() );
    for ( auto it=sides.begin(); it != last; ++it ) {
      f->tag(*it);
    }
  }


  //----------------------------------------------------------------------------
  // create the cells
  //----------------------------------------------------------------------------

  // get the list of vertices
  const auto & cell_to_vertices =
      mesh_def.entities_crs( cell_t::dimension, vertex_t::dimension );

  const auto & ghost_cells = extra_mesh_info.ghost_ids.at({0,cell_t::dimension});
  const auto & ghost_cell2vertex = extra_mesh_info.ghost_connectivity.
    at({0,cell_t::dimension}).at({0,vertex_t::dimension});
  
  const auto & cell_lid_to_mid = context.index_map( index_spaces::cells );
  
  const auto & cell_region_ids = mesh_def.region_ids();

  // create the cells
  for(auto & cm: cell_lid_to_mid) {

    auto lid = cm.first;
    auto mid = cm.second;

    // clear the list
    elem_vs.clear();

    // set default region id
    size_t region_id{0};

    // search this ranks mesh definition for the matching offset
    auto r = flecsi::coloring::rank_owner( extra_mesh_info.cell_distribution, mid );

    // if i am the owner
    if ( r == rank ) {
      // figure out the original local id (note flecsi will have renumberd this)
      auto id = mid - extra_mesh_info.cell_distribution[rank];
      // get the list of vertices
      const auto & vs = cell_to_vertices.at(id);
      // create a list of vertex pointers
      elem_vs.resize( vs.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        elem_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return vertices[ flecsi_id ];
        }
      );
      // get the region id
      region_id = cell_region_ids[id]-1;
    }
    // otherwise, it is a ghost cell
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_cells.begin(), ghost_cells.end(), cm.second );
      assert( it != ghost_cells.end() && "Could not find ghost cell id" );
      auto i = std::distance( ghost_cells.begin(), it );
      // reserve space
      const auto & vs = ghost_cell2vertex.at(i);
      elem_vs.reserve( vs.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        elem_vs.emplace_back( vertices[id] );
      }
    }

    // create the cell (only need the number of vertices rn)
    // if this is a ghost cell, it is junk for now
    auto c = mesh.create_cell( elem_vs );
    c->global_id().set_global(mid);
    c->region() = region_id;
  }

}


////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 1D and 2D version
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TYPE >
void create_extras(
  flecsi_sp::io::mesh_definition<2, burton_mesh_t::real_t> & mesh_def,
  extra_mesh_info_t & extra_mesh_info,
  MESH_TYPE & mesh )
{

  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------

  // some constant expressions
  constexpr auto num_dims = burton_mesh_t::num_dimensions;

  // Alias the index spaces type
  using index_spaces = burton_mesh_t::index_spaces_t;

  // alias some other types
  using point_t = typename burton_mesh_t::point_t;
  using vertex_t = typename burton_mesh_t::vertex_t;
  using edge_t = typename burton_mesh_t::edge_t;
  using cell_t = typename burton_mesh_t::cell_t;
  using corner_t = typename burton_mesh_t::corner_t;
  using wedge_t = typename burton_mesh_t::wedge_t;
  using side_t = typename burton_mesh_t::side_t;

  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  const auto & vert_local2global = mesh_def.local_to_global(vertex_t::dimension);
  const auto & edge_local2global = mesh_def.local_to_global(edge_t::dimension);
  const auto & wedge_local2global = extra_mesh_info.local_to_global.
    at({wedge_t::domain, wedge_t::dimension});
  const auto & corner_local2global = extra_mesh_info.local_to_global.
    at({corner_t::domain, corner_t::dimension});
  
  const auto & vertex_mid_to_lid = context.reverse_index_map(index_spaces::vertices);
  const auto & edge_mid_to_lid = context.reverse_index_map(index_spaces::edges);
  const auto & cell_mid_to_lid = context.reverse_index_map(index_spaces::cells);
  const auto & wedge_mid_to_lid = context.reverse_index_map(index_spaces::wedges);
  const auto & corner_mid_to_lid = context.reverse_index_map(index_spaces::corners);


  std::vector<vertex_t*> entity_vs;
  std::vector<edge_t*> entity_es;
  std::vector<cell_t*> entity_cs;
  std::vector<wedge_t*> entity_ws;
  std::vector<corner_t*> entity_ns;

  auto vertices = mesh.template get_entities<vertex_t::dimension, vertex_t::domain>();
  auto edges = mesh.template get_entities<edge_t::dimension, edge_t::domain>();
  auto cells = mesh.template get_entities<cell_t::dimension, cell_t::domain>();
  
  //----------------------------------------------------------------------------
  // create the wedges
  //----------------------------------------------------------------------------
  
  const auto & wedge_lid_to_mid = context.index_map( index_spaces::wedges );
  std::vector<wedge_t*> wedges; wedges.reserve(wedge_lid_to_mid.size());
  
  
  auto wedge_pair = std::make_pair(wedge_t::domain, wedge_t::dimension);
  const auto & wedge_global2local = extra_mesh_info.global_to_local.at(wedge_pair);
  const auto & wedge_conn = extra_mesh_info.connectivity.at(wedge_pair);
  const auto & wedge_to_vertices = wedge_conn.at({0,0});
  const auto & wedge_to_edges    = wedge_conn.at({0,1});
  const auto & wedge_to_cells    = wedge_conn.at({0,2});
  
  const auto & ghost_wedges = extra_mesh_info.ghost_ids.at(wedge_pair);
  const auto & ghost_wedge_conn = extra_mesh_info.ghost_connectivity.at(wedge_pair);
  const auto & ghost_wedge_to_vertices = ghost_wedge_conn.at({0,0});
  const auto & ghost_wedge_to_edges = ghost_wedge_conn.at({0,1});
  const auto & ghost_wedge_to_cells = ghost_wedge_conn.at({0,2});
  
 
  // create the wedges
  for(auto & em: wedge_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;
    
    // clear the lists
    entity_vs.clear();
    entity_es.clear();
    entity_cs.clear();

    
    // search this ranks mesh definition for the matching offset
    auto it = wedge_global2local.find( mid );
    
    // the wedge exists on this rank so create it
    if ( it != wedge_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = wedge_to_vertices.at(id);
      const auto & es = wedge_to_edges.at(id);
      const auto & cs = wedge_to_cells.at(id);
      // create a list of vertex pointers
      entity_vs.resize( vs.size() );
      entity_es.resize( es.size() );
      entity_cs.resize( cs.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        entity_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return &vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        entity_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return &edges[ flecsi_id ];
        }
      );
      std::transform(
        cs.begin(),
        cs.end(),
        entity_cs.begin(),
        [&](auto c) {
          auto global_id = extra_mesh_info.cell_distribution[rank] + c;
          auto flecsi_id = cell_mid_to_lid.at(global_id);
          return &cells[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost wedge
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_wedges.begin(), ghost_wedges.end(), mid );
      assert( it != ghost_wedges.end() && "Could not find ghost wedge id" );
      auto i = std::distance( ghost_wedges.begin(), it );
      // reserve space
      const auto & vs = ghost_wedge_to_vertices.at(i);
      const auto & es = ghost_wedge_to_edges.at(i);
      const auto & cs = ghost_wedge_to_cells.at(i);
      entity_vs.reserve( vs.size() );
      entity_es.reserve( es.size() );
      entity_cs.reserve( cs.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        entity_vs.emplace_back( &vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        entity_es.emplace_back( &edges[id] );
      }
      for ( auto c : cs ) {
        auto id = cell_mid_to_lid.at(c);
        entity_cs.emplace_back( &cells[id] );
      }
    }

    // make the wedge
    auto new_wedge = mesh.template make<wedge_t, wedge_t::domain>();
    new_wedge->global_id().set_global(mid);
    wedges.emplace_back( new_wedge );
    // add the connectivity
    mesh.template init_entity<wedge_t::domain, vertex_t::domain, 
      wedge_t::dimension, vertex_t::dimension>(new_wedge, entity_vs);
    mesh.template init_entity<wedge_t::domain, edge_t::domain,
      wedge_t::dimension, edge_t::dimension>(new_wedge, entity_es);
    mesh.template init_entity<wedge_t::domain, cell_t::domain,
      wedge_t::dimension, cell_t::dimension>(new_wedge, entity_cs);
  }

  //----------------------------------------------------------------------------
  // create the corners
  //----------------------------------------------------------------------------
  
  const auto & corner_lid_to_mid = context.index_map( index_spaces::corners );
  std::vector<corner_t*> corners; corners.reserve(corner_lid_to_mid.size());
  
  auto corner_pair = std::make_pair(corner_t::domain, corner_t::dimension);
  const auto & corner_global2local = extra_mesh_info.global_to_local.at(corner_pair);
  const auto & corner_conn = extra_mesh_info.connectivity.at(corner_pair);
  const auto & corner_to_vertices = corner_conn.at({0,0});
  const auto & corner_to_edges    = corner_conn.at({0,1});
  const auto & corner_to_cells    = corner_conn.at({0,2});
  const auto & corner_to_wedges   = corner_conn.at(wedge_pair);
  
  const auto & ghost_corners = extra_mesh_info.ghost_ids.at(corner_pair);
  const auto & ghost_corner_conn = extra_mesh_info.ghost_connectivity.at(corner_pair);
  const auto & ghost_corner_to_vertices = ghost_corner_conn.at({0,0});
  const auto & ghost_corner_to_edges = ghost_corner_conn.at({0,1});
  const auto & ghost_corner_to_cells = ghost_corner_conn.at({0,2});
  const auto & ghost_corner_to_wedges = ghost_corner_conn.at(wedge_pair);
  
 
  // create the corners
  for(auto & em: corner_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;
    
    // clear the lists
    entity_vs.clear();
    entity_es.clear();
    entity_cs.clear();
    entity_ws.clear();

    
    // search this ranks mesh definition for the matching offset
    auto it = corner_global2local.find( mid );
    
    // the corner exists on this rank so create it
    if ( it != corner_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = corner_to_vertices.at(id);
      const auto & es = corner_to_edges.at(id);
      const auto & cs = corner_to_cells.at(id);
      const auto & ws = corner_to_wedges.at(id);
      // create a list of vertex pointers
      entity_vs.resize( vs.size() );
      entity_es.resize( es.size() );
      entity_cs.resize( cs.size() );
      entity_ws.resize( ws.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        entity_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return &vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        entity_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return &edges[ flecsi_id ];
        }
      );
      std::transform(
        cs.begin(),
        cs.end(),
        entity_cs.begin(),
        [&](auto c) {
          auto global_id = extra_mesh_info.cell_distribution[rank] + c;
          auto flecsi_id = cell_mid_to_lid.at(global_id);
          return &cells[ flecsi_id ];
        }
      );
      std::transform(
        ws.begin(),
        ws.end(),
        entity_ws.begin(),
        [&](auto w) {
          auto global_id = wedge_local2global[w];
          auto flecsi_id = wedge_mid_to_lid.at(global_id);
          return wedges[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost corner
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_corners.begin(), ghost_corners.end(), mid );
      assert( it != ghost_corners.end() && "Could not find ghost corner id" );
      auto i = std::distance( ghost_corners.begin(), it );
      // reserve space
      const auto & vs = ghost_corner_to_vertices.at(i);
      const auto & es = ghost_corner_to_edges.at(i);
      const auto & cs = ghost_corner_to_cells.at(i);
      const auto & ws = ghost_corner_to_wedges.at(i);
      entity_vs.reserve( vs.size() );
      entity_es.reserve( es.size() );
      entity_cs.reserve( cs.size() );
      entity_ws.reserve( ws.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        entity_vs.emplace_back( &vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        entity_es.emplace_back( &edges[id] );
      }
      for ( auto c : cs ) {
        auto id = cell_mid_to_lid.at(c);
        entity_cs.emplace_back( &cells[id] );
      }
      for ( auto w : ws ) {
        auto id = wedge_mid_to_lid.at(w);
        entity_ws.emplace_back( wedges[id] );
      }
    }

    // make the corner
    auto new_corner = mesh.template make<corner_t, corner_t::domain>();
    new_corner->global_id().set_global(mid);
    corners.emplace_back( new_corner );
    // add the connectivity
    mesh.template init_entity<corner_t::domain, vertex_t::domain, 
      corner_t::dimension, vertex_t::dimension>(new_corner, entity_vs);
    mesh.template init_entity<corner_t::domain, edge_t::domain, 
      corner_t::dimension, edge_t::dimension>(new_corner, entity_es);
    mesh.template init_entity<corner_t::domain, cell_t::domain, 
      corner_t::dimension, cell_t::dimension>(new_corner, entity_cs);
    mesh.template init_entity<corner_t::domain, wedge_t::domain, 
      corner_t::dimension, wedge_t::dimension>(new_corner, entity_ws);
  }


  //----------------------------------------------------------------------------
  // create the sides
  //----------------------------------------------------------------------------
  
  const auto & side_lid_to_mid = context.index_map( index_spaces::sides );
  std::vector<side_t*> sides; sides.reserve(side_lid_to_mid.size());
  
  auto side_pair = std::make_pair(side_t::domain, side_t::dimension);
  const auto & side_global2local = extra_mesh_info.global_to_local.at(side_pair);
  const auto & side_conn = extra_mesh_info.connectivity.at(side_pair);
  const auto & side_to_vertices = side_conn.at({0,0});
  const auto & side_to_edges    = side_conn.at({0,1});
  const auto & side_to_cells    = side_conn.at({0,2});
  const auto & side_to_wedges   = side_conn.at(wedge_pair);
  const auto & side_to_corners  = side_conn.at(corner_pair);
  
  const auto & ghost_sides = extra_mesh_info.ghost_ids.at(side_pair);
  const auto & ghost_side_conn = extra_mesh_info.ghost_connectivity.at(side_pair);
  const auto & ghost_side_to_vertices = ghost_side_conn.at({0,0});
  const auto & ghost_side_to_edges = ghost_side_conn.at({0,1});
  const auto & ghost_side_to_cells = ghost_side_conn.at({0,2});
  const auto & ghost_side_to_wedges = ghost_side_conn.at(wedge_pair);
  const auto & ghost_side_to_corners = ghost_side_conn.at(corner_pair);
  
 
  // create the sides
  for(auto & em: side_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;
    
    // clear the lists
    entity_vs.clear();
    entity_es.clear();
    entity_cs.clear();
    entity_ws.clear();
    entity_ns.clear();

    
    // search this ranks mesh definition for the matching offset
    auto it = side_global2local.find( mid );
    
    // the side exists on this rank so create it
    if ( it != side_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = side_to_vertices.at(id);
      const auto & es = side_to_edges.at(id);
      const auto & cs = side_to_cells.at(id);
      const auto & ws = side_to_wedges.at(id);
      const auto & ns = side_to_corners.at(id);
      // create a list of vertex pointers
      entity_vs.resize( vs.size() );
      entity_es.resize( es.size() );
      entity_cs.resize( cs.size() );
      entity_ws.resize( ws.size() );
      entity_ns.resize( ns.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        entity_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return &vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        entity_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return &edges[ flecsi_id ];
        }
      );
      std::transform(
        cs.begin(),
        cs.end(),
        entity_cs.begin(),
        [&](auto c) {
          auto global_id = extra_mesh_info.cell_distribution[rank] + c;
          auto flecsi_id = cell_mid_to_lid.at(global_id);
          return &cells[ flecsi_id ];
        }
      );
      std::transform(
        ws.begin(),
        ws.end(),
        entity_ws.begin(),
        [&](auto w) {
          auto global_id = wedge_local2global[w];
          auto flecsi_id = wedge_mid_to_lid.at(global_id);
          return wedges[ flecsi_id ];
        }
      );
      std::transform(
        ns.begin(),
        ns.end(),
        entity_ns.begin(),
        [&](auto c) {
          auto global_id = corner_local2global[c];
          auto flecsi_id = corner_mid_to_lid.at(global_id);
          return corners[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost side
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_sides.begin(), ghost_sides.end(), mid );
      assert( it != ghost_sides.end() && "Could not find ghost side id" );
      auto i = std::distance( ghost_sides.begin(), it );
      // reserve space
      const auto & vs = ghost_side_to_vertices.at(i);
      const auto & es = ghost_side_to_edges.at(i);
      const auto & cs = ghost_side_to_cells.at(i);
      const auto & ws = ghost_side_to_wedges.at(i);
      const auto & ns = ghost_side_to_corners.at(i);
      entity_vs.reserve( vs.size() );
      entity_es.reserve( es.size() );
      entity_cs.reserve( cs.size() );
      entity_ws.reserve( ws.size() );
      entity_ns.reserve( ns.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        entity_vs.emplace_back( &vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        entity_es.emplace_back( &edges[id] );
      }
      for ( auto c : cs ) {
        auto id = cell_mid_to_lid.at(c);
        entity_cs.emplace_back( &cells[id] );
      }
      for ( auto w : ws ) {
        auto id = wedge_mid_to_lid.at(w);
        entity_ws.emplace_back( wedges[id] );
      }
      for ( auto c : ns ) {
        auto id = corner_mid_to_lid.at(c);
        entity_ns.emplace_back( corners[id] );
      }
    }

    // make the side
    auto new_side = mesh.template make<side_t, side_t::domain>();
    new_side->global_id().set_global(mid);
    // add the connectivity
    mesh.template init_entity<side_t::domain, vertex_t::domain, 
      side_t::dimension, vertex_t::dimension>(new_side, entity_vs);
    mesh.template init_entity<side_t::domain, edge_t::domain, 
      side_t::dimension, edge_t::dimension>(new_side, entity_es);
    mesh.template init_entity<side_t::domain, cell_t::domain, 
      side_t::dimension, cell_t::dimension>(new_side, entity_cs);
    mesh.template init_entity<side_t::domain, wedge_t::domain, 
      side_t::dimension, wedge_t::dimension>(new_side, entity_ws);
    mesh.template init_entity<side_t::domain, corner_t::domain, 
      side_t::dimension, corner_t::dimension>(new_side, entity_ns);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the cells
///
/// \remarks This is the 1D and 2D version
////////////////////////////////////////////////////////////////////////////////
template< typename MESH_TYPE >
void create_extras(
  flecsi_sp::io::mesh_definition<3, burton_mesh_t::real_t> & mesh_def,
  extra_mesh_info_t & extra_mesh_info,
  MESH_TYPE & mesh )
{

  //----------------------------------------------------------------------------
  // Initial compile time setup
  //----------------------------------------------------------------------------

  // some constant expressions
  constexpr auto num_dims = burton_mesh_t::num_dimensions;

  // Alias the index spaces type
  using index_spaces = burton_mesh_t::index_spaces_t;

  // alias some other types
  using point_t = burton_mesh_t::point_t;
  using vertex_t = burton_mesh_t::vertex_t;
  using edge_t = burton_mesh_t::edge_t;
  using face_t = burton_mesh_t::face_t;
  using cell_t = burton_mesh_t::cell_t;
  using corner_t = burton_mesh_t::corner_t;
  using wedge_t = burton_mesh_t::wedge_t;
  using side_t = burton_mesh_t::side_t;

  //----------------------------------------------------------------------------
  // Get context information
  //----------------------------------------------------------------------------

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  const auto & vert_local2global = mesh_def.local_to_global(vertex_t::dimension);
  const auto & edge_local2global = mesh_def.local_to_global(edge_t::dimension);
  const auto & face_local2global = mesh_def.local_to_global(face_t::dimension);
  const auto & wedge_local2global = extra_mesh_info.local_to_global.at({wedge_t::domain, wedge_t::dimension});
  const auto & corner_local2global = extra_mesh_info.local_to_global.at({corner_t::domain, corner_t::dimension});
  
  const auto & vertex_mid_to_lid = context.reverse_index_map(index_spaces::vertices);
  const auto & edge_mid_to_lid = context.reverse_index_map(index_spaces::edges);
  const auto & face_mid_to_lid = context.reverse_index_map(index_spaces::faces);
  const auto & cell_mid_to_lid = context.reverse_index_map(index_spaces::cells);
  const auto & wedge_mid_to_lid = context.reverse_index_map(index_spaces::wedges);
  const auto & corner_mid_to_lid = context.reverse_index_map(index_spaces::corners);


  std::vector<vertex_t*> entity_vs;
  std::vector<edge_t*> entity_es;
  std::vector<face_t*> entity_fs;
  std::vector<cell_t*> entity_cs;
  std::vector<wedge_t*> entity_ws;
  std::vector<corner_t*> entity_ns;

  auto vertices = mesh.template get_entities<vertex_t::dimension, vertex_t::domain>();
  auto edges = mesh.template get_entities<edge_t::dimension, edge_t::domain>();
  auto faces = mesh.template get_entities<face_t::dimension, face_t::domain>();
  auto cells = mesh.template get_entities<cell_t::dimension, cell_t::domain>();
  
  //----------------------------------------------------------------------------
  // create the wedges
  //----------------------------------------------------------------------------
  
  const auto & wedge_lid_to_mid = context.index_map( index_spaces::wedges );
  std::vector<wedge_t*> wedges; wedges.reserve(wedge_lid_to_mid.size());
  
  
  auto wedge_pair = std::make_pair(wedge_t::domain, wedge_t::dimension);
  const auto & wedge_global2local = extra_mesh_info.global_to_local.at(wedge_pair);
  const auto & wedge_conn = extra_mesh_info.connectivity.at(wedge_pair);
  const auto & wedge_to_vertices = wedge_conn.at({0,0});
  const auto & wedge_to_edges    = wedge_conn.at({0,1});
  const auto & wedge_to_faces    = wedge_conn.at({0,2});
  const auto & wedge_to_cells    = wedge_conn.at({0,3});
  
  const auto & ghost_wedges = extra_mesh_info.ghost_ids.at(wedge_pair);
  const auto & ghost_wedge_conn = extra_mesh_info.ghost_connectivity.at(wedge_pair);
  const auto & ghost_wedge_to_vertices = ghost_wedge_conn.at({0,0});
  const auto & ghost_wedge_to_edges = ghost_wedge_conn.at({0,1});
  const auto & ghost_wedge_to_faces = ghost_wedge_conn.at({0,2});
  const auto & ghost_wedge_to_cells = ghost_wedge_conn.at({0,3});
  
 
  // create the wedges
  for(auto & em: wedge_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;
    
    // clear the lists
    entity_vs.clear();
    entity_es.clear();
    entity_fs.clear();
    entity_cs.clear();

    
    // search this ranks mesh definition for the matching offset
    auto it = wedge_global2local.find( mid );
    
    // the wedge exists on this rank so create it
    if ( it != wedge_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = wedge_to_vertices.at(id);
      const auto & es = wedge_to_edges.at(id);
      const auto & fs = wedge_to_faces.at(id);
      const auto & cs = wedge_to_cells.at(id);
      // create a list of vertex pointers
      entity_vs.resize( vs.size() );
      entity_es.resize( es.size() );
      entity_fs.resize( fs.size() );
      entity_cs.resize( cs.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        entity_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return &vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        entity_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return &edges[ flecsi_id ];
        }
      );
      std::transform(
        fs.begin(),
        fs.end(),
        entity_fs.begin(),
        [&](auto f) {
          auto global_id = face_local2global[f];
          auto flecsi_id = face_mid_to_lid.at(global_id);
          return &faces[ flecsi_id ];
        }
      );
      std::transform(
        cs.begin(),
        cs.end(),
        entity_cs.begin(),
        [&](auto c) {
          auto global_id = extra_mesh_info.cell_distribution[rank] + c;
          auto flecsi_id = cell_mid_to_lid.at(global_id);
          return &cells[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost wedge
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_wedges.begin(), ghost_wedges.end(), mid );
      assert( it != ghost_wedges.end() && "Could not find ghost wedge id" );
      auto i = std::distance( ghost_wedges.begin(), it );
      // reserve space
      const auto & vs = ghost_wedge_to_vertices.at(i);
      const auto & es = ghost_wedge_to_edges.at(i);
      const auto & fs = ghost_wedge_to_faces.at(i);
      const auto & cs = ghost_wedge_to_cells.at(i);
      entity_vs.reserve( vs.size() );
      entity_es.reserve( es.size() );
      entity_fs.reserve( fs.size() );
      entity_cs.reserve( cs.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        entity_vs.emplace_back( &vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        entity_es.emplace_back( &edges[id] );
      }
      for ( auto f : fs ) {
        auto id = face_mid_to_lid.at(f);
        entity_fs.emplace_back( &faces[id] );
      }
      for ( auto c : cs ) {
        auto id = cell_mid_to_lid.at(c);
        entity_cs.emplace_back( &cells[id] );
      }
    }

    // make the wedge
    auto new_wedge = mesh.template make<wedge_t, wedge_t::domain>();
    new_wedge->global_id().set_global(mid);
    wedges.emplace_back( new_wedge );
    // add the connectivity
    mesh.template init_entity<wedge_t::domain, vertex_t::domain, 
      wedge_t::dimension, vertex_t::dimension>(new_wedge, entity_vs);
    mesh.template init_entity<wedge_t::domain, edge_t::domain,
      wedge_t::dimension, edge_t::dimension>(new_wedge, entity_es);
    mesh.template init_entity<wedge_t::domain, face_t::domain,
      wedge_t::dimension, face_t::dimension>(new_wedge, entity_fs);
    mesh.template init_entity<wedge_t::domain, cell_t::domain,
      wedge_t::dimension, cell_t::dimension>(new_wedge, entity_cs);
  }

  //----------------------------------------------------------------------------
  // create the corners
  //----------------------------------------------------------------------------
  
  const auto & corner_lid_to_mid = context.index_map( index_spaces::corners );
  std::vector<corner_t*> corners; corners.reserve(corner_lid_to_mid.size());
  
  auto corner_pair = std::make_pair(corner_t::domain, corner_t::dimension);
  const auto & corner_global2local = extra_mesh_info.global_to_local.at(corner_pair);
  const auto & corner_conn = extra_mesh_info.connectivity.at(corner_pair);
  const auto & corner_to_vertices = corner_conn.at({0,0});
  const auto & corner_to_edges    = corner_conn.at({0,1});
  const auto & corner_to_faces    = corner_conn.at({0,2});
  const auto & corner_to_cells    = corner_conn.at({0,3});
  const auto & corner_to_wedges   = corner_conn.at(wedge_pair);
  
  const auto & ghost_corners = extra_mesh_info.ghost_ids.at(corner_pair);
  const auto & ghost_corner_conn = extra_mesh_info.ghost_connectivity.at(corner_pair);
  const auto & ghost_corner_to_vertices = ghost_corner_conn.at({0,0});
  const auto & ghost_corner_to_edges = ghost_corner_conn.at({0,1});
  const auto & ghost_corner_to_faces = ghost_corner_conn.at({0,2});
  const auto & ghost_corner_to_cells = ghost_corner_conn.at({0,3});
  const auto & ghost_corner_to_wedges = ghost_corner_conn.at(wedge_pair);
  
 
  // create the corners
  for(auto & em: corner_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;
    
    // clear the lists
    entity_vs.clear();
    entity_es.clear();
    entity_fs.clear();
    entity_cs.clear();
    entity_ws.clear();

    
    // search this ranks mesh definition for the matching offset
    auto it = corner_global2local.find( mid );
    
    // the corner exists on this rank so create it
    if ( it != corner_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = corner_to_vertices.at(id);
      const auto & es = corner_to_edges.at(id);
      const auto & fs = corner_to_faces.at(id);
      const auto & cs = corner_to_cells.at(id);
      const auto & ws = corner_to_wedges.at(id);
      // create a list of vertex pointers
      entity_vs.resize( vs.size() );
      entity_es.resize( es.size() );
      entity_fs.resize( fs.size() );
      entity_cs.resize( cs.size() );
      entity_ws.resize( ws.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        entity_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return &vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        entity_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return &edges[ flecsi_id ];
        }
      );
      std::transform(
        fs.begin(),
        fs.end(),
        entity_fs.begin(),
        [&](auto f) {
          auto global_id = face_local2global[f];
          auto flecsi_id = face_mid_to_lid.at(global_id);
          return &faces[ flecsi_id ];
        }
      );
      std::transform(
        cs.begin(),
        cs.end(),
        entity_cs.begin(),
        [&](auto c) {
          auto global_id = extra_mesh_info.cell_distribution[rank] + c;
          auto flecsi_id = cell_mid_to_lid.at(global_id);
          return &cells[ flecsi_id ];
        }
      );
      std::transform(
        ws.begin(),
        ws.end(),
        entity_ws.begin(),
        [&](auto w) {
          auto global_id = wedge_local2global[w];
          auto flecsi_id = wedge_mid_to_lid.at(global_id);
          return wedges[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost corner
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_corners.begin(), ghost_corners.end(), mid );
      assert( it != ghost_corners.end() && "Could not find ghost corner id" );
      auto i = std::distance( ghost_corners.begin(), it );
      // reserve space
      const auto & vs = ghost_corner_to_vertices.at(i);
      const auto & es = ghost_corner_to_edges.at(i);
      const auto & fs = ghost_corner_to_faces.at(i);
      const auto & cs = ghost_corner_to_cells.at(i);
      const auto & ws = ghost_corner_to_wedges.at(i);
      entity_vs.reserve( vs.size() );
      entity_es.reserve( es.size() );
      entity_fs.reserve( fs.size() );
      entity_cs.reserve( cs.size() );
      entity_ws.reserve( ws.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        entity_vs.emplace_back( &vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        entity_es.emplace_back( &edges[id] );
      }
      for ( auto f : fs ) {
        auto id = face_mid_to_lid.at(f);
        entity_fs.emplace_back( &faces[id] );
      }
      for ( auto c : cs ) {
        auto id = cell_mid_to_lid.at(c);
        entity_cs.emplace_back( &cells[id] );
      }
      for ( auto w : ws ) {
        auto id = wedge_mid_to_lid.at(w);
        entity_ws.emplace_back( wedges[id] );
      }
    }

    // make the corner
    auto new_corner = mesh.template make<corner_t, corner_t::domain>();
    new_corner->global_id().set_global(mid);
    corners.emplace_back( new_corner );
    // add the connectivity
    mesh.template init_entity<corner_t::domain, vertex_t::domain, 
      corner_t::dimension, vertex_t::dimension>(new_corner, entity_vs);
    mesh.template init_entity<corner_t::domain, edge_t::domain, 
      corner_t::dimension, edge_t::dimension>(new_corner, entity_es);
    mesh.template init_entity<corner_t::domain, face_t::domain, 
      corner_t::dimension, face_t::dimension>(new_corner, entity_fs);
    mesh.template init_entity<corner_t::domain, cell_t::domain, 
      corner_t::dimension, cell_t::dimension>(new_corner, entity_cs);
    mesh.template init_entity<corner_t::domain, wedge_t::domain, 
      corner_t::dimension, wedge_t::dimension>(new_corner, entity_ws);
  }


  //----------------------------------------------------------------------------
  // create the sides
  //----------------------------------------------------------------------------
  
  const auto & side_lid_to_mid = context.index_map( index_spaces::sides );
  std::vector<side_t*> sides; sides.reserve(side_lid_to_mid.size());
  
  auto side_pair = std::make_pair(side_t::domain, side_t::dimension);
  const auto & side_global2local = extra_mesh_info.global_to_local.at(side_pair);
  const auto & side_conn = extra_mesh_info.connectivity.at(side_pair);
  const auto & side_to_vertices = side_conn.at({0,0});
  const auto & side_to_edges    = side_conn.at({0,1});
  const auto & side_to_faces    = side_conn.at({0,2});
  const auto & side_to_cells    = side_conn.at({0,3});
  const auto & side_to_wedges   = side_conn.at(wedge_pair);
  const auto & side_to_corners  = side_conn.at(corner_pair);
  
  const auto & ghost_sides = extra_mesh_info.ghost_ids.at(side_pair);
  const auto & ghost_side_conn = extra_mesh_info.ghost_connectivity.at(side_pair);
  const auto & ghost_side_to_vertices = ghost_side_conn.at({0,0});
  const auto & ghost_side_to_edges = ghost_side_conn.at({0,1});
  const auto & ghost_side_to_faces = ghost_side_conn.at({0,2});
  const auto & ghost_side_to_cells = ghost_side_conn.at({0,3});
  const auto & ghost_side_to_wedges = ghost_side_conn.at(wedge_pair);
  const auto & ghost_side_to_corners = ghost_side_conn.at(corner_pair);
  
 
  // create the sides
  for(auto & em: side_lid_to_mid) {

    auto lid = em.first;
    auto mid = em.second;
    
    // clear the lists
    entity_vs.clear();
    entity_es.clear();
    entity_fs.clear();
    entity_cs.clear();
    entity_ws.clear();
    entity_ns.clear();

    
    // search this ranks mesh definition for the matching offset
    auto it = side_global2local.find( mid );
    
    // the side exists on this rank so create it
    if ( it != side_global2local.end() ) {
      // get the list of vertices
      auto id = it->second;
      const auto & vs = side_to_vertices.at(id);
      const auto & es = side_to_edges.at(id);
      const auto & fs = side_to_faces.at(id);
      const auto & cs = side_to_cells.at(id);
      const auto & ws = side_to_wedges.at(id);
      const auto & ns = side_to_corners.at(id);
      // create a list of vertex pointers
      entity_vs.resize( vs.size() );
      entity_es.resize( es.size() );
      entity_fs.resize( fs.size() );
      entity_cs.resize( cs.size() );
      entity_ws.resize( ws.size() );
      entity_ns.resize( ns.size() );
      // transform the list of vertices to mesh ids
      std::transform(
        vs.begin(),
        vs.end(),
        entity_vs.begin(),
        [&](auto v) {
          auto global_id = vert_local2global[v];
          auto flecsi_id = vertex_mid_to_lid.at(global_id);
          return &vertices[ flecsi_id ];
        }
      );
      std::transform(
        es.begin(),
        es.end(),
        entity_es.begin(),
        [&](auto e) {
          auto global_id = edge_local2global[e];
          auto flecsi_id = edge_mid_to_lid.at(global_id);
          return &edges[ flecsi_id ];
        }
      );
      std::transform(
        fs.begin(),
        fs.end(),
        entity_fs.begin(),
        [&](auto f) {
          auto global_id = face_local2global[f];
          auto flecsi_id = face_mid_to_lid.at(global_id);
          return &faces[ flecsi_id ];
        }
      );
      std::transform(
        cs.begin(),
        cs.end(),
        entity_cs.begin(),
        [&](auto c) {
          auto global_id = extra_mesh_info.cell_distribution[rank] + c;
          auto flecsi_id = cell_mid_to_lid.at(global_id);
          return &cells[ flecsi_id ];
        }
      );
      std::transform(
        ws.begin(),
        ws.end(),
        entity_ws.begin(),
        [&](auto w) {
          auto global_id = wedge_local2global[w];
          auto flecsi_id = wedge_mid_to_lid.at(global_id);
          return wedges[ flecsi_id ];
        }
      );
      std::transform(
        ns.begin(),
        ns.end(),
        entity_ns.begin(),
        [&](auto c) {
          auto global_id = corner_local2global[c];
          auto flecsi_id = corner_mid_to_lid.at(global_id);
          return corners[ flecsi_id ];
        }
      );
    }
    
    // otherwise it is a ghost side
    else {
      // find out what its connectivity info is
      auto it = std::find( ghost_sides.begin(), ghost_sides.end(), mid );
      assert( it != ghost_sides.end() && "Could not find ghost side id" );
      auto i = std::distance( ghost_sides.begin(), it );
      // reserve space
      const auto & vs = ghost_side_to_vertices.at(i);
      const auto & es = ghost_side_to_edges.at(i);
      const auto & fs = ghost_side_to_faces.at(i);
      const auto & cs = ghost_side_to_cells.at(i);
      const auto & ws = ghost_side_to_wedges.at(i);
      const auto & ns = ghost_side_to_corners.at(i);
      entity_vs.reserve( vs.size() );
      entity_es.reserve( es.size() );
      entity_fs.reserve( fs.size() );
      entity_cs.reserve( cs.size() );
      entity_ws.reserve( ws.size() );
      entity_ns.reserve( ns.size() );
      // now convert my global ids to flecsi local ids
      for ( auto v : vs ) {
        auto id = vertex_mid_to_lid.at(v);
        entity_vs.emplace_back( &vertices[id] );
      }
      for ( auto e : es ) {
        auto id = edge_mid_to_lid.at(e);
        entity_es.emplace_back( &edges[id] );
      }
      for ( auto f : fs ) {
        auto id = face_mid_to_lid.at(f);
        entity_fs.emplace_back( &faces[id] );
      }
      for ( auto c : cs ) {
        auto id = cell_mid_to_lid.at(c);
        entity_cs.emplace_back( &cells[id] );
      }
      for ( auto w : ws ) {
        auto id = wedge_mid_to_lid.at(w);
        entity_ws.emplace_back( wedges[id] );
      }
      for ( auto c : ns ) {
        auto id = corner_mid_to_lid.at(c);
        entity_ns.emplace_back( corners[id] );
      }
    }

    // make the side
    auto new_side = mesh.template make<side_t, side_t::domain>();
    new_side->global_id().set_global(mid);
    // add the connectivity
    mesh.template init_entity<side_t::domain, vertex_t::domain, 
      side_t::dimension, vertex_t::dimension>(new_side, entity_vs);
    mesh.template init_entity<side_t::domain, edge_t::domain, 
      side_t::dimension, edge_t::dimension>(new_side, entity_es);
    mesh.template init_entity<side_t::domain, face_t::domain, 
      side_t::dimension, face_t::dimension>(new_side, entity_fs);
    mesh.template init_entity<side_t::domain, cell_t::domain, 
      side_t::dimension, cell_t::dimension>(new_side, entity_cs);
    mesh.template init_entity<side_t::domain, wedge_t::domain, 
      side_t::dimension, wedge_t::dimension>(new_side, entity_ws);
    mesh.template init_entity<side_t::domain, corner_t::domain, 
      side_t::dimension, corner_t::dimension>(new_side, entity_ns);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the subspaces
///
/// \remarks This is the 1D version
////////////////////////////////////////////////////////////////////////////////
template<
  typename MESH_DEFINITION,
  typename MESH_TYPE,
  bool Enabled = ( std::decay_t<MESH_DEFINITION>::dimension() == 1 ),
  typename std::enable_if_t< Enabled >** = nullptr
>
void create_subspaces( MESH_DEFINITION && mesh_def, MESH_TYPE && mesh )
{

  using mesh_t = typename std::decay_t< MESH_TYPE >;

  // Alias the index spaces type
  using index_spaces = typename mesh_t::index_spaces_t;
  using index_subspaces = typename mesh_t::index_subspaces_t;

  // get the type of storage
  using vertex_t = typename mesh_t::vertex_t;

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  // get the colorings
  const auto & cell_coloring = context.coloring(index_spaces::cells);
  // get the entity maps
  // - lid = local id - the id of the entity local to this processor
  // - mid = mesh id - the original id of the entity ( usually from file )
  const auto & cell_mid_to_lid =
     context.reverse_index_map( index_spaces::cells );

  // storage for the unique set of vertices
  std::vector<vertex_t*> my_verts;

  // determine the unique list of ids
  const auto & cells = mesh.cells();

  auto add_ents = [&](auto && c) {
    auto cell_lid = cell_mid_to_lid.at( c.id );
    const auto & cell = cells[cell_lid];
    for ( auto v : mesh.vertices(cell) ) my_verts.push_back( v );
  };

  for(auto c : cell_coloring.exclusive) add_ents(c);
  for(auto c : cell_coloring.shared) add_ents(c);

  // sort and remove any duplicates
  std::sort( my_verts.begin(), my_verts.end() );
  ristra::utils::remove_duplicates( my_verts );

  // add them to the index space
  auto & vert_subspace =
     mesh.template get_index_subspace< index_subspaces::overlapping_vertices >();
  for ( auto v : my_verts )
     vert_subspace.push_back( v->template global_id<0>() );
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the subspaces
///
/// \remarks This is the 2D version
////////////////////////////////////////////////////////////////////////////////
template<
  typename MESH_DEFINITION,
  typename MESH_TYPE,
  bool Enabled = ( std::decay_t<MESH_DEFINITION>::dimension() == 2 ),
  typename = std::enable_if_t< Enabled >
>
void create_subspaces( MESH_DEFINITION && mesh_def, MESH_TYPE && mesh )
{

  using mesh_t = typename std::decay_t< MESH_TYPE >;

  // Alias the index spaces type
  using index_spaces = typename mesh_t::index_spaces_t;
  using index_subspaces = typename mesh_t::index_subspaces_t;

  // get the type of storage
  using vertex_t = typename mesh_t::vertex_t;
  using edge_t = typename mesh_t::edge_t;

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  // get the colorings
  const auto & cell_coloring = context.coloring(index_spaces::cells);
  // get the entity maps
  // - lid = local id - the id of the entity local to this processor
  // - mid = mesh id - the original id of the entity ( usually from file )
  const auto & cell_mid_to_lid =
    context.reverse_index_map( index_spaces::cells );

  // storage for the unique set of vertices and edges
  std::vector<vertex_t*> my_verts;
  std::vector<edge_t*> my_edges;

  // determine the unique list of ids
  const auto & cells = mesh.cells();

  auto add_ents = [&](auto && c) {
    auto cell_lid = cell_mid_to_lid.at( c.id );
    const auto & cell = cells[cell_lid];
    for ( auto v : mesh.vertices(cell) ) my_verts.push_back( v );
    for ( auto e : mesh.edges(cell) ) my_edges.push_back( e );
  };

  for(auto c : cell_coloring.exclusive) add_ents(c);
  for(auto c : cell_coloring.shared) add_ents(c);

  // sort and remove any duplicates
  std::sort( my_verts.begin(), my_verts.end() );
  std::sort( my_edges.begin(), my_edges.end() );
  ristra::utils::remove_duplicates( my_verts );
  ristra::utils::remove_duplicates( my_edges );

  // add them to the index space
  auto & vert_subspace =
    mesh.template get_index_subspace< index_subspaces::overlapping_vertices >();
  for ( auto v : my_verts )
    vert_subspace.push_back( v->global_id() );

  auto & edge_subspace =
    mesh.template get_index_subspace< index_subspaces::overlapping_edges >();
  for ( auto e : my_edges )
    edge_subspace.push_back( e->global_id() );
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper function to create the subspaces
///
/// \remarks This is the 3D version
////////////////////////////////////////////////////////////////////////////////
template<
  typename MESH_DEFINITION,
  typename MESH_TYPE,
  bool Enabled = ( std::decay_t<MESH_DEFINITION>::dimension() == 3 ),
  typename std::enable_if_t< Enabled >* = nullptr
>
void create_subspaces( MESH_DEFINITION && mesh_def, MESH_TYPE && mesh )
{

  using mesh_t = typename std::decay_t< MESH_TYPE >;

  // Alias the index spaces type
  using index_spaces = typename mesh_t::index_spaces_t;
  using index_subspaces = typename mesh_t::index_subspaces_t;

  // get the type of storage
  using vertex_t = typename mesh_t::vertex_t;
  using edge_t = typename mesh_t::edge_t;
  using face_t = typename mesh_t::face_t;

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  // get the colorings
  const auto & cell_coloring = context.coloring(index_spaces::cells);
  // get the entity maps
  // - lid = local id - the id of the entity local to this processor
  // - mid = mesh id - the original id of the entity ( usually from file )
  const auto & cell_mid_to_lid =
    context.reverse_index_map( index_spaces::cells );

  // storage for the unique set of vertices and edges
  std::vector<vertex_t*> my_verts;
  std::vector<edge_t*> my_edges;
  std::vector<face_t*> my_faces;

  // determine the unique list of ids
  const auto & cells = mesh.cells();

  auto add_ents = [&](auto && c) {
    auto cell_lid = cell_mid_to_lid.at( c.id );
    const auto & cell = cells[cell_lid];
    for ( auto v : mesh.vertices(cell) ) my_verts.push_back( v );
    for ( auto e : mesh.edges(cell) ) my_edges.push_back( e );
    for ( auto f : mesh.faces(cell) ) my_faces.push_back( f );
  };

  for(auto c : cell_coloring.exclusive) add_ents(c);
  for(auto c : cell_coloring.shared) add_ents(c);

  // sort and remove any duplicates
  std::sort( my_verts.begin(), my_verts.end() );
  std::sort( my_edges.begin(), my_edges.end() );
  std::sort( my_faces.begin(), my_faces.end() );
  ristra::utils::remove_duplicates( my_verts );
  ristra::utils::remove_duplicates( my_edges );
  ristra::utils::remove_duplicates( my_faces );

  // add them to the index space
  auto & vert_subspace =
    mesh.template get_index_subspace< index_subspaces::overlapping_vertices >();
  for ( auto v : my_verts )
    vert_subspace.push_back( v->global_id() );

  auto & edge_subspace =
    mesh.template get_index_subspace< index_subspaces::overlapping_edges >();
  for ( auto e : my_edges )
    edge_subspace.push_back( e->global_id() );

  auto & face_subspace =
    mesh.template get_index_subspace< index_subspaces::overlapping_faces >();
  for ( auto f : my_faces )
    face_subspace.push_back( f->global_id() );
}


////////////////////////////////////////////////////////////////////////////////
/// \brief Helper to make corners
////////////////////////////////////////////////////////////////////////////////
void make_corners(
  flecsi_sp::io::mesh_definition<2, burton_mesh_t::real_t> & mesh_def,
  extra_mesh_info_t & extra_mesh_info,
  std::map<size_t, flecsi::coloring::index_coloring_t> & coloring,
  std::map<size_t, flecsi::coloring::coloring_info_t> & coloring_info,
  std::map< dom_dim_t, std::map<dom_dim_t, size_t> > counts )
{
  
  // get communication info
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // get the number of dimensions
  constexpr auto num_dims = 2;//mesh_def.dimension();

  // get the connectvitiy from the cells to the vertices.  there is a
  // corner per cell, per vertex
  const auto & cells_to_vertices = mesh_def.entities_crs(num_dims, 0);
  const auto & cells_to_edges = mesh_def.entities_crs(num_dims, 1);
  const auto & edges_to_vertices = mesh_def.entities_crs(1, 0);
  const auto & vert_local2global = mesh_def.local_to_global(0);
  const auto & edge_local2global = mesh_def.local_to_global(1);
  auto num_cells = cells_to_vertices.size();
  auto num_verts = mesh_def.num_entities(0);

  //----------------------------------------------------------------------------
  // Ascertain local corner connectviity
  //----------------------------------------------------------------------------
  
  using vertex_t = burton_mesh_t::vertex_t;
  using edge_t = burton_mesh_t::edge_t;
  using face_t = burton_mesh_t::face_t;
  using cell_t = burton_mesh_t::cell_t;
  using corner_t = burton_mesh_t::corner_t;
  using wedge_t = burton_mesh_t::wedge_t;
  using side_t = burton_mesh_t::side_t;

  // create storage
  using connectivity_t = std::decay_t< decltype(cells_to_vertices) >;
  connectivity_t cells_to_corners;
  connectivity_t cells_to_wedges;
  connectivity_t cells_to_sides;

  // create entries for the connectivities we care about
  auto corner_pair = std::make_pair(corner_t::domain, corner_t::dimension);
  auto & corner_conn = extra_mesh_info.connectivity[corner_pair];
  auto & corners_to_verts = corner_conn[{0,0}];
  auto & corners_to_edges = corner_conn[{0,1}];
  auto & corners_to_cells = corner_conn[{0,2}];
  
  auto wedge_pair = std::make_pair(wedge_t::domain, wedge_t::dimension);
  auto & wedge_conn = extra_mesh_info.connectivity[wedge_pair];
  auto & wedges_to_verts = wedge_conn[{0,0}];
  auto & wedges_to_edges = wedge_conn[{0,1}];
  auto & wedges_to_cells = wedge_conn[{0,2}];

  auto side_pair = std::make_pair(side_t::domain, side_t::dimension);
  auto & side_conn = extra_mesh_info.connectivity[side_pair];
  auto & sides_to_verts = side_conn[{0,0}];
  auto & sides_to_edges = side_conn[{0,1}];
  auto & sides_to_cells = side_conn[{0,2}];

  auto & corners_to_wedges = corner_conn[wedge_pair];
  auto & sides_to_wedges = side_conn[wedge_pair];
  auto & sides_to_corners = side_conn[corner_pair];

  // some temporary storage
  std::vector< size_t > corner_ids, wedge_ids, side_ids;

  size_t corner_cnt{0}, wedge_cnt{0}, side_cnt{0};
  for ( size_t cell_id=0; cell_id<num_cells; ++cell_id )
  {

    //-------------------------------------------------------------------------
    // First assign a corner to each cell vertex

    // clear any temporary storage
    corner_ids.clear();
    wedge_ids.clear();
    side_ids.clear();

    // get the list of cell vertices
    const auto & verts = cells_to_vertices.at(cell_id);
    const auto & edges = cells_to_edges.at(cell_id);


    // reserve some storage
    corner_ids.reserve( verts.size() );
    wedge_ids.reserve( 2*edges.size() );
    side_ids.reserve( verts.size() );
    
    //-------------------------------------------------------------------------
    //! create a temporary lambda function for searching for edges
    //! \param[in] pa,pb  the point indexes of the edge to match
    //! \return the edge id for that particular point pair
    auto _find_edge = [&]( const auto & pa, const auto & pb  )
    {
      auto edge = std::find_if(
        edges.begin(), edges.end(),
        [&]( const auto & e )
        {
          auto verts = edges_to_vertices.at(e);
          assert( verts.size() == 2 && "should be two vertices per edge" );
          return ( (verts[0] == pa && verts[1] == pb) ||
                   (verts[0] == pb && verts[1] == pa) );
        }
      );
      assert( edge != edges.end() && "should have found an edge");
      return edge;
    };
    
    //-------------------------------------------------------------------------
    // Now loop over each vertex-pair forming an edge and assign wedges
    
    // find the first edge that goes from the last vertex to the first
    auto edge0 = _find_edge( verts.front(), verts.back() );

    // doesnt matter what we set it to, just keep track of the old wedge
    size_t old_wedge = 0;
    size_t old_corner = 0;

    auto side_start = side_cnt;

    // loop over the vertices in pairs
    for (
        auto v1=verts.begin(), v0 = std::prev(verts.end());
        v1!=verts.end();
        ++v1
    ) {

      // get the next vertex, if it goes off the end of the array, then
      // make it start over
      auto v2 = std::next(v1);
      if ( v2==verts.end() ) v2 = verts.begin();

      // then find the next edge
      auto edge1 = _find_edge( *v1, *v2 );

      // there are two wedges attached to this vertex.  set the
      // connectivity for both
      auto wedge0 = wedge_cnt;
      wedge_ids.emplace_back( wedge_cnt );
      wedges_to_cells.append( cell_id );
      wedges_to_edges.append( *edge0 );
      wedges_to_verts.append( *v1 );
      ++wedge_cnt;

      auto wedge1 = wedge_cnt;
      wedge_ids.emplace_back( wedge_cnt );
      wedges_to_cells.append( cell_id );
      wedges_to_edges.append( *edge1 );
      wedges_to_verts.append( *v1 );
      ++wedge_cnt;
      
      // keep track of this cells corners
      auto corner0 = corner_cnt;
      corner_ids.emplace_back( corner_cnt );
      // store all the connectivity
      corners_to_cells.append( cell_id );
      corners_to_edges.push_back( {*edge0, *edge1} );
      corners_to_verts.append( *v1 );
      // bump the corner counter
      ++corner_cnt;
      
      // side consists of edge1 and centroid.
      // Set the connectivity for both
      side_ids.emplace_back( side_cnt );
      sides_to_cells.append( cell_id );
      sides_to_edges.append( *edge0 );
      sides_to_verts.push_back( {*v0, *v1} );
      ++side_cnt;

      corners_to_wedges.push_back( {wedge0, wedge1} );
      sides_to_wedges.push_back( {old_wedge, wedge0} );
      sides_to_corners.push_back( {old_corner, corner0} );

      // update the old vertex and edge (saves having to find it again)
      v0 = v1;
      edge0 = edge1;
      old_wedge = wedge1;
      old_corner = corner0;
    } // verts

    // now fix the first one
    auto offset = sides_to_wedges.offsets[side_start];
    sides_to_wedges.indices[offset] = old_wedge;

    offset = sides_to_corners.offsets[side_start];
    sides_to_corners.indices[offset] = old_corner;



    // add corner ids to the cell
    cells_to_corners.append(corner_ids.begin(), corner_ids.end());
    cells_to_wedges.append(wedge_ids.begin(), wedge_ids.end());
    cells_to_sides.append(side_ids.begin(), side_ids.end());

  }

  //----------------------------------------------------------------------------
  // Figure out global corner numbering
  //----------------------------------------------------------------------------
  

  // the mpi data type for size_t
  const auto mpi_size_t = flecsi::utils::mpi_typetraits_u<size_t>::type();

  // exchange total corners
  std::vector<size_t> dist(size+1);
  auto ret = MPI_Allgather( &corner_cnt, 1, mpi_size_t, dist.data()+1,
    1, mpi_size_t, MPI_COMM_WORLD);
  dist[0] = 0;
  for ( size_t i=0; i<size; ++i ) dist[i+1] += dist[i];

  // create global numbering
  auto & corner_local2global = extra_mesh_info.local_to_global[corner_pair];
  auto & corner_global2local = extra_mesh_info.global_to_local[corner_pair];
  corner_local2global.resize(corner_cnt);
  for ( size_t i=0; i<corner_cnt; ++i ) {
    auto global_id = dist[rank] + i;
    corner_local2global[i] = global_id;
    corner_global2local[global_id] = i;
  }
  
  // and for the wedges
  ret = MPI_Allgather( &wedge_cnt, 1, mpi_size_t, dist.data()+1,
    1, mpi_size_t, MPI_COMM_WORLD);
  for ( size_t i=0; i<size; ++i ) dist[i+1] += dist[i];

  auto & wedge_local2global = extra_mesh_info.local_to_global[wedge_pair];
  auto & wedge_global2local = extra_mesh_info.global_to_local[wedge_pair];
  wedge_local2global.resize(wedge_cnt);
  for ( size_t i=0; i<wedge_cnt; ++i ) {
    auto global_id = dist[rank] + i;
    wedge_local2global[i] = global_id;
    wedge_global2local[global_id] = i;
  }

  // and for the sides
  ret = MPI_Allgather( &side_cnt, 1, mpi_size_t, dist.data()+1,
    1, mpi_size_t, MPI_COMM_WORLD);
  for ( size_t i=0; i<size; ++i ) dist[i+1] += dist[i];

  auto & side_local2global = extra_mesh_info.local_to_global[side_pair];
  auto & side_global2local = extra_mesh_info.global_to_local[side_pair];
  side_local2global.resize(side_cnt);
  for ( size_t i=0; i<side_cnt; ++i ) {
    auto global_id = dist[rank] + i;
    side_local2global[i] = global_id;
    side_global2local[global_id] = i;
  }

  //----------------------------------------------------------------------------
  // Apply coloring
  //----------------------------------------------------------------------------
  size_t dummy_counts;

  // Alias the index spaces type
  using index_spaces = burton_mesh_t::index_spaces_t;
  const auto & cells = coloring.at( index_spaces::cells );

  using corner_t = burton_mesh_t::corner_t;
  auto & corners = coloring[ index_spaces::corners ];
  auto & corner_info = coloring_info[ index_spaces::corners ];
  flecsi::coloring::color_entities( cells_to_corners, corner_local2global,
      corner_global2local, cells, corners, corner_info, dummy_counts );

  using wedge_t = burton_mesh_t::wedge_t;
  auto & wedges = coloring[ index_spaces::wedges ];
  auto & wedge_info = coloring_info[ index_spaces::wedges ];
  flecsi::coloring::color_entities( cells_to_wedges, wedge_local2global,
      wedge_global2local, cells, wedges, wedge_info, dummy_counts );

  using side_t = burton_mesh_t::side_t;
  auto & sides = coloring[ index_spaces::sides ];
  auto & side_info = coloring_info[ index_spaces::sides ];
  flecsi::coloring::color_entities( cells_to_sides, side_local2global,
      side_global2local, cells, sides, side_info, dummy_counts );

  //----------------------------------------------------------------------------
  // Finish counts
  //----------------------------------------------------------------------------

  // the size everything is based off
  const auto & ghost_cells_to_vertices = 
    extra_mesh_info.ghost_connectivity[{0,num_dims}][{0,0}];
  auto n = cells_to_vertices.indices.size() + ghost_cells_to_vertices.indices.size();
  
  // figure out everyones counts
  std::vector<size_t> color_sizes(size);
  ret = MPI_Allgather( &n, 1, mpi_size_t, color_sizes.data(), 
      1, mpi_size_t, MPI_COMM_WORLD);
  auto color_sizes_times_2 = color_sizes;
  for ( auto & i : color_sizes_times_2 ) i *= 2;

  // create a adjacency object
  auto & context = flecsi::execution::context_t::instance();
  flecsi::coloring::adjacency_info_t ai;

  // corners <-> verts
  ai.index_space = index_spaces::corners_to_vertices;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.color_sizes = color_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::vertices_to_corners;
  ai.from_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> edges
  ai.index_space = index_spaces::corners_to_edges;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::edges_to_corners;
  ai.from_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> cells
  ai.index_space = index_spaces::corners_to_cells;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = color_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::cells_to_corners;
  ai.from_index_space = index_spaces::entity_map[cell_t::domain][cell_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> verts
  ai.index_space = index_spaces::wedges_to_vertices;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::vertices_to_wedges;
  ai.from_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> edges
  ai.index_space = index_spaces::wedges_to_edges;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::edges_to_wedges;
  ai.from_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> cells
  ai.index_space = index_spaces::wedges_to_cells;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::cells_to_wedges;
  ai.from_index_space = index_spaces::entity_map[cell_t::domain][cell_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // sides <-> verts
  ai.index_space = index_spaces::sides_to_vertices;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::vertices_to_sides;
  ai.from_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);

  // sides <-> edges
  ai.index_space = index_spaces::sides_to_edges;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = color_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::edges_to_sides;
  ai.from_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);

  // sides <-> cells
  ai.index_space = index_spaces::sides_to_cells;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = color_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::cells_to_sides;
  ai.from_index_space = index_spaces::entity_map[cell_t::domain][cell_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);
  
  // corners <-> wedges
  ai.index_space = index_spaces::corners_to_wedges;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::wedges_to_corners;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> sides
  ai.index_space = index_spaces::corners_to_sides;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::sides_to_corners;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // sides <-> wedges
  ai.index_space = index_spaces::sides_to_wedges;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.color_sizes = color_sizes_times_2;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::wedges_to_sides;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);


  // these counts are reduntant but I think would be better suited to collapsing the
  // above code
  auto & corner_counts = counts[corner_pair];
  corner_counts[std::make_pair(0,0)] = n;
  corner_counts[std::make_pair(0,1)] = n * 2;
  corner_counts[std::make_pair(0,2)] = n;

  auto & wedge_counts = counts[wedge_pair];
  wedge_counts[std::make_pair(0,0)] = n * 2;
  wedge_counts[std::make_pair(0,1)] = n * 2;
  wedge_counts[std::make_pair(0,2)] = n * 2;

  auto & side_counts = counts[side_pair];
  side_counts[std::make_pair(0,0)] = n * 2;
  side_counts[std::make_pair(0,1)] = n;
  side_counts[std::make_pair(0,2)] = n;


  auto num_corners =
    corners.exclusive.size() + corners.shared.size() + corners.ghost.size();
  corner_counts[wedge_pair] = num_corners * 2;
  corner_counts[side_pair] = num_corners * 2;

  auto num_wedges =
    wedges.exclusive.size() + wedges.shared.size() + wedges.ghost.size();
  wedge_counts[corner_pair] = num_wedges;
  wedge_counts[side_pair] = num_wedges;

  auto num_sides =
    sides.exclusive.size() + sides.shared.size() + sides.ghost.size();
  side_counts[corner_pair] = num_sides * 2;
  side_counts[wedge_pair] = num_sides * 2;

  //----------------------------------------------------------------------------
  // Exchange ghost connectivity
  //----------------------------------------------------------------------------

  // Since corners/wedges/sides do not overlap, and are associated with cells,
  // we can just pack everything together and send it at once

  // first count
  std::vector<size_t> sendcounts(size, 0);

  for ( auto c : cells.shared ) {

    // get the entities for this cell
    const auto & cs = cells_to_corners.at(c.offset);
    const auto & ws = cells_to_wedges.at(c.offset);
    const auto & ss = cells_to_sides.at(c.offset);

    // loop over shared ranks
    for ( auto r : c.shared ) {
      // we will be sending global id, number of vertices, plus vertices
      if ( r != rank ) {
        // send cell id, 3 sizes, and 3 lists (corners, wedges, sides)
        // ... for each corner: id, vertex_id, two edge ids, two wedges
        // ... for each wedge: id, vertex_id, edge_id
        // ... for each side: id, two vertex ids, edge id, two wedges, two corers
        sendcounts[r] += 4 + cs.size()*6 + ws.size()*3 + ss.size()*8;
      }
    }

  }
  
  // finish displacements
  std::vector<size_t> senddispls(size+1);
  senddispls[0] = 0;
  for ( size_t r=0; r<size; ++r ) senddispls[r+1] = senddispls[r] + sendcounts[r];

  // now fill
  std::vector<size_t> sendbuf( senddispls[size], 0);
  std::fill( sendcounts.begin(), sendcounts.end(), 0 );

  for ( auto c : cells.shared ) {

    // get the entities for this cell
    const auto & cs = cells_to_corners.at(c.offset);
    const auto & ws = cells_to_wedges.at(c.offset);
    const auto & ss = cells_to_sides.at(c.offset);

    // loop over shared ranks
    for ( auto r : c.shared ) {
      // we will be sending global id, number of vertices, plus vertices
      if ( r != rank ) {
        // Get beginning offset
        auto j = senddispls[r] + sendcounts[r];
        // send cell id, 3 sizes, and 3 lists (corners, wedges, sides)
        sendbuf[j++] = c.id;
        // ... for each corner: id, vertex_id, two edge ids
        sendbuf[j++] = cs.size();
        for ( auto i : cs ) {
          sendbuf[j++] = corner_local2global[i];
          for ( auto v : corners_to_verts.at(i) )
            sendbuf[j++] = vert_local2global[v];
          for ( auto e : corners_to_edges.at(i) )
            sendbuf[j++] = edge_local2global[e];
          for ( auto w : corners_to_wedges.at(i) )
            sendbuf[j++] = wedge_local2global[w];
        }
        // ... for each wedge: id, vertex_id, edge_id
        sendbuf[j++] = ws.size();
        for ( auto i : ws ) {
          sendbuf[j++] = wedge_local2global[i];
          for ( auto v : wedges_to_verts.at(i) )
            sendbuf[j++] = vert_local2global[v];
          for ( auto e : wedges_to_edges.at(i) )
            sendbuf[j++] = edge_local2global[e];
        }
        // ... for each side: id, two vertex ids, edge id
        sendbuf[j++] = ss.size();
        for ( auto i : ss ) {
          sendbuf[j++] = side_local2global[i];
          for ( auto v : sides_to_verts.at(i) )
            sendbuf[j++] = vert_local2global[v];
          for ( auto e : sides_to_edges.at(i) )
            sendbuf[j++] = edge_local2global[e];
          for ( auto w : sides_to_wedges.at(i) )
            sendbuf[j++] = wedge_local2global[w];
          for ( auto c : sides_to_corners.at(i) )
            sendbuf[j++] = corner_local2global[c];
        }
        // increment send counter
        sendcounts[r] = j - senddispls[r];
      }
    }

  }
  
  // a quick check
  for ( size_t r=0; r<size; ++r )
    assert( senddispls[r+1] == senddispls[r] + sendcounts[r] );

  std::vector<size_t> recvcounts(size);
  ret = MPI_Alltoall(sendcounts.data(), 1, mpi_size_t, recvcounts.data(),
      1, mpi_size_t, MPI_COMM_WORLD);
  if ( ret != MPI_SUCCESS ) clog_error( "Error communicating vertex counts" );

  // how much info will we be receiving
  std::vector<size_t> recvdispls(size+1);
  recvdispls[0] = 0;
  for ( size_t r=0; r<size; ++r )
    recvdispls[r+1] = recvdispls[r] + recvcounts[r];

  // now send the actual info
  std::vector<size_t> recvbuf( recvdispls[size] );
  ret = flecsi::coloring::alltoallv(sendbuf, sendcounts, senddispls, recvbuf,
      recvcounts, recvdispls, MPI_COMM_WORLD );
  if ( ret != MPI_SUCCESS ) clog_error( "Error communicating vertices" );

  // now upack connectivity
  auto & ghost_cell_ids = extra_mesh_info.ghost_ids[{0,2}];

  auto & ghost_corner_ids = extra_mesh_info.ghost_ids[corner_pair];
  auto & ghost_corner_conn = extra_mesh_info.ghost_connectivity[corner_pair];
  auto & ghost_corners_to_verts = ghost_corner_conn[{0,0}];
  auto & ghost_corners_to_edges = ghost_corner_conn[{0,1}];
  auto & ghost_corners_to_cells = ghost_corner_conn[{0,2}];

  auto & ghost_wedge_ids = extra_mesh_info.ghost_ids[wedge_pair];
  auto & ghost_wedge_conn = extra_mesh_info.ghost_connectivity[wedge_pair];
  auto & ghost_wedges_to_verts = ghost_wedge_conn[{0,0}];
  auto & ghost_wedges_to_edges = ghost_wedge_conn[{0,1}];
  auto & ghost_wedges_to_cells = ghost_wedge_conn[{0,2}];

  auto & ghost_side_ids = extra_mesh_info.ghost_ids[side_pair];
  auto & ghost_side_conn = extra_mesh_info.ghost_connectivity[side_pair];
  auto & ghost_sides_to_verts = ghost_side_conn[{0,0}];
  auto & ghost_sides_to_edges = ghost_side_conn[{0,1}];
  auto & ghost_sides_to_cells = ghost_side_conn[{0,2}];
  auto & ghost_sides_to_wedges = ghost_side_conn[wedge_pair];
  auto & ghost_sides_to_corners = ghost_side_conn[corner_pair];

  auto & ghost_corners_to_wedges = ghost_corner_conn[wedge_pair];

  auto add_entities = [&]( auto & new_ids, auto & global_ids ) {
    global_ids.reserve( global_ids.size() + new_ids.size() );
    for ( auto & i : new_ids ) {
      auto global_id = i;
      i = global_ids.size();
      global_ids.emplace_back( global_id );
    }
  };

  for ( size_t i, r=0, ghost_cell_id=0; r<size; ++r ) {
    for ( i=recvdispls[r]; i<recvdispls[r+1]; ) {
      // clear storage
      corner_ids.clear();
      wedge_ids.clear();
      side_ids.clear();
      // global id
      auto global_cell_id = recvbuf[i]; i++;
      // corners
      auto n = recvbuf[i]; i++;
      corner_ids.reserve(n);
      for ( size_t j=0; j<n; ++j ) {
        corner_ids.emplace_back(recvbuf[i]); i++;
        auto v0 = recvbuf[i]; i++;
        auto e0 = recvbuf[i]; i++;
        auto e1 = recvbuf[i]; i++;
        auto w0 = recvbuf[i]; i++;
        auto w1 = recvbuf[i]; i++;
        ghost_corners_to_cells.append( global_cell_id );
        ghost_corners_to_edges.push_back( {e0, e1} );
        ghost_corners_to_verts.append( v0 );
        ghost_corners_to_wedges.push_back( {w0, w1} );
      }
      // wedges
      n = recvbuf[i]; i++;
      wedge_ids.reserve(n);
      for ( size_t j=0; j<n; ++j ) {
        wedge_ids.emplace_back(recvbuf[i]); i++;
        auto v0 = recvbuf[i]; i++;
        auto e0 = recvbuf[i]; i++;
        ghost_wedges_to_cells.append( global_cell_id );
        ghost_wedges_to_edges.append( e0 );
        ghost_wedges_to_verts.append( v0 );
      }
      // sides
      n = recvbuf[i]; i++;
      side_ids.reserve(n);
      for ( size_t j=0; j<n; ++j ) {
        side_ids.emplace_back(recvbuf[i]); i++;
        auto v0 = recvbuf[i]; i++;
        auto v1 = recvbuf[i]; i++;
        auto e0 = recvbuf[i]; i++;
        auto w0 = recvbuf[i]; i++;
        auto w1 = recvbuf[i]; i++;
        auto c0 = recvbuf[i]; i++;
        auto c1 = recvbuf[i]; i++;
        ghost_sides_to_cells.append( global_cell_id );
        ghost_sides_to_edges.append( e0 );
        ghost_sides_to_verts.push_back( {v0, v1} );
        ghost_sides_to_wedges.push_back( {w0, w1} );
        ghost_sides_to_corners.push_back( {c0, c1} );
      }
      // now create the new entities (ids will change from global to local)
      add_entities( corner_ids, ghost_corner_ids );
      add_entities( wedge_ids, ghost_wedge_ids );
      add_entities( side_ids, ghost_side_ids );
      // global cell id should be in sync!
      assert( global_cell_id == ghost_cell_ids[ghost_cell_id] );
      // bump ghost cell id
      ++ghost_cell_id;
    } // recvdispls
    assert( i==recvdispls[r+1] && "comm messed up" );
  } // ranks

}

////////////////////////////////////////////////////////////////////////////////
/// \brief Helper to make corners
////////////////////////////////////////////////////////////////////////////////
void make_corners(
  flecsi_sp::io::mesh_definition<3, burton_mesh_t::real_t> & mesh_def,
  extra_mesh_info_t & extra_mesh_info,
  std::map<size_t, flecsi::coloring::index_coloring_t> & coloring,
  std::map<size_t, flecsi::coloring::coloring_info_t> & coloring_info,
  std::map< dom_dim_t, std::map<dom_dim_t, size_t> > counts )
{
  
  // get communication info
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // get the number of dimensions
  constexpr auto num_dims = 3;//mesh_def.dimension();

  // get the connectvitiy from the cells to the vertices.  there is a
  // corner per cell, per vertex
  const auto & cells_to_vertices = mesh_def.entities_crs(num_dims, 0);
  const auto & cells_to_edges = mesh_def.entities_crs(num_dims, 1);
  const auto & cells_to_faces = mesh_def.entities_crs(num_dims, 2);
  const auto & faces_to_edges = mesh_def.entities_crs(2, 1);
  const auto & faces_to_vertices = mesh_def.entities_crs(2, 0);
  const auto & edges_to_vertices = mesh_def.entities_crs(1, 0);
  const auto & vert_local2global = mesh_def.local_to_global(0);
  const auto & edge_local2global = mesh_def.local_to_global(1);
  const auto & face_local2global = mesh_def.local_to_global(2);
  const auto & cell_local2global = mesh_def.local_to_global(3);
  auto num_cells = cells_to_vertices.size();
  auto num_verts = mesh_def.num_entities(0);

  const auto & face_owners = mesh_def.face_owners();

  //----------------------------------------------------------------------------
  // Ascertain local corner connectviity
  //----------------------------------------------------------------------------
  
  using vertex_t = burton_mesh_t::vertex_t;
  using edge_t = burton_mesh_t::edge_t;
  using face_t = burton_mesh_t::face_t;
  using cell_t = burton_mesh_t::cell_t;
  using corner_t = burton_mesh_t::corner_t;
  using wedge_t = burton_mesh_t::wedge_t;
  using side_t = burton_mesh_t::side_t;

  // create storage
  using connectivity_t = std::decay_t< decltype(cells_to_vertices) >;
  connectivity_t cells_to_corners;
  connectivity_t cells_to_wedges;
  connectivity_t cells_to_sides;

  // create entries for the connectivities we care about
  auto corner_pair = std::make_pair(corner_t::domain, corner_t::dimension);
  auto & corner_conn = extra_mesh_info.connectivity[corner_pair];
  auto & corners_to_verts = corner_conn[{0,0}];
  auto & corners_to_edges = corner_conn[{0,1}];
  auto & corners_to_faces = corner_conn[{0,2}];
  auto & corners_to_cells = corner_conn[{0,3}];
  
  auto wedge_pair = std::make_pair(wedge_t::domain, wedge_t::dimension);
  auto & wedge_conn = extra_mesh_info.connectivity[wedge_pair];
  auto & wedges_to_verts = wedge_conn[{0,0}];
  auto & wedges_to_edges = wedge_conn[{0,1}];
  auto & wedges_to_faces = wedge_conn[{0,2}];
  auto & wedges_to_cells = wedge_conn[{0,3}];

  auto side_pair = std::make_pair(side_t::domain, side_t::dimension);
  auto & side_conn = extra_mesh_info.connectivity[side_pair];
  auto & sides_to_verts = side_conn[{0,0}];
  auto & sides_to_edges = side_conn[{0,1}];
  auto & sides_to_faces = side_conn[{0,2}];
  auto & sides_to_cells = side_conn[{0,3}];

  auto & corners_to_wedges = corner_conn[wedge_pair];
  auto & sides_to_wedges = side_conn[wedge_pair];
  auto & sides_to_corners = side_conn[corner_pair];

  // some temporary storage
  std::vector< size_t > corner_ids, wedge_ids, side_ids;

  size_t corner_cnt{0}, wedge_cnt{0}, side_cnt{0};
  for ( size_t cell_id=0; cell_id<num_cells; ++cell_id )
  {

    //-------------------------------------------------------------------------
    // setup some storage for the cells
    
    // get the list of cell vertices
    const auto & cell_faces = cells_to_faces.at(cell_id);
    const auto & cell_edges = cells_to_edges.at(cell_id);
    const auto & cell_verts = cells_to_vertices.at(cell_id);


    // clear any temporary storage
    corner_ids.clear();
    wedge_ids.clear();
    side_ids.clear();

    // reserve some storage
    corner_ids.reserve( cell_verts.size() );
    wedge_ids.reserve( 4*cell_edges.size() );
    side_ids.reserve( 2*cell_edges.size() );

    //-------------------------------------------------------------------------
    // First loop over vertices and count the corners
    std::map<size_t, size_t> vertex_to_corner;
    std::vector<size_t> corner_to_vertex;

    for ( auto vertex_id : cell_verts ) {
      // keep track of this cells corners
      corner_ids.emplace_back( corner_cnt );
      // store all the connectivity
      corners_to_cells.append( cell_id );
      corners_to_verts.append( vertex_id );
      // temporariliy store which corner belongs to which vertex
      vertex_to_corner[vertex_id] = corner_cnt;
      // bump the corner counter
      ++corner_cnt;
    }
    
    //-------------------------------------------------------------------------
    // Now loop over faces
    std::map<size_t, std::vector<size_t>>
      corner_edges, corner_faces, corner_wedges;

    for ( auto face_id : cell_faces ) {

      // copy the vertices/faces of this face
      auto face_verts = faces_to_vertices.at(face_id).vec();
      const auto & face_edges = faces_to_edges.at(face_id);

      // check if this cell owns this face, if it doesnt, rotate the
      // vertices
      if ( face_owners[face_id] != cell_local2global[cell_id] )
        std::reverse( face_verts.begin(), face_verts.end() );


      //! create a temporary lambda function for searching for edges
      //! \param[in] pa,pb  the point indexes of the edge to match
      //! \return the edge id for that particular point pair
      auto _find_edge = [&]( const auto & pa, const auto & pb  )
      {
        auto edge = std::find_if(
          face_edges.begin(), face_edges.end(),
          [&]( const auto & e )
          {
            auto verts = edges_to_vertices.at(e);
            assert( verts.size() == 2 && "should be two vertices per edge" );
            return ( (verts[0] == pa && verts[1] == pb) ||
                    (verts[0] == pb && verts[1] == pa) );
          }
        );
        assert( edge != face_edges.end() && "should have found an edge");
        return edge;
      };

      //-------------------------------------------------------------------------
      // Now loop over each vertex-pair forming an edge and assign wedges
    
      // find the first edge that goes from the last vertex to the first
      auto edge0 = _find_edge( face_verts.front(), face_verts.back() );

      // doesnt matter what we set it to, just keep track of the old wedge
      size_t old_wedge = 0;

      auto side_start = side_cnt;

      // loop over the vertices in pairs
      for (
          auto v1=face_verts.begin(), v0 = std::prev(face_verts.end());
          v1!=face_verts.end();
          ++v1
      ) {

        // get the next vertex, if it goes off the end of the array, then
        // make it start over
        auto v2 = std::next(v1);
        if ( v2==face_verts.end() ) v2 = face_verts.begin();

        // then find the next edge
        auto edge1 = _find_edge( *v1, *v2 );

        // there are two wedges attached to this vertex.  set the
        // connectivity for both
        auto wedge0 = wedge_cnt;
        wedge_ids.emplace_back( wedge_cnt );
        wedges_to_cells.append( cell_id );
        wedges_to_faces.append( face_id );
        wedges_to_edges.append( *edge0 );
        wedges_to_verts.append( *v1 );
        ++wedge_cnt;

        auto wedge1 = wedge_cnt;
        wedge_ids.emplace_back( wedge_cnt );
        wedges_to_cells.append( cell_id );
        wedges_to_faces.append( face_id );
        wedges_to_edges.append( *edge1 );
        wedges_to_verts.append( *v1 );
        ++wedge_cnt;
        
        // There are two conres attached to this vertex.  set the
        // connectivity for both
        auto corner0 = vertex_to_corner.at(*v0);
        auto corner1 = vertex_to_corner.at(*v1);
        
        corner_faces[corner1].emplace_back( face_id );
        auto & ce = corner_edges[corner1];
        ce.emplace_back( *edge0 );
        ce.emplace_back( *edge1 );
        
        // side consists of edge1 and centroid.
        // Set the connectivity for both
        side_ids.emplace_back( side_cnt );
        sides_to_cells.append( cell_id );
        sides_to_faces.append( face_id );
        sides_to_edges.append( *edge0 );
        sides_to_verts.push_back( {*v0, *v1} );
        ++side_cnt;

        auto & cw = corner_wedges[corner1];
        cw.emplace_back( wedge0 );
        cw.emplace_back( wedge1 );
        sides_to_wedges.push_back( {old_wedge, wedge0} );
        sides_to_corners.push_back( {corner0, corner1} );

        // update the old vertex and edge (saves having to find it again)
        v0 = v1;
        edge0 = edge1;
        old_wedge = wedge1;
      } // verts

    // now fix the first one
    auto offset = sides_to_wedges.offsets[side_start];
    sides_to_wedges.indices[offset] = old_wedge;


    } // faces

    // finish the corners
    for ( size_t i=0; i<corner_ids.size(); ++i ) {
      auto cn = corner_ids[i];
      const auto & cw = corner_wedges.at(cn);
      corners_to_wedges.append( cw.begin(), cw.end() );
      const auto & cf = corner_faces.at(cn);
      corners_to_faces.append( cf.begin(), cf.end() );
      auto & ce = corner_edges.at(cn);
      std::sort( ce.begin(), ce.end() );
      auto last = std::unique( ce.begin(), ce.end() );
      corners_to_edges.append( ce.begin(), last );
    }


    // add corner ids to the cell
    cells_to_corners.append(corner_ids.begin(), corner_ids.end());
    cells_to_wedges.append(wedge_ids.begin(), wedge_ids.end());
    cells_to_sides.append(side_ids.begin(), side_ids.end());
  }

  //----------------------------------------------------------------------------
  // Figure out global corner numbering
  //----------------------------------------------------------------------------
  

  // the mpi data type for size_t
  const auto mpi_size_t = flecsi::utils::mpi_typetraits_u<size_t>::type();

  // exchange total corners
  std::vector<size_t> dist(size+1);
  auto ret = MPI_Allgather( &corner_cnt, 1, mpi_size_t, dist.data()+1,
    1, mpi_size_t, MPI_COMM_WORLD);
  dist[0] = 0;
  for ( size_t i=0; i<size; ++i ) dist[i+1] += dist[i];

  // create global numbering
  auto & corner_local2global = extra_mesh_info.local_to_global[corner_pair];
  auto & corner_global2local = extra_mesh_info.global_to_local[corner_pair];
  corner_local2global.resize(corner_cnt);
  for ( size_t i=0; i<corner_cnt; ++i ) {
    auto global_id = dist[rank] + i;
    corner_local2global[i] = global_id;
    corner_global2local[global_id] = i;
  }
  
  // and for the wedges
  ret = MPI_Allgather( &wedge_cnt, 1, mpi_size_t, dist.data()+1,
    1, mpi_size_t, MPI_COMM_WORLD);
  for ( size_t i=0; i<size; ++i ) dist[i+1] += dist[i];

  auto & wedge_local2global = extra_mesh_info.local_to_global[wedge_pair];
  auto & wedge_global2local = extra_mesh_info.global_to_local[wedge_pair];
  wedge_local2global.resize(wedge_cnt);
  for ( size_t i=0; i<wedge_cnt; ++i ) {
    auto global_id = dist[rank] + i;
    wedge_local2global[i] = global_id;
    wedge_global2local[global_id] = i;
  }

  // and for the sides
  ret = MPI_Allgather( &side_cnt, 1, mpi_size_t, dist.data()+1,
    1, mpi_size_t, MPI_COMM_WORLD);
  for ( size_t i=0; i<size; ++i ) dist[i+1] += dist[i];

  auto & side_local2global = extra_mesh_info.local_to_global[side_pair];
  auto & side_global2local = extra_mesh_info.global_to_local[side_pair];
  side_local2global.resize(side_cnt);
  for ( size_t i=0; i<side_cnt; ++i ) {
    auto global_id = dist[rank] + i;
    side_local2global[i] = global_id;
    side_global2local[global_id] = i;
  }

  //----------------------------------------------------------------------------
  // Apply coloring
  //----------------------------------------------------------------------------
  size_t dummy_counts;
  
  using index_spaces = burton_mesh_t::index_spaces_t;

  const auto & cells = coloring.at( index_spaces::cells );

  using corner_t = burton_mesh_t::corner_t;
  auto & corners = coloring[ index_spaces::corners ];
  auto & corner_info = coloring_info[ index_spaces::corners ];
  flecsi::coloring::color_entities( cells_to_corners, corner_local2global,
      corner_global2local, cells, corners, corner_info, dummy_counts );

  using wedge_t = burton_mesh_t::wedge_t;
  auto & wedges = coloring[ index_spaces::wedges ];
  auto & wedge_info = coloring_info[ index_spaces::wedges ];
  flecsi::coloring::color_entities( cells_to_wedges, wedge_local2global,
      wedge_global2local, cells, wedges, wedge_info, dummy_counts );

  using side_t = burton_mesh_t::side_t;
  auto & sides = coloring[ index_spaces::sides ];
  auto & side_info = coloring_info[ index_spaces::sides ];
  flecsi::coloring::color_entities( cells_to_sides, side_local2global,
      side_global2local, cells, sides, side_info, dummy_counts );

  //----------------------------------------------------------------------------
  // Finish counts
  //----------------------------------------------------------------------------

  // the size everything is based off
  auto num_corners = corners.exclusive.size() + corners.shared.size() + corners.ghost.size();
  auto num_sides = sides.exclusive.size() + sides.shared.size() + sides.ghost.size();
  auto num_wedges = wedges.exclusive.size() + wedges.shared.size() + wedges.ghost.size();
  
  // figure out everyones counts
  std::vector<size_t> corner_sizes(size), side_sizes(size), wedge_sizes(size);
  ret = MPI_Allgather( &num_corners, 1, mpi_size_t, corner_sizes.data(), 
      1, mpi_size_t, MPI_COMM_WORLD);
  ret = MPI_Allgather( &num_sides, 1, mpi_size_t, side_sizes.data(), 
      1, mpi_size_t, MPI_COMM_WORLD);
  ret = MPI_Allgather( &num_wedges, 1, mpi_size_t, wedge_sizes.data(), 
      1, mpi_size_t, MPI_COMM_WORLD);

  // create a adjacency object
  auto & context = flecsi::execution::context_t::instance();
  flecsi::coloring::adjacency_info_t ai;

  // corners <-> verts
  ai.index_space = index_spaces::corners_to_vertices;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.color_sizes = corner_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::vertices_to_corners;
  ai.from_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> edges
  ai.index_space = index_spaces::corners_to_edges;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = side_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::edges_to_corners;
  ai.from_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> faces
  ai.index_space = index_spaces::corners_to_faces;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[face_t::domain][face_t::dimension];
  ai.color_sizes = side_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::faces_to_corners;
  ai.from_index_space = index_spaces::entity_map[face_t::domain][face_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> cells
  ai.index_space = index_spaces::corners_to_cells;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = corner_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::cells_to_corners;
  ai.from_index_space = index_spaces::entity_map[cell_t::domain][cell_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> verts
  ai.index_space = index_spaces::wedges_to_vertices;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::vertices_to_wedges;
  ai.from_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> edges
  ai.index_space = index_spaces::wedges_to_edges;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::edges_to_wedges;
  ai.from_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> faces
  ai.index_space = index_spaces::wedges_to_faces;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[face_t::domain][face_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::faces_to_wedges;
  ai.from_index_space = index_spaces::entity_map[face_t::domain][face_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // wedges <-> cells
  ai.index_space = index_spaces::wedges_to_cells;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::cells_to_wedges;
  ai.from_index_space = index_spaces::entity_map[cell_t::domain][cell_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  context.add_adjacency(ai);

  // sides <-> verts
  ai.index_space = index_spaces::sides_to_vertices;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::vertices_to_sides;
  ai.from_index_space = index_spaces::entity_map[vertex_t::domain][vertex_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);

  // sides <-> edges
  ai.index_space = index_spaces::sides_to_edges;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = side_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::edges_to_sides;
  ai.from_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);

  // sides <-> faces
  ai.index_space = index_spaces::sides_to_faces;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[face_t::domain][face_t::dimension];
  ai.color_sizes = side_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::faces_to_sides;
  ai.from_index_space = index_spaces::entity_map[face_t::domain][face_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);

  // sides <-> cells
  ai.index_space = index_spaces::sides_to_cells;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[edge_t::domain][edge_t::dimension];
  ai.color_sizes = side_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::cells_to_sides;
  ai.from_index_space = index_spaces::entity_map[cell_t::domain][cell_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);
  
  // corners <-> wedges
  ai.index_space = index_spaces::corners_to_wedges;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::wedges_to_corners;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // corners <-> sides
  ai.index_space = index_spaces::corners_to_sides;
  ai.from_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::sides_to_corners;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[corner_t::domain][corner_t::dimension];
  context.add_adjacency(ai);

  // sides <-> wedges
  ai.index_space = index_spaces::sides_to_wedges;
  ai.from_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  ai.to_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.color_sizes = wedge_sizes;
  context.add_adjacency(ai);

  ai.index_space = index_spaces::wedges_to_sides;
  ai.from_index_space = index_spaces::entity_map[wedge_t::domain][wedge_t::dimension];
  ai.to_index_space = index_spaces::entity_map[side_t::domain][side_t::dimension];
  context.add_adjacency(ai);


  // these counts are reduntant but I think would be better suited to collapsing the
  // above code
  auto & corner_counts = counts[corner_pair];
  corner_counts[std::make_pair(0,0)] = num_corners;
  corner_counts[std::make_pair(0,1)] = num_sides;
  corner_counts[std::make_pair(0,2)] = num_sides;
  corner_counts[std::make_pair(0,3)] = num_corners;

  auto & wedge_counts = counts[wedge_pair];
  wedge_counts[std::make_pair(0,0)] = num_wedges;
  wedge_counts[std::make_pair(0,1)] = num_wedges;
  wedge_counts[std::make_pair(0,2)] = num_wedges;
  wedge_counts[std::make_pair(0,3)] = num_wedges;

  auto & side_counts = counts[side_pair];
  side_counts[std::make_pair(0,0)] = num_wedges;
  side_counts[std::make_pair(0,1)] = num_sides;
  side_counts[std::make_pair(0,2)] = num_sides;
  side_counts[std::make_pair(0,3)] = num_sides;


  corner_counts[wedge_pair] = num_wedges;
  corner_counts[side_pair] = num_wedges;

  wedge_counts[corner_pair] = num_wedges;
  wedge_counts[side_pair] = num_wedges;

  side_counts[corner_pair] = num_wedges;
  side_counts[wedge_pair] = num_wedges;

  //----------------------------------------------------------------------------
  // Exchange ghost connectivity
  //----------------------------------------------------------------------------

  // Since corners/wedges/sides do not overlap, and are associated with cells,
  // we can just pack everything together and send it at once

  // first count
  std::vector<size_t> sendcounts(size, 0);

  for ( auto c : cells.shared ) {

    // get the entities for this cell
    const auto & cs = cells_to_corners.at(c.offset);
    const auto & ws = cells_to_wedges.at(c.offset);
    const auto & ss = cells_to_sides.at(c.offset);


    size_t edge_cnt{0}, face_cnt{0}, wedge_cnt{0};
    for ( auto i : cs ) {
      const auto & ce = corners_to_edges.at(i);
      const auto & cf = corners_to_faces.at(i);
      const auto & cw = corners_to_wedges.at(i);
      edge_cnt += 1 + ce.size();
      face_cnt += 1 + cf.size();
      wedge_cnt += 1 + cw.size();
    }
    auto corner_cnt = edge_cnt + face_cnt + wedge_cnt;
    

    // loop over shared ranks
    for ( auto r : c.shared ) {
      // we will be sending global id, number of vertices, plus vertices
      if ( r != rank ) {
        // send cell id, 3 sizes, and 3 lists (corners, wedges, sides)
        // ... for each corner: id, vertex_id, n1, n1 edge ids, n2, n2 faces, 
        //     n3, n3 wedges
        // ... for each wedge: id, vertex_id, edge_id, face_id
        // ... for each side: id, two vertex ids, edge id, face id, two wedges,
        //     two corers
        sendcounts[r] += 4 + (cs.size()*2+corner_cnt) + ws.size()*4
          + ss.size()*9;
      }
    }

  }
  
  // finish displacements
  std::vector<size_t> senddispls(size+1);
  senddispls[0] = 0;
  for ( size_t r=0; r<size; ++r ) senddispls[r+1] = senddispls[r] + sendcounts[r];

  // now fill
  std::vector<size_t> sendbuf( senddispls[size], 0);
  std::fill( sendcounts.begin(), sendcounts.end(), 0 );

  for ( auto c : cells.shared ) {

    // get the entities for this cell
    const auto & cs = cells_to_corners.at(c.offset);
    const auto & ws = cells_to_wedges.at(c.offset);
    const auto & ss = cells_to_sides.at(c.offset);

    // loop over shared ranks
    for ( auto r : c.shared ) {
      // we will be sending global id, number of vertices, plus vertices
      if ( r != rank ) {
        // Get beginning offset
        auto j = senddispls[r] + sendcounts[r];
        // send cell id, 3 sizes, and 3 lists (corners, wedges, sides)
        sendbuf[j++] = c.id;
        // ... for each corner: id, vertex_id, two edge ids
        sendbuf[j++] = cs.size();
        for ( auto i : cs ) {
          sendbuf[j++] = corner_local2global[i];
          for ( auto v : corners_to_verts.at(i) )
            sendbuf[j++] = vert_local2global[v];
          const auto & ce = corners_to_edges.at(i);
          sendbuf[j++] = ce.size();
          for ( auto e : ce )
            sendbuf[j++] = edge_local2global[e];
          const auto & cf = corners_to_faces.at(i);
          sendbuf[j++] = cf.size();
          for ( auto f : cf )
            sendbuf[j++] = face_local2global[f];
          const auto & cw = corners_to_wedges.at(i);
          sendbuf[j++] = cw.size();
          for ( auto w : cw )
            sendbuf[j++] = wedge_local2global[w];
        }
        // ... for each wedge: id, vertex_id, edge_id, face_id
        sendbuf[j++] = ws.size();
        for ( auto i : ws ) {
          sendbuf[j++] = wedge_local2global[i];
          for ( auto v : wedges_to_verts.at(i) )
            sendbuf[j++] = vert_local2global[v];
          for ( auto e : wedges_to_edges.at(i) )
            sendbuf[j++] = edge_local2global[e];
          for ( auto f : wedges_to_faces.at(i) )
            sendbuf[j++] = face_local2global[f];
        }
        // ... for each side: id, two vertex ids, edge id
        sendbuf[j++] = ss.size();
        for ( auto i : ss ) {
          sendbuf[j++] = side_local2global[i];
          for ( auto v : sides_to_verts.at(i) )
            sendbuf[j++] = vert_local2global[v];
          for ( auto e : sides_to_edges.at(i) )
            sendbuf[j++] = edge_local2global[e];
          for ( auto f : sides_to_faces.at(i) )
            sendbuf[j++] = face_local2global[f];
          for ( auto w : sides_to_wedges.at(i) )
            sendbuf[j++] = wedge_local2global[w];
          for ( auto c : sides_to_corners.at(i) )
            sendbuf[j++] = corner_local2global[c];
        }
        // increment send counter
        sendcounts[r] = j - senddispls[r];
      }
    }

  }
  
  // a quick check
  for ( size_t r=0; r<size; ++r )
    assert( senddispls[r+1] == senddispls[r] + sendcounts[r] );

  std::vector<size_t> recvcounts(size);
  ret = MPI_Alltoall(sendcounts.data(), 1, mpi_size_t, recvcounts.data(),
      1, mpi_size_t, MPI_COMM_WORLD);
  if ( ret != MPI_SUCCESS ) clog_error( "Error communicating vertex counts" );

  // how much info will we be receiving
  std::vector<size_t> recvdispls(size+1);
  recvdispls[0] = 0;
  for ( size_t r=0; r<size; ++r )
    recvdispls[r+1] = recvdispls[r] + recvcounts[r];

  // now send the actual info
  std::vector<size_t> recvbuf( recvdispls[size] );
  ret = flecsi::coloring::alltoallv(sendbuf, sendcounts, senddispls, recvbuf,
      recvcounts, recvdispls, MPI_COMM_WORLD );
  if ( ret != MPI_SUCCESS ) clog_error( "Error communicating vertices" );

  // now upack connectivity
  auto & ghost_cell_ids = extra_mesh_info.ghost_ids[{0,3}];

  auto & ghost_corner_ids = extra_mesh_info.ghost_ids[corner_pair];
  auto & ghost_corner_conn = extra_mesh_info.ghost_connectivity[corner_pair];
  auto & ghost_corners_to_verts = ghost_corner_conn[{0,0}];
  auto & ghost_corners_to_edges = ghost_corner_conn[{0,1}];
  auto & ghost_corners_to_faces = ghost_corner_conn[{0,2}];
  auto & ghost_corners_to_cells = ghost_corner_conn[{0,3}];

  auto & ghost_wedge_ids = extra_mesh_info.ghost_ids[wedge_pair];
  auto & ghost_wedge_conn = extra_mesh_info.ghost_connectivity[wedge_pair];
  auto & ghost_wedges_to_verts = ghost_wedge_conn[{0,0}];
  auto & ghost_wedges_to_edges = ghost_wedge_conn[{0,1}];
  auto & ghost_wedges_to_faces = ghost_wedge_conn[{0,2}];
  auto & ghost_wedges_to_cells = ghost_wedge_conn[{0,3}];

  auto & ghost_side_ids = extra_mesh_info.ghost_ids[side_pair];
  auto & ghost_side_conn = extra_mesh_info.ghost_connectivity[side_pair];
  auto & ghost_sides_to_verts = ghost_side_conn[{0,0}];
  auto & ghost_sides_to_edges = ghost_side_conn[{0,1}];
  auto & ghost_sides_to_faces = ghost_side_conn[{0,2}];
  auto & ghost_sides_to_cells = ghost_side_conn[{0,3}];
  auto & ghost_sides_to_wedges = ghost_side_conn[wedge_pair];
  auto & ghost_sides_to_corners = ghost_side_conn[corner_pair];

  auto & ghost_corners_to_wedges = ghost_corner_conn[wedge_pair];

  auto add_entities = [&]( auto & new_ids, auto & global_ids ) {
    global_ids.reserve( global_ids.size() + new_ids.size() );
    for ( auto & i : new_ids ) {
      auto global_id = i;
      i = global_ids.size();
      global_ids.emplace_back( global_id );
    }
  };

  std::vector<size_t> tmp;

  for ( size_t i, r=0, ghost_cell_id=0; r<size; ++r ) {
    for ( i=recvdispls[r]; i<recvdispls[r+1]; ) {
      // clear storage
      corner_ids.clear();
      wedge_ids.clear();
      side_ids.clear();
      // global id
      auto global_cell_id = recvbuf[i]; i++;
      // corners
      auto n = recvbuf[i]; i++;
      corner_ids.reserve(n);
      for ( size_t j=0; j<n; ++j ) {
        corner_ids.emplace_back(recvbuf[i]); i++;
        auto v0 = recvbuf[i]; i++;
        auto nn = recvbuf[i]; i++;
        tmp.clear(); tmp.reserve(nn);
        for ( size_t k=0; k<nn; ++k ) { tmp.emplace_back(recvbuf[i]); i++; }
        ghost_corners_to_edges.append( tmp.begin(), tmp.end() );
        nn = recvbuf[i]; i++;
        tmp.clear(); tmp.reserve(nn);
        for ( size_t k=0; k<nn; ++k ) { tmp.emplace_back(recvbuf[i]); i++; }
        ghost_corners_to_faces.append( tmp.begin(), tmp.end() );
        nn = recvbuf[i]; i++;
        tmp.clear(); tmp.reserve(nn);
        for ( size_t k=0; k<nn; ++k ) { tmp.emplace_back(recvbuf[i]); i++; }
        ghost_corners_to_wedges.append( tmp.begin(), tmp.end() );
        ghost_corners_to_cells.append( global_cell_id );
        ghost_corners_to_verts.append( v0 );
      }
      // wedges
      n = recvbuf[i]; i++;
      wedge_ids.reserve(n);
      for ( size_t j=0; j<n; ++j ) {
        wedge_ids.emplace_back(recvbuf[i]); i++;
        auto v0 = recvbuf[i]; i++;
        auto e0 = recvbuf[i]; i++;
        auto f0 = recvbuf[i]; i++;
        ghost_wedges_to_cells.append( global_cell_id );
        ghost_wedges_to_faces.append( f0 );
        ghost_wedges_to_edges.append( e0 );
        ghost_wedges_to_verts.append( v0 );
      }
      // sides
      n = recvbuf[i]; i++;
      side_ids.reserve(n);
      for ( size_t j=0; j<n; ++j ) {
        side_ids.emplace_back(recvbuf[i]); i++;
        auto v0 = recvbuf[i]; i++;
        auto v1 = recvbuf[i]; i++;
        auto e0 = recvbuf[i]; i++;
        auto f0 = recvbuf[i]; i++;
        auto w0 = recvbuf[i]; i++;
        auto w1 = recvbuf[i]; i++;
        auto c0 = recvbuf[i]; i++;
        auto c1 = recvbuf[i]; i++;
        ghost_sides_to_cells.append( global_cell_id );
        ghost_sides_to_faces.append( f0 );
        ghost_sides_to_edges.append( e0 );
        ghost_sides_to_verts.push_back( {v0, v1} );
        ghost_sides_to_wedges.push_back( {w0, w1} );
        ghost_sides_to_corners.push_back( {c0, c1} );
      }
      // now create the new entities (ids will change from global to local)
      add_entities( corner_ids, ghost_corner_ids );
      add_entities( wedge_ids, ghost_wedge_ids );
      add_entities( side_ids, ghost_side_ids );
      // global cell id should be in sync!
      assert( global_cell_id == ghost_cell_ids[ghost_cell_id] );
      // bump ghost cell id
      ++ghost_cell_id;
    } // recvdispls
    assert( i==recvdispls[r+1] && "comm messed up" );
  } // ranks

}

////////////////////////////////////////////////////////////////////////////////
/// \brief the main cell coloring driver
////////////////////////////////////////////////////////////////////////////////
void partition_mesh( utils::char_array_t filename, std::size_t max_entries )
{
  // set some compile time constants
  constexpr auto num_dims = burton_mesh_t::num_dimensions;
  constexpr auto num_domains = burton_mesh_t::num_domains;
  constexpr auto cell_dim = num_dims;
  constexpr auto thru_dim = 0;

  // make some type aliases
  using real_t = burton_mesh_t::real_t;
  using size_t = burton_mesh_t::size_t;
  using exodus_definition_t = flecsi_sp::io::exodus_definition<num_dims, real_t>;
  using entity_info_t = flecsi::coloring::entity_info_t;
  using vertex_t = burton_mesh_t::vertex_t;
  using edge_t = burton_mesh_t::edge_t;
  using face_t = burton_mesh_t::face_t;
  using cell_t = burton_mesh_t::cell_t;
  using corner_t = burton_mesh_t::corner_t;
  using wedge_t = burton_mesh_t::wedge_t;
  using side_t = burton_mesh_t::side_t;

  // Alias the index spaces type
  using index_spaces = burton_mesh_t::index_spaces_t;

  // Get the context instance.
  auto & context = flecsi::execution::context_t::instance();

  // Create a communicator instance to get neighbor information.
  auto communicator = std::make_unique<flecsi::coloring::mpi_communicator_t>();
  auto comm_size = communicator->size();
  auto rank = communicator->rank();

  enum class file_type_t {
    exodus,
    partitioned_exodus,
    unknown
  };

  // try to figure out what kind of file it is
  auto get_file_type = []( auto str ) {
    auto base = ristra::utils::basename(str);
    auto i = base.rfind( '.', base.length() );
    size_t cnt = 0;
    while ( i != std::string::npos ) {
      auto ext = base.substr(i+1, base.length()-1);
      if ( ext == "g" || ext == "exo" ) {
        if ( cnt != 0 ) return file_type_t::partitioned_exodus;
        else            return file_type_t::exodus;
      }
      base = base.substr(0, i);
      i = base.rfind( '.', base.length() );
      ++cnt;
    }
    return file_type_t::unknown;
  };

  // load the mesh
  auto filename_string = filename.str();
  auto file_type = get_file_type( filename.str() );

  // Here we pull the mesh definition out of the 
  using globals::mesh_def;
  using globals::extra_mesh_info;
  
  extra_mesh_info = std::make_unique<extra_mesh_info_t>();
  bool needs_partitioning{false};

  if ( file_type == file_type_t::exodus ) {
    mesh_def = std::make_unique<exodus_definition_t>( filename_string );
    needs_partitioning = comm_size > 1;
  }
  else if ( file_type == file_type_t::partitioned_exodus ) {
    // get extension, which should be number of ranks
    auto ext = ristra::utils::file_extension(filename_string);
    // how many digits in padded number
    auto n = ext.size();
    // figure out this processors filename 
    auto my_filename = filename_string + "." + ristra::utils::zero_padded(rank, n);
    mesh_def = std::make_unique<exodus_definition_t>( my_filename, false );
  }
  else
    THROW_IMPLEMENTED_ERROR( "Unknown mesh file type" );

  // create a vector of colorings and color info for each dimensional entity
  auto & entities = context.coloring_map();
  auto & gathered_color_info = context.coloring_info_map();

  std::map< size_t,  flecsi::coloring::coloring_info_t > entity_color_info;

  //----------------------------------------------------------------------------
  // Cell Coloring
  //----------------------------------------------------------------------------

  // Cells index coloring.
  auto & cells = entities[index_spaces::cells];
  auto & cell_color_info = entity_color_info[index_spaces::cells];
  
  // distributed compressed row storage
  flecsi::coloring::dcrs_t dcrs;

  if ( needs_partitioning ) {

    if ( rank == 0 ) std::cout << "Partitioning mesh..." << std::flush;
    // Create the dCRS representation for the distributed colorer.
    // This essentialy makes the graph of the dual mesh.
    mesh_def->create_graph( num_dims, 0, num_dims, dcrs );

    // Create a colorer instance to generate the primary coloring.
    auto colorer = std::make_unique<flecsi::coloring::parmetis_colorer_t>();

    // Create the primary coloring and partition the mesh.
    auto partitioning = colorer->new_color(dcrs);
    if ( rank == 0 ) std::cout << "done." << std::endl;

    // now migrate the entities to their respective ranks
    if ( rank == 0 ) std::cout << "Migrating mesh..." << std::flush;
    flecsi::coloring::migrate( num_dims, partitioning, dcrs, *mesh_def );
    if ( rank == 0 ) std::cout << "done." << std::endl;
  
  } // needs_partitioning
  
  //----------------------------------------------------------------------------
  // Cell Closure.  However many layers of ghost cells are needed are found
  // here.
  //----------------------------------------------------------------------------
  
  // now that partitioning is over, recreate the dual graph so that we can
  // determine ghost cells.  This time we want corner cells, so 1 shared 
  // vertex is eneough to consider a cell as a neighbor.
  // FIXME: Here is where you would specify the depth of ghost cells.
  mesh_def->create_graph( num_dims, 0, 1, dcrs );
  extra_mesh_info->cell_distribution = dcrs.distribution;

  // now we can build all the other connectivity
  mesh_def->build_connectivity();
  
  //----------------------------------------------------------------------------
  // Dump the partitioned mesh if requested
  //----------------------------------------------------------------------------

  // figure out this ranks file name
  if ( file_type == file_type_t::exodus )
  {
    auto basename = ristra::utils::basename( filename_string );
    auto output_prefix = ristra::utils::remove_extension( basename );
    auto output_extension = ristra::utils::file_extension(basename); 
    auto digits = ristra::utils::num_digits(comm_size);
    auto output_filename = output_prefix + "-partitioned." + output_extension + "." +
    ristra::utils::zero_padded(comm_size, digits) + "." +
    ristra::utils::zero_padded(rank, digits);
    auto exo_def = dynamic_cast<exodus_definition_t*>(mesh_def.get());
    if (rank == 0 && comm_size > 1)
      std::cout << "Writing partitioned mesh to " << output_filename << std::endl;
    exo_def->write( output_filename );
  }


  //----------------------------------------------------------------------------
  // Identify exclusive, shared, and ghost for other entities
  //----------------------------------------------------------------------------
  
  // Find exclusive, shared, and ghost cells..
  // Not a fan of sets, but alot of machinery is built up on it.
  flecsi::coloring::get_owner_info( dcrs, cells, cell_color_info );

  // Find exclusive, shared, and ghost for the auxiliary index spaces (i.e.
  // edges, vertices, faces)
  std::map< dom_dim_t, std::map<dom_dim_t, size_t> > counts;

  for ( int i=0; i<num_dims; ++i ) {
    const auto & cells2entity = mesh_def->entities_crs(num_dims, i);
    const auto & local2global = mesh_def->local_to_global(i);
    const auto & global2local = mesh_def->global_to_local(i);
    auto index_space_id = index_spaces::entity_map[0][i];
    flecsi::coloring::color_entities( cells2entity, local2global, global2local,
        cells, entities[index_space_id], entity_color_info[index_space_id],
        counts[{0,num_dims}][{0,i}] );
  }
  
 
  //----------------------------------------------------------------------------
  // Count connectivity sizes
  //----------------------------------------------------------------------------
 
  // already have cell to everthing connectivity, but still need ghost
  // connectivitiy
  for ( int to_dim=0; to_dim<num_dims; ++to_dim ) {
      flecsi::coloring::ghost_connectivity( *mesh_def, num_dims, to_dim,
          cells, extra_mesh_info->ghost_ids[{0,num_dims}],
          extra_mesh_info->ghost_connectivity[{0,num_dims}][{0,to_dim}] );
  }

  // Need rest of counts
  for ( int from_dim=1; from_dim<num_dims; ++from_dim ) {

    // get the local-to-global mapping for the from entity
    const auto & local2global = mesh_def->local_to_global(from_dim);

    // loop over the "to" dimension
    for ( int to_dim=0; to_dim<from_dim; ++to_dim ) {
      
      // storage for connectivity counts.  This is needed because 
      // there is some overlap beween the mesh definition and the 
      // exclusive/shared/ghost entities
      std::map<size_t, size_t> my_counts;

      // get local connecitivy and count it
      const auto & conn = mesh_def->entities_crs(from_dim, to_dim);
      for ( size_t i=0; i<conn.size(); ++i ) {
        auto global_id = local2global[i];
        my_counts.emplace(global_id, conn.at(i).size());
      }

      // get the ghost connectivity and count it
      auto & ghost_ids = extra_mesh_info->ghost_ids[{0,from_dim}];
      auto & ghost_conn = extra_mesh_info->ghost_connectivity[{0,from_dim}][{0,to_dim}];
      auto from_is = index_spaces::entity_map[0][from_dim];
      flecsi::coloring::ghost_connectivity( *mesh_def, from_dim, to_dim,
          entities.at(from_is), ghost_ids, ghost_conn );

      // get the ghost entitiy connecitivity and count it.
      for ( size_t i=0; i<ghost_ids.size(); ++i ) {
        auto global_id = ghost_ids[i];
        my_counts.emplace(global_id, ghost_conn.at(i).size());
      }

      // now sum the counts
      size_t cnt{0};
      for ( auto i : my_counts ) cnt += i.second;
      counts[std::make_pair(0,from_dim)][std::make_pair(0,to_dim)] = cnt;
    }
  }
  
  
  //----------------------------------------------------------------------------
  // Identify exclusive, shared and ghost for corners and wedges and side
  //----------------------------------------------------------------------------

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS

  make_corners( *mesh_def, *extra_mesh_info, entities, entity_color_info, counts);


#endif // FLECSI_SP_BURTON_MESH_EXTRAS

  //----------------------------------------------------------------------------
  // Add the results to the context
  //----------------------------------------------------------------------------

  // keep track of index spaces added
  std::vector<size_t> registered_index_spaces;

  // Gather the coloring info from all colors
  for ( int dom=0; dom<num_domains; ++dom ) {
    for ( int dim=0; dim<num_dims+1; ++dim ) {
      // skip empty slots
      auto index_space_id = index_spaces::entity_map[dom][dim];
      if ( entities.count(index_space_id) == 0 ||
           entity_color_info.count(index_space_id) == 0 )
        continue;
      // gather and set to context directly
      gathered_color_info[index_space_id] = 
        communicator->gather_coloring_info(entity_color_info[index_space_id]);
      // add index space to list
      registered_index_spaces.push_back( index_space_id );
    }
  }

  //----------------------------------------------------------------------------
  // add adjacency information
  //----------------------------------------------------------------------------
  
  // loop over each dimension and determine the adjacency sizes
  for ( int from_dim = 0; from_dim<=num_dims; ++from_dim ) {
    for (int to_dim = num_dims+1; to_dim-- > 0; ) {

      // skip the case where both dimensions are the same
      if ( from_dim == to_dim ) continue;
      
      // populate the adjacency information
      flecsi::coloring::adjacency_info_t ai;
      ai.index_space =
        index_spaces::connectivity_map[ from_dim ][ to_dim ];
      ai.from_index_space = index_spaces::entity_map[0][ from_dim ];
      ai.to_index_space = index_spaces::entity_map[0][to_dim ];
      ai.color_sizes.resize(comm_size);

      // loop over all cells and count the number of adjacencies
      size_t cnt = (to_dim>from_dim) ?
        counts.at(std::make_pair(0,to_dim)).at(std::make_pair(0,from_dim)) :
        counts.at(std::make_pair(0,from_dim)).at(std::make_pair(0,to_dim));

      // gather the results
      ai.color_sizes = communicator->gather_sizes( cnt );

      // add the result to the context
      context.add_adjacency(ai);

    }

  }


  //----------------------------------------------------------------------------
  // add index subspace mappings
  //----------------------------------------------------------------------------

  // loop over all dimensions, getting the list of ids that are
  // connected to owned cells
  for ( int to_dim = 0; to_dim<num_dims; ++to_dim ) {
    // connectivity
    const auto & ents = mesh_def->entities_crs(num_dims, to_dim);
    // storage for subspace ids
    std::vector<size_t> ids;
    // get all exclusive and shared ids.
    for ( auto c : cells.exclusive ) {
      auto start = ents.offsets[c.offset];
      auto end = ents.offsets[c.offset+1];
      ids.insert( ids.end(), &ents.indices[start], &ents.indices[end] );
    }
    for ( auto c : cells.shared ) {
      auto start = ents.offsets[c.offset];
      auto end = ents.offsets[c.offset+1];
      ids.insert( ids.end(), &ents.indices[start], &ents.indices[end] );
    }
    // get number of unique entries
    std::sort( ids.begin(), ids.end() );
    auto last = std::unique( ids.begin(), ids.end() );
    auto count = std::distance( ids.begin(), last );
    // subspace id happens to be the same as the dim
    context.add_index_subspace(to_dim, count);
  }


  //----------------------------------------------------------------------------
  // Allow sparse index spaces for any of the main index spaces
  //----------------------------------------------------------------------------

  for ( auto i : registered_index_spaces ) {
    flecsi::execution::context_t::sparse_index_space_info_t isi;
    isi.max_entries_per_index = max_entries;
    // figure out the maximum number of entities
    const auto & coloring = entities.at(i);
    auto num_ents =
      coloring.exclusive.size() +
      coloring.shared.size() +
      coloring.ghost.size();

    // worst case scenario (not sure what we get by allocating all this)
    isi.exclusive_reserve = communicator->get_max_request_size(
      isi.max_entries_per_index*num_ents);
    isi.index_space = i;
    context.set_sparse_index_space_info(isi);
  }
  

  clog(info) << "Finished mesh partitioning." << std::endl;
  if (rank == 0) std::cout << "Finished mesh partitioning." << std::endl;

} // partition_mesh


////////////////////////////////////////////////////////////////////////////////
/// \brief the main mesh initialization driver
////////////////////////////////////////////////////////////////////////////////
void initialize_mesh(
  utils::client_handle_w<burton_mesh_t> mesh
) {

  //----------------------------------------------------------------------------
  // Fill the mesh information
  //----------------------------------------------------------------------------

  // some constant expressions
  constexpr auto num_dims = burton_mesh_t::num_dimensions;

  // alias some types
  using real_t = burton_mesh_t::real_t;

  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  
  // Here we pull the mesh definition out of the 
  using globals::mesh_def;
  using globals::extra_mesh_info;
  
  if (rank == 0) std::cout << "Initializing mesh..." << std::flush;

  // fill the mesh
  create_cells( *mesh_def, *extra_mesh_info, mesh );
#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
  create_extras( *mesh_def, *extra_mesh_info, mesh );
#endif

  // initialize the mesh
  mesh.init();

  // create the subspaces
  create_subspaces( *mesh_def, mesh );

  // delete any leftover global memory
  mesh_def.reset();
  extra_mesh_info.reset();
  
  if (rank == 0) std::cout << "done" << std::endl;

} // initialize_mesh

///////////////////////////////////////////////////////////////////////////////
// Task Registration
///////////////////////////////////////////////////////////////////////////////
flecsi_register_mpi_task(partition_mesh, flecsi_sp::burton);
flecsi_register_task(initialize_mesh, flecsi_sp::burton, loc, index|flecsi::leaf);

#ifdef BURTON_ENABLE_APPLICATION_TLT_INIT
void application_tlt_init(int argc, char **argv);
#endif
} // namespace burton
} // namespace flecsi_sp



