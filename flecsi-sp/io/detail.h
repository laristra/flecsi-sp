/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <utility>
#include <vector>

//! \file
//!
//! some technical details used in IO of meshes

namespace flecsi_sp {
namespace io {
namespace detail {

using size_t = std::size_t;

//==============================================================================
//! \brief Transpose a connectivity array.
//! \tparam CONNECTIVITY_TYPE Something that behaves like a std::vector<> of
//! some integer type.
//! \param in The connectivity array to transpose.
//! \param out The transposed connectivity array.
//==============================================================================
template<typename CONNECTIVITY_TYPE>
void
transpose(const CONNECTIVITY_TYPE & in, CONNECTIVITY_TYPE & out) {
  // invert the list
  auto num_from = in.size();

  for (size_t from = 0; from < num_from; ++from) {
    for (auto to : in[from]) {
      if (to >= out.size())
        out.resize(to + 1);
      out[to].push_back(from);
    }
  }
}

//==============================================================================
//! \brief Create connectivity among cells, edges, and vertices from known
//! cell-to-vertex connectivity.
//! \tparam CONNECTIVITY_TYPE Something that behaves like a
//! std::vector<std::vector<>> of some integer type.
//! \tparam FUNCTION Callable that takes a const reference to something that
//! behaves like a std::vector<> of some integer type and a
//! \em CONNECTIVITY_TYPE.  The intention is to build edges from vertices.
//==============================================================================
template<typename CONNECTIVITY_TYPE, typename FUNCTION>
void
build_connectivity(
    const CONNECTIVITY_TYPE & cell_to_vertex,
    CONNECTIVITY_TYPE & cell_to_edge,
    CONNECTIVITY_TYPE & edge_to_vertex,
    CONNECTIVITY_TYPE & sorted_edge_to_vertex,
    FUNCTION && build_edges_from_vertices) {

  // resize the new connectivity
  cell_to_edge.reserve(cell_to_edge.size() + cell_to_vertex.size());

  // invert the face id to vertices map (the one with sorted vertices)
  size_t edgeid = 0;  // starts from zero because cumulative list of edges
  std::map<std::vector<size_t>, size_t> edges;
  for ( const auto & vs : sorted_edge_to_vertex )
      edges[vs] = edgeid++;

  // starting id is not always zero, assume we are building in blocks
  // of assending edge ids
  edgeid = edge_to_vertex.size();

  // loop over cells, adding all of their edges to the table
  for (const auto & these_verts : cell_to_vertex) {
    // reference the storage for this cell's edges
    cell_to_edge.resize(cell_to_edge.size() + 1);
    auto & these_edges = cell_to_edge.back();

    // now build the edges for the cell
    CONNECTIVITY_TYPE new_edges;
    std::forward<FUNCTION>(build_edges_from_vertices)(these_verts, new_edges);

    // now look for exsiting vertex pairs in the edge-to-vertex master list
    // or add a new edge to the list.  add the matched edge id to the
    // cell-to-edge list
    for (auto && vs : new_edges) {
      // sort the vertices
      auto sorted_vs = vs;
      std::sort(sorted_vs.begin(), sorted_vs.end());
      // if we dont find the edge
      if (edges.find(sorted_vs) == edges.end()) {
        // add to the local reverse map
        edges.insert({sorted_vs, edgeid});
        // add to the original sorted and unsorted maps
        sorted_edge_to_vertex.emplace_back(std::move(sorted_vs));
        edge_to_vertex.emplace_back(std::move(vs));
        // add to the list of edges
        these_edges.push_back(edgeid);
        // bump counter
        edgeid++;
      } else {
        // if we do find the edge
        // just add the id to the list of edges
        these_edges.push_back(edges[sorted_vs]);
      }
    }
  }  // for
}  // build_connectivity

template<
  typename CONNECTIVITY_TYPE,
  typename SORTED_CONNECTIVITY_TYPE,
  typename FUNCTION>
void
new_build_connectivity(
    const CONNECTIVITY_TYPE & cell_to_vertex,
    CONNECTIVITY_TYPE & cell_to_edge,
    CONNECTIVITY_TYPE & edge_to_vertex,
    SORTED_CONNECTIVITY_TYPE & sorted_vertices_to_edges,
    FUNCTION && build_edges_from_vertices) {

  using index_t = typename CONNECTIVITY_TYPE::value_type;

  // some temporary storage
  CONNECTIVITY_TYPE new_edges;
  std::vector<index_t> sorted_vs;
  std::vector<index_t> these_edges;
  
  // loop over cells, adding all of their edges to the table
  for (const auto & these_verts : cell_to_vertex) {

    // clear this cells edges
    these_edges.clear();

    // now build the edges for the cell
    new_edges.clear();
    std::forward<FUNCTION>(build_edges_from_vertices)(these_verts, new_edges);

    // now look for exsiting vertex pairs in the edge-to-vertex master list
    // or add a new edge to the list.  add the matched edge id to the
    // cell-to-edge list
    for (const auto & vs : new_edges) {
      // sort the vertices
      sorted_vs.assign(vs.begin(), vs.end());
      std::sort(sorted_vs.begin(), sorted_vs.end());
      // if we dont find the edge
      auto it = sorted_vertices_to_edges.find(sorted_vs);
      if (it == sorted_vertices_to_edges.end())
      {
        // add to the local reverse map
        auto edgeid = sorted_vertices_to_edges.size();
        sorted_vertices_to_edges.emplace(sorted_vs, edgeid);
        // add to the original sorted and unsorted maps
        edge_to_vertex.push_back(vs);
        // add to the list of edges
        these_edges.push_back(edgeid);
      } else {
        // if we do find the edge
        // just add the id to the list of edges
        these_edges.push_back(it->second);
      }
    }

    // now add this cells edges
    cell_to_edge.push_back( these_edges );

  }  // for
}  // build_connectivity

//==============================================================================
//! \brief Intersect two connectivity arrays; that is, given X-to-Y and Y-to-Z
//! connectivity, build X-to-Z connectivity.
//! \tparam CONNECTIVITY_TYPE Something that behaves like a
//! std::vector<std::vector<>> of some integer type.
//==============================================================================
template<typename CONNECTIVITY_TYPE>
void
intersect(
    const CONNECTIVITY_TYPE & cell_to_face,
    const CONNECTIVITY_TYPE & face_to_edge,
    CONNECTIVITY_TYPE & cell_to_edge) {

  // resize the result array to size
  cell_to_edge.reserve(cell_to_edge.size() + cell_to_face.size());

  // loop over cells, collecting their face edges
  for (const auto & these_faces : cell_to_face) {
    // reference the storage for this cells edges
    cell_to_edge.resize(cell_to_edge.size() + 1);
    auto & these_edges = cell_to_edge.back();

    // now add the edges of this cell
    for (auto f : these_faces)
      for (auto e : face_to_edge.at(f))
        these_edges.push_back(e);

    // sort them and remove the non-unique ones
    std::sort(these_edges.begin(), these_edges.end());
    auto last = std::unique(these_edges.begin(), these_edges.end());
    these_edges.erase(last, these_edges.end());
  }  // cells
}  // intersect

template<typename CONNECTIVITY_TYPE>
void
new_intersect(
    const CONNECTIVITY_TYPE & cell_to_face,
    const CONNECTIVITY_TYPE & face_to_edge,
    CONNECTIVITY_TYPE & cell_to_edge) {

  // some temporary storage
  using index_t = typename CONNECTIVITY_TYPE::value_type;
  std::vector<index_t> these_edges;

  // loop over cells, collecting their face edges
  for (const auto & these_faces : cell_to_face) {
    // reference the storage for this cells edges
    these_edges.clear();

    // now add the edges of this cell
    for (auto f : these_faces)
      for (auto e : face_to_edge.at(f)) {
        // set insertion is not the quickest, but maintains order
        auto it = std::find( these_edges.begin(), these_edges.end(), e );
        if ( it == these_edges.end() ) these_edges.push_back(e);
      }

    // sort them and remove the non-unique ones
    // ... This mightbe quicker but does not maintain order
    // std::sort(these_edges.begin(), these_edges.end());
    // auto last = std::unique(these_edges.begin(), these_edges.end());
    // these_edges.erase(last, these_edges.end());

    // now add this cells edges
    cell_to_edge.push_back( these_edges );

  }  // cells
}  // intersect


//==============================================================================
// Lambda function to process a block

template< typename T, typename U, typename V, typename W>
void filter_block(
    const T & counts_in,
    const U & indices_in,
    size_t min,
    size_t max,
    size_t & counter,
    V & offsets_out,
    W & indices_out)
{
  auto num_offsets = counts_in.size();
  auto num_indices = indices_in.size();

  // if this is outside the range, just increment counter and return
  if ( counter + num_offsets <= min || counter > max )
  {
    counter += num_offsets;
    return;
  }

  // storage for element verts
  offsets_out.reserve( offsets_out.size() + num_offsets );
  indices_out.reserve( indices_out.size() + num_indices );

  if (offsets_out.empty()) offsets_out.emplace_back(0);
  
  // create cells in mesh
  size_t base = 0;
  for (size_t e = 0; e < num_offsets; ++e, ++counter) {
    // get the number of nodes
    auto cnt = counts_in[e];
    // the global id is within the range
    if ( counter >= min && counter <= max )
    {
      // copy local vertices into vector ( exodus uses 1 indexed arrays )
      for (int v = 0; v < cnt; v++)
        indices_out.emplace_back(indices_in[base + v] - 1);
      // add the row
      offsets_out.emplace_back( offsets_out.back() + cnt );
    }
    // base offset into elt_conn
    base += cnt;
  }
}

template< typename T, typename U, typename V, typename W, typename X >
void filter_sides(
    size_t ss_id,
    const T & side_set_node_cnt_list,
    const T & side_set_node_list,
    const T & side_set_elem_list,
    size_t cell_min,
    size_t cell_max,
    U & side_id_,
    V & element_to_sides_,
    W & side_to_vertices_,
    X & side_sets_ )
{

  auto num_side_in_set = side_set_elem_list.size();
  
  // filter sides for elementes i own
  std::vector<size_t> vs;
  for ( size_t j=0; j<num_side_in_set; ++j ) {
    auto global_id = side_set_elem_list[j] - 1;
    if (global_id >= cell_min && global_id <= cell_max )
    {
      auto & this_ss = side_sets_[ss_id];
      // get the vertices
      auto side_start = side_set_node_cnt_list[j];
      auto side_end = side_set_node_cnt_list[j+1];
      vs.clear();
      vs.reserve(side_end - side_start);
      for ( size_t k=side_start; k<side_end; ++k )
        vs.emplace_back( side_set_node_list[k] - 1 );
      // add side info
      auto side_id = side_to_vertices_.size();
      side_id_.emplace_back( ss_id-1 );
      element_to_sides_[global_id].push_back(side_id);
      side_to_vertices_.push_back(vs);
    }
  }

}

}  // namespace detail
}  // namespace io
}  // namespace flecsi_sp
