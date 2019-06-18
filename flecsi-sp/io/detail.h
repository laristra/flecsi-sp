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

template<typename CONNECTIVITY_TYPE, typename FUNCTION>
void
new_build_connectivity(
    const CONNECTIVITY_TYPE & cell_to_vertex,
    CONNECTIVITY_TYPE & cell_to_edge,
    CONNECTIVITY_TYPE & edge_to_vertex,
    FUNCTION && build_edges_from_vertices) {

  using index_t = typename CONNECTIVITY_TYPE::value_type;

  // some temporary storage
  CONNECTIVITY_TYPE new_edges;
  std::vector<index_t> sorted_vs;
  std::vector<index_t> these_edges;
  
  // a map for searching: vertices -> edge id
  std::map<std::vector<index_t>, index_t> edges;

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
      if (edges.find(sorted_vs) == edges.end()) {
        // add to the local reverse map
        auto edgeid = edges.size();
        edges.emplace(sorted_vs, edgeid);
        // add to the original sorted and unsorted maps
        edge_to_vertex.push_back(vs);
        // add to the list of edges
        these_edges.push_back(edgeid);
      } else {
        // if we do find the edge
        // just add the id to the list of edges
        these_edges.push_back(edges[sorted_vs]);
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

}  // namespace detail
}  // namespace io
}  // namespace flecsi_sp
