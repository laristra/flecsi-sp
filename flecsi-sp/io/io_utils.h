/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

/// \file

// user includes
#include <flecsi/utils/logging.h>

#include "H5Cpp.h"

// system includes
#include <algorithm>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace flecsi_sp {
namespace io {

namespace detail {
//==============================================================================
//! \brief Transpose a connectivity array
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
//! \brief Create connectivity through
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
  auto cellid = cell_to_edge.size();

  // invert the face id to vertices map (the one with sorted vertices)
  size_t edgeid = 0; // starts from zero because cumulative list of edges
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
        sorted_edge_to_vertex.emplace_back( std::move(sorted_vs) );
        edge_to_vertex.emplace_back(std::move(vs));
        // add to the list of edges
        these_edges.push_back(edgeid);
        // bump counter
        edgeid++;
      }
      // if we do find the edge
      else {
        // just add the id to the list of edges  
        these_edges.push_back(edges[sorted_vs]);
      }
    }

    cellid++;
    
  } // for
}

//==============================================================================
//! \brief Intersect two connectivity arrays array
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

  } // cells
}

} // namespace detail

} // namespace io
} // namespace flecsi
