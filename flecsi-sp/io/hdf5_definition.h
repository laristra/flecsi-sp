/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

/// \file

// user includes
#include <flecsi/topology/mesh_definition.h>
#include <flecsi/utils/logging.h>

// thirdparty includes
#include <hdf5II.h>

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
#if 0
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
#endif
////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class hdf5_base__ {

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the size type
  using size_t = std::size_t;
  //! \brief the counter type
  using counter_t = flecsi::utils::counter_t;
  //! \brief the floating point type
  using real_t = T;
  //! \brief the type used for indexing arrays
  using index_t = std::size_t;

  //! \brief an alias for the vector class
  template<typename U>
  using vector = typename std::vector<U>;

  //! \brief an alias for the matrix class
  template<typename U>
  using sparse_matrix = vector<vector<U>>;

  //! \brief the data type for an index vector
  using index_vector_t = vector<index_t>;

  //! \brief the data type for connectivity
  using connectivity_t = sparse_matrix<index_t>;

  //! the number of dimensions
  static constexpr size_t num_dims = D;

  enum class block_t { tri, quad, polygon, tet, hex, polyhedron, unknown };


  //============================================================================
  //! \brief open the file for reading or writing
  //! \param [in] name  The name of the file to open.
  //! \param [in] mode  The mode to open the file in.
  //! \return The hdf5 handle for the open file.
  //============================================================================
  static auto open(const std::string & name, std::ios_base::openmode mode) {
    //FIXME
  }

  //============================================================================
  //! \brief close the file once completed reading or writing
  //! \param [in] hdf5_id  The hdf5 file id.
  //============================================================================
  static void close(int hdf5_id) {
    //FIXME
  }

#if 0
  //============================================================================
  //! \brief check the integer status
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return true if hdf5 is using 64 bit integers
  //============================================================================
  static auto is_int64(int hdf5_id) {
    return (hdf5_int64_status(hdf5_id) & EX_IDS_INT64_API);
  }
#endif

  //============================================================================
  //! \brief Helper function to make and initialize a set of hdf5 parameters.
  //! \return the hdf5 parameters
  //============================================================================
  static auto make_params() {
    hdf5_init_params par;
    strcpy(par.title, "Hdf5 output from flecsi.");
    par.num_dim = num_dims;
    par.num_nodes = 0;
    par.num_edge = 0;
    par.num_edge_blk = 0;
    par.num_face = 0;
    par.num_face_blk = 0;
    par.num_elem = 0;
    par.num_elem_blk = 0;
    par.num_node_sets = 0;
    par.num_edge_sets = 0;
    par.num_face_sets = 0;
    par.num_side_sets = 0;
    par.num_elem_sets = 0;
    par.num_node_maps = 0;
    par.num_edge_maps = 0;
    par.num_face_maps = 0;
    par.num_elem_maps = 0;
    return par;
  }

  //============================================================================
  //! \brief read the hdf5 parameters from a file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the hdf5 parameters
  //============================================================================
  static auto read_params(int hdf5_id) {
    hdf5_init_params hdf5_params;
    //FIXME
    return hdf5_params;
  }

  //============================================================================
  //! \brief write the hdf5_ parameters to a file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \param [in] hdf5_params  The hdf5 parameters
  //============================================================================
  static void write_params(int hdf5_id, const hdf5_init_params & hdf5_params) {
    // verify mesh dimension
    if (num_dims != hdf5_params.num_dim)
      clog_fatal(
          "HDF5 dimension mismatch: Expected dimension ("
          << num_dims << ") /= HDF5 dimension (" << hdf5_params.num_dim
          << ")");

    // put the initialization parameters
    auto status = hdf5_put_init_ext(hdf5_id, &hdf5_params);
    if (status)
      clog_fatal(
          "Problem putting hdf5 file parameters, hdf5_put_init_ext() returned "
          << status);
  }

  //============================================================================
  //! \brief read the coordinates of the mesh from a file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the vertex coordinates
  //============================================================================
  static auto read_point_coords(int hdf5_id, size_t num_nodes) {

    // get the number of nodes
    if (num_nodes <= 0)
      clog_fatal(
          "Exodus file has zero nodes, or parmeters haven't been read yet.");

    // read nodes
    vector<real_t> vertex_coord(num_dims * num_nodes);

    // hdf5 is kind enough to fetch the data in the real type we ask for
    auto status = hdf5_get_coord(
        hdf5_id, vertex_coord.data(), vertex_coord.data() + num_nodes,
        vertex_coord.data() + 2 * num_nodes);

    if (status)
      clog_fatal(
          "Problem getting vertex coordinates from hdf5 file, "
          << " hdf5_get_coord() returned " << status);

    return vertex_coord;
  }

  //============================================================================
  //! \brief write the coordinates of the mesh from a file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \param [in] vertex_coord  the vertex coordinates
  //============================================================================
  template<typename V>
  static void write_point_coords(int hdf5_id, const V & vertex_coord) {

    if (vertex_coord.empty())
      return;

    auto num_nodes = vertex_coord.size() / num_dims;

    // hdf5 is kind enough to fetch the data in the real type we ask for
    auto status = hdf5_put_coord(
        hdf5_id, vertex_coord.data(), vertex_coord.data() + num_nodes,
        vertex_coord.data() + 2 * num_nodes);

    if (status)
      clog_fatal(
          "Problem putting vertex coordinates to hdf5 file, "
          << " hdf5_put_coord() returned " << status);
  }

  //============================================================================
  //! \brief read the block ids from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename V>
  static void write_node_set(
      int hdf5_id,
      size_t node_set_id,
      const std::string & name,
      const vector<V> & vertex_list) {

    // some type aliases
    using hdf5_index_t = U;

    if (vertex_list.empty())
      return;

    // set the node set parameters
    hdf5_index_t num_dist_in_set = 0;
    hdf5_index_t num_nodes_this_set = vertex_list.size();
    auto status = hdf5_put_node_set_param(
        hdf5_id, node_set_id, num_nodes_this_set, num_dist_in_set);
    if (status)
      clog_fatal(
          "Problem writing node set param to hdf5 file, "
          << " hdf5_put_node_set_param() returned " << status);

    // copy the vertex ids
    vector<hdf5_index_t> node_set;
    node_set.reserve(vertex_list.size());

    for (auto v : vertex_list)
      node_set.push_back(v + 1); // hdf5 uses 1-based ids

    // write the node set
    status = hdf5_put_node_set(hdf5_id, node_set_id, node_set.data());
    if (status)
      clog_fatal(
          "Problem writing node set to hdf5 file, "
          << " hdf5_put_node_set() returned " << status);

    // write the set name
    // FIXME following
    /*
    status = hdf5_put_name(hdf5_id, EX_NODE_SET, node_set_id, name.c_str());
    if (status)
      clog_fatal(
          "Problem writing node set name to hdf5 file, "
          << " hdf5_put_name() returned " << status);
    */
  }

  //============================================================================
  //! \brief read the block ids from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto
  read_block_ids(int hdf5_id, hdf5_entity_type obj_type, size_t num_blocks) {
    // some type aliases
    using hdf5_index_t = U;

    // final storage
    vector<index_t> ids(num_blocks);

    if (num_blocks > 0) {

      // get the ids first
      vector<hdf5_index_t> block_ids(num_blocks);
      auto status = hdf5_get_ids(hdf5_id, obj_type, block_ids.data());
      if (status)
        clog_fatal(
            "Problem reading block ids, hdf5_get_ids() returned " << status);

      // now convert them
      std::transform(
          block_ids.begin(), block_ids.end(), ids.begin(),
          [](auto id) { return id; });
    }

    // now return them
    return ids;
  }

  //============================================================================
  //! \brief read the element blocks from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename ENTITY_CONN>
  static void write_block(
      int hdf5_id,
      size_t blk_id,
      const std::string & name,
      hdf5_entity_type entity_type,
      const char * entity_description,
      size_t num_entities,
      ENTITY_CONN && entity_conn) {

    // some type aliases
    using hdf5_index_t = U;

    // check if face data is provided instead of node data
    auto is_face_data = (num_dims == 3 && entity_type == EX_ELEM_BLOCK);

    // guess how many elements are in each connectivity slot
    auto num_conn_guess = num_dims * num_dims;

    // build the connectivitiy list for the block
    vector<hdf5_index_t> entity_nodes;
    vector<int> entity_node_counts;
    entity_nodes.reserve(num_entities * num_conn_guess);
    entity_node_counts.reserve(num_entities);

    // temporary storage for each connecitivity slot
    vector<hdf5_index_t> temp_entity_nodes;
    temp_entity_nodes.reserve(num_conn_guess);

    for (size_t e = 0; e < num_entities; ++e) {
      // get the connectivity from the user-provided function
      std::forward<ENTITY_CONN>(entity_conn)(e, temp_entity_nodes);
      // store them in the master list ( hdf5 wants 1-index arrays )
      for (auto i : temp_entity_nodes)
        entity_nodes.emplace_back(i + 1);
      entity_node_counts.push_back(temp_entity_nodes.size());
      // reset temp storage
      temp_entity_nodes.clear();
    }

    // the total size needed to hold the element connectivity
    hdf5_index_t num_nodes_this_blk = entity_nodes.size();
    hdf5_index_t num_entries_this_blk = entity_node_counts.size();

    // set the block header
    hdf5_index_t num_attr_per_entry = 0;
    hdf5_index_t num_nodes_per_entry = is_face_data ? 0 : num_nodes_this_blk;
    hdf5_index_t num_edges_per_entry = 0;
    hdf5_index_t num_faces_per_entry = is_face_data ? num_nodes_this_blk : 0;
    auto status = hdf5_put_block(
        hdf5_id, entity_type, blk_id, entity_description, num_entries_this_blk,
        num_nodes_per_entry, num_edges_per_entry, num_faces_per_entry,
        num_attr_per_entry);
    if (status)
      clog_fatal(
          "Problem writing block to hdf5 file, "
          << " hdf5_put_block() returned " << status);

    // write the block name
    status = hdf5_put_name(hdf5_id, entity_type, blk_id, name.c_str());
    if (status)
      clog_fatal(
          "Problem writing block name to hdf5 file, "
          << " hdf5_put_name() returned " << status);

    // write connectivity
    auto node_conn = is_face_data ? nullptr : entity_nodes.data();
    auto elem_edge_conn = nullptr;
    auto elem_face_conn = is_face_data ? entity_nodes.data() : nullptr;
    status = hdf5_put_conn(
        hdf5_id, entity_type, blk_id, node_conn, elem_edge_conn, elem_face_conn);
    if (status)
      clog_fatal(
          "Problem writing block connectivity to hdf5 file, "
          << " hdf5_put_conn() returned " << status);

    // write counts
    status = hdf5_put_entity_count_per_polyhedra(
        hdf5_id, entity_type, blk_id, entity_node_counts.data());
    if (status)
      clog_fatal(
          "Problem writing block counts to hdf5_id file, "
          << " hdf5_put_entity_count_per_polyhedra() returned " << status);

  }; // write block

  //============================================================================
  //! \brief read the element blocks from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto read_block(
      int hdf5_id,
      hdf5_entity_id blk_id,
      hdf5_entity_type entity_type,
      connectivity_t & entities) {
    // some type aliases
    using hdf5_index_t = U;

    // get the info about this block
    hdf5_index_t num_elem_this_blk = 0;
    hdf5_index_t num_faces_per_elem = 0;
    hdf5_index_t num_edges_per_elem = 0;
    hdf5_index_t num_nodes_per_elem = 0;
    hdf5_index_t num_attr = 0;
    char elem_type[MAX_STR_LENGTH];
    auto status = hdf5_get_block(
        hdf5_id, entity_type, blk_id, elem_type, &num_elem_this_blk,
        &num_nodes_per_elem, &num_edges_per_elem, &num_faces_per_elem,
        &num_attr);
    if (status)
      clog_fatal("Problem reading block, hdf5_get_block() returned " << status);

    //------------------------------------------------------------------------
    // polygon data
    if (strcasecmp("nsided", elem_type) == 0) {

      // the number of nodes per element is really the number of nodes
      // in the whole block
      auto num_nodes_this_blk = num_nodes_per_elem;

      // get the number of nodes per element
      vector<int> elem_node_counts(num_elem_this_blk);
      status = hdf5_get_entity_count_per_polyhedra(
          hdf5_id, entity_type, blk_id, elem_node_counts.data());
      if (status)
        clog_fatal(
            "Problem getting element node numbers, "
            << "hdf5_get_entity_count_per_polyhedra() returned " << status);

      // read element definitions
      vector<hdf5_index_t> elem_nodes(num_nodes_this_blk);
      status = hdf5_get_conn(
          hdf5_id, entity_type, blk_id, elem_nodes.data(), nullptr, nullptr);
      if (status)
        clog_fatal(
            "Problem getting element connectivity, hdf5_get_elem_conn() "
            << "returned " << status);

      // storage for element verts
      vector<index_t> elem_vs;
      elem_vs.reserve(num_dims * elem_node_counts[0]);

      // create cells in mesh
      size_t base = 0;
      for (counter_t e = 0; e < num_elem_this_blk; ++e) {
        elem_vs.clear();
        // get the number of nodes
        num_nodes_per_elem = elem_node_counts[e];
        // copy local vertices into vector ( hdf5 uses 1 indexed arrays )
        for (int v = 0; v < num_nodes_per_elem; v++)
          elem_vs.emplace_back(elem_nodes[base + v] - 1);
        // add the row
        entities.push_back(elem_vs);
        // base offset into elt_conn
        base += num_nodes_per_elem;
      }

      return block_t::polygon;

    }
    //------------------------------------------------------------------------
    // polygon data
    else if (strcasecmp("nfaced", elem_type) == 0) {

      // the number of faces per element is really the number of
      // faces in the whole block ( includes duplicate / overlapping
      // nodes )
      auto num_face_this_blk = num_faces_per_elem;

      // get the number of nodes per element
      vector<int> elem_face_counts(num_face_this_blk);
      status = hdf5_get_entity_count_per_polyhedra(
          hdf5_id, entity_type, blk_id, elem_face_counts.data());
      if (status)
        clog_fatal(
            "Problem reading element node info, "
            << "hdf5_get_entity_count_per_polyhedra() returned " << status);

      // read element definitions
      vector<hdf5_index_t> elem_faces(num_face_this_blk);
      status = hdf5_get_conn(
          hdf5_id, entity_type, blk_id, nullptr, nullptr, elem_faces.data());
      if (status)
        clog_fatal(
            "Problem getting element connectivity, hdf5_get_conn() "
            << "returned " << status);

      // storage for element faces
      vector<index_t> elem_fs;
      elem_fs.reserve(elem_face_counts[0]);

      // create cells in mesh
      for (counter_t e = 0, base = 0; e < num_elem_this_blk; ++e) {
        // reset storage
        elem_fs.clear();
        // get the number of faces
        num_faces_per_elem = elem_face_counts[e];
        // copy local vertices into vector ( hdf5 uses 1 indexed arrays )
        for (int v = 0; v < num_faces_per_elem; v++)
          elem_fs.emplace_back(elem_faces[base + v] - 1);
        // base offset into elt_conn
        base += num_faces_per_elem;
        // add vertex list to master list
        entities.push_back(elem_fs);
      }

      return block_t::polyhedron;

    }
    //------------------------------------------------------------------------
    // fixed element size
    else {

      // read element definitions
      vector<hdf5_index_t> elt_conn(num_elem_this_blk * num_nodes_per_elem);
      status = hdf5_get_elem_conn(hdf5_id, blk_id, elt_conn.data());
      if (status)
        clog_fatal(
            "Problem getting element connectivity, hdf5_get_elem_conn() "
            << "returned " << status);

      // storage for element verts
      vector<index_t> elem_vs;
      elem_vs.reserve(num_nodes_per_elem);

      // create cells in mesh
      for (counter_t e = 0; e < num_elem_this_blk; ++e) {
        elem_vs.clear();
        // base offset into elt_conn
        auto b = e * num_nodes_per_elem;
        // copy local vertices into vector ( hdf5 uses 1 indexed arrays )
        for (int v = 0; v < num_nodes_per_elem; v++)
          elem_vs.emplace_back(elt_conn[b + v] - 1);
        // add the row
        entities.push_back(elem_vs);
      }

      // return element type
      if (strcasecmp("tri3", elem_type) == 0)
        return block_t::tri;
      else if (
          strcasecmp("quad4", elem_type) == 0 ||
          strcasecmp("shell4", elem_type) == 0)
        return block_t::quad;
      else if (
          strcasecmp("tet4", elem_type) == 0 ||
          strcasecmp("tetra", elem_type) == 0)
        return block_t::tet;
      else if (strcasecmp("hex8", elem_type) == 0)
        return block_t::hex;
      else {
        clog_fatal("Unknown block type, " << elem_type);
        return block_t::unknown;
      }

    } // element type
    //------------------------------------------------------------------------
  }

  //============================================================================
  //! \brief read the element blocks from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto
  read_face_block(int hdf5_id, hdf5_entity_id blk_id, connectivity_t & faces) {

    return read_block<U>(hdf5_id, blk_id, EX_FACE_BLOCK, faces);
  }

  //============================================================================
  //! \brief read the element blocks from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename CONN_TYPE>
  static void write_face_block(
      int hdf5_id,
      size_t blk_id,
      const std::string & name,
      size_t num_faces,
      CONN_TYPE && face_conn) {

  /* FIXME
    write_block<U>(
        hdf5_id, blk_id, name, EX_FACE_BLOCK, "nsided", num_faces,
        std::forward<CONN_TYPE>(face_conn));
   */

  }

  //============================================================================
  //! \brief read the element blocks from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto read_element_block(
      int hdf5_id,
      hdf5_entity_id elem_blk_id,
      connectivity_t & elements) {

    return read_block<U>(hdf5_id, elem_blk_id, EX_ELEM_BLOCK, elements);
  }

  //============================================================================
  //! \brief read the element blocks from an hdf5 file.
  //! \param [in] hdf5_id  The hdf5 file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename CONN_TYPE>
  static void write_element_block(
      int hdf5_id,
      size_t blk_id,
      const std::string & name,
      size_t num_elems,
      CONN_TYPE && element_conn) {

    auto entity_desc = (num_dims == 3) ? "nfaced" : "nsided";
    write_block<U>(
        hdf5_id, blk_id, name, EX_ELEM_BLOCK, entity_desc, num_elems,
        std::forward<CONN_TYPE>(element_conn));
  }
};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class hdf5_definition__ {};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class hdf5_definition__<2, T> : public flecsi::topology::mesh_definition__<2>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = hdf5_base__<2, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition__<2>;

  //! the number of dimensions
  using mesh_definition_t::dimension;

  //! the floating point type
  using real_t = typename base_t::real_t;
  //! the index type
  using index_t = typename base_t::index_t;

  //! the vector type
  template<typename U>
  using vector = typename base_t::template vector<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Default constructor
  hdf5_definition__() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  hdf5_definition__(const std::string & filename) {
    read(filename);
  }

  /// Copy constructor (disabled)
  hdf5_definition__(const hdf5_definition__ &) = delete;

  /// Assignment operator (disabled)
  hdf5_definition__ & operator=(const hdf5_definition__ &) = delete;

  /// Destructor
  ~hdf5_definition__() = default;

  //============================================================================
  //! \brief Implementation of hdf5 mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void read(const std::string & name) {

    clog(info) << "Reading mesh from: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file

    // open the hdf5 file
    auto hdf5_id = base_t::open(name, std::ios_base::in);
    if (hdf5_id < 0)
      clog_fatal("Problem reading hdf5 file");

    // get the initialization parameters
    auto exo_params = base_t::read_params(hdf5_id);

    // check the integer type used in the hdf5 file
    auto int64 = base_t::is_int64(hdf5_id);

    //--------------------------------------------------------------------------
    // read coordinates

    vertices_ = base_t::read_point_coords(hdf5_id, exo_params.num_nodes);
    clog_assert(
        vertices_.size() == dimension() * exo_params.num_nodes,
        "Mismatch in read vertices");

    auto num_vertices = vertices_.size() / dimension();

    //--------------------------------------------------------------------------
    // element blocks

    auto num_elem_blk = exo_params.num_elem_blk;
    vector<index_t> elem_blk_ids;

    auto & cell_vertices_ref = entities_[2][0];

    // get the element block ids
    if (int64)
      elem_blk_ids = base_t::template read_block_ids<long long>(
          hdf5_id, EX_ELEM_BLOCK, num_elem_blk);
    else
      elem_blk_ids = base_t::template read_block_ids<int>(
          hdf5_id, EX_ELEM_BLOCK, num_elem_blk);

    // read each block
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      if (int64)
        base_t::template read_element_block<long long>(
            hdf5_id, elem_blk_ids[iblk], cell_vertices_ref);
      else
        base_t::template read_element_block<int>(
            hdf5_id, elem_blk_ids[iblk], cell_vertices_ref);
    }

    // check some assertions
    clog_assert(
        cell_vertices_ref.size() == exo_params.num_elem,
        "Mismatch in read blocks");

    //--------------------------------------------------------------------------
    // build the edges

    // reference storage for the cell edges and edge vertices
    auto & cell_edges_ref = entities_[2][1];
    auto & edge_vertices_ref = entities_[1][0];

    // temprary storage for matching edges
    auto edge_vertices_sorted = std::make_unique<connectivity_t>();

    // build the connecitivity array
    detail::build_connectivity(
        cell_vertices_ref, cell_edges_ref, edge_vertices_ref,
        *edge_vertices_sorted, [](const auto & vs, auto & edge_vs) {
          using connectivity_type = std::decay_t<decltype(edge_vs)>;
          using list_type = std::decay_t<decltype(*edge_vs.begin())>;
          for (auto v0 = std::prev(vs.end()), v1 = vs.begin(); v1 != vs.end();
               v0 = v1, ++v1)
            edge_vs.emplace_back(list_type{*v0, *v1});
        });

    // clear temporary storage
    edge_vertices_sorted.reset();

    // make storage for the edge vertices
    auto num_edges = edge_vertices_ref.size();

    //--------------------------------------------------------------------------
    // Create the remainder of the connectivities

    entities_[1][2].reserve(num_edges);
    entities_[0][2].reserve(num_vertices);
    entities_[0][1].reserve(num_vertices);

    detail::transpose(entities_[2][1], entities_[1][2]);
    detail::transpose(entities_[2][0], entities_[0][2]);
    detail::transpose(entities_[1][0], entities_[0][1]);

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(hdf5_id);
  }

  //============================================================================
  //! \brief Implementation of hdf5 mesh write for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  template<typename U = int>
  void write(
      const std::string & name,
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
          element_sets = {},
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
          node_sets = {}) const {

    clog(info) << "Writing mesh to: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file

    // open the hdf5 file
    auto hdf5_id = base_t::open(name, std::ios_base::out);

    // write the initialization parameters
    auto exo_params = base_t::make_params();
    auto num_cells = num_entities(dimension());
    exo_params.num_nodes = num_entities(0);
    exo_params.num_node_sets = node_sets.size();
    exo_params.num_elem_blk = element_sets.size() ? element_sets.size() : 1;

    if (element_sets.size()) {
      for (const auto & set : element_sets)
        exo_params.num_elem += set.second.size();
    } else
      exo_params.num_elem = num_cells;

    base_t::write_params(hdf5_id, exo_params);

    // check the integer type used in the hdf5 file
    auto int64 = base_t::is_int64(hdf5_id);

    //--------------------------------------------------------------------------
    // write coordinates

    base_t::write_point_coords(hdf5_id, vertices_);

    //--------------------------------------------------------------------------
    // element blocks

    const auto & cell_vertices = entities_.at(2).at(0);
    int elem_blk_id = 1;

    // add the whole element block
    if (!element_sets.size()) {

      auto cell_vertices_func = [&](auto c, auto & vert_list) {
        const auto & vs = cell_vertices[c];
        vert_list.insert(vert_list.end(), vs.begin(), vs.end());
      };

      if (int64)
        base_t::template write_element_block<long long>(
            hdf5_id, elem_blk_id++, "cells", num_cells, cell_vertices_func);
      else
        base_t::template write_element_block<int>(
            hdf5_id, elem_blk_id++, "cells", num_cells, cell_vertices_func);

    } // element block

    // or add the element sets
    for (const auto & set : element_sets) {

      const auto & set_name = set.first;
      const auto & set_elements = set.second;
      auto num_cells_set = set_elements.size();

      if (num_cells_set == 0)
        continue;

      auto cell_vertices_func = [&](auto c, auto & vert_list) {
        const auto & vs = cell_vertices[set_elements[c]];
        vert_list.insert(vert_list.end(), vs.begin(), vs.end());
      };

      if (int64)
        base_t::template write_element_block<long long>(
            hdf5_id, elem_blk_id++, set_name, num_cells_set, cell_vertices_func);
      else
        base_t::template write_element_block<int>(
            hdf5_id, elem_blk_id++, set_name, num_cells_set, cell_vertices_func);

    } // sets

    //--------------------------------------------------------------------------
    // Node sets

    int node_set_id = 1;

    for (const auto & set : node_sets) {

      const auto & set_name = set.first;
      const auto & set_nodes = set.second;
      auto num_nodes_set = set_nodes.size();

      if (num_nodes_set == 0)
        continue;

      if (int64)
        base_t::template write_node_set<long long>(
            hdf5_id, ++node_set_id, set_name, set_nodes);
      else
        base_t::template write_node_set<int>(
            hdf5_id, ++node_set_id, set_name, set_nodes);

    } // sets

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(hdf5_id);
  }

  //============================================================================
  // Required Overrides
  //============================================================================

  /// Return the number of entities of a particular dimension
  /// \param [in] dim  The entity dimension to query.
  size_t num_entities(size_t dim) const override {
    switch (dim) {
      case 0:
        return vertices_.size() / dimension();
      case 1:
      case 2:
        return entities_.at(dim).at(0).size();
      default:
        clog_fatal(
            "Dimension out of range: 0 < " << dim << " </ " << dimension());
        return 0;
    }
  }

  /// Return the set of vertices of a particular entity.
  /// \param [in] dimension  The entity dimension to query.
  /// \param [in] entity_id  The id of the entity in question.
  const auto & entities(size_t from_dim, size_t to_dim) const {
    return entities_.at(from_dim).at(to_dim);
  } // vertices

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in question.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override {
    return entities_.at(from_dim).at(to_dim).at(from_id);
  } // vertices

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  template<typename POINT_TYPE>
  auto vertex(size_t vertex_id) const {
    auto num_vertices = vertices_.size() / dimension();
    POINT_TYPE p;
    for (int i = 0; i < dimension(); ++i)
      p[i] = vertices_[i * num_vertices + vertex_id];
    return p;
  } // vertex

private:
  //============================================================================
  // Private data
  //============================================================================

  //! \brief storage for element verts
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;
};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the hdf5 file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class hdf5_definition__<3, T> : public flecsi::topology::mesh_definition__<3>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = hdf5_base__<3, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition__<3>;

  //! the number of dimensions
  using mesh_definition_t::dimension;

  //! the floating point type
  using real_t = typename base_t::real_t;
  //! the index type
  using index_t = typename base_t::index_t;

  //! the vector type
  template<typename U>
  using vector = typename base_t::template vector<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Default constructor
  hdf5_definition__() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  hdf5_definition__(const std::string & filename) {
    read(filename);
  }

  /// Copy constructor (disabled)
  hdf5_definition__(const hdf5_definition__ &) = delete;

  /// Assignment operator (disabled)
  hdf5_definition__ & operator=(const hdf5_definition__ &) = delete;

  /// Destructor
  ~hdf5_definition__() = default;

  //============================================================================
  //! \brief Implementation of hdf5 mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void read(const std::string & name) {

    clog(info) << "Reading mesh from: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file

    // open the hdf5 file
    auto hdf5_id = base_t::open(name, std::ios_base::in);

    // get the initialization parameters
    auto exo_params = base_t::read_params(hdf5_id);

    // check the integer type used in the hdf5 file
    auto int64 = base_t::is_int64(hdf5_id);

    //--------------------------------------------------------------------------
    // read coordinates

    vertices_ = base_t::read_point_coords(hdf5_id, exo_params.num_nodes);

    //--------------------------------------------------------------------------
    // read face blocks

    auto num_face_blk = exo_params.num_face_blk;
    vector<index_t> face_blk_ids;

    auto & face_vertices_ref = entities_[2][0];

    // get the face block ids
    if (int64)
      face_blk_ids = base_t::template read_block_ids<long long>(
          hdf5_id, EX_FACE_BLOCK, num_face_blk);
    else
      face_blk_ids = base_t::template read_block_ids<int>(
          hdf5_id, EX_FACE_BLOCK, num_face_blk);

    // read each block
    for (int iblk = 0; iblk < num_face_blk; iblk++) {
      // first get the vertices
      if (int64)
        base_t::template read_face_block<long long>(
            hdf5_id, face_blk_ids[iblk], face_vertices_ref);
      else
        base_t::template read_face_block<int>(
            hdf5_id, face_blk_ids[iblk], face_vertices_ref);
      // rotate the vertices so the lowest one is first
      for (auto & vs : face_vertices_ref) {
        auto lowest = std::min_element(vs.begin(), vs.end());
        std::rotate(vs.begin(), lowest, vs.end());
      }
    }

    //--------------------------------------------------------------------------
    // element blocks

    auto num_elem_blk = exo_params.num_elem_blk;
    vector<index_t> elem_blk_ids;

    auto & cell_faces_ref = entities_[3][2];
    auto & cell_vertices_ref = entities_[3][0];

    // temprary storage for matching faces
    auto face_vertices_sorted = std::make_unique<connectivity_t>();

    // get the element block ids
    if (int64)
      elem_blk_ids = base_t::template read_block_ids<long long>(
          hdf5_id, EX_ELEM_BLOCK, num_elem_blk);
    else
      elem_blk_ids = base_t::template read_block_ids<int>(
          hdf5_id, EX_ELEM_BLOCK, num_elem_blk);

    // read each block
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {

      // read the block info, the actual type may change based on
      // the type of the block
      connectivity_t results;
      typename base_t::block_t block_type;
      if (int64)
        block_type = base_t::template read_element_block<long long>(
            hdf5_id, elem_blk_ids[iblk], results);
      else
        block_type = base_t::template read_element_block<int>(
            hdf5_id, elem_blk_ids[iblk], results);

      //--------------------------------
      // make sure the block type isnt unknown
      if (block_type == base_t::block_t::unknown) {

        clog_fatal("Unknown block type");

      }

      //--------------------------------
      // if the block type is a polyhedra, we get a list of faces, so these
      // need to be transcribed into vertices
      else if (block_type == base_t::block_t::polyhedron) {

        // insert the faces into the cell face list
        cell_faces_ref.insert(
            cell_faces_ref.end(), results.begin(), results.end());

        // collect the cell vertices
        detail::intersect(
            results, // i.e. cell_faces for this block
            face_vertices_ref, cell_vertices_ref);

      }

      //--------------------------------
      // if the block type is fixed element size, we need to build the faces
      else {

        // insert the faces into the cell vertex list
        cell_vertices_ref.insert(
            cell_vertices_ref.end(), results.begin(), results.end());

        // build the connecitivity array
        detail::build_connectivity(
            results, // i.e. cell_vertices for each block
            cell_faces_ref, face_vertices_ref, *face_vertices_sorted,
            [=](const auto & vs, auto & face_vs) {
              using connectivity_type = std::decay_t<decltype(face_vs)>;
              using list_type = std::decay_t<decltype(*face_vs.begin())>;
              // hardcoded for hex
              if (block_type == base_t::block_t::hex) {
                face_vs.emplace_back(list_type{vs[0], vs[1], vs[5], vs[4]});
                face_vs.emplace_back(list_type{vs[1], vs[2], vs[6], vs[5]});
                face_vs.emplace_back(list_type{vs[2], vs[3], vs[7], vs[6]});
                face_vs.emplace_back(list_type{vs[3], vs[0], vs[4], vs[7]});
                face_vs.emplace_back(list_type{vs[0], vs[3], vs[2], vs[1]});
                face_vs.emplace_back(list_type{vs[4], vs[5], vs[6], vs[7]});
              }
              // this is for a tet
              else if (block_type == base_t::block_t::tet) {
                face_vs.emplace_back(list_type{vs[0], vs[1], vs[3]});
                face_vs.emplace_back(list_type{vs[1], vs[2], vs[3]});
                face_vs.emplace_back(list_type{vs[2], vs[0], vs[3]});
                face_vs.emplace_back(list_type{vs[0], vs[2], vs[1]});
              } // block type
            });

      } // block type

    } // blocks

    // clear temporary storage
    face_vertices_sorted.reset();

    // check some assertions
    clog_assert(
        cell_vertices_ref.size() == exo_params.num_elem,
        "Mismatch in read blocks");

    //--------------------------------------------------------------------------
    // build the edges

    // make storage for the face edges
    auto & face_edges_ref = entities_[2][1];
    auto & edge_vertices_ref = entities_[1][0];

    // temprary storage for matching edges
    auto edge_vertices_sorted = std::make_unique<connectivity_t>();

    // build the connecitivity array
    detail::build_connectivity(
        face_vertices_ref, face_edges_ref, edge_vertices_ref,
        *edge_vertices_sorted, [](const auto & vs, auto & edge_vs) {
          using connectivity_type = std::decay_t<decltype(edge_vs)>;
          using list_type = std::decay_t<decltype(*edge_vs.begin())>;
          for (auto v0 = std::prev(vs.end()), v1 = vs.begin(); v1 != vs.end();
               v0 = v1, ++v1)
            edge_vs.emplace_back(list_type{*v0, *v1});
        });

    // clear temporary storage
    edge_vertices_sorted.reset();

    // make storage for the edge vertices
    auto num_edges = edge_vertices_ref.size();

    // Determine cell edges
    auto & cell_edges_ref = entities_[3][1];
    detail::intersect(cell_faces_ref, face_edges_ref, cell_edges_ref);

    //--------------------------------------------------------------------------
    // Create the remainder of the connectivities

    auto num_vertices = vertices_.size() / dimension();
    auto num_faces = face_vertices_ref.size();

    entities_[0][1].reserve(num_vertices);
    entities_[0][2].reserve(num_vertices);
    entities_[0][3].reserve(num_vertices);
    entities_[1][2].reserve(num_edges);
    entities_[1][3].reserve(num_edges);
    entities_[2][3].reserve(num_faces);

    detail::transpose(entities_[1][0], entities_[0][1]);
    detail::transpose(entities_[2][0], entities_[0][2]);
    detail::transpose(entities_[3][0], entities_[0][3]);
    detail::transpose(entities_[2][1], entities_[1][2]);
    detail::transpose(entities_[3][1], entities_[1][3]);
    detail::transpose(entities_[3][2], entities_[2][3]);

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(hdf5_id);
  }

  //============================================================================
  //! \brief Implementation of hdf5 mesh write for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  template<typename U = int>
  void write(
      const std::string & name,
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
          element_sets = {},
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
          node_sets = {}) const {

    clog(info) << "Writing mesh to: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file

    // open the hdf5 file
    auto hdf5_id = base_t::open(name, std::ios_base::out);

    // write the initialization parameters
    auto exo_params = base_t::make_params();
    auto num_faces = num_entities(dimension() - 1);
    auto num_cells = num_entities(dimension());
    exo_params.num_nodes = num_entities(0);
    exo_params.num_face = num_faces;
    exo_params.num_face_blk = 1;
    exo_params.num_node_sets = node_sets.size();
    exo_params.num_elem_blk = element_sets.size() ? element_sets.size() : 1;

    if (element_sets.size()) {
      for (const auto & set : element_sets)
        exo_params.num_elem += set.second.size();
    } else
      exo_params.num_elem = num_cells;

    base_t::write_params(hdf5_id, exo_params);

    // check the integer type used in the hdf5 file
    auto int64 = base_t::is_int64(hdf5_id);

    //--------------------------------------------------------------------------
    // write coordinates

    base_t::write_point_coords(hdf5_id, vertices_);

    //--------------------------------------------------------------------------
    // face blocks

    const auto & face_vertices = entities_.at(2).at(0);

    auto face_vertices_func = [&](auto f, auto & vert_list) {
      const auto & vs = face_vertices[f];
      vert_list.insert(vert_list.end(), vs.begin(), vs.end());
    };

    if (int64)
      base_t::template write_face_block<long long>(
          hdf5_id, 1, "faces", num_faces, face_vertices_func);
    else
      base_t::template write_face_block<int>(
          hdf5_id, 1, "faces", num_faces, face_vertices_func);

    //--------------------------------------------------------------------------
    // element blocks

    const auto & cell_faces = entities_.at(3).at(2);
    int elem_blk_id = 1;

    // add the whole element block
    if (!element_sets.size()) {

      auto cell_faces_func = [&](auto c, auto & face_list) {
        const auto & fs = cell_faces[c];
        face_list.insert(face_list.end(), fs.begin(), fs.end());
      };

      if (int64)
        base_t::template write_element_block<long long>(
            hdf5_id, elem_blk_id++, "cells", num_cells, cell_faces_func);
      else
        base_t::template write_element_block<int>(
            hdf5_id, elem_blk_id++, "cells", num_cells, cell_faces_func);

    } // element block

    // or add the element sets
    for (const auto & set : element_sets) {

      const auto & set_name = set.first;
      const auto & set_elements = set.second;
      auto num_cells_set = set_elements.size();

      if (num_cells_set == 0)
        continue;

      auto cell_faces_func = [&](auto c, auto & face_list) {
        const auto & fs = cell_faces[set_elements[c]];
        face_list.insert(face_list.end(), fs.begin(), fs.end());
      };

      if (int64)
        base_t::template write_element_block<long long>(
            hdf5_id, elem_blk_id++, set_name, num_cells_set, cell_faces_func);
      else
        base_t::template write_element_block<int>(
            hdf5_id, elem_blk_id++, set_name, num_cells_set, cell_faces_func);

    } // sets

    //--------------------------------------------------------------------------
    // Node sets

    int node_set_id = 1;

    for (const auto & set : node_sets) {

      const auto & set_name = set.first;
      const auto & set_nodes = set.second;
      auto num_nodes_set = set_nodes.size();

      if (num_nodes_set == 0)
        continue;

      if (int64)
        base_t::template write_node_set<long long>(
            hdf5_id, ++node_set_id, set_name, set_nodes);
      else
        base_t::template write_node_set<int>(
            hdf5_id, ++node_set_id, set_name, set_nodes);

    } // sets

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(hdf5_id);
  }

  //============================================================================
  // Required Overrides
  //============================================================================

  /// Return the number of entities of a particular dimension
  /// \param [in] dim  The entity dimension to query.
  size_t num_entities(size_t dim) const override {
    switch (dim) {
      case 0:
        return vertices_.size() / dimension();
      default:
        return entities_.at(dim).at(0).size();
    }
  }

  /// Return the set of vertices of a particular entity.
  /// \param [in] dimension  The entity dimension to query.
  /// \param [in] entity_id  The id of the entity in question.
  const auto & entities(size_t from_dim, size_t to_dim) const {
    return entities_.at(from_dim).at(to_dim);
  } // vertices

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in question.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override {
    return entities_.at(from_dim).at(to_dim).at(from_id);
  } // vertices

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  template<typename POINT_TYPE>
  auto vertex(size_t vertex_id) const {
    auto num_vertices = vertices_.size() / dimension();
    POINT_TYPE p;
    for (int i = 0; i < dimension(); ++i)
      p[i] = vertices_[i * num_vertices + vertex_id];
    return p;
  } // vertex

private:
  //============================================================================
  // Private data
  //============================================================================

  //! \brief storage for element verts
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;
};

} // namespace io
} // namespace flecsi
