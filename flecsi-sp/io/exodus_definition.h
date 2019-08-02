/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

/// \file

// user includes
#include <flecsi-sp/io/detail.h>

#include <flecsi/coloring/dcrs_utils.h>
#include <flecsi/topology/mesh_definition.h>
#include <flecsi/utils/logging.h>
#include <flecsi/utils/logging.h>

#include <ristra/assertions/errors.h>

// thirdparty includes
#include <exodusII.h>

// system includes
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <ios>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace flecsi_sp {
namespace io {

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class exodus_base {

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

  //! An enumeration to keep track of element types
  enum class block_t { tri, quad, polygon, tet, hex, polyhedron, unknown };

  //============================================================================
  //! \brief open the file for reading or writing
  //! \param [in] name  The name of the file to open.
  //! \param [in] mode  The mode to open the file in.
  //! \return The exodus handle for the open file.
  //============================================================================
  static auto open(const std::string & name, std::ios_base::openmode mode) {

#ifdef DEBUG
    // useful for debug
    ex_opts(EX_ABORT | EX_VERBOSE);
#endif

    // size of floating point variables used in app.
    int app_word_size = sizeof(real_t);

    if ((mode & std::ios_base::in) == std::ios_base::in) {

      // size of floating point stored in name.
      int exo_word_size = 0;
      // the version number
      float version;

      // open the file
      auto exo_id = ex_open(
          name.c_str(), EX_READ, &app_word_size, &exo_word_size, &version);
      if (exo_id < 0)
        clog_fatal(
            "Problem opening exodus file, ex_open() returned " << exo_id);

      return exo_id;

    } else if ((mode & std::ios_base::out) == std::ios_base::out) {

      // size of floating point to be stored in file.
      // change to float to save space
      int exo_word_size = sizeof(real_t);

      // determine the file creation mode
      int cmode = (mode & std::ios_base::app) == std::ios_base::app
                      ? EX_NOCLOBBER
                      : EX_CLOBBER;

      // create file
      auto exo_id =
          ex_create(name.c_str(), cmode, &app_word_size, &exo_word_size);
      if (exo_id < 0)
        clog_fatal(
            "Problem writing exodus file, ex_create() returned " << exo_id);

      return exo_id;

    } else {

      clog_fatal("Unknown file mode");
      return -1;
    }
  }

  //============================================================================
  //! \brief close the file once completed reading or writing
  //! \param [in] exo_id  The exodus file id.
  //============================================================================
  static void close(int exo_id) {
    auto status = ex_close(exo_id);
    if (status)
      clog_fatal("Problem closing exodus file, ex_close() returned " << exo_id);
  }

  //============================================================================
  //! \brief check the integer status
  //! \param [in] exo_id  The exodus file id.
  //! \return true if exodus is using 64 bit integers
  //============================================================================
  static auto is_int64(int exo_id) {
    return (ex_int64_status(exo_id) & EX_IDS_INT64_API);
  }

  //============================================================================
  //! \brief Helper function to make and initialize a set of exodus parameters.
  //! \return the exodus parameters
  //============================================================================
  static auto make_params() {
    ex_init_params exopar;
    std::strcpy(exopar.title, "Exodus II output from flecsi.");
    exopar.num_dim = num_dims;
    exopar.num_nodes = 0;
    exopar.num_edge = 0;
    exopar.num_edge_blk = 0;
    exopar.num_face = 0;
    exopar.num_face_blk = 0;
    exopar.num_elem = 0;
    exopar.num_elem_blk = 0;
    exopar.num_node_sets = 0;
    exopar.num_edge_sets = 0;
    exopar.num_face_sets = 0;
    exopar.num_side_sets = 0;
    exopar.num_elem_sets = 0;
    exopar.num_node_maps = 0;
    exopar.num_edge_maps = 0;
    exopar.num_face_maps = 0;
    exopar.num_elem_maps = 0;
    return exopar;
  }

  //============================================================================
  //! \brief read the exodus parameters from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the exodus parameters
  //============================================================================
  static auto read_params(int exo_id) {
    ex_init_params exo_params;

    // get the initialization parameters
    auto status = ex_get_init_ext(exo_id, &exo_params);
    if (status)
      clog_fatal(
          "Problem getting exodus file parameters, ex_get_init_ext() returned "
          << status);

    // verify mesh dimension
    if (num_dims != exo_params.num_dim)
      clog_fatal(
          "Exodus dimension mismatch: Expected dimension ("
          << num_dims << ") /= Exodus dimension (" << exo_params.num_dim
          << ")");

    return exo_params;
  }

  //============================================================================
  //! \brief write the exodus parameters from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] exo_params  The  exodus parameters
  //============================================================================
  static void write_params(int exo_id, const ex_init_params & exo_params) {
    // verify mesh dimension
    if (num_dims != exo_params.num_dim)
      clog_fatal(
          "Exodus dimension mismatch: Expected dimension ("
          << num_dims << ") /= Exodus dimension (" << exo_params.num_dim
          << ")");

    // put the initialization parameters
    auto status = ex_put_init_ext(exo_id, &exo_params);
    if (status)
      clog_fatal(
          "Problem putting exodus file parameters, ex_put_init_ext() returned "
          << status);
  }
    
  //============================================================================
  //! \brief write the global mapping a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<typename U, typename V>
  static void write_node_map(int exo_id, const V & ids) {
  
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(ids.size());
    std::transform( ids.begin(), ids.end(), mapping.begin(),
        [](auto v){ return v+1; } );
    auto status = ex_put_node_num_map(exo_id, mapping.data());
    if (status)
      clog_fatal(
          "Problem putting node number map, ex_put_node_num_map() returned "
          << status);

  }
  
  //============================================================================
  //! \brief read the global mapping from file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<typename U>
  static auto read_node_map(int exo_id, size_t num_nodes) {
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(num_nodes);
    auto status = ex_get_node_num_map(exo_id, mapping.data());
    if (status)
      clog_fatal(
          "Problem getting node number map, ex_get_elem_num_map() returned "
          << status);
    if constexpr (std::is_same_v<ex_index_t, index_t>) {
      for ( auto & v : mapping ) v--;
      return mapping;
    }
    else {
      std::vector<index_t> ids;
      ids.reserve(num_nodes);
      for ( auto v : mapping ) ids.emplace_back(v-1);
      return ids;
    }
  }
  
  //============================================================================
  //! \brief write the global mapping a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<typename U, typename V>
  static void write_element_map(int exo_id, const V & ids) {
  
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(ids.size());
    std::transform( ids.begin(), ids.end(), mapping.begin(),
        [](auto v){ return v+1; } );
    auto status = ex_put_elem_num_map(exo_id, mapping.data());
    if (status)
      clog_fatal(
          "Problem putting node number map, ex_put_elem_num_map() returned "
          << status);

  }

  
  //============================================================================
  //! \brief read the global mapping from file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<typename U>
  static auto read_element_map(int exo_id, size_t num_elem) {
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(num_elem);
    auto status = ex_get_elem_num_map(exo_id, mapping.data());
    if (status)
      clog_fatal(
          "Problem getting node number map, ex_get_elem_num_map() returned "
          << status);
    if constexpr (std::is_same_v<ex_index_t, index_t>) {
      for ( auto & v : mapping ) v--;
      return mapping;
    }
    else {
      std::vector<index_t> ids;
      ids.reserve(num_elem);
      for ( auto v : mapping ) ids.emplace_back(v-1);
      return ids;
    }
  }

  //============================================================================
  //! \brief read the coordinates of the mesh from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the vertex coordinates
  //============================================================================
  static auto read_point_coords(int exo_id, size_t num_nodes) {

    // get the number of nodes
    if (num_nodes <= 0)
      clog_fatal(
          "Exodus file has zero nodes, or parmeters haven't been read yet.");

    // read nodes
    vector<real_t> vertex_coord(num_dims * num_nodes);

    // exodus is kind enough to fetch the data in the real type we ask for
    auto status = ex_get_coord(
        exo_id, vertex_coord.data(), vertex_coord.data() + num_nodes,
        vertex_coord.data() + 2 * num_nodes);

    if (status)
      clog_fatal(
          "Problem getting vertex coordinates from exodus file, "
          << " ex_get_coord() returned " << status);

    return vertex_coord;
  }

  //============================================================================
  //! \brief write the coordinates of the mesh from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] vertex_coord  the vertex coordinates
  //============================================================================
  template<typename V>
  static void write_point_coords(int exo_id, const V & vertex_coord) {

    if (vertex_coord.empty())
      return;

    auto num_nodes = vertex_coord.size() / num_dims;

    // exodus is kind enough to fetch the data in the real type we ask for
    auto status = ex_put_coord(
        exo_id, vertex_coord.data(), vertex_coord.data() + num_nodes,
        vertex_coord.data() + 2 * num_nodes);

    if (status)
      clog_fatal(
          "Problem putting vertex coordinates to exodus file, "
          << " ex_put_coord() returned " << status);
  }

  //============================================================================
  //! \brief read the block ids from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename V>
  static void write_node_set(
      int exoid,
      size_t node_set_id,
      const std::string & name,
      const vector<V> & vertex_list) {

    // some type aliases
    using ex_index_t = U;

    if (vertex_list.empty())
      return;

    // set the node set parameters
    ex_index_t num_dist_in_set = 0;
    ex_index_t num_nodes_this_set = vertex_list.size();
    auto status = ex_put_node_set_param(
        exoid, node_set_id, num_nodes_this_set, num_dist_in_set);
    if (status)
      clog_fatal(
          "Problem writing node set param to exodus file, "
          << " ex_put_node_set_param() returned " << status);

    // copy the vertex ids
    vector<ex_index_t> node_set;
    node_set.reserve(vertex_list.size());

    for (auto v : vertex_list)
      node_set.push_back(v + 1); // exodus uses 1-based ids

    // write the node set
    status = ex_put_node_set(exoid, node_set_id, node_set.data());
    if (status)
      clog_fatal(
          "Problem writing node set to exodus file, "
          << " ex_put_node_set() returned " << status);

    // write the set name
    status = ex_put_name(exoid, EX_NODE_SET, node_set_id, name.c_str());
    if (status)
      clog_fatal(
          "Problem writing node set name to exodus file, "
          << " ex_put_name() returned " << status);
  }

  //============================================================================
  //! \brief read the block ids from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto
  read_block_ids(int exoid, ex_entity_type obj_type, size_t num_blocks) {
    // some type aliases
    using ex_index_t = U;

    // final storage
    vector<index_t> ids(num_blocks);

    if (num_blocks > 0) {

      // get the ids first
      vector<ex_index_t> block_ids(num_blocks);
      auto status = ex_get_ids(exoid, obj_type, block_ids.data());
      if (status)
        clog_fatal(
            "Problem reading block ids, ex_get_ids() returned " << status);

      // now convert them
      std::transform(
          block_ids.begin(), block_ids.end(), ids.begin(),
          [](auto id) { return id; });
    }

    // now return them
    return ids;
  }

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename ENTITY_CONN>
  static void write_block(
      int exoid,
      size_t blk_id,
      const std::string & name,
      ex_entity_type entity_type,
      const std::string & entity_description,
      size_t num_entities,
      ENTITY_CONN && entity_conn) {

    // some type aliases
    using ex_index_t = U;

    //-------------------------------------------------------------------------
    if ( entity_description == "nfaced" || entity_description == "nsided" )
    {

      // check if face data is provided instead of node data
      auto is_face_data = (num_dims == 3 && entity_type == EX_ELEM_BLOCK);

      // guess how many elements are in each connectivity slot
      auto num_conn_guess = num_dims * num_dims;

      // build the connectivitiy list for the block
      vector<ex_index_t> entity_nodes;
      vector<int> entity_node_counts;
      entity_nodes.reserve(num_entities * num_conn_guess);
      entity_node_counts.reserve(num_entities);

      // temporary storage for each connecitivity slot
      vector<ex_index_t> temp_entity_nodes;
      temp_entity_nodes.reserve(num_conn_guess);

      for (size_t e = 0; e < num_entities; ++e) {
        // get the connectivity from the user-provided function
        std::forward<ENTITY_CONN>(entity_conn)(e, temp_entity_nodes);
        // store them in the master list ( exodus wants 1-index arrays )
        for (auto i : temp_entity_nodes)
          entity_nodes.emplace_back(i + 1);
        entity_node_counts.push_back(temp_entity_nodes.size());
        // reset temp storage
        temp_entity_nodes.clear();
      }

      // the total size needed to hold the element connectivity
      ex_index_t num_nodes_this_blk = entity_nodes.size();
      ex_index_t num_entries_this_blk = entity_node_counts.size();

      // set the block header
      ex_index_t num_attr_per_entry = 0;
      ex_index_t num_nodes_per_entry = is_face_data ? 0 : num_nodes_this_blk;
      ex_index_t num_edges_per_entry = 0;
      ex_index_t num_faces_per_entry = is_face_data ? num_nodes_this_blk : 0;
      auto status = ex_put_block(
          exoid, entity_type, blk_id, entity_description.c_str(), num_entries_this_blk,
          num_nodes_per_entry, num_edges_per_entry, num_faces_per_entry,
          num_attr_per_entry);
      if (status)
        clog_fatal(
            "Problem writing block to exodus file, "
            << " ex_put_block() returned " << status);

      // write the block name
      status = ex_put_name(exoid, entity_type, blk_id, name.c_str());
      if (status)
        clog_fatal(
            "Problem writing block name to exodus file, "
            << " ex_put_name() returned " << status);

      // write connectivity
      auto node_conn = is_face_data ? nullptr : entity_nodes.data();
      auto elem_edge_conn = nullptr;
      auto elem_face_conn = is_face_data ? entity_nodes.data() : nullptr;
      status = ex_put_conn(
          exoid, entity_type, blk_id, node_conn, elem_edge_conn, elem_face_conn);
      if (status)
        clog_fatal(
            "Problem writing block connectivity to exodus file, "
            << " ex_put_conn() returned " << status);

      // write counts
      status = ex_put_entity_count_per_polyhedra(
          exoid, entity_type, blk_id, entity_node_counts.data());
      if (status)
        clog_fatal(
            "Problem writing block counts to exodus file, "
            << " ex_put_entity_count_per_polyhedra() returned " << status);

    }

    //-------------------------------------------------------------------------
    else {
      
      vector<ex_index_t> entity_nodes, temp_entity_nodes;
      
      for (size_t e = 0; e < num_entities; ++e) {
        temp_entity_nodes.clear();
        // get the connectivity from the user-provided function
        std::forward<ENTITY_CONN>(entity_conn)(e, temp_entity_nodes);
        // store them in the master list ( exodus wants 1-index arrays )
        for (auto i : temp_entity_nodes)
          entity_nodes.emplace_back(i + 1);
      }

      // write the block header
      ex_index_t num_entries_this_blk = num_entities;
      ex_index_t num_nodes_per_entry = temp_entity_nodes.size();
      ex_index_t num_edges_per_entry = 0;
      ex_index_t num_faces_per_entry = 0;
      ex_index_t num_attr_per_entry = 0;
      
      auto status = ex_put_block(
          exoid, entity_type, blk_id, entity_description.c_str(),
          num_entries_this_blk, num_nodes_per_entry, num_edges_per_entry,
          num_faces_per_entry, num_attr_per_entry);
      if (status)
        clog_fatal(
            "Problem writing block to exodus file, "
            << " ex_put_block() returned " << status);

      // write the block name
      status = ex_put_name(exoid, entity_type, blk_id, name.c_str());
      if (status)
        clog_fatal(
            "Problem writing block name to exodus file, "
            << " ex_put_name() returned " << status);

      // write connectivity
      auto node_conn = entity_nodes.data();
      auto elem_edge_conn = nullptr;
      auto elem_face_conn = nullptr;
      status = ex_put_conn(
          exoid, entity_type, blk_id, node_conn, elem_edge_conn, elem_face_conn);
      if (status)
        clog_fatal(
            "Problem writing block connectivity to exodus file, "
            << " ex_put_conn() returned " << status);
    } // type

  }; // write block

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto read_block(
      int exoid,
      ex_entity_id blk_id,
      ex_entity_type entity_type,
      vector<int> & counts,
      vector<U> & indices )
  {

    // some type aliases
    using ex_index_t = U;

    // get the info about this block
    ex_index_t num_elem_this_blk = 0;
    ex_index_t num_faces_per_elem = 0;
    ex_index_t num_edges_per_elem = 0;
    ex_index_t num_nodes_per_elem = 0;
    ex_index_t num_attr = 0;
    char elem_type[MAX_STR_LENGTH];
    auto status = ex_get_block(
        exoid, entity_type, blk_id, elem_type, &num_elem_this_blk,
        &num_nodes_per_elem, &num_edges_per_elem, &num_faces_per_elem,
        &num_attr);
    if (status)
      clog_fatal("Problem reading block, ex_get_block() returned " << status);

    //------------------------------------------------------------------------
    // polygon data
    if (strcasecmp("nsided", elem_type) == 0) {

      // the number of nodes per element is really the number of nodes
      // in the whole block
      auto num_nodes_this_blk = num_nodes_per_elem;

      // get the number of nodes per element
      counts.resize(num_elem_this_blk);
      status = ex_get_entity_count_per_polyhedra(
          exoid, entity_type, blk_id, counts.data());
      if (status)
        clog_fatal(
            "Problem getting element node numbers, "
            << "ex_get_entity_count_per_polyhedra() returned " << status);

      // read element definitions
      indices.resize(num_nodes_this_blk);
      status = ex_get_conn(
          exoid, entity_type, blk_id, indices.data(), nullptr, nullptr);
      if (status)
        clog_fatal(
            "Problem getting element connectivity, ex_get_elem_conn() "
            << "returned " << status);
      
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
      counts.resize(num_elem_this_blk);
      status = ex_get_entity_count_per_polyhedra(
          exoid, entity_type, blk_id, counts.data());
      if (status)
        clog_fatal(
            "Problem reading element node info, "
            << "ex_get_entity_count_per_polyhedra() returned " << status);

      // read element definitions
      indices.resize(num_face_this_blk);
      status = ex_get_conn(
          exoid, entity_type, blk_id, nullptr, nullptr, indices.data());
      if (status)
        clog_fatal(
            "Problem getting element connectivity, ex_get_conn() "
            << "returned " << status);

      return block_t::polyhedron;

    }
    //------------------------------------------------------------------------
    // fixed element size
    else {

      // set the counts
      counts.resize(num_elem_this_blk);
      std::fill( counts.begin(), counts.end(), num_nodes_per_elem );

      // read element definitions
      indices.resize(num_elem_this_blk * num_nodes_per_elem);
      status = ex_get_elem_conn(exoid, blk_id, indices.data());
      if (status)
        clog_fatal(
            "Problem getting element connectivity, ex_get_elem_conn() "
            << "returned " << status);

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
      else if (strcasecmp("hex8", elem_type) == 0 ||
               strcasecmp("hex", elem_type) == 0)
        return block_t::hex;
      else {
        clog_fatal("Unknown block type, " << elem_type);
        return block_t::unknown;
      }

    } // element type
    //------------------------------------------------------------------------
  }

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto
  read_face_block(
    int exoid,
    ex_entity_id blk_id, 
    vector<int> & counts,
    vector<U> & indices ) {
    return read_block<U>(exoid, blk_id, EX_FACE_BLOCK, counts, indices);
  }

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename CONN_TYPE>
  static void write_face_block(
      int exoid,
      size_t blk_id,
      const std::string & name,
      size_t num_faces,
      CONN_TYPE && face_conn) {

    write_block<U>(
        exoid, blk_id, name, EX_FACE_BLOCK, "nsided", num_faces,
        std::forward<CONN_TYPE>(face_conn));
  }

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  static auto read_element_block(
      int exoid,
      ex_entity_id elem_blk_id,
      vector<int> & counts,
      vector<U> & indices )
  {
    return read_block<U>(exoid, elem_blk_id, EX_ELEM_BLOCK, counts, indices);
  }

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename CONN_TYPE>
  static void write_element_block(
      int exoid,
      size_t blk_id,
      const std::string & name,
      const std::string & entity_desc,
      size_t num_elems,
      CONN_TYPE && element_conn) {

    write_block<U>(
        exoid, blk_id, name, EX_ELEM_BLOCK, entity_desc, num_elems,
        std::forward<CONN_TYPE>(element_conn));
  }
};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the multi-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class exodus_definition {};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the one-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class exodus_definition<1, T> : public flecsi::topology::mesh_definition_u<1>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = exodus_base<1, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition_u<1>;

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
  exodus_definition() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  exodus_definition(const std::string & filename) {
    read(filename);
  }

  /// Copy constructor (disabled)
  exodus_definition(const exodus_definition &) = delete;

  /// Assignment operator (disabled)
  exodus_definition & operator=(const exodus_definition &) = delete;

  /// Destructor
  ~exodus_definition() = default;

  //============================================================================
  //! \brief Implementation of exodus mesh read for burton specialization.
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

    // open the exodus file
    auto exoid = base_t::open(name, std::ios_base::in);
    if (exoid < 0)
      clog_fatal("Problem reading exodus file");

    // get the initialization parameters
    auto exo_params = base_t::read_params(exoid);

    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    //--------------------------------------------------------------------------
    // read coordinates

    vertices_ = base_t::read_point_coords(exoid, exo_params.num_nodes);
    clog_assert(
        vertices_.size() == dimension() * exo_params.num_nodes,
        "Mismatch in read vertices");

    auto num_vertices = vertices_.size() / dimension();

    //--------------------------------------------------------------------------
    // element blocks

    auto num_elem_blk = exo_params.num_elem_blk;
    vector<index_t> elem_blk_ids;

    auto & cell_vertices_ref = entities_[1][0];

    // get the element block ids
    if (int64)
      elem_blk_ids = base_t::template read_block_ids<long long>(
          exoid, EX_ELEM_BLOCK, num_elem_blk);
    else
      elem_blk_ids = base_t::template read_block_ids<int>(
          exoid, EX_ELEM_BLOCK, num_elem_blk);

    // read each block
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {
      if (int64)
        base_t::template read_element_block<long long>(
            exoid, elem_blk_ids[iblk], cell_vertices_ref);
      else
        base_t::template read_element_block<int>(
            exoid, elem_blk_ids[iblk], cell_vertices_ref);
    }

    // check some assertions
    clog_assert(
        cell_vertices_ref.size() == exo_params.num_elem,
        "Mismatch in read blocks");

    //--------------------------------------------------------------------------
    // Create the remainder of the connectivities

    entities_[0][1].reserve(num_vertices);

    detail::transpose(entities_[1][0], entities_[0][1]);

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(exoid);
  }

  //============================================================================
  //! \brief Implementation of exodus mesh write for burton specialization.
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

    // open the exodus file
    auto exoid = base_t::open(name, std::ios_base::out);

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

    base_t::write_params(exoid, exo_params);

    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    //--------------------------------------------------------------------------
    // write coordinates

    base_t::write_point_coords(exoid, vertices_);

    //--------------------------------------------------------------------------
    // element blocks

    const auto & cell_vertices = entities_.at(1).at(0);
    int elem_blk_id = 1;

    // add the whole element block
    if (!element_sets.size()) {

      auto cell_vertices_func = [&](auto c, auto & vert_list) {
        const auto & vs = cell_vertices[c];
        vert_list.insert(vert_list.end(), vs.begin(), vs.end());
      };

      if (int64)
        base_t::template write_element_block<long long>(
            exoid, elem_blk_id++, "cells", num_cells, cell_vertices_func);
      else
        base_t::template write_element_block<int>(
            exoid, elem_blk_id++, "cells", num_cells, cell_vertices_func);

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
            exoid, elem_blk_id++, set_name, num_cells_set, cell_vertices_func);
      else
        base_t::template write_element_block<int>(
            exoid, elem_blk_id++, set_name, num_cells_set, cell_vertices_func);

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
            exoid, ++node_set_id, set_name, set_nodes);
      else
        base_t::template write_node_set<int>(
            exoid, ++node_set_id, set_name, set_nodes);

    } // sets

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(exoid);
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
  const std::vector<std::vector<size_t>> & entities(size_t from_dim, size_t to_dim) const override {
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

  void create_graph(
      size_t from_dimension,
      size_t to_dimension,
      size_t min_connections,
      flecsi::coloring::dcrs_t & dcrs ) const override
  {
    constexpr auto num_dims = dimension();

    if ( from_dimension != num_dims || to_dimension != 0 )
      THROW_RUNTIME_ERROR( "Incorrect dimensions provided to create_graph" );

    flecsi::coloring::make_dcrs_distributed<num_dims>(
      *this, from_dimension, to_dimension, min_connections, dcrs);
  }

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
/// \brief This is the two-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class exodus_definition<2, T> : public flecsi::topology::mesh_definition_u<2, T>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = exodus_base<2, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition_u<2>;

  //! the number of dimensions
  using mesh_definition_t::dimension;

  //! the floating point type
  using real_t = typename base_t::real_t;

  //! the index type
  using index_t = typename base_t::index_t;

  using typename flecsi::topology::mesh_definition_u<2, T>::byte_t;

  //! the vector type
  template<typename U>
  using vector = typename base_t::template vector<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;
  using crs_t = flecsi::coloring::crs_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Default constructor
  exodus_definition() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  exodus_definition(const std::string & filename, bool serial=true) {
    read(filename, serial);
  }

  /// Copy constructor (disabled)
  exodus_definition(const exodus_definition &) = delete;

  /// Assignment operator (disabled)
  exodus_definition & operator=(const exodus_definition &) = delete;

  /// Destructor
  ~exodus_definition() = default;

  //============================================================================
  //! \brief Implementation of exodus mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void read(const std::string & name, bool serial) {

    clog(info) << "Reading mesh from: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file

    // open the exodus file
    auto exoid = base_t::open(name, std::ios_base::in);
    if (exoid < 0)
      clog_fatal("Problem reading exodus file");

    // get the initialization parameters
    auto exo_params = base_t::read_params(exoid);

    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    //--------------------------------------------------------------------------
    // Figure out partitioning 
    
    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    // figure out the number of cells for this rank.
    // if this is reading a bunch of exodus files, there is
    // nothing to do
    size_t num_cells = exo_params.num_elem;
    size_t cell_min{0}, cell_max{num_cells-1};

    if (serial) {
      vector<index_t> cell_partitioning;
      flecsi::coloring::subdivide( num_cells, comm_size, cell_partitioning );

      cell_min = cell_partitioning[comm_rank];
      cell_max = cell_partitioning[comm_rank+1] - 1;
      num_cells = cell_max - cell_min + 1;
    }

    constexpr auto num_dims = base_t::num_dims;
    auto & cell2vertices_ = local_connectivity_[num_dims][0];
    cell2vertices_.clear();

    //--------------------------------------------------------------------------
    // element blocks
    
    // offsets always have initial zero
    cell2vertices_.offsets.emplace_back( 0 );
    size_t cell_counter{0};

    // the number of blocks and some storage for block ids
    auto num_elem_blk = exo_params.num_elem_blk;
    vector<index_t> elem_blk_ids;

    // get the element block ids
    if (int64)
      elem_blk_ids = base_t::template read_block_ids<long long>(
          exoid, EX_ELEM_BLOCK, num_elem_blk);
    else
      elem_blk_ids = base_t::template read_block_ids<int>(
          exoid, EX_ELEM_BLOCK, num_elem_blk);

    // read each block
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {

      typename base_t::block_t block_type;
      auto cell_start = cell2vertices_.size();

      vector<int> counts;
      if (int64) {
        vector<long long> indices;
        block_type = base_t::template read_element_block<long long>(
            exoid, elem_blk_ids[iblk], counts, indices );
        detail::filter_block( counts, indices, cell_min, cell_max, cell_counter,
            cell2vertices_.offsets, cell2vertices_.indices );
      }
      else {
        vector<int> indices;
        block_type = base_t::template read_element_block<int>(
            exoid, elem_blk_ids[iblk], counts, indices );
        detail::filter_block( counts, indices, cell_min, cell_max, cell_counter,
            cell2vertices_.offsets, cell2vertices_.indices );
      }
      
      // add the block type
      auto num_added = cell2vertices_.size() - cell_start;
      for (size_t i=0; i<num_added; ++i) {
        cell_block_id_.emplace_back(elem_blk_ids[iblk]);
        cell_type_.emplace_back(block_type);
      }

    }

    // check some assertions
    clog_assert(
        cell2vertices_.size() == num_cells,
        "Mismatch in read blocks");
    
    // create the local to global cell mapping
    auto & cell_local2global_ = local_to_global_[num_dims];
    
    // for serial files, make up the indices
    if ( serial ) {
      cell_local2global_.clear();
      cell_local2global_.resize(num_cells);
      for ( size_t i=0; i<num_cells; ++i ) {
        cell_local2global_[i] = cell_min + i;
      }
    }
    // read them for parallel files
    else if (int64) {
      cell_local2global_ =
        base_t::template read_element_map<long long>(exoid, num_cells);
    }
    else {
      cell_local2global_ =
        base_t::template read_element_map<int>(exoid, num_cells);
    }
    
    // invert the global id to local id map
    auto & cell_global2local_ = global_to_local_[num_dims];
    for ( size_t i=0; i<num_cells; ++i ) {
      auto global_id = cell_local2global_[i];
      cell_global2local_.emplace( std::make_pair( global_id, i ) );
    }
    
    //--------------------------------------------------------------------------
    // read coordinates

    auto coordinates = base_t::read_point_coords(exoid, exo_params.num_nodes);
    clog_assert(
        coordinates.size() == dimension() * exo_params.num_nodes,
        "Mismatch in read vertices");

    auto num_vertices = coordinates.size() / dimension();
    
    //--------------------------------------------------------------------------
    // Renumber vertices

    // get the vertex maps
    auto & vertex_global2local_ = global_to_local_[0];
    auto & vertex_local2global_ = local_to_global_[0];

    if ( serial ) {

      // Throw away vertices that are not mine
      size_t local_vertices{0};
      vertex_global2local_.clear();

      for ( size_t i=0; i<num_cells; ++i ) {
        for ( auto j=cell2vertices_.offsets[i]; j<cell2vertices_.offsets[i+1]; ++j ) {
          auto global_id = cell2vertices_.indices[j];
          auto res = vertex_global2local_.emplace(
              std::make_pair(global_id, local_vertices) );
          if ( res.second ) local_vertices++;
        }
      }

      // invert the global id to local id map
      vertex_local2global_.clear();
      vertex_local2global_.resize(local_vertices);

      for ( const auto & global_local : vertex_global2local_ ) {
        vertex_local2global_[global_local.second] = global_local.first;
      }
      
      // convert element conectivity to local ids and extract only coordinates
      // that are needed
      vertices_.clear();
      vertices_.resize(local_vertices*num_dims);

      for ( auto & v : cell2vertices_.indices ) {
        auto global_id = v;
        v = vertex_global2local_.at(global_id);
        for ( int d=0; d<num_dims; ++d ) {
          vertices_[ v*num_dims + d ] = 
            coordinates[d * num_vertices + global_id];
        }
      }
    } // serial
    // parallel
    else {
    
      if (int64) {
        vertex_local2global_ =
          base_t::template read_node_map<long long>(exoid, num_vertices);
      }
      else {
        vertex_local2global_ =
          base_t::template read_node_map<int>(exoid, num_vertices);
      }
    
      vertices_.clear();
      vertices_.resize(num_vertices*num_dims);

      for ( size_t i=0; i<num_vertices; ++i ) {
        auto global_id = vertex_local2global_[i];
        vertex_global2local_.emplace( std::make_pair( global_id, i ) );
        for (int d=0; d<num_dims; ++d)
          vertices_[i*num_dims + d] = coordinates[d*num_vertices + i];
      }

    } // parallel

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(exoid);

  }

  //============================================================================
  //! \brief Implementation of exodus mesh write for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void write( const std::string & name ) const {
    
    // some mesh params
    constexpr auto num_dims = dimension();
    auto num_verts = num_entities(0);
    auto num_cells = num_entities(num_dims);

    //--------------------------------------------------------------------------
    // Open file

    // open the exodus file
    auto exoid = base_t::open(name, std::ios_base::out);
    
    // figure out the unique block ids/types
    std::map<index_t, typename base_t::block_t> blocks;
    for ( size_t i=0; i<num_cells; ++i )
      blocks[ cell_block_id_[i] ] = cell_type_[i];


    // write the initialization parameters
    auto exo_params = base_t::make_params();
    exo_params.num_nodes = num_verts;
    exo_params.num_node_sets = 0;
    exo_params.num_elem_blk = blocks.size();
    exo_params.num_elem = num_cells;

    base_t::write_params(exoid, exo_params);

    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    //--------------------------------------------------------------------------
    // write coordinates
    std::vector<real_t> coordinates(vertices_.size());
    for ( int d=0; d<num_dims; ++d ){
      for ( size_t i=0; i<num_verts; ++i ) {
        coordinates[i + d*num_verts] = vertices_[d + i*num_dims];
      }
    }
    base_t::write_point_coords(exoid, coordinates);

    if (int64) { 
      base_t::template write_node_map<long long>(exoid, local_to_global_.at(0));
    }
    else {
      base_t::template write_node_map<int>(exoid, local_to_global_.at(0));
    }

    //--------------------------------------------------------------------------
    // element blocks

    const auto & cell_vertices = local_connectivity_.at(2).at(0);
    const auto & cells_local2global = local_to_global_.at(num_dims);
    
    // need to keep track of global ids, since the order might change
    std::vector<index_t> cell_global_ids;
    cell_global_ids.reserve(num_cells);
    
    for ( auto [bid, btype] : blocks ) { 
      
      // figure the block string
      std::string type_str;
      switch (btype) {
      case (base_t::block_t::quad):
        type_str = "quad4";
        break;
      case (base_t::block_t::tri):
        type_str = "tri3";
        break;
      default:
        type_str = "nsided";
      };
      
      // collect all cells with this block type
      std::vector<index_t> cells_this_blk;
      for ( size_t i=0; i<num_cells; ++i )
        if (cell_block_id_[i] == bid) {
          cells_this_blk.emplace_back(i);
          cell_global_ids.emplace_back( cells_local2global[i] );
        }
      auto num_cells_this_blk = cells_this_blk.size();

      // add the whole element block
      auto cell_vertices_func = [&](auto c, auto & vert_list) {
        auto local_id = cells_this_blk[c];
        const auto & vs = cell_vertices.at(local_id);
        vert_list.insert(vert_list.end(), vs.begin(), vs.end());
      };

      if (int64) {
        base_t::template write_element_block<long long>(
            exoid, bid, "cells", type_str, num_cells_this_blk, cell_vertices_func);
      }
      else {
        base_t::template write_element_block<int>(
            exoid, bid, "cells", type_str, num_cells_this_blk, cell_vertices_func);
      }

    }
    
    // write final element mapping
    if (int64) {
      base_t::template write_element_map<long long>(exoid, cell_global_ids);
    }
    else {
      base_t::template write_element_map<int>(exoid, cell_global_ids);
    }

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(exoid);
  }


  void build_connectivity() override {
    
    //--------------------------------------------------------------------------
    // build the edges

    // reference storage for the cell edges and edge vertices
    const auto & cells2vertices = local_connectivity_.at(2).at(0);
    auto & cells2edges = local_connectivity_[2][1];
    auto & edges2vertices = local_connectivity_[1][0];

    // build the connecitivity array
    std::map<std::vector<index_t>, index_t> sorted_vertices_to_edges;
    detail::new_build_connectivity(
        cells2vertices, cells2edges, edges2vertices,
        sorted_vertices_to_edges,
        [](const auto & vs, auto & edge_vs) {
          for (
            auto v0 = std::prev(vs.end()), v1 = vs.begin();
            v1 != vs.end();
            v0 = v1, ++v1
          )
            edge_vs.push_back({*v0, *v1});
        });


    // now figure out global ids
    auto & edge_local2global = local_to_global_[1];
    auto & edge_global2local = global_to_local_[1];
    flecsi::coloring::match_ids( *this, 1, edge_local2global, edge_global2local );
    
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
        return local_connectivity_.at(dim).at(0).size();
      default:
        clog_fatal(
            "Dimension out of range: 0 < " << dim << " </ " << dimension());
        return 0;
    }
  }

  const std::vector<size_t> &
  local_to_global(size_t dim) const override
  {
    return local_to_global_.at(dim);
  }
  
  const std::map<size_t, size_t> &
  global_to_local(size_t dim) const override
  {
    return global_to_local_.at(dim);
  }


  const crs_t &
  entities_crs(size_t from_dim, size_t to_dim) const override {
    return local_connectivity_.at(from_dim).at(to_dim);
  }

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in question.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override
  {
    const auto & connectivity = local_connectivity_.at(from_dim).at(to_dim);
    auto start = connectivity.offsets[from_id];
    auto end = connectivity.offsets[from_id+1];
    std::vector<size_t> result( &connectivity.indices[start],
        &connectivity.indices[end] );
    return result;
  }

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  void vertex(size_t vertex_id, real_t * coords) const override {
    constexpr auto num_dims = dimension();
    for ( int i=0; i<num_dims; ++i ) {
      coords[i] =
        vertices_[ vertex_id*num_dims + i ];
    }
  } // vertex
  
  void create_graph(
      size_t from_dimension,
      size_t to_dimension,
      size_t min_connections,
      flecsi::coloring::dcrs_t & dcrs ) const override
  {
    constexpr auto num_dims = dimension();

    if ( from_dimension != num_dims || to_dimension != 0 )
      THROW_RUNTIME_ERROR( "Incorrect dimensions provided to create_graph" );

    flecsi::coloring::make_dcrs_distributed<num_dims>(
      *this, from_dimension, to_dimension, min_connections, dcrs);
  }

  void erase(
    size_t dimension,
    const vector<size_t> & local_ids ) override
  {
    
    if (local_ids.empty()) return;

    // assume sorted
    assert( std::is_sorted( local_ids.begin(), local_ids.end() )
        && "entries to delete are not sorted" );

    constexpr auto num_dims = base_t::num_dims;

    //--------------------------------------------------------------------------
    // Erase elements

    // erase any vertices that are no longer used.
    auto & cells2verts = local_connectivity_.at(num_dims).at(0);
    std::vector<index_t> check_vertices;
    
    // erase the local mapping
    auto & cell_local2global = local_to_global_.at(num_dims);
    auto & cell_global2local = global_to_local_.at(num_dims);

    auto num_remove = local_ids.size();
    auto num_cells = cell_local2global.size();

    for (auto i = num_remove; i-- > 0; ) {

      // get global and local id
      auto local_id = local_ids[i];
      auto global_id = cell_local2global[local_id];

      // add vertices to check to the list
      auto start = cells2verts.offsets[local_id];
      auto end = cells2verts.offsets[local_id+1];
      for ( auto j=start; j<end; ++j )
        check_vertices.emplace_back( cells2verts.indices[j] );
      

      // shift all local to global ids
      for ( auto j=local_id; j<num_cells; ++j ) {
        cell_local2global[j] = cell_local2global[j+1];
        cell_type_[j] = cell_type_[j+1];
        cell_block_id_[j] = cell_block_id_[j];
      }
    
      // and the global mapping ( need to decrement entries )
      auto count = cell_global2local.erase( global_id );
      assert( count > 0 && "global id not found!" );
      for ( auto & pair : cell_global2local )
        if ( pair.second>local_id )
          pair.second--;
      
      // decrement number of cells
      num_cells--;
    }

    // resize
    cell_local2global.resize(num_cells);
    cell_block_id_.resize(num_cells);
    cell_type_.resize(num_cells);
    
    // erase cells to vertices info
    cells2verts.erase(local_ids);

    //--------------------------------------------------------------------------
    // Determine unused vertices

    // sort vertices to check and remove duplicates
    std::sort( check_vertices.begin(), check_vertices.end() );
    auto last = std::unique( check_vertices.begin(), check_vertices.end() );

    // find vertices that are no longer used
    std::vector<index_t> delete_vertices;
    delete_vertices.reserve(
        std::distance(check_vertices.begin(), check_vertices.end()) );

    for ( auto it = check_vertices.begin(); it != last; ++it ) {
      // mark unused first
      bool used = false;
      // iterate over all cell-vertices, and bounce if its used
      for ( auto i : cells2verts.indices ) {
        if ( i == *it ) {
          used = true;
          break;
        }
      }
      // if vertex is unused, add to list to delete
      if ( !used ) delete_vertices.emplace_back( *it );
    }
    
    //--------------------------------------------------------------------------
    // Delete vertices

    num_remove = delete_vertices.size();
    if ( num_remove ) {
   
      //------------
      // Delete

      // erase the local mapping
      auto num_vertices = vertices_.size() / num_dims;
      auto & vert_local2global = local_to_global_.at(0);
      auto & vert_global2local = global_to_local_.at(0);

      auto vertex_count = num_vertices;
      for (auto i = num_remove; i-- > 0; ) {

        // get global and local id
        auto local_id = delete_vertices[i];
        auto global_id = vert_local2global[local_id];

        // shift all local to global ids, and vertices
        for ( auto j=local_id; j<vertex_count; ++j ) {
          vert_local2global[j] = vert_local2global[j+1];
          for ( int d=0; d<num_dims; ++d )
            vertices_[ j*num_dims + d ] = vertices_[ (j+1)*num_dims + d ];
        }
      
        // and delete the global mapping
        auto count = vert_global2local.erase( global_id );
        assert( count > 0 && "global id not found!" );
        for ( auto & pair : vert_global2local )
          if ( pair.second>local_id )
            pair.second--;
        
        // decrement number of vertices
        vertex_count--;
      }

      // resize
      vert_local2global.resize(vertex_count);
      vertices_.resize(vertex_count*num_dims);
      
      //------------
      // Renumber
      
      // storage for renumberings
      vector<index_t>
        old2new( num_vertices, std::numeric_limits<size_t>::max() );
      
      // add the number of vertices to the oend
      delete_vertices.emplace_back( num_vertices );

      // figure new numbering
      size_t i{0};
      vertex_count = 0;

      for (auto it=delete_vertices.begin(); it != delete_vertices.end(); ++it)
      {
        for ( ; i<*it; ++i ) {
          old2new[i] = vertex_count;
          vertex_count++;
        }
        i++;
      }
 
      assert( vertex_count == num_vertices - num_remove &&
          "Vertex deletion messed up" );

      // renumber cell connectivity
      for ( auto & i : cells2verts.indices ) {
        i = old2new[i];
        assert( i< vertex_count && "Messed up renumbering #1" );
      }

    } // delete vertices
    
  }

  void pack(
    size_t dimension,
    size_t local_id,
    std::vector<byte_t> & buffer ) const override
  {
    // get number of vertices
    constexpr auto num_dims = base_t::num_dims;
    const auto & cells2verts = local_connectivity_.at(num_dims).at(0);
    const auto & start = cells2verts.offsets[local_id];
    const auto & end   = cells2verts.offsets[local_id+1];
    auto num_verts = end - start;
    // get numbering maps
    const auto & cells_local2global = local_to_global_.at(num_dims);
    const auto & verts_local2global = local_to_global_.at(0);
    // add mesh global id to buffer
    flecsi::topology::cast_insert( &cells_local2global[local_id], 1, buffer );
    // add cell block info
    flecsi::topology::cast_insert( &cell_block_id_[local_id], 1, buffer );
    flecsi::topology::cast_insert( &cell_type_[local_id], 1, buffer );
    // add num verts to buffer
    flecsi::topology::cast_insert( &num_verts, 1, buffer );
    // add global vertex indices to buffer
    for (auto i=start; i<end; ++i) {
      auto lid = cells2verts.indices[i];
      auto gid = verts_local2global[lid];
      flecsi::topology::cast_insert( &gid, 1, buffer );
    }
    // also pack coordinates too, not ideal but easiest
    for (auto i=start; i<end; ++i) {
      auto lid = cells2verts.indices[i];
      auto coord = vertices_.data() + lid*num_dims;
      flecsi::topology::cast_insert( coord, num_dims, buffer );
    }
  }

  void unpack(
    size_t dimension,
    size_t local_id,
    byte_t const * & buffer ) override
  {

    // get the numbering and connectivity maps
    constexpr auto num_dims = base_t::num_dims;
    auto & cells2verts = local_connectivity_.at(num_dims).at(0);
    auto & cells_local2global = local_to_global_.at(num_dims);
    auto & cells_global2local = global_to_local_.at(num_dims);
    
    // compute the new local id
    local_id = cells_local2global.size();
    // get global mesh id
    size_t global_id;
    flecsi::topology::uncast( buffer, 1, &global_id );
    // add mapping
    cells_local2global.emplace_back( global_id );
    auto res = cells_global2local.emplace( std::make_pair(global_id, local_id) );
    assert( res.second && "Global id already exists but shouldn't" );
    
    // get cell info
    size_t cell_block_id;
    flecsi::topology::uncast( buffer, 1, &cell_block_id );
    cell_block_id_.emplace_back(cell_block_id);

    typename base_t::block_t cell_type;
    flecsi::topology::uncast( buffer, 1, &cell_type );
    cell_type_.emplace_back(cell_type);
    
    // get the number of vertices
    size_t num_verts;
    flecsi::topology::uncast( buffer, 1, &num_verts );
   
    // retrieve global vertex ids
    auto start = cells2verts.offsets.back();
    auto end = start + num_verts;
    cells2verts.offsets.emplace_back( end );
    cells2verts.indices.resize( end );
    flecsi::topology::uncast( buffer, num_verts, &cells2verts.indices[start] );

    // unpack vertex coordinates
    vector<real_t> coords(num_dims*num_verts);
    flecsi::topology::uncast( buffer, coords.size(), coords.data() );

    // need to convert global vertex ids to local ids
    auto & verts_local2global = local_to_global_.at(0);
    auto & verts_global2local = global_to_local_.at(0);
    for ( size_t i=0; i<num_verts; ++i ) {
      auto iv = start + i;
      auto gid = cells2verts.indices[iv];
      auto lid = verts_global2local.size();
      auto res = verts_global2local.emplace( std::make_pair(gid,lid) );
      // inserted, so add entries to maps
      if ( res.second ) {
        verts_local2global.emplace_back( gid );
        vertices_.reserve( vertices_.size() + num_dims );
        for ( int dim=0; dim<num_dims; ++dim )
          vertices_.emplace_back( coords[i*num_dims + dim] );
      }
      // not inserted, so get local id
      else {
        lid = res.first->second;
      }
      // now remap indices
      cells2verts.indices[iv] = lid;
    }

  }


private:
  //============================================================================
  // Private data
  //============================================================================

  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;

  //! \brief storage for cell to vertex connectivity
  std::map<index_t, std::map<index_t, crs_t>> local_connectivity_;
        
  //! \brief global/local id maps
  std::map<index_t, std::map<index_t, index_t>> global_to_local_;
  std::map<index_t, vector<index_t>> local_to_global_;
  
  //! regions
  std::vector<index_t> cell_block_id_;
  std::vector<typename base_t::block_t> cell_type_;

};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the three-dimensional mesh reader and writer based on the
///        Exodus format.
///
/// io_base_t provides registrations of the exodus file extensions.
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class exodus_definition<3, T> : public flecsi::topology::mesh_definition_u<3>
{

public:
  //============================================================================
  // Typedefs
  //============================================================================

  //! the instantiated base type
  using base_t = exodus_base<3, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition_u<3>;

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
  using crs_t = flecsi::coloring::crs_t;

  //============================================================================
  // Constructors
  //============================================================================

  //! \brief Default constructor
  exodus_definition() = default;

  //! \brief Constructor with filename
  //! \param [in] filename  The name of the file to load
  exodus_definition(const std::string & filename, bool serial=true) {
    read(filename, serial);
  }

  /// Copy constructor (disabled)
  exodus_definition(const exodus_definition &) = delete;

  /// Assignment operator (disabled)
  exodus_definition & operator=(const exodus_definition &) = delete;

  /// Destructor
  ~exodus_definition() = default;

  //============================================================================
  //! \brief Implementation of exodus mesh read for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void read(const std::string & name, bool serial) {

    clog(info) << "Reading mesh from: " << name << std::endl;

    //--------------------------------------------------------------------------
    // Open file

    // open the exodus file
    auto exoid = base_t::open(name, std::ios_base::in);
    if (exoid < 0)
      clog_fatal("Problem reading exodus file");

    // get the initialization parameters
    auto exo_params = base_t::read_params(exoid);

    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    //--------------------------------------------------------------------------
    // Figure out partitioning 
    
    int comm_size, comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    // figure out the number of cells for this rank.
    // if this is reading a bunch of exodus files, there is
    // nothing to do
    size_t num_cells = exo_params.num_elem;
    size_t cell_min{0}, cell_max{num_cells-1};

    if (serial) {
      vector<index_t> cell_partitioning;
      flecsi::coloring::subdivide( num_cells, comm_size, cell_partitioning );

      cell_min = cell_partitioning[comm_rank];
      cell_max = cell_partitioning[comm_rank+1] - 1;
      num_cells = cell_max - cell_min + 1;
    }

    constexpr auto num_dims = base_t::num_dims;
    
    //--------------------------------------------------------------------------
    // read face blocks

    auto num_face_blk = exo_params.num_face_blk;
    vector<index_t> face_blk_ids;

    crs_t temp_face2vertices;

    // get the face block ids
    if (int64)
      face_blk_ids = base_t::template read_block_ids<long long>(
          exoid, EX_FACE_BLOCK, num_face_blk);
    else
      face_blk_ids = base_t::template read_block_ids<int>(
          exoid, EX_FACE_BLOCK, num_face_blk);

    // read each block (add all faces right now.  we will
    // filter out the unused ones later
    for (int iblk = 0; iblk < num_face_blk; iblk++) {
      vector<int> counts;
      // first get the vertices
      if (int64) {
        vector<long long> indices;
        base_t::template read_face_block<long long>(
            exoid, face_blk_ids[iblk], counts, indices);
        for (auto & i : indices) i--;
        for (size_t i=0, j=0; i<counts.size(); ++i) {
          auto n = counts[i];
          temp_face2vertices.append(&indices[j], &indices[j+n]);
          j += n;
        }
      }
      else {
        vector<int> indices;
        base_t::template read_face_block<int>(
            exoid, face_blk_ids[iblk], counts, indices);
        for (auto & i : indices) i--;
        for (size_t i=0, j=0; i<counts.size(); ++i) {
          auto n = counts[i];
          temp_face2vertices.append(&indices[j], &indices[j+n]);
          j += n;
        }
      }
    }
    
    // sort the face vertices and add them to the search map
    std::vector<index_t> sorted_vs;
    for ( auto vs : temp_face2vertices ) {
      sorted_vs.assign( vs.begin(), vs.end() );
      auto face_id = sorted_vertices_to_faces_.size();
      sorted_vertices_to_faces_.emplace( sorted_vs, face_id );
    }
    
    //--------------------------------------------------------------------------
    // read coordinates

    auto coordinates = base_t::read_point_coords(exoid, exo_params.num_nodes);
    clog_assert(
        coordinates.size() == dimension() * exo_params.num_nodes,
        "Mismatch in read vertices");

    auto num_vertices = coordinates.size() / dimension();
    
    
    //--------------------------------------------------------------------------
    // element blocks
    
    // offsets always have initial zero
    auto & cell2vertices_ = local_connectivity_[num_dims][0];
    auto & cell2faces_ = local_connectivity_[num_dims][2];

    cell2vertices_.clear();
    cell2faces_.clear();
    
    size_t cell_counter{0};

    // need to keep track of full polygon connectivitiy
    crs_t temp_cell2faces;
    std::vector<index_t> temp_cellids;

    // the number of blocks and some storage for block ids
    auto num_elem_blk = exo_params.num_elem_blk;
    vector<index_t> elem_blk_ids;

    // get the element block ids
    if (int64)
      elem_blk_ids = base_t::template read_block_ids<long long>(
          exoid, EX_ELEM_BLOCK, num_elem_blk);
    else
      elem_blk_ids = base_t::template read_block_ids<int>(
          exoid, EX_ELEM_BLOCK, num_elem_blk);

    // a lambda to add face info
    auto track_faces = [](auto counter_start, const auto & counts,
        auto & indices, auto & cell2faces, auto & cellids )
    {
      for ( auto & i : indices ) i--;
      for (size_t i=0, j=0; i<counts.size(); ++i)
      {
        auto n = counts[i];
        cell2faces.append( &indices[j], &indices[j+n] );
        cellids.emplace_back(counter_start + i);
        j += n;
      }
    };

    // read each block
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {

      typename base_t::block_t block_type;
      crs_t new_entities;
      auto counter_start = cell_counter;

      if (int64) {
        vector<int> temp_counts;
        vector<long long> temp_indices;
        block_type = base_t::template read_element_block<long long>(
            exoid, elem_blk_ids[iblk], temp_counts, temp_indices );
        detail::filter_block( temp_counts, temp_indices, cell_min, cell_max,
            cell_counter, new_entities.offsets, new_entities.indices );
        if (block_type == base_t::block_t::polyhedron)
          track_faces( counter_start, temp_counts, temp_indices,
              temp_cell2faces, temp_cellids );
      }
      else {
        vector<int> temp_counts, temp_indices;
        block_type = base_t::template read_element_block<int>(
            exoid, elem_blk_ids[iblk], temp_counts, temp_indices );
        detail::filter_block( temp_counts, temp_indices, cell_min, cell_max,
            cell_counter, new_entities.offsets, new_entities.indices );
        if (block_type == base_t::block_t::polyhedron)
          track_faces( counter_start, temp_counts, temp_indices,
              temp_cell2faces, temp_cellids );
      }
      
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
        for ( auto fs : new_entities )
          cell2faces_.append(fs.begin(), fs.end());

        // collect the cell vertices
        detail::new_intersect(new_entities, temp_face2vertices, cell2vertices_);

      }

      //--------------------------------
      // if the block type is fixed element size, we need to build the faces
      else {
        
        // insert the faces into the cell face list
        for ( auto vs : new_entities )
          cell2vertices_.append(vs.begin(), vs.end());

        // build the connecitivity array
        detail::new_build_connectivity(
            new_entities, // i.e. cell_vertices for each block
            cell2faces_, temp_face2vertices,
            sorted_vertices_to_faces_,
            [=](const auto & vs, auto & face_vs) {
              // hardcoded for hex
              if (block_type == base_t::block_t::hex) {
                face_vs.push_back({vs[0], vs[1], vs[5], vs[4]});
                face_vs.push_back({vs[1], vs[2], vs[6], vs[5]});
                face_vs.push_back({vs[2], vs[3], vs[7], vs[6]});
                face_vs.push_back({vs[3], vs[0], vs[4], vs[7]});
                face_vs.push_back({vs[0], vs[3], vs[2], vs[1]});
                face_vs.push_back({vs[4], vs[5], vs[6], vs[7]});
              }
              // this is for a tet
              else if (block_type == base_t::block_t::tet) {
                face_vs.push_back({vs[0], vs[1], vs[3]});
                face_vs.push_back({vs[1], vs[2], vs[3]});
                face_vs.push_back({vs[2], vs[0], vs[3]});
                face_vs.push_back({vs[0], vs[2], vs[1]});
              } // block type
            });

      } // block type


      // add the block type
      for (size_t i=0; i<new_entities.size(); ++i) {
        cell_block_id_.emplace_back(elem_blk_ids[iblk]);
        cell_type_.emplace_back(block_type);
      }

    } // blocks
    
    // check some assertions
    clog_assert(
        cell2vertices_.size() == num_cells,
        "Mismatch in read blocks");

    clog_assert( cell_block_id_.size() == num_cells, "Mismatch in block types" );
    
   
    //--------------------------------------------------------------------------
    // If there were face blocks, i.e. polyhedra, then we need to pick out the
    // actual faces that were used.  We also need to flip any faces whose owner
    // cell was discarded (.i.e. belongs to another proc).  Exodus does not
    // tell you how a cell uses a face, so the conventioned used here is that
    // the cell with the lowest id owns the face.

    auto & face2vertices_ = local_connectivity_[2][0];

    //--------------------------------
    if (num_face_blk > 0) {

      // reset storage
      face2vertices_.clear();
      sorted_vertices_to_faces_.clear(); // dont need anymore

      // figure out which faces are actually indexed
      std::vector<index_t> used_faces;
      for ( auto fs : cell2faces_ )
        used_faces.insert( used_faces.end(), fs.begin(), fs.end() );

      std::sort(used_faces.begin(), used_faces.end());
      auto last = std::unique(used_faces.begin(), used_faces.end());

      // figure out the global owner of all the faces
      std::map<index_t, index_t> temp_face2cell;
      for ( size_t i=0; i<temp_cell2faces.size(); ++i )
        for ( auto f : temp_cell2faces.at(i) ) {
          temp_face2cell.try_emplace(f, temp_cellids[i] );
        }

      // flip any faces whose owner cells were discarded
      for ( auto [fid, cid] : temp_face2cell ) {
        if (cid < cell_min || cid > cell_max ) {
          auto start = temp_face2vertices.offsets[fid];
          auto end = temp_face2vertices.offsets[fid+1];
          std::reverse( &temp_face2vertices.indices[start],
              &temp_face2vertices.indices[end] );
        }
      }

      // now remap the faces
      std::map<size_t, size_t> old2new;
      for ( auto fit=used_faces.begin(); fit != last; ++fit ) {
        auto old_id = *fit;
        auto new_id = face2vertices_.size();
        old2new.emplace(old_id, new_id);
        const auto & vs = temp_face2vertices[old_id];
        face2vertices_.append(vs.begin(), vs.end());
        sorted_vs.assign(vs.begin(), vs.end());
        std::sort( sorted_vs.begin(), sorted_vs.end() );
        sorted_vertices_to_faces_.emplace( sorted_vs, new_id );
      }

      // remap cell to faces
      for ( auto & i : cell2faces_.indices ) i = old2new.at(i);
      
    } // num_face_blk

    //--------------------------------
    // if there wer no face blocks, then the created faces are correct
    else {
      std::swap( face2vertices_, temp_face2vertices );
    }

    //--------------------------------------------------------------------------
    // Figure global cell ids
    
    // create the local to global cell mapping
    auto & cell_local2global_ = local_to_global_[num_dims];
    
    // for serial files, make up the indices
    if ( serial ) {
      cell_local2global_.clear();
      cell_local2global_.resize(num_cells);
      for ( size_t i=0; i<num_cells; ++i ) {
        cell_local2global_[i] = cell_min + i;
      }
    }
    // read them for parallel files
    else {
      
      if (int64) {
        cell_local2global_ =
          base_t::template read_element_map<long long>(exoid, num_cells);
      }
      else {
        cell_local2global_ =
          base_t::template read_element_map<int>(exoid, num_cells);
      }

    }
    
    // invert the global id to local id map
    auto & cell_global2local_ = global_to_local_[num_dims];
    for ( size_t i=0; i<num_cells; ++i ) {
      auto global_id = cell_local2global_[i];
      cell_global2local_.emplace( std::make_pair( global_id, i ) );
    }

    //--------------------------------------------------------------------------
    // Renumber vertices

    // get the vertex maps
    auto & vertex_global2local_ = global_to_local_[0];
    auto & vertex_local2global_ = local_to_global_[0];

    if ( serial ) {

      // Throw away vertices that are not mine
      size_t local_vertices{0};
      vertex_global2local_.clear();

      for ( size_t i=0; i<num_cells; ++i ) {
        for ( auto j=cell2vertices_.offsets[i]; j<cell2vertices_.offsets[i+1]; ++j ) {
          auto global_id = cell2vertices_.indices[j];
          auto res = vertex_global2local_.emplace(
              std::make_pair(global_id, local_vertices) );
          if ( res.second ) local_vertices++;
        }
      }

      // invert the global id to local id map
      auto & vertex_local2global_ = local_to_global_[0];
      vertex_local2global_.clear();
      vertex_local2global_.resize(local_vertices);

      for ( const auto & global_local : vertex_global2local_ ) {
        vertex_local2global_[global_local.second] = global_local.first;
      }

      
      // convert element conectivity to local ids and extract only coordinates
      // that are needed
      vertices_.clear();
      vertices_.resize(local_vertices*num_dims);

      for ( auto & v : cell2vertices_.indices ) {
        auto global_id = v;
        v = vertex_global2local_.at(global_id);
        for ( int d=0; d<num_dims; ++d ) {
          vertices_[ v*num_dims + d ] = 
            coordinates[d * num_vertices + global_id];
        }
      }
    
      for ( auto & v : face2vertices_.indices ) 
        v = vertex_global2local_.at(v);

    } // serial
    // parallel
    else {

      if (int64) {
        vertex_local2global_ =
          base_t::template read_node_map<long long>(exoid, num_vertices);
      }
      else {
        vertex_local2global_ =
          base_t::template read_node_map<int>(exoid, num_vertices);
      }
    
      vertices_.clear();
      vertices_.resize(num_vertices*num_dims);

      for ( size_t i=0; i<num_vertices; ++i ) {
        auto global_id = vertex_local2global_[i];
        vertex_global2local_.emplace( std::make_pair( global_id, i ) );
        for (int d=0; d<num_dims; ++d)
          vertices_[i*num_dims + d] = coordinates[d*num_vertices + i];
      }

    } // parallel

    // determine face owners
    // The convention is the first cell to use the face is the owner.
    auto num_faces = face2vertices_.size();
    constexpr auto max_id = std::numeric_limits<size_t>::max();
    //face_owner_.clear();
    face_owner_.resize( num_faces, max_id );

    for ( size_t c=0; c<num_cells; ++c ){
      for ( auto f : cell2faces_.at(c) ) {
        if ( face_owner_[f] == max_id ) {
          face_owner_[f] = cell_local2global_[c];
        }
      }
    }

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(exoid);
  }
  
  //============================================================================
  //! \brief Implementation of exodus mesh write for burton specialization.
  //!
  //! \param[in] name Read burton mesh \e m from \e name.
  //! \param[out] m Populate burton mesh \e m with contents of \e name.
  //!
  //! \return Exodus error code. 0 on success.
  //============================================================================
  void write( const std::string & name ) const {
    
    //--------------------------------------------------------------------------
    // Open file

    // open the exodus file
    auto exoid = base_t::open(name, std::ios_base::out);
    
    // check the integer type used in the exodus file
    auto int64 = base_t::is_int64(exoid);

    
    // some mesh params
    constexpr auto num_dims = dimension();
    auto num_verts = num_entities(0);
    auto num_cells = num_entities(num_dims);
    

    //--------------------------------------------------------------------------
    // Check blocks

    // figure out the unique block ids/types
    size_t max_block_id{0};
    std::map<index_t, typename base_t::block_t> blocks;
    for ( size_t i=0; i<num_cells; ++i ) {
      max_block_id = std::max( max_block_id, cell_block_id_[i] );
      blocks[ cell_block_id_[i] ] = cell_type_[i];
    }
    
    // check if there are polyhedra, they need face blocks
    bool has_polyhedra = false;
    for ( auto [bid, btype] : blocks )
      if ( btype == base_t::block_t::polyhedron ) {
        has_polyhedra = true;
        break;
      }
    
    //--------------------------------------------------------------------------
    // setup the initialization parameters

    auto exo_params = base_t::make_params();
    exo_params.num_nodes = num_verts;
    exo_params.num_node_sets = 0;
    exo_params.num_elem = num_cells;
    exo_params.num_face_blk = has_polyhedra ? 1 : 0;
    exo_params.num_elem_blk = blocks.size();


    //--------------------------------------------------------------------------
    // face blocks
      
    const auto & cell_faces = local_connectivity_.at(3).at(2);
    std::map< index_t, index_t > face_renumbering;

    if ( has_polyhedra ) {
        
      auto face_blk_id = max_block_id + 1;
      std::vector<index_t> faces_this_blk;
    
      const auto & face_vertices = local_connectivity_.at(2).at(0);
      const auto & cell_faces = local_connectivity_.at(3).at(2);
      const auto & cell_local2global = local_to_global_.at(num_dims);
    
      // Figure out first cell to use a face      
      auto num_faces = face_vertices.size();
      std::vector<bool> face_visited( num_faces, false );
      std::vector<bool> face_flipped( num_faces, false );

    for ( auto [bid, btype] : blocks ) { 
      for ( size_t c=0; c<num_cells; ++c ){
        if (cell_block_id_[c] != bid) continue;
        for ( auto f : cell_faces.at(c) ) {
          if ( face_visited[f] ) continue;
          face_visited[f] = true;
          face_flipped[f] = ( face_owner_[f] != cell_local2global[c] );
        }
      }
    }

      // collect all faces
      for ( auto [bid, btype] : blocks ) { 
        if ( btype != base_t::block_t::polyhedron ) continue;
      
        // collect all cells with this block type
        for ( size_t i=0; i<num_cells; ++i )
          if (cell_block_id_[i] == bid) {
            for ( auto f : cell_faces.at(i) ) {
              faces_this_blk.emplace_back( f );
            }
          }
      }

      // sort the faces
      std::sort( faces_this_blk.begin(), faces_this_blk.end() );
      auto last = std::unique( faces_this_blk.begin(), faces_this_blk.end() );

      // write the block
      auto num_faces_this_blk = std::distance(faces_this_blk.begin(), last);
      for ( auto i=0; i<num_faces_this_blk; ++i )
        face_renumbering[ faces_this_blk[i] ] = i;
    
      // add the whole face block
      auto face_vertices_func = [&](auto f, auto & vert_list) {
        auto local_id = faces_this_blk[f];
        const auto & vs = face_vertices.at(local_id);
        // flip if someone else owns
        if ( face_flipped[local_id] ) {
          vert_list.insert(vert_list.end(), vs.rbegin(), vs.rend());
        }
        else {
          vert_list.insert(vert_list.end(), vs.begin(), vs.end());
        }
      };
    
      // write the params before we write a face block
      exo_params.num_face = num_faces_this_blk;
      base_t::write_params(exoid, exo_params);

      if (int64)
        base_t::template write_face_block<long long>(
            exoid, face_blk_id, "faces", num_faces_this_blk, face_vertices_func);
      else
        base_t::template write_face_block<int>(
          exoid, face_blk_id, "faces", num_faces_this_blk, face_vertices_func);

    } // has_polyhedra

    else {

      // just write the params
      base_t::write_params(exoid, exo_params);

    }

    //--------------------------------------------------------------------------
    // write coordinates
    std::vector<real_t> coordinates(vertices_.size());
    for ( int d=0; d<num_dims; ++d ){
      for ( size_t i=0; i<num_verts; ++i ) {
        coordinates[i + d*num_verts] = vertices_[d + i*num_dims];
      }
    }
    base_t::write_point_coords(exoid, coordinates);

    if (int64) { 
      base_t::template write_node_map<long long>(exoid, local_to_global_.at(0));
    }
    else {
      base_t::template write_node_map<int>(exoid, local_to_global_.at(0));
    }
    

    //--------------------------------------------------------------------------
    // element blocks

    const auto & cell_vertices = local_connectivity_.at(3).at(0);
    const auto & cells_local2global = local_to_global_.at(num_dims);

    // need to keep track of global ids, since the order might change
    std::vector<index_t> cell_global_ids;
    cell_global_ids.reserve(num_cells);

    size_t block_counter{1};
    
    for ( auto [bid, btype] : blocks ) { 

      auto this_blk_id = bid;
      
      // collect all cells with this block type
      std::vector<index_t> cells_this_blk;
      for ( size_t i=0; i<num_cells; ++i )
        if (cell_block_id_[i] == bid) {
          cells_this_blk.emplace_back(i);
          cell_global_ids.emplace_back( cells_local2global[i] );
        }
      auto num_cells_this_blk = cells_this_blk.size();
      
      //-----------------------------------------------------------------------
      // polyhedra
      if ( btype == base_t::block_t::polyhedron ) {

        auto cell_faces_func = [&](auto c, auto & face_list) {
          auto local_id = cells_this_blk[c];
          const auto & fs = cell_faces[local_id];
          for ( auto f : fs ) face_list.emplace_back( face_renumbering.at(f) );
        };

        if (int64)
          base_t::template write_element_block<long long>(
              exoid, this_blk_id, "cells", "nfaced", num_cells_this_blk, cell_faces_func);
        else
          base_t::template write_element_block<int>(
            exoid, this_blk_id, "cells", "nfaced", num_cells_this_blk, cell_faces_func);

      }

      //-----------------------------------------------------------------------
      // hex/tet
      else {

        std::string type_str = (btype == base_t::block_t::hex) ? "hex8" : "tet4"; 

        // add the whole element block
        auto cell_vertices_func = [&](auto c, auto & vert_list) {
          auto local_id = cells_this_blk[c];
          const auto & vs = cell_vertices.at(local_id);
          vert_list.insert(vert_list.end(), vs.begin(), vs.end());
        };

        if (int64) {
          base_t::template write_element_block<long long>(
              exoid, this_blk_id, "cells", type_str, num_cells_this_blk, cell_vertices_func);
        }
        else {
          base_t::template write_element_block<int>(
              exoid, this_blk_id, "cells", type_str, num_cells_this_blk, cell_vertices_func);
        }

      } // type

      // bump block counter
      block_counter++;

    }
        
    // write final element mapping
    if (int64) {
      base_t::template write_element_map<long long>(exoid, cell_global_ids);
    }
    else {
      base_t::template write_element_map<int>(exoid, cell_global_ids);
    }

    //--------------------------------------------------------------------------
    // close the file
    base_t::close(exoid);
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
        return local_connectivity_.at(dim).at(0).size();
    }
  }

  const std::vector<size_t> &
  local_to_global(size_t dim) const override
  {
    return local_to_global_.at(dim);
  }
  
  const std::map<size_t, size_t> &
  global_to_local(size_t dim) const override
  {
    return global_to_local_.at(dim);
  }


  const crs_t &
  entities_crs(size_t from_dim, size_t to_dim) const override {
    return local_connectivity_.at(from_dim).at(to_dim);
  }

  /// return the set of vertices of a particular entity.
  /// \param [in] dimension  the entity dimension to query.
  /// \param [in] entity_id  the id of the entity in question.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override
  {
    const auto & connectivity = local_connectivity_.at(from_dim).at(to_dim);
    auto start = connectivity.offsets[from_id];
    auto end = connectivity.offsets[from_id+1];
    std::vector<size_t> result( &connectivity.indices[start],
        &connectivity.indices[end] );
    return result;
  }

  /// Return the vertex coordinates for a certain id.
  /// \param [in] vertex_id  The id of the vertex to query.
  void vertex(size_t vertex_id, real_t * coords) const override {
    constexpr auto num_dims = dimension();
    for ( int i=0; i<num_dims; ++i ) {
      coords[i] =
        vertices_[ vertex_id*num_dims + i ];
    }
  } // vertex
  
  void create_graph(
      size_t from_dimension,
      size_t to_dimension,
      size_t min_connections,
      flecsi::coloring::dcrs_t & dcrs ) const override
  {
    constexpr auto num_dims = dimension();

    if ( from_dimension != num_dims || to_dimension != 0 )
      THROW_RUNTIME_ERROR( "Incorrect dimensions provided to create_graph" );

    flecsi::coloring::make_dcrs_distributed<num_dims>(
      *this, from_dimension, to_dimension, min_connections, dcrs);
  }
  
  void pack(
    size_t dimension,
    size_t local_id,
    std::vector<byte_t> & buffer ) const override
  {
    constexpr auto num_dims = base_t::num_dims;
    // get connectivity maps
    const auto & cells2verts = local_connectivity_.at(num_dims).at(0);
    const auto & cells2faces = local_connectivity_.at(num_dims).at(2);
    const auto & faces2verts = local_connectivity_.at(2).at(0);
    // get numbering maps
    const auto & cells_local2global = local_to_global_.at(num_dims);
    const auto & verts_local2global = local_to_global_.at(0);

    // add mesh global id to buffer
    flecsi::topology::cast_insert( &cells_local2global[local_id], 1, buffer );

    flecsi::topology::cast_insert( &cell_block_id_[local_id], 1, buffer );
    flecsi::topology::cast_insert( &cell_type_[local_id], 1, buffer );

    //--------------------------------------------------------------------------
    // pack vertices
    {
      // get number of vertices
      const auto & vs = cells2verts.at(local_id);
      auto num_verts = vs.size();
      // add num verts to buffer
      flecsi::topology::cast_insert( &num_verts, 1, buffer );
      // add global vertex indices to buffer
      for (auto lid : vs) {
        auto gid = verts_local2global[lid];
        flecsi::topology::cast_insert( &gid, 1, buffer );
      }
      // also pack coordinates too, not ideal but easiest
      for (auto lid : vs) {
        auto coord = vertices_.data() + lid*num_dims;
        flecsi::topology::cast_insert( coord, num_dims, buffer );
      }
    } // verts

    //--------------------------------------------------------------------------
    // now pack faces
    {
      // get number of faces
      const auto & fs = cells2faces.at(local_id);
      auto num_faces = fs.size();
      // add num faces to buffer
      flecsi::topology::cast_insert( &num_faces, 1, buffer );
      // loop over faces
      for ( auto fid : fs ) {
        // get face verts
        const auto & vs = faces2verts.at(fid);
        auto num_verts = vs.size();
        std::vector<size_t> vs_copy( vs.begin(), vs.end() );
        // if this cell is not the owner, then reverse the face verts
        if ( face_owner_[fid] != cells_local2global[local_id] )
          std::reverse( vs_copy.begin(), vs_copy.end() );
        // add number of verts to buffer
        flecsi::topology::cast_insert( &num_verts, 1, buffer );
        // add face verts to buffer
        for ( auto vid : vs_copy ) {
          auto gid = verts_local2global[vid];
          flecsi::topology::cast_insert( &gid, 1, buffer );
        }
      } // faces

    }
  }

  void unpack(
    size_t dimension,
    size_t local_id,
    byte_t const * & buffer ) override
  {

    // get the numbering and connectivity maps
    constexpr auto num_dims = base_t::num_dims;

    auto & cells2verts = local_connectivity_.at(num_dims).at(0);
    auto & cells2faces = local_connectivity_.at(num_dims).at(2);
    auto & faces2verts = local_connectivity_.at(2).at(0);

    auto & cells_local2global = local_to_global_.at(num_dims);
    auto & cells_global2local = global_to_local_.at(num_dims);
    auto & verts_local2global = local_to_global_.at(0);
    auto & verts_global2local = global_to_local_.at(0);
    
    // compute the new local id
    local_id = cells_local2global.size();
    // get global mesh id
    size_t global_id;
    flecsi::topology::uncast( buffer, 1, &global_id );
    // add mapping
    cells_local2global.emplace_back( global_id );
    auto res = cells_global2local.emplace( std::make_pair(global_id, local_id) );
    assert( res.second && "Global id already exists but shouldn't" );

    // get cell info
    size_t cell_block_id;
    flecsi::topology::uncast( buffer, 1, &cell_block_id );
    cell_block_id_.emplace_back(cell_block_id);

    typename base_t::block_t cell_type;
    flecsi::topology::uncast( buffer, 1, &cell_type );
    cell_type_.emplace_back(cell_type);
   
    //--------------------------------------------------------------------------
    // verts
    {
      // get the number of vertices
      size_t num_verts;
      flecsi::topology::uncast( buffer, 1, &num_verts );
   
      // retrieve global vertex ids
      std::vector<size_t> vs(num_verts);
      flecsi::topology::uncast( buffer, num_verts, vs.data() );

      // unpack vertex coordinates
      vector<real_t> coords(num_dims*num_verts);
      flecsi::topology::uncast( buffer, coords.size(), coords.data() );

      // need to convert global vertex ids to local ids
      for ( size_t i=0; i<num_verts; ++i ) {
        auto gid = vs[i];
        auto lid = verts_global2local.size();
        auto res = verts_global2local.emplace( std::make_pair(gid,lid) );
        // inserted, so add entries to maps
        if ( res.second ) {
          verts_local2global.emplace_back( gid );
          vertices_.reserve( vertices_.size() + num_dims );
          for ( int dim=0; dim<num_dims; ++dim )
            vertices_.emplace_back( coords[i*num_dims + dim] );
        }
        // not inserted, so get local id
        else {
          lid = res.first->second;
        }
        // now remap indices
        vs[i] = lid;
      }

      // now append to the connectivity list
      cells2verts.append( vs.begin(), vs.end() );
    } // verts

    //--------------------------------------------------------------------------
    // faces
    {
      // get the number of faces
      size_t num_faces;
      flecsi::topology::uncast( buffer, 1, &num_faces );

      // storage to keep track of new cell faces
      std::vector<size_t> fs;
      fs.reserve( num_faces );
        
      // storage for vertices
      std::vector<size_t> vs, sorted_vs;

      // loop over faces and unpack verts
      for ( size_t i=0; i<num_faces; ++i ) {

        // the number of vertices
        size_t num_verts;
        flecsi::topology::uncast( buffer, 1, &num_verts );

        // the ids of the vertices
        vs.clear(); vs.resize(num_verts);
        flecsi::topology::uncast( buffer, num_verts, vs.data() );
        // convert vertex global ids to local ids
        for ( auto & v : vs ) v = verts_global2local.at(v);

        // sort vertices for search
        sorted_vs.assign(vs.begin(), vs.end());
        std::sort( sorted_vs.begin(), sorted_vs.end() );

        // search for the face to see if it exists
        auto it = sorted_vertices_to_faces_.find( sorted_vs );
        // no match, add it to the connectivity list
        if ( it == sorted_vertices_to_faces_.end() ) {
          auto new_face_id = faces2verts.size();
          fs.emplace_back( new_face_id );
          faces2verts.append( vs.begin(), vs.end() );
          face_owner_.emplace_back( global_id );
          sorted_vertices_to_faces_.emplace( sorted_vs, new_face_id );
          assert( faces2verts.size() == sorted_vertices_to_faces_.size() );
        }
        // if there was a match, use the id
        else {
          fs.emplace_back(it->second);
        }

      }
      
      // now append to the connectivity list
      cells2faces.append( fs.begin(), fs.end() );

    } // faces
    
  }

  void erase(
    size_t dimension,
    const vector<size_t> & local_ids ) override
  {
    
    if (local_ids.empty()) return;

    // assume sorted
    assert( std::is_sorted( local_ids.begin(), local_ids.end() )
        && "entries to delete are not sorted" );

    constexpr auto num_dims = base_t::num_dims;

    //--------------------------------------------------------------------------
    // Erase elements

    // erase any vertices and faces that are no longer used.
    auto & cells2verts = local_connectivity_.at(num_dims).at(0);
    auto & cells2faces = local_connectivity_.at(num_dims).at(2);
    auto & faces2verts = local_connectivity_.at(2).at(0);
    std::vector<index_t> check_vertices, check_faces;
    
    // erase the local mapping
    auto & cell_local2global = local_to_global_.at(num_dims);
    auto & cell_global2local = global_to_local_.at(num_dims);

    auto num_remove = local_ids.size();
    auto num_cells = cell_local2global.size();

    for (auto i = num_remove; i-- > 0; ) {

      // get global and local id
      auto local_id = local_ids[i];
      auto global_id = cell_local2global[local_id];

      // add vertices to check to the list
      for ( auto v : cells2verts.at(local_id) )
        check_vertices.emplace_back(v);
      for ( auto f : cells2faces.at(local_id) )
        check_faces.emplace_back(f);
      

      // shift all local to global ids
      for ( auto j=local_id; j<num_cells; ++j ) {
        cell_local2global[j] = cell_local2global[j+1];
        cell_type_[j] = cell_type_[j+1];
        cell_block_id_[j] = cell_block_id_[j];
      }
    
      // and the global mapping ( need to decrement entries )
      auto count = cell_global2local.erase( global_id );
      assert( count > 0 && "global id not found!" );
      for ( auto & pair : cell_global2local )
        if ( pair.second>local_id )
          pair.second--;
      
      // decrement number of cells
      num_cells--;
    }

    // resize
    cell_local2global.resize(num_cells);
    cell_block_id_.resize(num_cells);
    cell_type_.resize(num_cells);
    
    // erase cells to vertices info
    cells2verts.erase(local_ids);
    cells2faces.erase(local_ids);

    //--------------------------------------------------------------------------
    // Determine unused vertices/faces
    
    auto remove_duplicates = [](auto & list) {
      std::sort( list.begin(), list.end() );
      auto last = std::unique( list.begin(), list.end() );
      list.erase( last, list.end() );
    };

    // sort vertices to check and remove duplicates
    remove_duplicates( check_vertices );
    remove_duplicates( check_faces );
    
    //--------------------------------------------------------------------------
    // Delete faces
    
    // find faces that are no longer used
    std::vector<index_t> delete_items;
    delete_items.reserve( check_faces.size() );
    
    for ( auto f_rm : check_faces ) {
      // mark unused first
      bool used = false;
      // iterate over all cell-vertices, and bounce if its used
      for ( auto f : cells2faces.indices ) {
        if ( f == f_rm ) {
          used = true;
          break;
        }
      }
      // if face is unused, add to list to delete
      if ( !used ) delete_items.emplace_back( f_rm );
    }
    
    // Now actually delete
    if ( !delete_items.empty() ) {
      
      // store the number of faces
      auto num_faces = faces2verts.size();
      num_remove = delete_items.size();
    
      // delete face to vertex connectivity
      faces2verts.erase(delete_items);
      
      // renumber cell to face connectivity
      // storage for renumberings
      constexpr auto max_id = std::numeric_limits<size_t>::max();
      vector<index_t> old2new( num_faces, max_id );
      
      // add the number of faces to the end for the logic below to work
      delete_items.emplace_back( num_faces );
    
      // figure new numbering
      size_t i{0}, face_count{0};

      for (auto it=delete_items.begin(); it != delete_items.end(); ++it)
      {
        for ( ; i<*it; ++i ) {
          old2new[i] = face_count;
          face_count++;
        }
        i++;
      }
      
      assert( face_count == num_faces - num_remove &&
          "Face deletion messed up" );

      // renumber cell connectivity
      for ( auto & i : cells2faces.indices ) {
        i = old2new[i];
        assert( i<face_count && "Messed up renumbering #1" );
      }


      // renumber face owners
      vector<index_t> new_face_owner(face_count);
      for ( size_t old_id=0; old_id<num_faces; ++old_id ) {
        auto new_id = old2new[old_id];
        if ( new_id != max_id )
          new_face_owner[new_id] = face_owner_[old_id];
      }
      std::swap( new_face_owner, face_owner_ );

      // reverse any vertices whose owner was deleted
      for ( size_t f=0; f<face_count; ++f ) {
        // get the old owner
        auto owner_id = face_owner_[f];
        // if it isn't the global list, reverse it
        if ( !cell_global2local.count(owner_id) ) {
          auto start = faces2verts.offsets[f];
          auto end = faces2verts.offsets[f+1];
          std::reverse( &faces2verts.indices[start], &faces2verts.indices[end] );
          // reset the id
          face_owner_[f] = max_id;
        }
      }
   
      // figure out any of the reset face owners
      for ( size_t c=0; c<cells2verts.size(); ++c ){
        for ( auto f : cells2faces.at(c) ) {
          if ( face_owner_[f] == max_id ) {
            face_owner_[f] = cell_local2global[c];
          }
        }
      }

    } // remove faces
    

    //--------------------------------------------------------------------------
    // Delete vertices

    // find vertices that are no longer used
    delete_items.reserve( check_vertices.size() );
    delete_items.clear();

    for ( auto v_rm : check_vertices ) {
      // mark unused first
      bool used = false;
      // iterate over all cell-vertices, and bounce if its used
      for ( auto v : cells2verts.indices ) {
        if ( v == v_rm ) {
          used = true;
          break;
        }
      }
      // if vertex is unused, add to list to delete
      if ( !used ) delete_items.emplace_back( v_rm );
    }
    
    if ( !delete_items.empty() ) {
   
      //------------
      // Delete
    
      num_remove = delete_items.size();

      // erase the local mapping
      auto num_vertices = vertices_.size() / num_dims;
      auto & vert_local2global = local_to_global_.at(0);
      auto & vert_global2local = global_to_local_.at(0);

      auto vertex_count = num_vertices;
      for (auto i = num_remove; i-- > 0; ) {

        // get global and local id
        auto local_id = delete_items[i];
        auto global_id = vert_local2global[local_id];

        // shift all local to global ids, and vertices
        for ( auto j=local_id; j<vertex_count; ++j ) {
          vert_local2global[j] = vert_local2global[j+1];
          for ( int d=0; d<num_dims; ++d )
            vertices_[ j*num_dims + d ] = vertices_[ (j+1)*num_dims + d ];
        }
      
        // and delete the global mapping
        auto count = vert_global2local.erase( global_id );
        assert( count > 0 && "global id not found!" );
        for ( auto & pair : vert_global2local )
          if ( pair.second>local_id )
            pair.second--;
        
        // decrement number of vertices
        vertex_count--;
      }

      // resize
      vert_local2global.resize(vertex_count);
      vertices_.resize(vertex_count*num_dims);
      
      //------------
      // Renumber
      
      // storage for renumberings
      vector<index_t>
        old2new( num_vertices, std::numeric_limits<size_t>::max() );
      
      // add the number of vertices to the end for the logic below to work
      delete_items.emplace_back( num_vertices );

      // figure new numbering
      size_t i{0};
      vertex_count = 0;

      for (auto it=delete_items.begin(); it != delete_items.end(); ++it)
      {
        for ( ; i<*it; ++i ) {
          old2new[i] = vertex_count;
          vertex_count++;
        }
        i++;
      }
 
      assert( vertex_count == num_vertices - num_remove &&
          "Vertex deletion messed up" );

      // renumber cell connectivity
      for ( auto & i : cells2verts.indices ) {
        i = old2new[i];
        assert( i< vertex_count && "Messed up renumbering #1" );
      }
      // renumber face connectivity
      for ( auto & i : faces2verts.indices ) {
        i = old2new[i];
        assert( i< vertex_count && "Messed up renumbering #1" );
      }

    } // delete vertices
      
    // delete the sorted face vertices.  you cant change the key of a map, so
    // you cant renumber the vertices.  just blow away the whole map and
    // rebuild it.
    sorted_vertices_to_faces_.clear();

    size_t face_cnt{0};
    std::vector<size_t> sorted_vs;
    for ( const auto & vs : faces2verts ) {
      sorted_vs.assign( vs.begin(), vs.end() );
      std::sort( sorted_vs.begin(), sorted_vs.end() );
      sorted_vertices_to_faces_.emplace(sorted_vs, face_cnt);
      ++face_cnt;
    }
    
  }

  void build_connectivity() override {

    //--------------------------------------------------------------------------
    // build the edges
      
    // reference storage for the cell edges and edge vertices
    const auto & faces2vertices = local_connectivity_.at(2).at(0);
    auto & faces2edges = local_connectivity_[2][1];
    auto & edges2vertices = local_connectivity_[1][0];

    // build the connecitivity array
    std::map<std::vector<index_t>, index_t> sorted_vertices_to_edges;
    detail::new_build_connectivity(
        faces2vertices, faces2edges, edges2vertices,
        sorted_vertices_to_edges,
        [](const auto & vs, auto & edge_vs) {
          for (
            auto v0 = std::prev(vs.end()), v1 = vs.begin();
            v1 != vs.end();
            v0 = v1, ++v1
          )
            edge_vs.push_back({*v0, *v1});
        });
    
    // now figure out cell to edge connectivity through the faces
    const auto & cells2faces = local_connectivity_.at(3).at(2);
    auto & cells2edges = local_connectivity_[3][1];
    detail::new_intersect(cells2faces, faces2edges, cells2edges);

    // now figure out global ids of edges
    auto & edge_local2global = local_to_global_[1];
    auto & edge_global2local = global_to_local_[1];
    flecsi::coloring::match_ids( *this, 1, edge_local2global, edge_global2local );
    
    // now figure out global ids of faces
    auto & face_local2global = local_to_global_[2];
    auto & face_global2local = global_to_local_[2];
    flecsi::coloring::match_ids( *this, 2, face_local2global, face_global2local );
    
  }

  const std::vector<size_t> & face_owners() const override {
    return face_owner_;
  }

private:
  //============================================================================
  // Private data
  //============================================================================
  
  //! \brief storage for vertex coordinates
  vector<real_t> vertices_;

  //! \brief storage for cell to vertex connectivity
  std::map<index_t, std::map<index_t, crs_t>> local_connectivity_;

  std::vector<index_t> face_owner_;
  std::map<std::vector<index_t>, index_t> sorted_vertices_to_faces_;

  //! \brief global/local id maps
  std::map<index_t, std::map<index_t, index_t>> global_to_local_;
  std::map<index_t, vector<index_t>> local_to_global_;

  //! regions
  std::vector<index_t> cell_block_id_;
  std::vector<typename base_t::block_t> cell_type_;

};

} // namespace io
} // namespace flecsi
