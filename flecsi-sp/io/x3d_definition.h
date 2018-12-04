/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#pragma once

/** \file
 * Definition of the X3D mesh format.
 * The implementation here is based on the LANL report "Summary of the FLAG X3D
 * Format" (LA-UR-04-9033).
 *
 * The X3D format provides for geometry, topology, and field data
 * on multimaterial, distributed, general unstructured polytopal meshes in two
 * and three dimensions.  Multiple materials can exist within the mesh, but each
 * cell within the mesh is only given a single material under the format
 * specification.
 *
 * Multiprocessor distribution is handled by having multiple X3D files describe
 * the complete mesh.  Each X3D file has a header that describes entities for
 * the processor to which said X3D file corresponds.  Mesh entities are listed
 * with local indices, but ghost information includes both the processor
 * identifier and that processor's local index for said ghost entity.
 *
 * NOTE: The current implementation here assumes a single X3D file.  Ghost
 *       information is currently ignored, as we assume we will be
 *       repartitioning the mesh anyway.  This will be sufficient for meshes
 *       small enough to be loaded entirely by a single processor, but may not
 *       be efficient.  This will need to be relaxed in the future for larger
 *       meshes.
 *
 * The X3D format specification allows for providing names, EOS identifiers, and
 * opacity identifiers for each material in the problem.  In addition,
 * cell-centered and node-centered field data can be stored within the X3D file.
 *
 * NOTE: The current implementation here will read the above material
 *       information, but we don't actually do anything with it.  In
 *       addition, we do NOT provide a method for reading of field
 *       data at all.  The reason for this is that a typical FleCSI
 *       application will need to register its own field data to live
 *       on the data client.  This includes materials and material
 *       properties.  If it is desired in the future to be able to use
 *       the data directly from the X3D file, we will need to find a
 *       good way to store such data so that it can be easily passed
 *       back to the application.
 *
 * The X3D format specification allows for defining nodes in a slave-master
 * configuration, whereby a single node (the slaved node) is constrained by two
 * or more other nodes (the master nodes).  This is useful for creating meshes
 * with dendrites (A.K.A. hanging nodes, T-junctions, etc.).
 *
 * NOTE: We currently do nothing with this constraint information; we read it
 *       simply to advance the file.  If a future application wishes to use such
 *       constraints, we will need to implement a good way to store this
 *       information such that it can be passed back to the host code.
 *
 * The X3D format specification treats a physical mesh face that is shared
 * between two mesh zones as two distinct faces.  For example, in the image
 * below, the physical mesh face between the two zones is listed twice with
 * ids 2 and 8.  Both of faces 2 and 8 will be composed of the two nodes with
 * ids 2 and 3, but the ordering will be different to indicate an outward
 * normal direction.  The X3D face connectivity section indicates which face
 * corresponds to a given "neighbor face"; it will tell us that face 8 is
 * a neighbor face of face 2, and vice versa.  Furthermore, since the neighbor
 * face could be on a different processor, we are given the processor ID of
 * the processor that owns the neighbor face.
 *
 *     4     3     6
 *     +-----+-----+
 *     |  3  |  7  |
 *     |4   2|8   6|
 *     |  1  |  5  |
 *     +-----+-----+
 *     1     2     5
 *
 * NOTE: The FleCSI connectivity we build below will collapse this duality of
 *       face ids to a single face.  One of the advanatages of having the dual
 *       specification is for constructing connectivity with parallelism.  Since
 *       we are not allowing reading of distributed X3D files at this time, we
 *       parse the neighbor face information simply to perform bookkeeping in
 *       the collapse process.
 */

// user includes
#include <flecsi-sp/io/detail.h>

#include <flecsi/topology/mesh_definition.h>
#include <flecsi/utils/logging.h>

// system includes
#include <cstddef>
#include <fstream>
#include <initializer_list>
#include <ios>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace flecsi_sp {
namespace io {

////////////////////////////////////////////////////////////////////////////////
/// \brief Datastructure to hold X3D file header parameters
////////////////////////////////////////////////////////////////////////////////
struct x3d_header {
  int process;  // processor for which this file is targetted (1-based)
  int numdim;  // spatial dimension
  int materials;  // total number of materials in the mesh
  // the remaining entries are all for the current processor
  std::size_t nodes;  // number of logically distinct nodes
  std::size_t faces;  // number of logically distinct faces
  std::size_t elements;  // number of cells
  std::size_t ghost_nodes;  // number of ghost nodes
  std::size_t slaved_nodes;  // number of slaved nodes
  // number of master nodes for most constrained slave
  // this will typically be 2 or 4 (depending on numdim), but should not be
  // zero even if there are no slaved nodes.
  int nodes_per_slave;
  int nodes_per_face;  // number of nodes needed for the most complex face
  int faces_per_cell;  // number of nodes needed for the most complex cell
  int node_data_fields;  // number of node-centered data fields for simulation
  int cell_data_fields;  // number of cell-centered data fields for simulation
};

////////////////////////////////////////////////////////////////////////////////
/// \brief This is the X3D mesh reader.
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class x3d_base__ {
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

  //! \brief an alias for the std::vector class
  template<typename U>
  using vector = typename std::vector<U>;

  //! \brief an alias for the matrix class
  template<typename U>
  using sparse_matrix = vector<vector<U>>;

  //! \brief the data type for an index vector
  using index_vector_t = vector<index_t>;

  //! \brief the data type for connectivty arrays
  using connectivity_t = sparse_matrix<index_t>;

  //! the number of dimensions
  static constexpr size_t num_dims = D;

  //! \brief Given an open file, read the next several lines assuming it is an
  //! X3D header, parse the information, and check for proper formatting.
  //! \param [in] f The open file object.
  //! \return an \e x3d_header object populated with the header data
  static auto verify_and_read_header(std::fstream & f) {
    std::string token;
    x3d_header header;
    // check for proper format
    std::getline(f, token);
    verify_string(token, X3D_TOKEN);
    // start of header block
    std::getline(f, token);
    verify_string(token, "header");
    // Get the contents
    f >> token >> header.process;
    f >> token >> header.numdim;
    f >> token >> header.materials;
    f >> token >> header.nodes;
    f >> token >> header.faces;
    f >> token >> header.elements;
    f >> token >> header.ghost_nodes;
    f >> token >> header.slaved_nodes;
    f >> token >> header.nodes_per_slave;
    f >> token >> header.nodes_per_face;
    f >> token >> header.faces_per_cell;
    f >> token >> header.node_data_fields;
    f >> token >> header.cell_data_fields;
    // end of header_block
    std::getline(f, token);  // to pick up \n from previous line
    std::getline(f, token);
    verify_string(token, "end_header");

    return header;
  }  // verify_and_read_header

  //! \brief Read material names from an X3D file, assuming the open file is
  //! at the correct location.
  //! \param [in] f The open file object.
  //! \param [in] nummats The number of materials to read.
  //! \return vector<std::string> of material names; order here corresponds
  //! to material ID's
  static auto get_material_names(std::fstream & f, int nummats) {
    vector<std::string> matnames(nummats);
    std::string token;
    // start of matnames block
    std::getline(f, token);
    verify_string(token, "matnames");
    // content
    for (int i = 0; i < nummats; ++i)
      f >> token >> matnames[i];
    // end of matnames block
    std::getline(f, token);  // to pick up \n from previous line
    std::getline(f, token);
    verify_string(token, "end_matnames");

    return matnames;
  }  // get_material_names

  //! \brief Read material EOS IDs from an X3D file, assuming the open file is
  //! at the correct location.
  //! \param [in] f The open file object.
  //! \param [in] nummats The number of materials to read.
  //! \return vector<int> of material EOS IDs; order here corresponds
  //! to material ID's
  static auto get_material_eosid(std::fstream & f, int nummats) {
    vector<int> mateosid(nummats);
    std::string token;
    // start of mateos block
    std::getline(f, token);
    verify_string(token, "mateos");
    // content
    for (int i = 0; i < nummats; ++i)
      f >> token >> mateosid[i];
    // end of mateos block
    std::getline(f, token);  // to pick up \n from previous line
    std::getline(f, token);
    verify_string(token, "end_mateos");

    return mateosid;
  }  // get_material_eosid

  //! \brief Read material opacity IDs from an X3D file, assuming the open file
  //! is at the correct location.
  //! \param [in] f The open file object.
  //! \param [in] nummats The number of materials to read.
  //! \return vector<int> of material opacity IDs; order here corresponds
  //! to material ID's
  static auto get_material_opacid(std::fstream & f, int nummats) {
    vector<int> matopacid(nummats);
    std::string token;
    // start of matopc block
    std::getline(f, token);
    verify_string(token, "matopc");
    // content
    for (int i = 0; i < nummats; ++i)
      f >> token >> matopacid[i];
    // end of matopc block
    std::getline(f, token);  // to pick up \n from previous line
    std::getline(f, token);
    verify_string(token, "end_matopc");

    return matopacid;
  }  // get_material_opacid

  //! \brief Read the node coordinate data from an X3D file, assuming the open
  //! file is at the correct location.
  //! \param [in] f The open file object.
  //! \param [in] numnodes The number of nodes to read.
  //! \param [in] dim The spatial dimension of the problem.
  //! \return sparse_matrix<real_t> with shape (numnodes, dim) of coordinates.
  //! Node order corresponds to node ID.
  //!
  //! The X3D format specification states that these coordinates should be
  //! specified as fully three-dimensional, with zeros for the entries not
  //! needed if the problem has less than three spatial dimensions.  Here we
  //! are bypassing this by only processing what we need.
  static auto get_node_coords(std::fstream & f, size_t numnodes, int dim) {
    sparse_matrix<real_t> coords;
    coords.reserve(numnodes);
    std::string token;
    // start of nodes block
    std::getline(f, token);
    verify_string(token, "nodes");
    // content
    for (size_t inode = 0; inode < numnodes; ++inode) {
      coords.resize(coords.size() + 1);
      auto & thisNode = coords.back();
      thisNode.resize(dim);
      f >> token;  // node ID
      for (int idim = 0; idim < dim; ++idim)
        f >> thisNode[idim];
      std::getline(f, token);  // read remainder of data
    }
    // end of nodes block
    std::getline(f, token);
    verify_string(token, "end_nodes");

    return coords;
  }  // get_node_coords

  //! \brief Read face-node connectivity from an X3D file, assuming the file is
  //! open to the correct location.
  //! \param [in] f The open file object.
  //! \param [in] numfaces The number of faces to read.
  //! \return connectivity_t of the nodes that make up each face.  See details.
  //!
  //! Faces that are not on a physical domain boundary are listed twice in the
  //! X3D format specification - once for each of the two cells on either side
  //! of the face.  The returned connectivity_t object will contain the node
  //! lists for all of the faces specified.  Furthermore, there is a final entry
  //! added to the returned object that will contain, for each face, the
  //! corresponding neighbor face ID.
  //!
  //! NOTE: This routine will need to be modified for distributed mesh reading.
  //!       This is because the face neighbor ID is given with local indexing to
  //!       the processor that owns that face.  That processor's ID is also
  //!       given in the X3D face specification, but we skip it in the current
  //!       implementation.
  static auto get_face_connectivity(std::fstream & f, size_t numfaces) {
    connectivity_t face_connectivity;
    // +1 in size here because of the added list of neighbor faces
    face_connectivity.reserve(numfaces + 1);
    index_vector_t neighbor_list;
    neighbor_list.reserve(numfaces);
    std::string token;
    size_t nnodes, node;
    // start of faces block
    std::getline(f, token);
    verify_string(token, "faces");
    // content
    for (size_t iface = 0; iface < numfaces; ++iface) {
      face_connectivity.resize(face_connectivity.size() + 1);
      auto & thisFace = face_connectivity.back();
      f >> token;  // face ID
      f >> nnodes;  // number of nodes for this face
      thisFace.resize(nnodes);
      for (size_t inode = 0; inode < nnodes; ++inode) {
        f >> node;
        thisFace[inode] = node - 1;  // X3D is 1-based
      }
      f >> token;  // owner processor ID
      f >> token;  // neighbor processor ID
      f >> node;  // ID of neighbor face
      // domain boundaries have a 0 in this field, so this will wrap around to
      // max(size_t)
      neighbor_list.emplace_back(node-1);
      // The remainder of the line contains unused fields, so skip it.
      // std::getline only will not work here in 3D where one can have more than
      // three nodes per face, and the line is too long, and gets split.  There
      // should be 5 more unused columns.
      f >> token >> token >> token >> token >> token;
      std::getline(f, token);  // to get \n from previous line
    }
    // end of faces block
    std::getline(f, token);
    verify_string(token, "end_faces");

    // add the neighbor information
    face_connectivity.emplace_back(neighbor_list);

    return face_connectivity;
  }  // get_face_connectivity

  //! \brief Read cell-face connectivity from an X3D file, assuming the file is
  //! open to the correct location.
  //! \param [in] f The open file object.
  //! \param [in] numcells The number of cells to read.
  //! \return connectivity_t of the faces that make up cells.
  static auto get_cell_connectivity(std::fstream & f, size_t numcells) {
    connectivity_t cell_connectivity;
    cell_connectivity.reserve(numcells);
    std::string token;
    size_t nfaces, face;
    // start of cells block
    std::getline(f, token);
    verify_string(token, "cells");
    // content
    for (size_t icell = 0; icell < numcells; ++icell) {
      cell_connectivity.resize(cell_connectivity.size() + 1);
      auto & thisCell = cell_connectivity.back();
      f >> token;  // cell ID
      f >> nfaces;  // number of faces for this cell
      thisCell.resize(nfaces);
      for (int iface = 0; iface < nfaces; ++iface) {
        f >> face;
        thisCell[iface] = face - 1;  // X3D is 1-based
      }
    }
    // end of cells block
    std::getline(f, token);  // to pick up \n from previous line
    std::getline(f, token);
    verify_string(token, "end_cells");

    return cell_connectivity;
  }  // get_cell_connectivity

  //! \brief Read slaved node data from an X3D file, assuming the file is open
  //! to the correct location.
  //! \param [in] f Open file object.
  //! \param [in] expected Number of expected slave nodes (for verification)
  //! \return connectivity_t which has an entry for each slaved node.
  //!
  //! Each slaved node has an entry within the returned connectivity_t object.
  //! For each slaved node, the entry is a vector whose first entry is the
  //! node ID of the slaved node.  The remaining entries (2 or 4 in number) are
  //! the node ID's of the master nodes for this slaved node.
  static auto get_slaved_node_connectivity(std::fstream & f, size_t expected) {
    connectivity_t slaved_nodes;
    slaved_nodes.reserve(expected);
    std::string token;
    int nodeID, nmasters;
    // start of slaved_nodes block
    // this block is different in that its start contains an entry for the
    // number of slaves to read
    f >> token;
    verify_string(token, "slaved_nodes");
    f >> token;  // number of slaved nodes listed next
    verify_string(token, std::to_string(expected));
    // content
    for (size_t islave = 0; islave < expected; ++islave) {
      slaved_nodes.resize(slaved_nodes.size() + 1);
      auto & thisSlave = slaved_nodes.back();
      f >> token;  // slaveID
      f >> nodeID;  // nodeID for this slave
      f >> nmasters;  // number of master nodes for this slave
      thisSlave.resize(nmasters + 1);  // also store the nodeID
      thisSlave[0] = nodeID - 1;  // X3D is 1-based
      for (size_t imaster = 1; imaster <= nmasters; ++imaster) {
        f >> nodeID;  // nodeID for the master
        thisSlave[imaster] = nodeID - 1;  // X3D is 1-based
      }
    }
    // end of slaved_nodes block
    std::getline(f, token);  // to pick up \n from previous line
    std::getline(f, token);
    verify_string(token, "end_slaved_nodes");

    return slaved_nodes;
  }  // get_slaved_node_connectivity

  //! \brief Read ghost node information from an X3D file, assuming the file is
  //! open to the correct location.
  //! \param [in] f Open file object.
  //! \param [in] expected Number of expected ghosts (for verification)
  //! \return connectivity_t which has an entry for each ghost
  //!
  //! Each ghost has an entry in the returned connectivit_t object.  For each
  //! ghost node, the entry is a vector whose first entry is the node ID for the
  //! ghost on the current processor.  The remaining two entries are the
  //! processor ID of the ghost node's owner, and the node ID according to the
  //! owner.
  static auto get_ghost_nodes(std::fstream & f, size_t expected) {
    connectivity_t ghost_nodes;
    ghost_nodes.reserve(expected);
    std::string token;
    size_t ghostID, proc, nodeID;
    // start of ghost_nodes block
    // this block is different in that its start contains an entry for the
    // number of ghosts to read
    f >> token;
    verify_string(token, "ghost_nodes");
    f >> token;  // number of ghost nodes listed next
    verify_string(token, std::to_string(expected));
    // content
    for (size_t ighost = 0; ighost < expected; ++ighost) {
      ghost_nodes.resize(ghost_nodes.size() + 1);
      auto & thisGhost = ghost_nodes.back();
      // this proc's ID for the node, the owner's proc ID, node ID on owner proc
      // last entry unused
      f >> ghostID >> proc >> nodeID >> token;
      // X3D is 1-based
      thisGhost.emplace_back(ghostID -1);
      thisGhost.emplace_back(proc - 1);
      thisGhost.emplace_back(nodeID - 1);
    }
    // end of ghost_nodes block
    std::getline(f, token);  // to get \n from previous line
    std::getline(f, token);
    verify_string(token, "end_ghost_nodes");

    return ghost_nodes;
  }  // get_ghost_nodes

 private:
  //! Indicator that this is an ASCII X3D file
  static constexpr auto X3D_TOKEN = "x3dtoflag ascii";

  //! Helper to check strings
  static void verify_string(const std::string given, const std::string wanted) {
    if (given != wanted)
      clog_fatal("Error parsing X3D file: expected '" << wanted <<
                 "' but received '" << given << "'");
  }
};  // class x3d_base__

////////////////////////////////////////////////////////////////////////////////
/// \brief The mesh definition for X3D
////////////////////////////////////////////////////////////////////////////////
template<int D, typename T>
class x3d_definition__;

////////////////////////////////////////////////////////////////////////////////
/// \brief The two-dimensional mesh definition based on the X3D format
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class x3d_definition__<2, T> : public flecsi::topology::mesh_definition_u<2> {
 public:
  //============================================================================
  // Typedefs
  //============================================================================
  //! the instantiated base type
  using base_t = x3d_base__<2, T>;

  //! the instantiated mesh definition type
  using mesh_definition_t = flecsi::topology::mesh_definition_u<2>;

  //! the number of dimensions
  using mesh_definition_t::dimension;

  //! the floating point type
  using real_t = typename base_t::real_t;

  //! the index type
  using index_t = typename base_t::index_t;

  //! the vector type
  template<typename U>
  using vector = typename base_t::template vector<U>;

  //! \brief the data type for an index vector
  using index_vector_t = typename base_t::index_vector_t;

  //! matrix type
  template<typename U>
  using matrix = typename base_t::template sparse_matrix<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;

  //============================================================================
  // Constructors
  //============================================================================
  //! \brief Default constructor
  x3d_definition__() = default;

  //! \brief Constructor with filename
  //! \param [in] filename The name of the file to load
  explicit x3d_definition__(const std::string & filename) {
    read(filename);
  }

  /// Copy constructor (disabled)
  x3d_definition__(const x3d_definition__ &) = delete;

  /// Assignment opeartor (disabled)
  x3d_definition__ & operator=(const x3d_definition__ &) = delete;

  /// Destructor
  ~x3d_definition__() = default;

  //============================================================================
  //! \brief Implementation of the X3D mesh reader.  Populates private data
  //! describing connectivity and vertex positions.
  //! \param [in] name Read mesh from \e name
  //============================================================================
  void read(const std::string & name) {
    clog(info) << "Reading X3D mesh from: " << name << std::endl;

    // Open the X3D file
    std::fstream file(name, std::ios_base::in);
    if (!file.is_open()) {
      clog_fatal("Error opening file: " << name);
    }

    // Verify and read header
    auto params = base_t::verify_and_read_header(file);

    // sanity checks
    if (dimension() != params.numdim)
      clog_fatal("Expected dimension " << dimension() <<
                 " in file " << name << " but received " <<
                 params.numdim);
    if (params.process != 1)
      clog_fatal("X3D reader does not support distributed X3D files.");

    // material information -- currently we do nothing with this
    base_t::get_material_names(file, params.materials);
    base_t::get_material_eosid(file, params.materials);
    base_t::get_material_opacid(file, params.materials);

    // connectivities
    auto & cell_vertices_ref = entities_[2][0];
    auto & edge_vertices_ref = entities_[1][0];
    auto & cell_edges_ref = entities_[2][1];

    // read coordinates
    vertices_ = base_t::get_node_coords(file, params.nodes, params.numdim);

    // read faces/edges to vertices connectivity
    // NOTE: each internal face is listed twice in this list from the X3D format
    //       and the last entry contains the items to deduplicate
    edge_vertices_ref = base_t::get_face_connectivity(file, params.faces);

    // read cells to faces/edges connectivity
    cell_edges_ref = base_t::get_cell_connectivity(file, params.elements);

    // Now we fix up the duplicated faces/edges and correct the IDs in the
    // cell connectivity
    dedupe_and_fixup(edge_vertices_ref, cell_edges_ref);

    // build cell to node connectivity from cell-to-edge and edge-to-node
    detail::intersect(cell_edges_ref, edge_vertices_ref, cell_vertices_ref);
    // and all the inversions
    // params.faces no longer the correct size
    entities_[1][2].reserve(edge_vertices_ref.size());
    entities_[0][1].reserve(params.nodes);
    entities_[0][2].reserve(params.nodes);
    detail::transpose(cell_edges_ref, entities_[1][2]);
    detail::transpose(edge_vertices_ref, entities_[0][1]);
    detail::transpose(cell_vertices_ref, entities_[0][2]);

    // cleanup
    file.close();
  }  // read

  //////////////////////////////////////////////////////////////////////////////
  // Required overrides from mesh_definition_u
  //////////////////////////////////////////////////////////////////////////////
  //! \brief Return the number of entities of a particular dimension
  //! \param [in] dim The entity dimension to query.
  size_t num_entities(size_t dim) const override {
    switch (dim) {
      case 0:
        return vertices_.size();
      case 1:
      case 2:
        return entities_.at(dim).at(0).size();
      default:
        clog_fatal(
            "Dimension out of range: 0 < " << dim << " </ " << dimension());
        return 0;
    }
  }  // num_entities

  //! \brief Return the set of entities of dimension \em to_dim
  //! that define entity \em id of dimension \em from_dim.
  //! \param [in] from_dim The dimension of the entity for which the definition
  //! is being requested.
  //! \param [in] to_dim The dimension of the entities of the definition.
  //! \param [in] from_id The id of the entitiy for which the definition is
  //! being requested.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override {
    return entities_.at(from_dim).at(to_dim).at(from_id);
  }  // entities

  //! \brief Return the set of entities of dimension \em to_dim that define
  //! entities of dimension \em from_dim.
  //! \param [in] from_dim The dimension of the entities for which the
  //! definition is being requested.
  //! \param [in] to_dim The dimension of the entities of the definition.
  const auto & entities(size_t from_dim, size_t to_dim) const {
    return entities_.at(from_dim).at(to_dim);
  }  // entities

  //! \brief Coordinates of a particular vertex.
  //! \tparam POINT_TYPE The data type that holds the coordinates; uses []
  //! \param [in] vertex_id The id of the vertex to query.
  //! \return POINT_TYPE object filled with coordinates.
  template<typename POINT_TYPE>
  auto vertex(const size_t vertex_id) const {
    POINT_TYPE p;
    for (int d = 0; d < dimension(); ++d)
      p[d] = vertices_[vertex_id][d];
    return p;
  }  // vertex

  //! \brief THIS METHOD REQUIRED BY partition_mesh BUT IS A NO-OP
  template<typename U = int>
  void write(
      const std::string & name,
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
      element_sets = {},
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
      node_sets = {}) const {
    clog(info) << "X3D WON'T BE WRITING A MESH TO: " << name << std::endl;
  }

 private:
  //! \brief Remove duplicated physical faces and correct cell-face IDs.
  //! \param [inout] edge_vertices The list of vertices for each edge, including
  //! duplicates.
  //! \param [inout] cell_edges The list of edges for each cell.
  //!
  //! On return, edge_vertices is stripped up duplicate edges and the list of
  //! duplicate edges.  Similarly, cell_edges has had all of its entries updated
  //! to mirror the ID's of only the remaining edges after deduplication.
  void dedupe_and_fixup(connectivity_t & edge_vertices,
                        connectivity_t & cell_edges) {
    // last entry contains the duplicates
    // a value in this vector of max(size_t) indicates a domain boundary
    auto dupes = edge_vertices.back();
    edge_vertices.pop_back();
    auto boundary = std::numeric_limits<size_t>::max();
    size_t pos(0), offset(0);
    // Lookup table for mapping original edge IDs to the new IDs created after
    // removal and accounting for duplicates.
    index_vector_t edge_LUT;
    edge_LUT.reserve(edge_vertices.size());
    // This is pretty brute force and expensive...
    for (auto d : dupes) {
      // a non-boundary face that has already been encountered
      if ((d < boundary) && (d < pos)) {
        // we've encountered this before, so point to it
        // note this might have been offset previously
        edge_LUT.push_back(edge_LUT[d]);
        // get rid of the duplicate
        edge_vertices.erase(edge_vertices.begin() + pos - offset);
        // keep track of how many times we have done this
        ++offset;
      } else {
        // valid face
        edge_LUT.push_back(pos - offset);
      }
      ++pos;
    }  // dupes
    // We've eliminated duplicate faces, so now fix the cell connectivity.
    // Again, brute force checking everything...
    for (auto & edges : cell_edges) {
      for (auto & e : edges)
        e = edge_LUT[e];
    }
  }  // dedupe_and_fixup

  //! \brief storage for element vertices
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for vertex coordinates
  matrix<real_t> vertices_;
};  // x3d_definition__<2, T>

////////////////////////////////////////////////////////////////////////////////
/// \brief The three-dimensional mesh definition based on the X3D format
////////////////////////////////////////////////////////////////////////////////
template<typename T>
class x3d_definition__<3, T> : public flecsi::topology::mesh_definition_u<3> {
 public:
  //============================================================================
  // Typedefs
  //============================================================================
  //! the instantiated base type
  using base_t = x3d_base__<3, T>;

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

  //! \brief the data type for an index vector
  using index_vector_t = typename base_t::index_vector_t;

  //! matrix type
  template<typename U>
  using matrix = typename base_t::template sparse_matrix<U>;

  //! the connectivity type
  using connectivity_t = typename base_t::connectivity_t;

  //============================================================================
  // Constructors
  //============================================================================
  //! \brief Default constructor
  x3d_definition__() = default;

  //! \brief Constructor with filename
  //! \param [in] filename The name of the file to load
  explicit x3d_definition__(const std::string & filename) {
    read(filename);
  }

  /// Copy constructor (disabled)
  x3d_definition__(const x3d_definition__ &) = delete;

  /// Assignment opeartor (disabled)
  x3d_definition__ & operator=(const x3d_definition__ &) = delete;

  /// Destructor
  ~x3d_definition__() = default;

  //============================================================================
  //! \brief Implementation of the X3D mesh reader.  Populates private data
  //! describing connectivity and vertex positions.
  //! \param [in] name Read mesh from \e name
  //============================================================================
  void read(const std::string & name) {
    clog(info) << "Reading X3D mesh from: " << name << std::endl;

    // Open the X3D file
    std::fstream file(name, std::ios_base::in);
    if (!file.is_open()) {
      clog_fatal("Error opening file: " << name);
    }

    // Verify and read header
    auto params = base_t::verify_and_read_header(file);

    // sanity checks
    if (dimension() != params.numdim)
      clog_fatal("Expected dimension " << dimension() <<
                 " in file " << name << " but received " <<
                 params.numdim);
    if (params.process != 1)
      clog_fatal("X3D reader does not support distributed X3D files.");

    // material information -- currently we do nothing with this
    base_t::get_material_names(file, params.materials);
    base_t::get_material_eosid(file, params.materials);
    base_t::get_material_opacid(file, params.materials);

    // connectivities
    auto & cell_faces_ref = entities_[3][2];
    auto & cell_edges_ref = entities_[3][1];
    auto & cell_vertices_ref = entities_[3][0];
    auto & face_edges_ref = entities_[2][1];
    auto & face_vertices_ref = entities_[2][0];
    auto & edge_vertices_ref = entities_[1][0];

    // read coordinates
    vertices_ = base_t::get_node_coords(file, params.nodes, params.numdim);

    // read faces to vertices connectivity
    // NOTE: each internal face is listed twice in this list from the X3D format
    //       and the last entry contains the items to deduplicate
    face_vertices_ref = base_t::get_face_connectivity(file, params.faces);

    // read cells to face connectivity
    cell_faces_ref = base_t::get_cell_connectivity(file, params.elements);

    // Now we fix up the duplicated faces and correct the IDs in the
    // cell connectivity
    dedupe_and_fixup(face_vertices_ref, cell_faces_ref);

    // build cell to node connectivity from cell-to-face and face-to-node
    detail::intersect(cell_faces_ref, face_vertices_ref, cell_vertices_ref);

    // build face to edge and edge to node connectivity from face to node
    auto edge_vertices_sorted = std::make_unique<connectivity_t>();
    detail::build_connectivity(
        face_vertices_ref, face_edges_ref, edge_vertices_ref,
        *edge_vertices_sorted, [](const auto & vs, auto & edge_vs) {
          using list_type = std::decay_t<decltype(*edge_vs.begin())>;
          for (auto v0 = std::prev(vs.end()), v1 = vs.begin(); v1 != vs.end();
               v0 = v1, ++v1)
            edge_vs.emplace_back(list_type{*v0, *v1});
        });

    // build cell to edge connectivity from cell to face and face to edge
    detail::intersect(cell_faces_ref, face_edges_ref, cell_edges_ref);

    // and all the inversions
    // params.faces no longer the correct size
    entities_[0][1].reserve(params.nodes);
    entities_[0][2].reserve(params.nodes);
    entities_[0][3].reserve(params.nodes);
    entities_[1][2].reserve(edge_vertices_ref.size());
    entities_[1][3].reserve(edge_vertices_ref.size());
    entities_[2][3].reserve(face_vertices_ref.size());

    detail::transpose(edge_vertices_ref, entities_[0][1]);
    detail::transpose(face_vertices_ref, entities_[0][2]);
    detail::transpose(cell_vertices_ref, entities_[0][3]);
    detail::transpose(face_edges_ref, entities_[1][2]);
    detail::transpose(cell_edges_ref, entities_[1][3]);
    detail::transpose(cell_faces_ref, entities_[2][3]);

    // cleanup
    file.close();
  }  // read

  //////////////////////////////////////////////////////////////////////////////
  // Required overrides from mesh_definition_u
  //////////////////////////////////////////////////////////////////////////////
  //! \brief Return the number of entities of a particular dimension
  //! \param [in] dim The entity dimension to query.
  size_t num_entities(size_t dim) const override {
    switch (dim) {
      case 0:
        return vertices_.size();
      default:
        return entities_.at(dim).at(0).size();
    }
  }  // num_entities

  //! \brief Return the set of entities of dimension \em to_dim
  //! that define entity \em id of dimension \em from_dim.
  //! \param [in] from_dim The dimension of the entity for which the definition
  //! is being requested.
  //! \param [in] to_dim The dimension of the entities of the definition.
  //! \param [in] from_id The id of the entitiy for which the definition is
  //! being requested.
  std::vector<size_t>
  entities(size_t from_dim, size_t to_dim, size_t from_id) const override {
    return entities_.at(from_dim).at(to_dim).at(from_id);
  }  // entities

  //! \brief Return the set of entities of dimension \em to_dim that define
  //! entities of dimension \em from_dim.
  //! \param [in] from_dim The dimension of the entities for which the
  //! definition is being requested.
  //! \param [in] to_dim The dimension of the entities of the definition.
  const auto & entities(size_t from_dim, size_t to_dim) const {
    return entities_.at(from_dim).at(to_dim);
  }  // entities

  //! \brief Coordinates of a particular vertex.
  //! \tparam POINT_TYPE The data type that holds the coordinates; uses []
  //! \param [in] vertex_id The id of the vertex to query.
  //! \return POINT_TYPE object filled with coordinates.
  template<typename POINT_TYPE>
  auto vertex(const size_t vertex_id) const {
    POINT_TYPE p;
    for (int d = 0; d < dimension(); ++d)
      p[d] = vertices_[vertex_id][d];
    return p;
  }  // vertex

  //! \brief THIS METHOD REQUIRED BY partition_mesh BUT IS A NO-OP
  template<typename U = int>
  void write(
      const std::string & name,
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
      element_sets = {},
      const std::initializer_list<std::pair<const char *, std::vector<U>>> &
      node_sets = {}) const {
    clog(info) << "X3D WON'T BE WRITING A MESH TO: " << name << std::endl;
  }

 private:
  //! \brief Remove duplicated physical faces and correct cell-face IDs.
  //! \param [inout] edge_vertices The list of vertices for each edge, including
  //! duplicates.
  //! \param [inout] cell_edges The list of edges for each cell.
  //!
  //! On return, edge_vertices is stripped up duplicate edges and the list of
  //! duplicate edges.  Similarly, cell_edges has had all of its entries updated
  //! to mirror the ID's of only the remaining edges after deduplication.
  void dedupe_and_fixup(connectivity_t & edge_vertices,
                        connectivity_t & cell_edges) {
    // last entry contains the duplicates
    // a value in this vector of max(size_t) indicates a domain boundary
    auto dupes = edge_vertices.back();
    edge_vertices.pop_back();
    auto boundary = std::numeric_limits<size_t>::max();
    size_t pos(0), offset(0);
    // Lookup table for mapping original edge IDs to the new IDs created after
    // removal and accounting for duplicates.
    index_vector_t edge_LUT;
    edge_LUT.reserve(edge_vertices.size());
    // This is pretty brute force and expensive...
    for (auto d : dupes) {
      // a non-boundary face that has already been encountered
      if ((d < boundary) && (d < pos)) {
        // we've encountered this before, so point to it
        // note this might have been offset previously
        edge_LUT.push_back(edge_LUT[d]);
        // get rid of the duplicate
        edge_vertices.erase(edge_vertices.begin() + pos - offset);
        // keep track of how many times we have done this
        ++offset;
      } else {
        // valid face
        edge_LUT.push_back(pos - offset);
      }
      ++pos;
    }  // dupes
    // We've eliminated duplicate faces, so now fix the cell connectivity.
    // Again, brute force checking everything...
    for (auto & edges : cell_edges) {
      for (auto & e : edges)
        e = edge_LUT[e];
    }
  }  // dedupe_and_fixup

  //! \brief storage for element vertices
  std::map<index_t, std::map<index_t, connectivity_t>> entities_;

  //! \brief storage for vertex coordinates
  matrix<real_t> vertices_;
};  // x3d_definition__<3, T>

}  // namespace io
}  // namespace flecsi_sp
