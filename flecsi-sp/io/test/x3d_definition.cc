/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2018 Triad National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

// the main test include
#include <cinchtest.h>

// user includes
#include <flecsi/topology/closure_utils.h>
#include <flecsi-sp/io/x3d_definition.h>

// system includes
#include <stdexcept>
#include <fstream>
#include <iostream>

// some type aliases confined to this test
template<int D>
using x3d_base_t =
    flecsi_sp::io::x3d_base__<D, double>;

template< int D >
using x3d_definition_t =
  flecsi_sp::io::x3d_definition__<D, double>;

// The flecsi library has undefined symbols in it.  It calls
// flecsi::execution::driver even though we are not using the execution model.
// Define an empty stub for linkage.
namespace flecsi {
namespace execution {
void driver(int argc, char ** argv) {}
}
}

////////////////////////////////////////////////////////////////////////////////
// x3d_base__ tests
////////////////////////////////////////////////////////////////////////////////
TEST(x3d_base__, verify_header) {
  std::fstream file("uniform_32.x3d.00001");
  auto header_params = x3d_base_t<2>::verify_and_read_header(file);

  ASSERT_EQ(header_params.process, 1);
  ASSERT_EQ(header_params.numdim, 2);
  ASSERT_EQ(header_params.materials, 1);
  ASSERT_EQ(header_params.nodes, 1089);
  ASSERT_EQ(header_params.faces, 4096);
  ASSERT_EQ(header_params.elements, 1024);
  ASSERT_EQ(header_params.ghost_nodes, 0);
  ASSERT_EQ(header_params.slaved_nodes, 0);
  ASSERT_EQ(header_params.nodes_per_slave, 2);
  ASSERT_EQ(header_params.nodes_per_face, 2);
  ASSERT_EQ(header_params.faces_per_cell, 4);
  ASSERT_EQ(header_params.node_data_fields, 0);
  ASSERT_EQ(header_params.cell_data_fields, 2);

  file.close();
}  // TEST

// The rest depend on the proper order in which they are called, so we'll do it
// all
TEST(x3d_base__, verify_file) {
  std::fstream file("uniform_32.x3d.00001");
  auto header_params = x3d_base_t<2>::verify_and_read_header(file);
  // material names
  auto matnames = x3d_base_t<2>::get_material_names(file,
                                                    header_params.materials);
  ASSERT_EQ(matnames.size(), 1);
  ASSERT_EQ(matnames[0], "1");
  // material EOS IDs
  auto mateosids = x3d_base_t<2>::get_material_eosid(file,
                                                     header_params.materials);
  ASSERT_EQ(mateosids.size(), 1);
  ASSERT_EQ(mateosids[0], -1);
  // material opacity IDs
  auto matopacids = x3d_base_t<2>::get_material_opacid(file,
                                                       header_params.materials);
  ASSERT_EQ(matopacids.size(), 1);
  ASSERT_EQ(matopacids[0], -1);
  // node coordinates
  auto nodeCoords = x3d_base_t<2>::get_node_coords(file,
                                                   header_params.nodes,
                                                   header_params.numdim);
  // correct number of nodes
  ASSERT_EQ(nodeCoords.size(), header_params.nodes);
  // correct spatial size
  ASSERT_EQ(nodeCoords[0].size(), header_params.numdim);
  // Some random points
  ASSERT_DOUBLE_EQ(nodeCoords[0][1], 0.0);
  ASSERT_DOUBLE_EQ(nodeCoords[2][0], 3.125e-2);
  ASSERT_DOUBLE_EQ(nodeCoords[113][1], 4.375e-1);
  ASSERT_DOUBLE_EQ(nodeCoords[594][1], 0.0);
  ASSERT_DOUBLE_EQ(nodeCoords[1080][0], 1.0);
  // face connectivity
  auto face_conn = x3d_base_t<2>::get_face_connectivity(file,
                                                        header_params.faces);
  // correct number of faces
  // +1 here because get_face_connectivity also appends a vector containing
  // face neighbor information
  ASSERT_EQ(face_conn.size(), header_params.faces+1);
  // Some random connectivities
  ASSERT_EQ(face_conn[0][0], 0);
  ASSERT_EQ(face_conn[0][1], 1);
  ASSERT_EQ(face_conn[80][0], 41);
  ASSERT_EQ(face_conn[80][1], 40);
  ASSERT_EQ(face_conn[286][0], 107);
  ASSERT_EQ(face_conn[286][1], 74);
  // cell connectivity
  auto cell_conn = x3d_base_t<2>::get_cell_connectivity(file,
                                                        header_params.elements);
  // correct number of cells
  ASSERT_EQ(cell_conn.size(), header_params.elements);
  // Random cell connection
  ASSERT_EQ(cell_conn[643].size(), 4);
  ASSERT_EQ(cell_conn[643][0], 2572);
  ASSERT_EQ(cell_conn[643][1], 2573);
  ASSERT_EQ(cell_conn[643][2], 2574);
  ASSERT_EQ(cell_conn[643][3], 2575);
  // slaved nodes - this file has none, so we just ensure correct size and make
  // sure we can properly read the block
  auto slave_conn =
      x3d_base_t<2>::get_slaved_node_connectivity(file,
                                                  header_params.slaved_nodes);
  ASSERT_EQ(slave_conn.size(), 0);
  // ghost nodes - this file has none, so we just ensure correct size and make
  // sure we can properly read the block
  auto ghost_nodes = x3d_base_t<2>::get_ghost_nodes(file,
                                                    header_params.ghost_nodes);
  ASSERT_EQ(ghost_nodes.size(), 0);

  file.close();
}  // TEST

////////////////////////////////////////////////////////////////////////////////
// x3d_definition__ tests
////////////////////////////////////////////////////////////////////////////////
//! \brief Read a simple uniform 2d mesh with 32x32 cells
TEST(x3d_definition_2d, simple) {
  x3d_definition_t<2> simple("uniform_32.x3d.00001");

  // num cells
  ASSERT_EQ(simple.num_entities(2), 32*32);
  // num faces
  ASSERT_EQ(simple.num_entities(1), 2*32*33);
  // num nodes
  ASSERT_EQ(simple.num_entities(0), 33*33);

  std::cout << "done" << std::endl;
}  // TEST
