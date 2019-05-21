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
    flecsi_sp::io::x3d_base<D, double>;

template< int D >
using x3d_definition_t =
  flecsi_sp::io::x3d_definition<D, double>;

// The flecsi library has undefined symbols in it.  It calls
// flecsi::execution::driver even though we are not using the execution model.
// Define an empty stub for linkage.
namespace flecsi {
namespace execution {
void driver(int argc, char ** argv) {}
}
}

////////////////////////////////////////////////////////////////////////////////
// x3d_base tests
////////////////////////////////////////////////////////////////////////////////
TEST(x3d_base, verify_header2d) {
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

TEST(x3d_base, verify_header3d) {
  std::fstream file("spun_uniform_32.x3d");
  auto header_params = x3d_base_t<3>::verify_and_read_header(file);

  ASSERT_EQ(header_params.process, 1);
  ASSERT_EQ(header_params.numdim, 3);
  ASSERT_EQ(header_params.materials, 1);
  ASSERT_EQ(header_params.nodes, 5478);
  ASSERT_EQ(header_params.faces, 24384);
  ASSERT_EQ(header_params.elements, 4032);
  ASSERT_EQ(header_params.ghost_nodes, 0);
  ASSERT_EQ(header_params.slaved_nodes, 231);
  ASSERT_EQ(header_params.nodes_per_slave, 2);
  ASSERT_EQ(header_params.nodes_per_face, 5);
  ASSERT_EQ(header_params.faces_per_cell, 7);
  ASSERT_EQ(header_params.node_data_fields, 0);
  ASSERT_EQ(header_params.cell_data_fields, 2);

  file.close();
}  // TEST

// The rest depend on the proper order in which they are called, so we'll do it
// all
TEST(x3d_base, verify_file2d) {
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

TEST(x3d_base, verify_file3d) {
  std::fstream file("spun_uniform_32.x3d");
  auto header_params = x3d_base_t<3>::verify_and_read_header(file);
  // material names
  auto matnames = x3d_base_t<3>::get_material_names(file,
                                                    header_params.materials);
  ASSERT_EQ(matnames.size(), 1);
  ASSERT_EQ(matnames[0], "1");
  // material EOS IDs
  auto mateosids = x3d_base_t<3>::get_material_eosid(file,
                                                     header_params.materials);
  ASSERT_EQ(mateosids.size(), 1);
  ASSERT_EQ(mateosids[0], -1);
  // material opacity IDs
  auto matopacids = x3d_base_t<3>::get_material_opacid(file,
                                                       header_params.materials);
  ASSERT_EQ(matopacids.size(), 1);
  ASSERT_EQ(matopacids[0], -1);
  // node coordinates
  auto nodeCoords = x3d_base_t<3>::get_node_coords(file,
                                                   header_params.nodes,
                                                   header_params.numdim);
  // correct number of nodes
  ASSERT_EQ(nodeCoords.size(), header_params.nodes);
  // correct spatial size
  ASSERT_EQ(nodeCoords[0].size(), header_params.numdim);
  // Some random points
  ASSERT_DOUBLE_EQ(nodeCoords[0][1], 0.0);
  ASSERT_DOUBLE_EQ(nodeCoords[2][0], 3.125e-2);
  ASSERT_DOUBLE_EQ(nodeCoords[113][2], 6.25e-2);
  ASSERT_DOUBLE_EQ(nodeCoords[594][1], 3.4375e-1);
  ASSERT_DOUBLE_EQ(nodeCoords[1080][0], 3.200825214724776e-1);
  ASSERT_DOUBLE_EQ(nodeCoords[5476][2], 1.950903220161283e-1);
  // face connectivity
  auto face_conn = x3d_base_t<3>::get_face_connectivity(file,
                                                        header_params.faces);
  // correct number of faces
  // +1 here because get_face_connectivity also appends a vector containing
  // face neighbor information
  ASSERT_EQ(face_conn.size(), header_params.faces+1);
  // Some random connectivities
  ASSERT_EQ(face_conn[0][0], 0);
  ASSERT_EQ(face_conn[0][1], 2);
  ASSERT_EQ(face_conn[80][0], 81);
  ASSERT_EQ(face_conn[80][1], 82);
  ASSERT_EQ(face_conn[286][0], 197);
  ASSERT_EQ(face_conn[286][1], 198);
  // face with dendrite
  ASSERT_EQ(face_conn[618].size(), 5);
  ASSERT_EQ(face_conn[618][0], 301);
  ASSERT_EQ(face_conn[618][1], 302);
  ASSERT_EQ(face_conn[618][2], 371);
  ASSERT_EQ(face_conn[618][3], 370);
  ASSERT_EQ(face_conn[618][4], 369);
  // cell connectivity
  auto cell_conn = x3d_base_t<3>::get_cell_connectivity(file,
                                                        header_params.elements);
  // correct number of cells
  ASSERT_EQ(cell_conn.size(), header_params.elements);
  // Random cell connection
  ASSERT_EQ(cell_conn[2143].size(), 7);
  ASSERT_EQ(cell_conn[2143][0], 8655);
  ASSERT_EQ(cell_conn[2143][1], 8662);
  ASSERT_EQ(cell_conn[2143][2], 8663);
  ASSERT_EQ(cell_conn[2143][3], 8667);
  ASSERT_EQ(cell_conn[2143][4], 8671);
  ASSERT_EQ(cell_conn[2143][5], 20606);
  ASSERT_EQ(cell_conn[2143][6], 20601);
  // slaved nodes - we don't do anything with this information, but just read it
  auto slave_conn =
      x3d_base_t<3>::get_slaved_node_connectivity(file,
                                                  header_params.slaved_nodes);
  ASSERT_EQ(slave_conn.size(), header_params.slaved_nodes);
  // ghost nodes - this file has none, so we just ensure correct size and make
  // sure we can properly read the block
  auto ghost_nodes = x3d_base_t<3>::get_ghost_nodes(file,
                                                    header_params.ghost_nodes);
  ASSERT_EQ(ghost_nodes.size(), header_params.ghost_nodes);

  file.close();
}  // TEST

////////////////////////////////////////////////////////////////////////////////
// x3d_definition tests
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

  // original face with ID 4 was a duplicate (with face ID 2)
  // check to ensure it was removed and that the new face ID 4 corresponds to
  // the original face with ID 5
  auto face4_vertices = simple.entities(1, 0, 4);
  // vertex IDs making up this 2d face
  ASSERT_EQ(face4_vertices[0], 2);
  ASSERT_EQ(face4_vertices[1], 4);

  // similarly, ensure that that the cells have the correct numbering
  auto cell0_faces = simple.entities(2, 1, 0);
  auto cell1_faces = simple.entities(2, 1, 1);
  ASSERT_EQ(cell0_faces[0], 0);
  ASSERT_EQ(cell0_faces[1], 1);
  ASSERT_EQ(cell0_faces[2], 2);
  ASSERT_EQ(cell0_faces[3], 3);
  ASSERT_EQ(cell1_faces[0], 2);
  ASSERT_EQ(cell1_faces[1], 4);
  ASSERT_EQ(cell1_faces[2], 5);
  ASSERT_EQ(cell1_faces[3], 6);
}  // TEST

//! \brief Read a 3d spun version of the simple uniform mesh
TEST(x3d_definition_3d, simple) {
  x3d_definition_t<3> simple("spun_uniform_32.x3d");

  // num cells
  ASSERT_EQ(simple.num_entities(3), 4032);
  // num faces
  // 24384 original X3D faces.  2556 of them are on the boundary.  The rest are
  // double counted.  So, account for the duplicates that we removed.
  auto numfaces = (24384 - 2556)/2 + 2556;
  ASSERT_EQ(simple.num_entities(2), numfaces);
  // num edges
  // the spun mesh is made up of 33 identical planes
  // each plane has 291 edges
  // each plane is connected by 166 edges
  auto numedges = 33*291 + 32*166;
  ASSERT_EQ(simple.num_entities(1), numedges);
  // num nodes
  ASSERT_EQ(simple.num_entities(0), 5478);

  // original face with ID 3 was a duplicate (with face ID 2)
  // check to ensure it was removed and that the new face ID 3 corresponds to
  // the original face with ID 4
  auto face3_vertices = simple.entities(2, 0, 3);
  // vertex IDs making up this 3d face
  ASSERT_EQ(face3_vertices[0], 3);
  ASSERT_EQ(face3_vertices[1], 4);
  ASSERT_EQ(face3_vertices[2], 7);
  ASSERT_EQ(face3_vertices[3], 6);

  // similarly, ensure that that the cells have the correct numbering
  // no easy way to check ID's of higher numbered faces because there is no easy
  // way to know how many we have removed up to that number.
  auto cell0_faces = simple.entities(3, 2, 0);
  auto cell1_faces = simple.entities(3, 2, 1);
  ASSERT_EQ(cell0_faces[2], 2);
  ASSERT_EQ(cell1_faces[0], 2);
}  // TEST
