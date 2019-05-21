/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

// the main test include
#include <cinchtest.h>

// user includes
#include <flecsi/topology/closure_utils.h>
#include <flecsi-sp/io/exodus_definition.h>

// some type aliases confined to this test
template< int D >
using exodus_definition_t = 
  flecsi_sp::io::exodus_definition<D, double>;

// The flecsi library has undefined symbols in it.  It calls 
// flecsi::execution::driver even though we are not using the execution model.
// Define an empty stub for linkage.
namespace flecsi {
namespace execution {
void driver(int argc, char ** argv) {}
}
}

////////////////////////////////////////////////////////////////////////////////
/// \brief This tests the 1d version of the exodus mesh definition
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_1d, simple) {

  exodus_definition_t<1> simple("uniform_32.exo");

} // TEST

////////////////////////////////////////////////////////////////////////////////
/// \brief This tests the 2d version of the exodus mesh definition
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_2d, simple) {

  exodus_definition_t<2> mixed("mixed.exo");
  exodus_definition_t<2> nfaced("voro.g");

} // TEST

////////////////////////////////////////////////////////////////////////////////
/// \brief This tests the 3d version of the exodus mesh definition
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_3d, simple) {

  exodus_definition_t<3> hex("box-hex.exo");
  exodus_definition_t<3> tet("box-tet.exo");
  exodus_definition_t<3> nfaced("box-nfaced.exo");

} // TEST


////////////////////////////////////////////////////////////////////////////////
/// \brief This excersizes some of the partitioning functionality in 1d
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_1d, neighbors) {

  exodus_definition_t<1> mesh("uniform_32.exo");

  // Primary partititon
  std::vector<size_t> selected_cells = { 8 };

  // The closure captures any cell that is adjacent to a cell in the
  // set of indices passed to the method. The closure includes the
  // initial set of indices.
  
  // here we include all cells that use a vertex of the original cell
  auto vertex_neighbors = 
    flecsi::topology::entity_neighbors<1,1,0>(mesh, selected_cells);

  ASSERT_EQ(
    vertex_neighbors, 
    std::set<size_t>({7, 8, 9})
  );
} // TEST

////////////////////////////////////////////////////////////////////////////////
/// \brief This excersizes some of the partitioning functionality in 2d
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_2d, neighbors) {

  exodus_definition_t<2> mesh("mixed.exo");

  // Primary partititon
  std::vector<size_t> selected_cells = { 8 };

  // The closure captures any cell that is adjacent to a cell in the
  // set of indices passed to the method. The closure includes the
  // initial set of indices.
  
  // here we include all cells that use a vertex of the original cell
  auto vertex_neighbors = 
    flecsi::topology::entity_neighbors<2,2,0>(mesh, selected_cells);

  ASSERT_EQ(
    vertex_neighbors, 
    std::set<size_t>({7, 8, 9, 57, 58, 59, 60, 61, 169, 170})
  );

  // now we include all neighbors that share an edge
  auto edge_neighbors = 
    flecsi::topology::entity_neighbors<2,2,1>(mesh, selected_cells);

  ASSERT_EQ(
    edge_neighbors, 
    std::set<size_t>({7, 8, 9, 58, 60})
  );
} // TEST

////////////////////////////////////////////////////////////////////////////////
/// \brief This excersizes some of the partitioning functionality in 3d
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_3d, neighbors_hex) {

  exodus_definition_t<3> mesh("box-hex.exo");

  // Primary partititon
  std::vector<size_t> selected_cells = { 21 };

  // The closure captures any cell that is adjacent to a cell in the
  // set of indices passed to the method. The closure includes the
  // initial set of indices.
  
  // here we include all cells that use a vertex of the original cell
  auto vertex_neighbors = 
    flecsi::topology::entity_neighbors<3,3,0>(mesh, selected_cells);

  ASSERT_EQ( vertex_neighbors.size(), 27 );
  ASSERT_EQ(
    vertex_neighbors, 
    std::set<size_t>({
      0, 1, 2, 4, 5, 6, 8, 9, 10, 16, 17, 18, 20, 21, 22, 24, 25, 26, 32, 33, 
      34, 36, 37, 38, 40, 41, 42
    })
  );

  // now we include all neighbors that share an edge
  auto edge_neighbors = 
    flecsi::topology::entity_neighbors<3,3,1>(mesh, selected_cells);

  ASSERT_EQ( edge_neighbors.size(), 19 );
  ASSERT_EQ(
    edge_neighbors, 
    std::set<size_t>({
      1, 4, 5, 6, 9, 16, 17, 18, 20, 21, 22, 24, 25, 26, 33, 36, 37, 38, 41
    })
  );
  
  // now we include all neighbors that share a face
  auto face_neighbors = 
    flecsi::topology::entity_neighbors<3,3,2>(mesh, selected_cells);

  ASSERT_EQ( face_neighbors.size(), 7 );
  ASSERT_EQ(
    face_neighbors, 
    std::set<size_t>({5, 17, 20, 21, 22, 25, 37})
  );

} // TEST

////////////////////////////////////////////////////////////////////////////////
/// \brief This excersizes some of the partitioning functionality in 3d
////////////////////////////////////////////////////////////////////////////////
TEST(exodus_definition_3d, neighbors_nfaced) {

  exodus_definition_t<3> mesh("box-nfaced.exo");

  // Primary partititon
  std::vector<size_t> selected_cells = { 21 };

  // The closure captures any cell that is adjacent to a cell in the
  // set of indices passed to the method. The closure includes the
  // initial set of indices.
  
  // here we include all cells that use a vertex of the original cell
  auto vertex_neighbors = 
    flecsi::topology::entity_neighbors<3,3,0>(mesh, selected_cells);

  ASSERT_EQ( vertex_neighbors.size(), 27 );
  ASSERT_EQ(
    vertex_neighbors, 
    std::set<size_t>({
      0, 1, 2, 4, 5, 6, 8, 9, 10, 16, 17, 18, 20, 21, 22, 24, 25, 26, 32, 33, 
      34, 36, 37, 38, 40, 41, 42
    })
  );

  // now we include all neighbors that share an edge
  auto edge_neighbors = 
    flecsi::topology::entity_neighbors<3,3,1>(mesh, selected_cells);

  ASSERT_EQ( edge_neighbors.size(), 19 );
  ASSERT_EQ(
    edge_neighbors, 
    std::set<size_t>({
      1, 4, 5, 6, 9, 16, 17, 18, 20, 21, 22, 24, 25, 26, 33, 36, 37, 38, 41
    })
  );
  
  // now we include all neighbors that share a face
  auto face_neighbors = 
    flecsi::topology::entity_neighbors<3,3,2>(mesh, selected_cells);

  ASSERT_EQ( face_neighbors.size(), 7 );
  ASSERT_EQ(
    face_neighbors, 
    std::set<size_t>({5, 17, 20, 21, 22, 25, 37})
  );

} // TEST
