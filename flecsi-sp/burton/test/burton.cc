/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
// \file
// \brief Tests general features of the burton mesh.
////////////////////////////////////////////////////////////////////////////////

// user includes
#include <cinchtest.h>
#include <flecsi/execution/execution.h>
#include <flecsi-sp/burton/burton_mesh.h>
#include <flecsi-sp/burton/mesh_interface.h>
#include <flecsi-sp/utils/types.h>

// system includes
#include <array>

// using statements
using std::cout;
using std::endl;

using mesh_t = flecsi_sp::burton::burton_mesh_t;
using index_spaces_t = mesh_t::index_spaces_t;

using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;
using array_t = std::array<int, 9>;

// a temporary data type
struct data_t
{
  int i;
  double x;
};

namespace flecsi_sp {
namespace burton {
namespace test {


////////////////////////////////////////////////////////////////////////////////
//! \brief Construct the test prefix string.
////////////////////////////////////////////////////////////////////////////////
auto prefix() 
{
  std::stringstream ss;
  ss << "burton_";
#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
  ss << "extras_";
#endif
  ss << FLECSI_SP_BURTON_MESH_DIMENSION << "d";
  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Update the goemetry
////////////////////////////////////////////////////////////////////////////////
void update_geometry( utils::client_handle_r<mesh_t> mesh ) {
  mesh.update_geometry( );
}

flecsi_register_task(update_geometry, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);


////////////////////////////////////////////////////////////////////////////////
//! \brief Test some initial connectivity
////////////////////////////////////////////////////////////////////////////////
void dump_test( utils::client_handle_r<mesh_t> mesh ) {

  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  
  // create an output file name
  std::stringstream ss;
  ss << prefix() << "_dump_rank" << rank << ".out";

  // dump the mesh to a file
  std::ofstream file( ss.str() );
  mesh.dump( file );
  file.close();

  // now verify its correct
  std::ifstream std_file( ss.str() + ".std" );
  std::ifstream res_file( ss.str() );
  std::istream_iterator<std::string> actual_begin(res_file);
  std::istream_iterator<std::string> expected_begin(std_file);
  std::istream_iterator<std::string> eos;

  CINCH_EXPECT_EQUAL_COLLECTIONS( actual_begin, eos, expected_begin, eos )
    << "DUMP TEST FAILED"; 

}

flecsi_register_task(dump_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief Test validity of mesh
////////////////////////////////////////////////////////////////////////////////
void validity_test( utils::client_handle_r<mesh_t> mesh ) {
  auto ret = mesh.is_valid( false );
  EXPECT_TRUE( ret ) << "VALIDITY TEST FAILED";
}

flecsi_register_task(validity_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);


////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh connectivity functions
////////////////////////////////////////////////////////////////////////////////
void connectivity_test( utils::client_handle_r<mesh_t> mesh ) {
  
  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // create an output file name
  std::stringstream ss;
  ss << prefix() << "_connectivity_rank" << rank << ".out";

  // dump the mesh to a file
  std::ofstream file( ss.str() );

  file << "Vertices in mesh:" << endl;

  for(auto v : mesh.vertices()) {
    file << "----------- vertex id: " << v.id()
      << " with coordinates " << v->coordinates() << endl;
  } // for

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
  file << "Edges in mesh:" << endl;

  for(auto e : mesh.edges()) {
    file << "----------- edge id: " << e.id()
      << " with midpoint " << e->midpoint() << endl;
  } // for

  file << "Faces in mesh:" << endl;

  for(auto f : mesh.faces()) 
    file << "----------- faces id: " << f.id()
       << " with centroid " << f->centroid() << endl;
#endif

  file << "Cells in mesh:" << endl;

  for(auto c : mesh.cells()) {
    file << "----------- cell id: " << c.id()
      << " with centroid " << c->centroid() << endl;
  } // for

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
  file << "Corners in mesh:" << endl;

  for(auto c : mesh.corners()) {
    file << "----------- corner id: " << c.id() << endl;
  } // for

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
  file << "Wedges in mesh:" << endl;

  for(auto c : mesh.wedges()) {
    file << "----------- wedge id: " << c.id() << endl;
  } // for
#endif
#endif

  file << "For each vertex:" << endl;

  for(auto v: mesh.vertices()) {
    file << "^^^^^^^^Vertex id: " << v.id() << endl;

    file << "    ----Cells:" << endl;
    for(auto c: mesh.cells(v))
      file << "    ++++ cell id: " << c.id() << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
    file << "    ----Faces:" << endl;
    for(auto f: mesh.faces(v))
      file << "    ++++ face id: " << f.id() << endl;

    file << "    ----Edges:" << endl;
    for(auto e: mesh.edges(v))
      file << "    ++++ edge id: " << e.id() << endl;
#endif

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    file << "    ----Corners:" << endl;
    for(auto c: mesh.corners(v))
      file << "    ++++ corner id: " << c.id() << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
    file << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(v))
      file << "    ++++ wedge id: " << w.id() << endl;
#endif
#endif

  } // for

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
  file << "For each edge:" << endl;

  for(auto e : mesh.edges()) {
    file << "^^^^^^^^Edge id: " << e.id() << endl;

    file << "    ----Cells:" << endl;
    for(auto c : mesh.cells(e))
      file << "    ++++ cell id: " << c.id() << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 2
    file << "    ----Faces:" << endl;
    for(auto f: mesh.faces(e))
      file << "    ++++ face id: " << f.id() << endl;
#endif

    file << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(e))
      file << "    ++++ vertex id: " << v.id() << endl;

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    file << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(e))
      file << "    ++++ corner id: " << cnr.id() << endl;

    file << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(e))
      file << "    ++++ wedge id: " << w.id() << endl;
#endif

  } // for
#endif  // DIMENSION > 1

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
  file << "For each face:" << endl;

  for(auto f : mesh.faces()) {
    file << "^^^^^^^^Face id: " << f.id() << endl;

    file << "    ----Cells:" << endl;
    for(auto c: mesh.cells(f))
      file << "    ++++ cell id: " << c.id() << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 2
    file << "    ----Edges:" << endl;
    for(auto e : mesh.edges(f))
      file << "    ++++ edge id: " << e.id() << endl;
#endif

    file << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(f))
      file << "    ++++ vertex id: " << v.id() << endl;

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    file << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(f))
      file << "    ++++ corner id: " << cnr.id() << endl;

    file << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(f))
      file << "    ++++ wedge id: " << w.id() << endl;
#endif

  } // for
#endif  // DIMENSION > 1

  file << "For each cell:" << endl;

  for(auto c : mesh.cells()) {
    file << "^^^^^^^^Cell id: " << c.id() << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
    file << "    ----Faces:" << endl;
    for(auto f: mesh.faces(c))
      file << "    ++++ face id: " << f.id() << endl;

    file << "    ----Edges:" << endl;
    for(auto e : mesh.edges(c))
      file << "    ++++ edge id: " << e.id() << endl;
#endif

    file << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(c))
      file << "    ++++ vertex id: " << v.id() << endl;

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    file << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(c))
      file << "    ++++ corner id: " << cnr.id() << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
    file << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(c))
      file << "    ++++ wedge id: " << w.id() << endl;
#endif
#endif

  } // for

  // close file before comparison
  file.close();

  // now verify its correct
  std::ifstream std_file( ss.str() + ".std" );
  std::ifstream res_file( ss.str() );
  std::istream_iterator<std::string> actual_begin(res_file);
  std::istream_iterator<std::string> expected_begin(std_file);
  std::istream_iterator<std::string> eos;

  CINCH_EXPECT_EQUAL_COLLECTIONS( actual_begin, eos, expected_begin, eos )
    << "CONNECTIVITY TEST FAILED";

} // TEST_F

flecsi_register_task(connectivity_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh geometry functions
////////////////////////////////////////////////////////////////////////////////
void geometry_test( utils::client_handle_r<mesh_t> mesh ) {
  
  // get the context
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();

  // create an output file name
  std::stringstream ss;
  ss << prefix() << "_geometry_rank" << rank << ".out";

  // dump the mesh to a file
  std::ofstream file( ss.str() );

  file << "For each cell:" << endl;

  for(auto c: mesh.cells()) {
    auto xc = c->centroid();
    auto vol = c->volume();

    file << "---- cell id: " << c.id()
      << " with centroid " << xc << " and volume " << vol << endl;

    for(auto v : mesh.vertices(c)){
      auto xv = v->coordinates();
      file << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for
  
#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
  file << "For each face:" << endl;

  for(auto f: mesh.faces()) {
    auto xc = f->centroid();
    auto a = f->area();
    auto n = f->normal();

    file << "---- face id: " << f.id()
      << " with midpoint " << xc << ", area " << a
      << " and normal " << n << endl;

    for(auto v : mesh.vertices(f)){
      auto xv = v->coordinates();
      file << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for
#endif  // DIMENSION > 1

  file << "For each edge:" << endl;

  for(auto e: mesh.edges()) {
    auto xc = e->midpoint();
    auto l = e->length();

    file << "---- edge id: " << e.id()
      << " with midpoint " << xc << ", length " << l;

#if FLECSI_SP_BURTON_MESH_DIMENSION < 3
    auto n = e->normal();
    file << " and normal " << n;
#endif
    
    file << endl;

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
    for(auto v : mesh.vertices(e)){
      auto xv = v->coordinates();
      file << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for
#endif

  } // for



#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
  
  file << "For each wedge:" << endl;

  for(auto w: mesh.wedges()) {

    file << "---- wedge id: " << w.id()
      << " with facet normal " << w->facet_normal() 
      << ", facet area " << w->facet_area()
      << endl;
    
  }

#if FLECSI_SP_BURTON_MESH_DIMENSION > 1
  file << "For each corner:" << endl;

  for (auto c: mesh.corners()) {
    
    file << "---- corner id: " << c.id() << endl;

    for (auto w: mesh.wedges(c))  
      file << "++++ wedge id: " << w.id()
          << " with facet normal " << w->facet_normal() 
          << ", facet area " << w->facet_area()
          << endl;

  }

#endif // DIMENSION > 1
#endif // FLECSI_SP_BURTON_MESH_EXTRAS


  // close file before comparison
  file.close();

  // now verify its correct
  std::ifstream std_file( ss.str() + ".std" );
  std::ifstream res_file( ss.str() );
  std::istream_iterator<std::string> actual_begin(res_file);
  std::istream_iterator<std::string> expected_begin(std_file);
  std::istream_iterator<std::string> eos;

  CINCH_EXPECT_EQUAL_COLLECTIONS( actual_begin, eos, expected_begin, eos )
    << "GEOMETRY TEST FAILED";

} // TEST_F

flecsi_register_task(geometry_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh normal functions
////////////////////////////////////////////////////////////////////////////////
void normals_test( utils::client_handle_r<mesh_t> mesh ) {
  
  for(auto f : mesh.faces(flecsi::owned)) {
    auto n = f->normal();
    auto fx = f->centroid();
    auto c = mesh.cells(f).front();
    auto cx = c->midpoint();
    auto delta = fx - cx;
    auto dot = dot_product( n, delta );
    EXPECT_GT( dot, 0 );
  } // for

} // TEST_F

flecsi_register_task(normals_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh subsets
////////////////////////////////////////////////////////////////////////////////
void subset_test( utils::client_handle_r<mesh_t> mesh ) {
  
  std::set< mesh_t::vertex_t* > overlapping_verts;
  const auto vs = mesh.vertices( mesh_t::subset_t::overlapping );

  for(auto c : mesh.cells(flecsi::owned))
  {
    for ( auto v : mesh.vertices(c) ) {
      overlapping_verts.emplace( v );
      auto found = false;
      for ( auto vb : vs ) {
        if ( vb.id() == v.id() ) {
          found = true;
          break;
        }
      }
      ASSERT_TRUE( found ) << "VERTEX NOT FOUND IN SUBSET";
    }
  }

  ASSERT_EQ(
      overlapping_verts.size(),
      mesh.num_vertices(mesh_t::subset_t::overlapping)
  ) << "VERTEX SUBSET SIZE MISMATCH";

} // TEST_F

flecsi_register_task(subset_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);


////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh state
////////////////////////////////////////////////////////////////////////////////
flecsi_register_field(mesh_t, hydro, cell_data, real_t, dense, 1, index_spaces_t::cells);
flecsi_register_field(mesh_t, hydro, face_data, data_t, dense, 1, index_spaces_t::faces);
flecsi_register_field(mesh_t, hydro, edge_data, vector_t, dense, 1, index_spaces_t::edges);
flecsi_register_field(mesh_t, hydro, vert_data, array_t, dense, 1, index_spaces_t::vertices);

void state_fill_test(
  utils::client_handle_r<mesh_t> mesh,
  utils::dense_handle_w<real_t> cell_data,
  utils::dense_handle_w<data_t> face_data,
  utils::dense_handle_w<vector_t> edge_data,
  utils::dense_handle_w<array_t> vert_data
) {
  
  const auto & context = flecsi::execution::context_t::instance();
  auto & verts_lid_to_mid = context.index_map( index_spaces_t::vertices );
  auto & edges_lid_to_mid = context.index_map( index_spaces_t::edges );
  auto & faces_lid_to_mid = context.index_map( index_spaces_t::faces );
  auto & cells_lid_to_mid = context.index_map( index_spaces_t::cells );

  // cells
  for(auto c: mesh.cells(flecsi::owned)) {
    auto id = cells_lid_to_mid.at( c.id() );
    cell_data(c) = id;
  }

  // faces
  for (auto f: mesh.faces(flecsi::owned)) {
    auto id = faces_lid_to_mid.at( f.id() );
    face_data(f).i = id;
    face_data(f).x = id;
  }

  // edges
  for (auto e: mesh.edges(flecsi::owned)) {
    auto id = edges_lid_to_mid.at( e.id() );
    for ( auto & x : edge_data(e) ) x = id;
  }

  // vertices
  for (auto v: mesh.vertices(flecsi::owned)) {
    auto id = verts_lid_to_mid.at( v.id() );
    for ( auto & x : vert_data(v) ) x = id;
  }

} // TEST_F

void state_check_test(
  utils::client_handle_r<mesh_t> mesh,
  utils::dense_handle_r<real_t> cell_data,
  utils::dense_handle_r<data_t> face_data,
  utils::dense_handle_r<vector_t> edge_data,
  utils::dense_handle_r<array_t> vert_data
) {
  
  const auto & context = flecsi::execution::context_t::instance();
  auto & verts_lid_to_mid = context.index_map( index_spaces_t::vertices );
  auto & edges_lid_to_mid = context.index_map( index_spaces_t::edges );
  auto & faces_lid_to_mid = context.index_map( index_spaces_t::faces );
  auto & cells_lid_to_mid = context.index_map( index_spaces_t::cells );

  // cells
  for(auto c: mesh.cells()) {
    auto id = cells_lid_to_mid.at( c.id() );
    ASSERT_EQ( cell_data(c), id ) << "CELL DATA NOT COMMUNICATED PROPERLY";
  }

  // faces
  for (auto f: mesh.faces()) {
    auto id = faces_lid_to_mid.at( f.id() );
    ASSERT_EQ( face_data(f).i, id ) << "FACE DATA NOT COMMUNICATED PROPERLY";
    ASSERT_EQ( face_data(f).x, id ) << "FACE DATA NOT COMMUNICATED PROPERLY";
  }

  // edges
  for (auto e: mesh.edges()) {
    auto id = edges_lid_to_mid.at( e.id() );
    for ( auto & x : edge_data(e) )
      ASSERT_EQ(x, id) << "EDGE DATA NOT COMMUNICATED PROPERLY";
  }

  // vertices
  for (auto v: mesh.vertices()) {
    auto id = verts_lid_to_mid.at( v.id() );
    for ( auto & x : vert_data(v) )
      ASSERT_EQ(x, id) << "VERTEX DATA NOT COMMUNICATED PROPERLY";
  }

} // TEST_F

flecsi_register_task(state_fill_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

flecsi_register_task(state_check_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS


////////////////////////////////////////////////////////////////////////////////
//! \brief test the classification of wedges and corners
////////////////////////////////////////////////////////////////////////////////
void extras_test( utils::client_handle_r<mesh_t> mesh ) {
  
  auto crns_excl = mesh.corners(flecsi::exclusive);

  for (auto cl : mesh.cells(flecsi::exclusive)) {
    for (auto cn : mesh.corners(cl) ) {
      auto found = false;
      for ( auto cn_excl : crns_excl )
        if ( cn_excl.id() == cn.id() ) {
          found = true;
          break;
        }
      ASSERT_TRUE( found ) << "CORNER NOT FOUND IN SUBSET";
    }
  } // for

} // TEST_F

flecsi_register_task(extras_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh state at corners and wedges
////////////////////////////////////////////////////////////////////////////////
flecsi_register_field(mesh_t, hydro, wedge_data, int, dense, 1, index_spaces_t::wedges);
flecsi_register_field(mesh_t, hydro, cornr_data, double, dense, 1, index_spaces_t::corners);

void extra_state_fill_test(
  utils::client_handle_r<mesh_t> mesh,
  utils::dense_handle_w<int> wedge_data,
  utils::dense_handle_w<double> corner_data
) {
  
  const auto & context = flecsi::execution::context_t::instance();
  auto & wedges_lid_to_mid = context.index_map( index_spaces_t::wedges );
  auto & corners_lid_to_mid = context.index_map( index_spaces_t::corners );

  // corners
  for(auto c: mesh.corners(flecsi::owned)) {
    auto id = corners_lid_to_mid.at( c.id() );
    corner_data(c) = id;
  }

  // wedges
  for (auto w: mesh.wedges(flecsi::owned)) {
    auto id = wedges_lid_to_mid.at( w.id() );
    wedge_data(w) = id;
  }

} // TEST_F

void extra_state_check_test(
  utils::client_handle_r<mesh_t> mesh,
  utils::dense_handle_r<int> wedge_data,
  utils::dense_handle_r<double> corner_data
) {
  
  const auto & context = flecsi::execution::context_t::instance();
  auto & wedges_lid_to_mid = context.index_map( index_spaces_t::wedges );
  auto & corners_lid_to_mid = context.index_map( index_spaces_t::corners );

  // corners
  for(auto c: mesh.corners()) {
    auto id = corners_lid_to_mid.at( c.id() );
    ASSERT_EQ( corner_data(c), id ) << "CORNER DATA NOT COMMUNICATED PROPERLY";
  }

  // wedges
  for (auto w: mesh.wedges()) {
    auto id = wedges_lid_to_mid.at( w.id() );
    ASSERT_EQ( wedge_data(w), id ) << "WEDGE DATA NOT COMMUNICATED PROPERLY";
  }

} // TEST_F

flecsi_register_task(extra_state_fill_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

flecsi_register_task(extra_state_check_test, flecsi_sp::burton::test, loc,
    index|flecsi::leaf);

#endif

} // namespace
} // namespace
} // namespace

namespace flecsi {
namespace execution {



////////////////////////////////////////////////////////////////////////////////
//! \brief the driver for all tests
////////////////////////////////////////////////////////////////////////////////
void driver(int argc, char ** argv)
{

  // get the mesh handle
  auto mesh_handle = flecsi_get_client_handle(mesh_t, meshes, mesh0);
  auto f = 
    flecsi_execute_task(update_geometry, flecsi_sp::burton::test, index, mesh_handle);
  f.wait(); // dont go forward until this is done!

  // launch the validity test task
  flecsi_execute_task(validity_test, flecsi_sp::burton::test, index, mesh_handle);
  
  // launch the dump test task
  flecsi_execute_task(dump_test, flecsi_sp::burton::test, index, mesh_handle);
  
  // launch the connectivity test task
  flecsi_execute_task(connectivity_test, flecsi_sp::burton::test, index, mesh_handle);

  // launch the geometry test task
  flecsi_execute_task(geometry_test, flecsi_sp::burton::test, index, mesh_handle);

  // launch the normals test task
  flecsi_execute_task(normals_test, flecsi_sp::burton::test, index, mesh_handle);
  
  // launch the connectivity test task
  flecsi_execute_task(subset_test, flecsi_sp::burton::test, index, mesh_handle);

  // launch the state test task
  auto cell_data_handle = flecsi_get_handle(mesh_handle, hydro, cell_data, real_t, dense, 0);
  auto face_data_handle = flecsi_get_handle(mesh_handle, hydro, face_data, data_t, dense, 0);
  auto edge_data_handle = flecsi_get_handle(mesh_handle, hydro, edge_data, vector_t, dense, 0);
  auto vert_data_handle = flecsi_get_handle(mesh_handle, hydro, vert_data, array_t, dense, 0);

  flecsi_execute_task(
      state_fill_test,
      flecsi_sp::burton::test,
      index,
      mesh_handle,
      cell_data_handle,
      face_data_handle,
      edge_data_handle,
      vert_data_handle );

  flecsi_execute_task(
      state_check_test,
      flecsi_sp::burton::test,
      index,
      mesh_handle,
      cell_data_handle,
      face_data_handle,
      edge_data_handle,
      vert_data_handle );


#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
  
  // launch the extras test task
  flecsi_execute_task(extras_test, flecsi_sp::burton::test, index, mesh_handle);

  // now do the same for the corners and wedges
  auto wedge_data_handle = flecsi_get_handle(mesh_handle, hydro, wedge_data, int, dense, 0);
  auto corner_data_handle = flecsi_get_handle(mesh_handle, hydro, cornr_data, double, dense, 0);

  flecsi_execute_task(
      extra_state_fill_test,
      flecsi_sp::burton::test,
      index,
      mesh_handle,
      wedge_data_handle,
      corner_data_handle );

  flecsi_execute_task(
      extra_state_check_test,
      flecsi_sp::burton::test,
      index,
      mesh_handle,
      wedge_data_handle,
      corner_data_handle );

#endif

} // driver

} // namespace execution
} // namespace flecsi

////////////////////////////////////////////////////////////////////////////////
//! \brief Only here so test runs?
////////////////////////////////////////////////////////////////////////////////
TEST(burton, simple) {}

