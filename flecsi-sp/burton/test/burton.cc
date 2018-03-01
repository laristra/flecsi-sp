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
#include <flecsi-sp/utils/types.h>

// using statements
using std::cout;
using std::endl;

using mesh_t = flecsi_sp::burton::burton_mesh_t;
using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;

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
//! \brief Test some initial connectivity
////////////////////////////////////////////////////////////////////////////////
void dump_test( utils::client_handle_r__<mesh_t> mesh ) {

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

  CINCH_EXPECT_EQUAL_COLLECTIONS( actual_begin, eos, expected_begin, eos );

}

flecsi_register_task(dump_test, flecsi_sp::burton::test, loc,
    single|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh connectivity functions
////////////////////////////////////////////////////////////////////////////////
void connectivity_test( utils::client_handle_r__<mesh_t> mesh ) {
  
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

  file << "Edges in mesh:" << endl;

  for(auto e : mesh.edges()) {
    file << "----------- edge id: " << e.id()
      << " with midpoint " << e->midpoint() << endl;
  } // for

  file << "Faces in mesh:" << endl;

  for(auto f : mesh.faces()) 
    file << "----------- faces id: " << f.id()
       << " with centroid " << f->centroid() << endl;

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

  file << "Wedges in mesh:" << endl;

  for(auto c : mesh.wedges()) {
    file << "----------- wedge id: " << c.id() << endl;
  } // for
#endif

  file << "For each vertex:" << endl;

  for(auto v: mesh.vertices()) {
    file << "^^^^^^^^Vertex id: " << v.id() << endl;

    file << "    ----Cells:" << endl;
    for(auto c: mesh.cells(v))
      file << "    ++++ cell id: " << c.id() << endl;

    file << "    ----Faces:" << endl;
    for(auto f: mesh.faces(v))
      file << "    ++++ face id: " << f.id() << endl;

    file << "    ----Edges:" << endl;
    for(auto e: mesh.edges(v))
      file << "    ++++ edge id: " << e.id() << endl;

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    file << "    ----Corners:" << endl;
    for(auto c: mesh.corners(v))
      file << "    ++++ corner id: " << c.id() << endl;

    file << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(v))
			file << "    ++++ wedge id: " << w.id() << endl;
#endif

  } // for

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

  file << "For each cell:" << endl;

  for(auto c : mesh.cells()) {
    file << "^^^^^^^^Cell id: " << c.id() << endl;

    file << "    ----Faces:" << endl;
    for(auto f: mesh.faces(c))
      file << "    ++++ face id: " << f.id() << endl;

    file << "    ----Edges:" << endl;
    for(auto e : mesh.edges(c))
      file << "    ++++ edge id: " << e.id() << endl;

    file << "    ----Vertices:" << endl;
    for(auto v : mesh.vertices(c))
      file << "    ++++ vertex id: " << v.id() << endl;

#ifdef FLECSI_SP_BURTON_MESH_EXTRAS
    file << "    ----Corners:" << endl;
    for(auto cnr : mesh.corners(c))
      file << "    ++++ corner id: " << cnr.id() << endl;

    file << "    ----Wedges:" << endl;
    for(auto w: mesh.wedges(c))
			file << "    ++++ wedge id: " << w.id() << endl;
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

  CINCH_EXPECT_EQUAL_COLLECTIONS( actual_begin, eos, expected_begin, eos );

} // TEST_F

flecsi_register_task(connectivity_test, flecsi_sp::burton::test, loc,
    single|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh geometry functions
////////////////////////////////////////////////////////////////////////////////
void geometry_test( utils::client_handle_r__<mesh_t> mesh ) {
  
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

  file << "For each edge:" << endl;

  for(auto e: mesh.edges()) {
    auto xc = e->midpoint();
    auto l = e->length();

    file << "---- edge id: " << e.id()
      << " with midpoint " << xc << ", length " << l;

#if FLECSI_SP_BURTON_MESH_DIMENSION == 2
    auto n = e->normal();
    file << " and normal " << n;
#endif
    
    file << endl;

    for(auto v : mesh.vertices(e)){
      auto xv = v->coordinates();
      file << "++++ vertex id: " << v.id()
        << " with coordinates " << xv << endl;
    } // for

  } // for

  // close file before comparison
  file.close();

  // now verify its correct
  std::ifstream std_file( ss.str() + ".std" );
  std::ifstream res_file( ss.str() );
  std::istream_iterator<std::string> actual_begin(res_file);
  std::istream_iterator<std::string> expected_begin(std_file);
  std::istream_iterator<std::string> eos;

  CINCH_EXPECT_EQUAL_COLLECTIONS( actual_begin, eos, expected_begin, eos );

} // TEST_F

flecsi_register_task(geometry_test, flecsi_sp::burton::test, loc,
    single|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh normal functions
////////////////////////////////////////////////////////////////////////////////
void normals_test( utils::client_handle_r__<mesh_t> mesh ) {
  
  for(auto f : mesh.faces(flecsi::owned)) {
    auto n = f->normal();
    auto fx = f->centroid();
    auto c = mesh.cells(f).front();
    auto cx = c->centroid();
    auto delta = fx - cx;
    auto dot = dot_product( n, delta );
    ASSERT_GT( dot, 0 );
  } // for

} // TEST_F

flecsi_register_task(normals_test, flecsi_sp::burton::test, loc,
    single|flecsi::leaf);

// Data registration for the tests
flecsi_register_field(mesh_t, hydro, pressure, real_t, dense, 1, mesh_t::index_spaces_t::cells);
flecsi_register_field(mesh_t, hydro, velocity, vector_t, dense, 1, mesh_t::index_spaces_t::vertices);
flecsi_register_field(mesh_t, hydro, H, vector_t, dense, 1, mesh_t::index_spaces_t::edges);

////////////////////////////////////////////////////////////////////////////////
//! \brief test the mesh state
////////////////////////////////////////////////////////////////////////////////
void state_test(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_r__<real_t> press,
  utils::dense_handle_rw__<vector_t> vel,
  utils::dense_handle_rw__<vector_t> H
) {

  // cells
  for(auto c: mesh.cells()) {
    press(c) = c.id();
  } // for

  for(auto c: mesh.cells()) {
    ASSERT_EQ(c.id(), press(c));
  } // for

  // vertices
  for (auto v: mesh.vertices()) {
    vel(v)[0] = v.id();
    vel(v)[1] = 2.0*v.id();
  } // for

  for (auto v: mesh.vertices()) {
    ASSERT_EQ(v.id(), vel(v)[0]);
    ASSERT_EQ(2.0*v.id(), vel(v)[1]);
  } // for

  // edges
  for (auto e: mesh.edges()) {
    H(e)[0] = e.id()*e.id();
    H(e)[1] = e.id()*e.id()*e.id();
  } // for

  for (auto e: mesh.edges()) {
    ASSERT_EQ(e.id()*e.id(), H(e)[0]);
    ASSERT_EQ(e.id()*e.id()*e.id(), H(e)[1]);
  } // for

} // TEST_F

flecsi_register_task(state_test, flecsi_sp::burton::test, loc,
    single|flecsi::leaf);

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

  // launch the dump test task
  flecsi_execute_task(dump_test, flecsi_sp::burton::test, single, mesh_handle);
  
  // launch the connectivity test task
  flecsi_execute_task(connectivity_test, flecsi_sp::burton::test, single, mesh_handle);

  // launch the geometry test task
  flecsi_execute_task(geometry_test, flecsi_sp::burton::test, single, mesh_handle);

  // launch the normals test task
  flecsi_execute_task(normals_test, flecsi_sp::burton::test, single, mesh_handle);

  // launch the state test task
  auto p_handle = flecsi_get_handle(mesh_handle, hydro, pressure, real_t, dense, 0);
  auto v_handle = flecsi_get_handle(mesh_handle, hydro, velocity, vector_t, dense, 0);
  auto H_handle = flecsi_get_handle(mesh_handle, hydro, H, vector_t, dense, 0);
  flecsi_execute_task(state_test, flecsi_sp::burton::test, single, mesh_handle,
    p_handle, v_handle, H_handle);

} // driver

} // namespace execution
} // namespace flecsi

////////////////////////////////////////////////////////////////////////////////
//! \brief Only here so test runs?
////////////////////////////////////////////////////////////////////////////////
TEST(burton, simple) {}

