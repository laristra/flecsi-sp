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
using index_spaces_t = mesh_t::index_spaces_t;
using integer_t = mesh_t::integer_t;

namespace flecsi_sp {
namespace burton {
namespace test {

// Data registration for the tests
flecsi_register_field(mesh_t, hydro, face_field, integer_t, dense, 1, index_spaces_t::faces);

////////////////////////////////////////////////////////////////////////////////
//! \brief write to the field
////////////////////////////////////////////////////////////////////////////////
void field_write_test(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_w__<integer_t> face_field
) {

  const auto & context = flecsi::execution::context_t::instance();
  auto & lid_to_mid = context.index_map( index_spaces_t::faces );

  for ( auto f : mesh.faces() )
  {
    face_field(f) = lid_to_mid.at( f.id() );
  }
    
} // TEST_F

flecsi_register_task(field_write_test, flecsi_sp::burton::test, loc,
    single|flecsi::leaf);

////////////////////////////////////////////////////////////////////////////////
//! \brief read from the field
////////////////////////////////////////////////////////////////////////////////
void field_read_test(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_r__<integer_t> face_field
) {
  
  const auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  auto & lid_to_mid = context.index_map( index_spaces_t::faces );

  for ( auto f : mesh.faces() ) 
  {
    std::cout << "[rank:" << rank << "] lid=" << f.id() << " -> mid=" << lid_to_mid.at( f.id() ) << std::endl;
    EXPECT_EQ( face_field(f), lid_to_mid.at( f.id() ) );
  }

} // TEST_F

flecsi_register_task(field_read_test, flecsi_sp::burton::test, loc,
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
  auto face_field_handle = flecsi_get_handle(mesh_handle, hydro, face_field, integer_t, dense, 0);

  // test the field communication
  flecsi_execute_task(field_write_test, flecsi_sp::burton::test, single, mesh_handle, face_field_handle);
  flecsi_execute_task(field_read_test, flecsi_sp::burton::test, single, mesh_handle, face_field_handle);

} // driver

} // namespace execution
} // namespace flecsi

////////////////////////////////////////////////////////////////////////////////
//! \brief Only here so test runs?
////////////////////////////////////////////////////////////////////////////////
TEST(burton, simple) {}

