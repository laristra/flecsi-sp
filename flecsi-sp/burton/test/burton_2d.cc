/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
// \file
// \brief Tests general features of the burton mesh.
////////////////////////////////////////////////////////////////////////////////

// user includes
#include "burton_2d_test.h"

#include <flecsi/execution/execution.h>
#include <flecsi-sp/utils/types.h>

// using statements
using std::cout;
using std::endl;

using mesh_t = flecsi_sp::burton::test::burton_2d::mesh_t;

namespace flecsi_sp {
namespace burton {
namespace test {

void dump_mesh( utils::client_handle_r__<mesh_t> mesh) {
  std::cout << "dumping mesh" << std::endl;
  mesh.dump( CINCH_CAPTURE() );
}


flecsi_register_task(dump_mesh, flecsi_sp::burton::test, loc,
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

  auto mesh_handle = flecsi_get_client_handle(mesh_t, meshes, mesh0);
  
  flecsi_execute_task(dump_mesh, flecsi_sp::burton::test, single, mesh_handle);

  //std::string outfile = output_prefix()+".blessed";
  //CINCH_ASSERT(TRUE, CINCH_EQUAL_BLESSED( outfile.c_str() ));
} // driver

} // namespace execution
} // namespace flecsi


using namespace flecsi_sp::burton::test;


////////////////////////////////////////////////////////////////////////////////
//! \brief Only here so test runs?
////////////////////////////////////////////////////////////////////////////////
TEST_F(burton_2d, simple) {}

