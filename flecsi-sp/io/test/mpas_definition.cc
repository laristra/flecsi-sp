/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

// the main test include
#include <cinchtest.h>

// user includes
#include <flecsi/topology/closure_utils.h>
#include <flecsi-sp/io/mpas_definition.h>

// the MPAS mesh files hold native doubles
using mpas_def_t = 
  flecsi_sp::io::mpas_definition__<double>;

// The flecsi library has undefined symbols in it.  It calls 
// flecsi::execution::driver even though we are not using the execution model.
// Define an empty stub for linkage.
namespace flecsi {
namespace execution {
void driver(int argc, char ** argv) {}
}
}


////////////////////////////////////////////////////////////////////////////////
/// \brief This tests simple reading of the MPAS mesh
////////////////////////////////////////////////////////////////////////////////
TEST(mpas_definition, simple) {
  mpas_def_t mesh("init.h5");
} // TEST

