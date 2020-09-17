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

#include <ristra/initialization/arguments.h>
#include <ristra/initialization/input.h>

// flecsi insitu
#ifdef FLECSI_SP_ENABLE_CATALYST
#include "ristra/io/catalyst/utils/catalystUtils.hpp"
#endif //FLECSI_SP_ENABLE_CATALYST

// system includes
#include <array>
#include <string>
#include <assert.h> 
#include <chrono>
#include <thread>

// using statements
using std::cout;
using std::endl;

namespace flecsi_sp {
namespace burton {
namespace test {


////////////////////////////////////////////////////////////////////////////////
//! \brief Insitu Catalyst
////////////////////////////////////////////////////////////////////////////////

#ifdef FLECSI_SP_ENABLE_CATALYST

bool catalyst_validation(std::string inputFilename)
{
  clog(info) << "Catalyst validation" << std::endl;
  std::cout << "Catalyst validation" << std::endl;

  // Allow catalyst to write the file to disk: sleep 2 seconds
  std::this_thread::sleep_for(std::chrono::milliseconds(2000));

  //
  // Do image validation
  bool imagesSimilar = false;
  {
    // Get image baseline
    std::string baseline = replaceInString(inputFilename, "catalyst_scripts/catalyst-scripts.json", "catalyst_scripts/CatalystBaseline.png");

    // Get the image result from catalyst
    std::string catalystOutput = "views/RenderView1_0.png";

    // Compare images
    imagesSimilar = imageCompare(baseline, catalystOutput);
  }

  if (!imagesSimilar)
    return false;


  //
  // Do file validation
  bool filesSimilar = false;
  {
    // pvtu file
    {
      std::string baseline = replaceInString(inputFilename, "catalyst_scripts/catalyst-scripts.json", "catalyst_scripts/simOutagewww_0000.pvtu");
      std::string catalystPVTUOutput = "simOutagewww_0000.pvtu";

      filesSimilar = compareFiles(baseline, catalystPVTUOutput);
      std::cout << "filesSimilar: " << filesSimilar << std::endl;

      if (!filesSimilar)
        return false;
    }

    /*
    // rank 0 vtu
    {
      std::string baseline = replaceInString(inputFilename, "catalyst_scripts/catalyst-scripts.json", "catalyst_scripts/simOutagewww_0000/simOutagewww_0000_0.vtu");
      //std::cout << "Baseline vtu1: " << baseline << std::endl;

      std::string catalystVTUOutput = "simOutagewww_0000/simOutagewww_0000_0.vtu";
      //std::cout << "Catalyst vtu file: " << catalystVTUOutput << std::endl;

      filesSimilar = compareFiles(baseline, catalystVTUOutput);
      std::cout << "filesSimilar: " << filesSimilar << std::endl;

      if (!filesSimilar)
        return false;
    }

    // rank 1 vtu
    {
      std::string baseline = replaceInString(inputFilename, "catalyst_scripts/catalyst-scripts.json", "catalyst_scripts/simOutagewww_0000/simOutagewww_0000_1.vtu");
      //std::cout << "Baseline vtu1: " << baseline << std::endl;

      std::string catalystVTUOutput = "simOutagewww_0000/simOutagewww_0000_1.vtu";
      //std::cout << "Catalyst vtu file: " << catalystVTUOutput << std::endl;

      filesSimilar = compareFiles(baseline, catalystVTUOutput);
      std::cout << "filesSimilar: " << filesSimilar << std::endl;

      if (!filesSimilar)
        return false;
    }
    */
  } 
  
  return filesSimilar;
}

flecsi_register_task(catalyst_validation, flecsi_sp::burton::test, loc, index|flecsi::leaf);

#endif //FLECSI_SP_ENABLE_CATALYST


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
	std::cout << "Run driver" << std::endl;
#ifdef FLECSI_SP_ENABLE_CATALYST
  auto & cmdl = ristra::initialization::command_line_arguments_t::instance();
  auto args = cmdl.process_arguments(argc, argv);
  const auto & variables = args.variables;
  auto inputFilename = variables.as<std::string>("catalyst-file");
  std::cout << "inputFilename: " << inputFilename << std::endl;


  // Do validation
  auto catalystValidate = flecsi_execute_task(catalyst_validation, flecsi_sp::burton::test, index,   inputFilename);
  catalystValidate.wait();  
  bool isValid = catalystValidate.get();

  assert(isValid);

#endif //FLECSI_SP_ENABLE_CATALYST


} // driver

} // namespace execution
} // namespace flecsi

////////////////////////////////////////////////////////////////////////////////
//! \brief Only here so test runs?
////////////////////////////////////////////////////////////////////////////////
TEST(burton, simple) {}