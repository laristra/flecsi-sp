/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include <cinchlog.h>
#include <flecsi-sp/burton/burton_specialization_init.h>
#include <ristra/initialization/arguments.h>

using distribution_alg_t = flecsi_sp::burton::distribution_alg_t;
using partition_alg_t = flecsi_sp::burton::partition_alg_t;

namespace flecsi {
namespace execution {
 
// register command line arguments
auto register_mesh_args = 
  ristra::initialization::command_line_arguments_t::instance().register_group(
    "mesh", "Mesh specialization parameters" );

auto register_mesh_args_file = 
  ristra::initialization::command_line_arguments_t::instance().
    register_argument<std::string>("mesh", "mesh-file,m", "Specify mesh file" );

auto register_mesh_args_entries = 
  ristra::initialization::command_line_arguments_t::instance().
    register_argument<int>( "mesh", "max-entries,e",
        "Specify the maximum number of sparse entries" ); 

auto register_part_args =
  ristra::initialization::command_line_arguments_t::instance().
    register_argument<int>( "mesh", "partition-only,p",
        "Partition mesh and exit" );

auto register_dist_args =
  ristra::initialization::command_line_arguments_t::instance().
    register_argument<std::string>( "mesh", "rank-distribution",
        "How to distribute mesh files among ranks when using M-to-N partitioning"
        "  from pre-partitioned files. Options include: sequential (default), balanced,"
        " hostname.");

auto register_partition_args =
  ristra::initialization::command_line_arguments_t::instance().
    register_argument<std::string>( "mesh", "partition-algorithm",
        "Partition algorithm.  Options include: kway (default), geom, geomkway, "
        "refinekway, metis, naive.");

auto register_repart_args =
  ristra::initialization::command_line_arguments_t::instance().
  register_argument( "mesh", "repartition",
      "Refine partitioned mesh.");

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver.
///////////////////////////////////////////////////////////////////////////////
void specialization_tlt_init(int argc, char** argv)
{
#ifdef BURTON_ENABLE_APPLICATION_TLT_INIT
  application_tlt_init(argc, argv);
#endif // BURTON_ENABLE_APPLICATION_TLT_INIT
  clog(info) << "In specialization top-level-task init" << std::endl;

  // get the color
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  auto size = context.colors();

  //===========================================================================
  // Parse arguments
  //===========================================================================

  // check command line arguments
  auto & cmdl = ristra::initialization::command_line_arguments_t::instance();
  auto args = cmdl.process_arguments(argc, argv);
  const auto & variables = args.variables;

  if ( variables.count("help") ) exit(0);
  if ( args.error ) exit(-1);

  // get the input file
  auto mesh_filename_string = variables.as<std::string>("mesh-file");

  // override any inputs if need be
  if ( !mesh_filename_string.empty() ) {
    if ( rank == 0 )
      std::cout << "Using mesh file \"" << mesh_filename_string << "\"."
                << std::endl;
  }
  else {
    THROW_RUNTIME_ERROR( "No mesh file provided" );
  }
  
  // get the maximum number of entries
  std::size_t max_entries = variables.as<int>("max-entries", 5);
  if ( variables.count("max-entries") && rank == 0 )
      std::cout << "Setting max_entries to \"" << max_entries << "\"." << std::endl;
  
  int partition_only = variables.as<int>("partition-only", 0);
  if ( variables.count("partition-only") && rank == 0 ) {
    if (partition_only < size) {
      THROW_RUNTIME_ERROR(
          "Number of mpi ranks must be less than or equal to number "
          "of partitions");
    }
    std::cout << "Partitioning mesh into \"" << partition_only << "\" pieces." << std::endl;
  }

  auto distribution_str = variables.as<std::string>("rank-distribution", "sequential");
  
  distribution_alg_t distribution_alg;
  if (distribution_str ==  "balanced")
    distribution_alg = distribution_alg_t::balanced;
  else if (distribution_str ==  "sequential")
    distribution_alg = distribution_alg_t::sequential;
  else if (distribution_str ==  "hostname")
    distribution_alg = distribution_alg_t::hostname;
  else
    THROW_RUNTIME_ERROR("Unknown partition distribution type '" << distribution_str << "'");

  auto partition_str = variables.as<std::string>("partition-algorithm", "kway");
  
  partition_alg_t partition_alg;
  if (partition_str ==  "kway")
    partition_alg = partition_alg_t::kway;
  else if (partition_str ==  "geom")
    partition_alg = partition_alg_t::geom;
  else if (partition_str ==  "geomkway")
    partition_alg = partition_alg_t::geomkway;
  else if (partition_str ==  "refinekway")
    partition_alg = partition_alg_t::refinekway;
  else if (partition_str ==  "metis")
    partition_alg = partition_alg_t::metis;
  else if (partition_str ==  "naive")
    partition_alg = partition_alg_t::naive;
  else
    THROW_RUNTIME_ERROR("Unknown partition algorithm type '" << partition_str << "'");

  auto do_repart = variables.count("repartition");

  //===========================================================================
  // Partition mesh
  //===========================================================================

  clog(info) << "Partitioning mesh" << std::endl;

  // need to put the filename into a statically sized character array
  auto mesh_filename = flecsi_sp::utils::to_char_array( mesh_filename_string );

  // execute the mpi task to partition the mesh
  flecsi_execute_mpi_task(partition_mesh, flecsi_sp::burton, mesh_filename,
    max_entries, partition_only, distribution_alg, partition_alg, do_repart);

  if (partition_only) exit(0);
}

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver.
///////////////////////////////////////////////////////////////////////////////
void specialization_spmd_init(int argc, char** argv)
{
  clog(info) << "In specialization spimd init" << std::endl;

  // get a mesh handle and call the initialization task
  auto mesh_handle = flecsi_get_client_handle(
      flecsi_sp::burton::burton_mesh_t, meshes, mesh0);
  
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  flecsi_execute_task(initialize_mesh, flecsi_sp::burton, index, mesh_handle);
}

} // namespace
} // namespace
