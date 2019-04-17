/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Simple tests related to solving full hydro solutions.
///////////////////////////////////////////////////////////////////////////////

// hydro includes
#include <cinch/logging/cinchlog.h>
#include <flecsi-sp/burton/burton_specialization_init.h>
#include <flecsi-sp/burton/burton_specialization_arguments.h>
#include <ristra/utils/string_utils.h>

namespace flecsi {
namespace execution {



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

  //===========================================================================
  // Parse arguments
  //===========================================================================

  auto args = flecsi_sp::burton::process_arguments( argc, argv );

  // process the simple ones
  if ( args.count("h") )
    return;

  // get the input file
  auto mesh_filename_string =
    args.count("m") ? args.at("m") : std::string();

  // override any inputs if need be
  if ( !mesh_filename_string.empty() ) {
    if ( rank == 0 )
      std::cout << "Using mesh file \"" << mesh_filename_string << "\"."
                << std::endl;
  }
  else {
    throw_runtime_error( "No mesh file provided" );
  }
  
  // get the maximum number of entries
  std::size_t max_entries =
    args.count("e") ? std::stoi(args.at("e")) : 5;
  if ( args.count("e") && rank == 0 )
      std::cout << "Setting max_entries to \"" << max_entries << "\"." << std::endl;

  //===========================================================================
  // Partition mesh
  //===========================================================================  

  clog(info) << "Partitioning mesh" << std::endl;

  // need to put the filename into a statically sized character array
  auto mesh_filename = flecsi_sp::utils::to_char_array( mesh_filename_string );

  // io definition params
  // using real_t = flecsi_sp::burton::burton_mesh_t::real_t;
  // constexpr auto num_dims = flecsi_sp::burton::burton_mesh_t::num_dimensions;
  // using exo_def_t = flecsi_sp::io::exodus_definition__<num_dims, real_t>;
  // using mpas_def_t = flecsi_sp::io::mpas_definition_u<real_t>;

  auto extension = ristra::utils::file_extension(mesh_filename_string);

  // execute the mpi task to partition the mesh
  if(extension == "exo" || extension == "g") {
    flecsi_execute_mpi_task(partition_exo_mesh, flecsi_sp::burton,
                            mesh_filename, max_entries);
  }
  else if(extension == "h5") {
    flecsi_execute_mpi_task(partition_mpas_mesh, flecsi_sp::burton,
                            mesh_filename, max_entries);
  }
  else {
    throw_runtime_error("unrecognized file extension: " << extension);
  } 

}

///////////////////////////////////////////////////////////////////////////////
//! \brief The specialization initialization driver.
///////////////////////////////////////////////////////////////////////////////
void specialization_spmd_init(int argc, char** argv)
{
  clog(info) << "In specialization spimd init" << std::endl;

  //===========================================================================
  // Parse arguments
  //===========================================================================

  auto args = flecsi_sp::burton::process_arguments( argc, argv );
  // Assume arguments are sanitized

  // get the input file
  auto mesh_filename_string =
    args.count("m") ? args.at("m") : std::string();

  // need to put the filename into a statically sized character array
  auto mesh_filename = flecsi_sp::utils::to_char_array( mesh_filename_string );

  //===========================================================================
  // Load the mesh
  //===========================================================================

  // get a mesh handle and call the initialization task
  auto mesh_handle = flecsi_get_client_handle(
      flecsi_sp::burton::burton_mesh_t, meshes, mesh0);

  auto extension = ristra::utils::file_extension(mesh_filename_string);

  // execute the mpi task to partition the mesh
  if(extension == "exo" || extension == "g") {
    flecsi_execute_task(initialize_exo_mesh, flecsi_sp::burton, index,
                        mesh_handle, mesh_filename);
  }
  else if(extension == "h5") {
    flecsi_execute_task(initialize_mpas_mesh, flecsi_sp::burton, index,
                        mesh_handle, mesh_filename);
  }
  else {
    throw_runtime_error("unrecognized file extension: " << extension);
  }
}

} // namespace
} // namespace
