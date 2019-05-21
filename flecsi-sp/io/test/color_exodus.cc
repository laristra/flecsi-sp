
/*~-------------------------------------------------------------------------~~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/

#include <cinchtest.h>
#include <flecsi-sp-config.h>

#include <mpi.h>

#include <array>
#include <fstream>

#include <flecsi/coloring/dcrs_utils.h>
#include <flecsi/coloring/mpi_communicator.h>
#include <flecsi/coloring/parmetis_colorer.h>
#include <flecsi/topology/closure_utils.h>
#include <flecsi/utils/set_utils.h>
#include <flecsi-sp/io/exodus_definition.h>

clog_register_tag(coloring);

// The flecsi library has undefined symbols in it.  It calls 
// flecsi::execution::driver even though we are not using the execution model.
// Define an empty stub for linkage.
namespace flecsi {
namespace execution {
void driver(int argc, char ** argv) {}
}
}


// some type aliases
using exodus_definition_2d_t = flecsi_sp::io::exodus_definition<2, double>;

using exodus_definition_3d_t = flecsi_sp::io::exodus_definition<3, double>;

using std::vector;
  
// test constant
constexpr auto thru_dim = 0;

////////////////////////////////////////////////////////////////////////////////
/// \brief the main cell coloring driver
////////////////////////////////////////////////////////////////////////////////
template<int THRU_DIM, typename MD>
void
color_cells(const MD & md, const std::string & output_prefix) {

  constexpr auto cell_dim = MD::dimension();
  constexpr auto vertex_dim = 0;

  // Set the output rank
  clog_set_output_rank(1);

  using entry_info_t = flecsi::coloring::entity_info_t;

  int comm_size;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create a communicator instance to get neighbor information.
  auto communicator = std::make_shared<flecsi::coloring::mpi_communicator_t>();

  //----------------------------------------------------------------------------
  // Cell Coloring
  //----------------------------------------------------------------------------

  // Create the dCRS representation for the distributed colorer.
  // This essentialy makes the graph of the dual mesh.
  auto dcrs = flecsi::coloring::make_dcrs(md);

  // Create a colorer instance to generate the primary coloring.
  auto colorer = std::make_shared<flecsi::coloring::parmetis_colorer_t>();

  // Create the primary coloring.
  auto primary = colorer->color(dcrs);

  {
    clog_tag_guard(coloring);
    clog_container_one(info, "primary coloring", primary, clog::space);
  } // guard

  //----------------------------------------------------------------------------
  // Cell Closure.  However many layers of ghost cells are needed are found
  // here.
  //----------------------------------------------------------------------------

  // Compute the dependency closure of the primary cell coloring
  // through vertex intersections (specified by last argument "1").
  // To specify edge or face intersections, use 2 (edges) or 3 (faces).
  auto closure =
      flecsi::topology::entity_neighbors<cell_dim, cell_dim, THRU_DIM>(
          md, primary);

  {
    clog_tag_guard(coloring);
    clog_container_one(info, "closure", closure, clog::space);
  } // guard

  // Subtracting out the initial set leaves just the nearest
  // neighbors. This is similar to the image of the adjacency
  // graph of the initial indices.
  auto nearest_neighbors = flecsi::utils::set_difference(closure, primary);

  {
    clog_tag_guard(coloring);
    clog_container_one(
        info, "nearest neighbors", nearest_neighbors, clog::space);
  } // guard

  //----------------------------------------------------------------------------
  // Find one more layer of ghost cells now, these are needed to get some corner
  // cases right.
  //----------------------------------------------------------------------------

  // We can iteratively add halos of nearest neighbors, e.g.,
  // here we add the next nearest neighbors. For most mesh types
  // we actually need information about the ownership of these indices
  // so that we can deterministically assign rank ownership to vertices.
  auto nearest_neighbor_closure =
      flecsi::topology::entity_neighbors<cell_dim, cell_dim, THRU_DIM>(
          md, nearest_neighbors);

  {
    clog_tag_guard(coloring);
    clog_container_one(
        info, "nearest neighbor closure", nearest_neighbor_closure,
        clog::space);
  } // guard

  // Subtracting out the closure leaves just the
  // next nearest neighbors.
  auto next_nearest_neighbors =
      flecsi::utils::set_difference(nearest_neighbor_closure, closure);

  {
    clog_tag_guard(coloring);
    clog_container_one(
        info, "next nearest neighbor", next_nearest_neighbors, clog::space);
  } // guard

  // The union of the nearest and next-nearest neighbors gives us all
  // of the cells that might reference a vertex that we need.
  auto all_neighbors =
      flecsi::utils::set_union(nearest_neighbors, next_nearest_neighbors);

  {
    clog_tag_guard(coloring);
    clog_container_one(info, "all neighbors", all_neighbors, clog::space);
  } // guard

  //----------------------------------------------------------------------------
  // Find exclusive, shared, and ghost cells..
  //----------------------------------------------------------------------------

  // Get the rank and offset information for our nearest neighbor
  // dependencies. This also gives information about the ranks
  // that access our shared cells.
  auto cell_nn_info =
      communicator->get_primary_info(primary, nearest_neighbors);

  // Create a map version of the local info for lookups below.
  std::unordered_map<size_t, size_t> primary_indices_map;
  {
    size_t offset(0);
    for (auto i : primary) {
      primary_indices_map[offset++] = i;
    } // for
  }   // scope

  std::set<entry_info_t> exclusive_cells;
  std::set<entry_info_t> shared_cells;
  std::set<entry_info_t> ghost_cells;

  // Populate exclusive and shared cell information.
  {
    size_t offset(0);
    for (auto i : std::get<0>(cell_nn_info)) {
      if (i.size()) {
        shared_cells.insert(
            entry_info_t(primary_indices_map[offset], rank, offset, i));
      } else {
        exclusive_cells.insert(
            entry_info_t(primary_indices_map[offset], rank, offset, i));
      } // if
      ++offset;
    } // for
  }   // scope

  // Populate ghost cell information.
  {
    size_t offset(0);
    for (auto i : std::get<1>(cell_nn_info)) {
      ghost_cells.insert(i);
    } // for
  }   // scope

  {
    clog_tag_guard(coloring);
    clog_container_one(
        info, "exclusive cells ", exclusive_cells, clog::newline);
    clog_container_one(info, "shared cells ", shared_cells, clog::newline);
    clog_container_one(info, "ghost cells ", ghost_cells, clog::newline);
  } // guard

  //----------------------------------------------------------------------------
  // Create some maps for easy lookups when determining the other dependent
  // closures.
  //----------------------------------------------------------------------------

  // Create a map version for lookups below.
  std::unordered_map<size_t, entry_info_t> shared_cells_map;
  {
    for (auto i : shared_cells) {
      shared_cells_map[i.id] = i;
    } // for
  }   // scope

  // Get the rank and offset information for all relevant neighbor
  // dependencies. This information will be necessary for determining
  // shared vertices.
  auto cell_all_info = communicator->get_primary_info(primary, all_neighbors);

  // Create a map version of the remote info for lookups below.
  std::unordered_map<size_t, entry_info_t> remote_info_map;
  for (auto i : std::get<1>(cell_all_info)) {
    remote_info_map[i.id] = i;
  } // for

  // Get the intersection of our nearest neighbors with the nearest
  // neighbors of other ranks. This map of sets will only be populated
  // with intersections that are non-empty
  auto closure_intersection_map =
      communicator->get_intersection_info(nearest_neighbors);

  //----------------------------------------------------------------------------
  // Vertex Closure
  //----------------------------------------------------------------------------

  // Form the vertex closure
  auto vertex_closure =
      flecsi::topology::entity_closure<cell_dim, vertex_dim>(md, closure);

  // Assign vertex ownership
  vector<std::set<size_t>> vertex_requests(comm_size);
  std::set<entry_info_t> vertex_info;

  size_t offset(0);
  for (auto i : vertex_closure) {

    // Get the set of cells that reference this vertex.
    auto referencers =
        flecsi::topology::entity_referencers<cell_dim, vertex_dim>(md, i);

    {
      clog_tag_guard(coloring);
      clog_container_one(info, i << " referencers", referencers, clog::space);
    } // guard

    size_t min_rank(std::numeric_limits<size_t>::max());
    std::set<size_t> shared_vertices;

    // Iterate the direct referencers to assign vertex ownership.
    for (auto c : referencers) {

      // Check the remote info map to see if this cell is
      // off-color. If it is, compare it's rank for
      // the ownership logic below.
      if (remote_info_map.find(c) != remote_info_map.end()) {
        min_rank = std::min(min_rank, remote_info_map[c].rank);
        shared_vertices.insert(remote_info_map[c].rank);
      } else {
        // If the referencing cell isn't in the remote info map
        // it is a local cell.

        // Add our rank to compare for ownership.
        min_rank = std::min(min_rank, size_t(rank));

        // If the local cell is shared, we need to add all of
        // the ranks that reference it.
        if (shared_cells_map.find(c) != shared_cells_map.end()) {
          shared_vertices.insert(
              shared_cells_map[c].shared.begin(),
              shared_cells_map[c].shared.end());
        } // if
      }   // if

      // Iterate through the closure intersection map to see if the
      // indirect reference is part of another rank's closure, i.e.,
      // that it is an indirect dependency.
      for (auto ci : closure_intersection_map) {
        if (ci.second.find(c) != ci.second.end()) {
          shared_vertices.insert(ci.first);
        } // if
      }   // for
    }     // for

    if (min_rank == rank) {
      // This is a vertex that belongs to our rank.
      auto entry = entry_info_t(i, rank, offset, shared_vertices);
      vertex_info.insert(entry_info_t(i, rank, offset++, shared_vertices));
    } else {
      // Add remote vertex to the request for offset information.
      vertex_requests[min_rank].insert(i);
    } // if
  }   // for

  auto vertex_offset_info =
      communicator->get_entity_info(vertex_info, vertex_requests);

  std::set<entry_info_t> exclusive_vertices;
  std::set<entry_info_t> shared_vertices;
  std::set<entry_info_t> ghost_vertices;

  for (auto i : vertex_info) {
    if (i.shared.size()) {
      shared_vertices.insert(i);
    } else {
      exclusive_vertices.insert(i);
    } // if
  }   // for

  {
    size_t r(0);
    for (auto i : vertex_requests) {

      auto offset(vertex_offset_info[r].begin());
      for (auto s : i) {
        ghost_vertices.insert(entry_info_t(s, r, *offset));
        ++offset;
      } // for

      ++r;
    } // for
  }   // scope

  {
    clog_tag_guard(coloring);
    clog_container_one(
        info, "exclusive vertices ", exclusive_vertices, clog::newline);
    clog_container_one(
        info, "shared vertices ", shared_vertices, clog::newline);
    clog_container_one(info, "ghost vertices ", ghost_vertices, clog::newline);
  } // guard

  //----------------------------------------------------------------------------
  // output the result
  //----------------------------------------------------------------------------

  // figure out this ranks file name
  std::stringstream output_filename;
  output_filename << output_prefix;
  output_filename << "_rank";
  output_filename << std::setfill('0') << std::setw(6) << rank;
  output_filename << ".exo";

  // lambda function to convert to integer lists
  auto to_vec = [](const auto & list_in) {
    using value_type = typename std::decay_t<decltype(list_in.begin()->id)>;
    std::vector<value_type> list_out;
    list_out.reserve(list_in.size());
    for (auto & e : list_in)
      list_out.push_back(e.id);
    return list_out;
  };

  // write the exodus file
  using std::make_pair;
  md.write(
      output_filename.str(),
      {make_pair("exclusive cells", to_vec(exclusive_cells)),
       make_pair("shared cells", to_vec(shared_cells)),
       make_pair("ghost cells", to_vec(ghost_cells))},
      {make_pair("exclusive vertices", to_vec(exclusive_vertices)),
       make_pair("shared vertices", to_vec(shared_vertices)),
       make_pair("ghost vertices", to_vec(ghost_vertices))});

} // somerhing

////////////////////////////////////////////////////////////////////////////////
/// \brief the main cell coloring test
////////////////////////////////////////////////////////////////////////////////
TEST(color_exodus, 2d) {

  auto prefix_2d = std::string("mixed");
  exodus_definition_2d_t md2d(prefix_2d + ".exo");
  md2d.write("test2.exo");
  color_cells<thru_dim>(md2d, prefix_2d + "-colored");

}
  
TEST(color_exodus, 3d) {

  auto prefix_3d = std::string("box-hex");
  exodus_definition_3d_t md3d(prefix_3d + ".exo");
  md3d.write("test3.exo");
  color_cells<thru_dim>(md3d, prefix_3d + "-colored");

}
