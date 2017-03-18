/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef mmstate_h
#define mmstate_h

#include <set>

#include "config.h"

///
// In this example, num_mats is the logical number of materials one
// would use for a dense storage type.
///
static constexpr size_t num_mats = 5;

///
// The driver task is the top-level user task for a FleCSI application. The
// driver task is executed after runtime and specailization initialization
// have been executed. Semantically, the driver task is sequential. However,
// when built against a parallel backend to FleCSI, it will execute in a
// distributed-memory mode. The fact of parallel execution is not exposed to
// the user.
///
void driver(int argc, char ** argv) {

  // Mesh object
	minimal_mesh_t mesh;

  // Initialize the mesh
	init_mesh(mesh);

  // Convenience struct to hold cell-wise data
  struct cell_data_t {
    double density;
    double volume;
  }; // struct cell_data_t

  ///
  // This example uses a mixture of sparse and dense storage types. The
  // material volume fractions are stored using the sparse storage type,
  // while the macroscopic quantities use a dense storage type.
  ///
	flecsi_register_data(mesh, solver, density, double, sparse,
    1, cells, num_mats);
	flecsi_register_data(mesh, solver, macro, cell_data_t, dense,
    1, cells);

  // Initialize material data and cell volumes
  {
  ///
  // A mutator is an accessor for the sparse storage type that allows
  // the user to mutate the sparsity pattern of the data. Users can add
  // non-zero entries via the () operator described below.
  ///
  auto mats = flecsi_get_mutator(mesh, solver, density, double, sparse, 0, 3);

  ///
  // An accessor allows read/write access to registered data.
  ///
  auto cell = flecsi_get_accessor(mesh, solver, macro, cell_data_t, dense, 0);

  ///
  // Iterate over the cells and randomly initialize non-zero
  // material quantities.
  ///
  std::set<size_t> nzs;
  for(auto c: mesh.cells()) {
    cell(c).volume = c->volume();

    ///
    // Create a random set of indices that will be non-zero.
    ///
    int nz = num_mats*(rand()+1)/RAND_MAX;
    nz = nz > 0 ? nz : 1;

    int indices = 0;
    nzs.clear();
    while(nzs.size()<nz) {
      nzs.insert(num_mats*(rand()+1)/RAND_MAX);
    } // while

    ///
    // Iterate over the indices in the random set and assign a value.
    ///
    for(auto e: nzs) {
      std::cout << e << " ";

      ///
      // Set non-zero material values.
      //
      // The parameters are:
      // c The cell index (from the outer 'for' loop iterator variable).
      // e The material index (from the inner 'for' loop iterator variable).
      ///
      mats(c, e) = 1.0;
    } // for

    std::cout << std::endl;
  } // for

  } // scope

  ///-------------------------------------------------------------------------//
  // Algorithm 1: Compute average density using a cell-centric approach
  // on a per-cell.
  ///-------------------------------------------------------------------------//
  {
  ///
  // Here we use an accessor for the sparse storage type because we only
  // need read/write access (as opposed to mutate access).
  ///
  auto mats = flecsi_get_accessor(mesh, solver, density, double, sparse, 0);
  auto cell = flecsi_get_accessor(mesh, solver, macro, cell_data_t, dense, 0);

  ///
  // Iterate over the cells.
  ///
  for(auto c: mesh.cells()) {
    cell(c).density = 0.0;

    ///
    // Iterate over the materials within this cell.
    ///
    for(auto m: mats.entries(c)) {
      cell(c).density += mats(c, m);
    } // for

    // Normalize.
    cell(c).density /= cell(c).volume;

    // Print the result.
    std::cout << "average " << c->id<0>() << ": " <<
      cell(c).density << std::endl;
  } // for

  } // scope

  ///-------------------------------------------------------------------------//
  // Algorithm 2: Compute average density using a material-centrix approach.
  ///-------------------------------------------------------------------------//
  {
  // See explanation for Algorithm 1.
  auto mats = flecsi_get_accessor(mesh, solver, density, double, sparse, 0);
  auto cell = flecsi_get_accessor(mesh, solver, macro, cell_data_t, dense, 0);

  ///
  // Iterate over the cells to re-zero the macroscopic density.
  ///
  for(auto c: mesh.cells()) {
    cell(c).density = 0.0;
  } // for

  ///
  // Iterate over the materials.
  ///
  for(auto m: mats.entries()) {

    ///
    // Iterate over the cells that contain this material.
    ///
    for(auto c: mats.indices(m)) {
      cell(c).density += mats(c, m);
    } // for
  } // for

  ///
  // Iterate over the cells again to normalize and print the result.
  ///
  for(auto c: mesh.cells()) {
    // Normalize.
    cell(c).density /= cell(c).volume;

    // Print the result.
    std::cout << "average " << c->id<0>() << ": " <<
      cell(c).density << std::endl;
  } // for

  } // scope

} // driver

#endif // mmstate_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
