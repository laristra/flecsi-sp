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
#include <flecsi-sp/burton/flecsi_mesh_wrapper.h>
#include <flecsi-sp/burton/flecsi_state_wrapper.h>
#include <flecsi-sp/utils/types.h>

#include <portage/driver/mmdriver.h>

// system includes
#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

// using statements
using std::cout;
using std::endl;

using mesh_t = flecsi_sp::burton::burton_mesh_t;
using index_spaces_t = mesh_t::index_spaces_t;

using real_t = mesh_t::real_t;
using vector_t = mesh_t::vector_t;

namespace flecsi_sp {
namespace burton {
namespace test {



////////////////////////////////////////////////////////////////////////////////
// Register the mesh state
////////////////////////////////////////////////////////////////////////////////
flecsi_register_field(
  mesh_t,      
  hydro,
  density,
  real_t,
  dense,
  1,
  index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,      
  hydro,
  velocity,
  vector_t,
  sparse,
  1,
  index_spaces_t::cells
);

flecsi_register_field(
  mesh_t,
  hydro,
  node_coordinates,
  mesh_t::vector_t,
  dense,
  1,
  index_spaces_t::vertices
);

///////////////////////////////////////////////////////////////////////////////
//! \brief Tack on an iteration number to a string
///////////////////////////////////////////////////////////////////////////////
static auto zero_padded( 
  std::size_t n, std::size_t padding = 6 
)
{
  std::stringstream ss;
  ss << std::setw( padding ) << std::setfill( '0' ) << n;
  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////
/// \brief output the solution
/// \param [in] mesh  the mesh object
/// \param [in] iteration  the iteration count
/// \param [in] d  the bulk density
////////////////////////////////////////////////////////////////////////////////
void output( 
  utils::client_handle_r__<mesh_t> mesh, 
  size_t iteration,
  utils::dense_handle_r__<real_t> d
) {
  clog(info) << "OUTPUT MESH TASK" << std::endl;

  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto rank = context.color();
  auto size = context.colors();
  auto time = 0.0;

  constexpr auto num_dims = mesh_t::num_dimensions;

  // figure out this ranks file name
  auto output_filename = "burton.remap_test"
    + zero_padded(iteration) + "_rank"
    + zero_padded(rank) + ".vtk";

  std::ofstream file(output_filename);

  file << "# vtk DataFile Version 3.0" << std::endl;
  file << "Remapping test" << std::endl;
  file << "ASCII" << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;
 
  file.precision(14);
  file.setf( std::ios::scientific );

  file << "POINTS " << mesh.num_vertices() << " double" << std::endl;
  for (auto v : mesh.vertices()) {
    for (int d=0; d<num_dims; d++ ) 
      file << v->coordinates()[d] << " ";
    for (int d=num_dims; d<3; d++ )
      file << 0 << " ";
    file << std::endl;
  }

  size_t vert_cnt = 0;
  const auto & cells = mesh.cells();
  for ( auto c : cells ) vert_cnt += mesh.vertices(c).size() + 1;

  file << "CELLS " << cells.size() << " " << vert_cnt << std::endl;    
  for ( auto c : cells ) {
    const auto & vs = mesh.vertices(c);
    file << vs.size() << " ";
    for ( auto v : vs )
      file << v.id() << " ";
    file << std::endl;
  }

  file << "CELL_TYPES " << cells.size() << std::endl;
  for ( auto c : cells ) file << "7" << std::endl;
    
  file << "CELL_DATA " << cells.size() << std::endl;    
  file << "SCALARS density double 1" << std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  
  for ( auto c : cells ) {
    file << d(c) << std::endl;
  }

  file.close();

}

////////////////////////////////////////////////////////////////////////////////
/// \brief make a remapper object
////////////////////////////////////////////////////////////////////////////////
template<
  typename mesh_wrapper_a_t,
  typename state_wrapper_a_t,
  typename mesh_wrapper_b_t,
  typename state_wrapper_b_t,
  typename var_name_t
  >
auto make_remapper(
         mesh_wrapper_a_t & mesh_wrapper_a,
         state_wrapper_a_t & state_wrapper_a,
         mesh_wrapper_b_t & mesh_wrapper_b,
         state_wrapper_b_t & state_wrapper_b,
         var_name_t & var_names
         ) {

  std::cout << " Inside make_remapper" << std::endl;
  for ( auto v: var_names) {
    state_wrapper_a.check_map(v);
    state_wrapper_b.check_map(v);
  }
  std::cout << " Inside make_remapper" << std::endl;

  constexpr int dim = mesh_wrapper_a_t::mesh_t::num_dimensions;

  if constexpr( dim == 2) {

    Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_2ndOrder,
      mesh_t::num_dimensions,
      flecsi_mesh_t<mesh_t>,
      flecsi_state_t<mesh_t>,
      flecsi_mesh_t<mesh_t>,
      flecsi_state_t<mesh_t> > remapper(
                mesh_wrapper_a,
                state_wrapper_a,
                mesh_wrapper_b,
                state_wrapper_b );
      return std::move(remapper);

  } else if (dim == 3) {

    Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR3D,
      Portage::Interpolate_2ndOrder,
      mesh_t::num_dimensions,
      flecsi_mesh_t<mesh_t>,
      flecsi_state_t<mesh_t>,
      flecsi_mesh_t<mesh_t>,
      flecsi_state_t<mesh_t> > remapper(
                mesh_wrapper_a,
                state_wrapper_a,
                mesh_wrapper_b,
                state_wrapper_b );
      return std::move(remapper);

  } else {
    static_assert(dim!=3 || dim!=2, "Make_remapper dimensions are out of range");
  }
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Test the remap capabilities of portage
/// \param [in] mesh       the mesh object
/// \param [in] mat_state  a densely populated set of data
/// \param [in] coord0     the set of coordinates to be applied
////////////////////////////////////////////////////////////////////////////////
void remap_test(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_rw__<real_t> density_handle,
  utils::sparse_handle_rw__<vector_t> velocity_handle,
  utils::dense_handle_r__<vector_t> new_vertex_coords
) {
  
  constexpr auto num_dims = mesh_t::num_dimensions;
  constexpr auto epsilon = config::test_tolerance;
  
  // get the context
  auto & context = flecsi::execution::context_t::instance();
  auto comm_rank = context.color();
  auto comm_size = context.colors();

  //---------------------------------------------------------------------------
  // Compute some source and target mesh quantities.

  // Apply the new coordinates to the mesh and update its geometry.
  // Here we save some info we need.
  std::vector< real_t > target_cell_volume(mesh.num_cells());
  std::vector< vector_t > target_cell_centroid(mesh.num_cells());
  // update coordinates
  for ( auto v : mesh.vertices() )
    std::swap(v->coordinates(), new_vertex_coords(v));
  // compute new target mesh quantities
  for ( auto c: mesh.cells() ) {
    c->update(&mesh);
    target_cell_volume[c] = c->volume();
    target_cell_centroid[c] = c->centroid();
  }

  // Store the cell volumes for the source mesh 
  std::vector< real_t > source_volume(mesh.num_cells());
  // update coordinates back to old mesh
  for ( auto v : mesh.vertices() )
    std::swap(v->coordinates(), new_vertex_coords(v));
  // re-compute source mesh quantities
  for ( auto c: mesh.cells() ) {
    c->update(&mesh);
    source_volume[c] = c->volume();
  }


  //---------------------------------------------------------------------------
  // Set up portage mesh/data wrappers

  // Create the mesh wrapper objects
  flecsi_mesh_t<mesh_t> mesh_wrapper_a(mesh);
  flecsi_mesh_t<mesh_t> mesh_wrapper_b(mesh);

  mesh_wrapper_b.set_new_coordinates(
    &new_vertex_coords(0),
    target_cell_volume.data(),
    target_cell_centroid.data() );

  // Create temporary vectors for the data to be remapped
  std::vector< double > remap_density(mesh.num_cells(), 0);

  std::vector< double > velocity(mesh.num_cells()*num_dims, 0);
  std::vector< double > remap_velocity(mesh.num_cells()*num_dims, 0);

  // Fill the temporary vectors with the data that needs to be remapped
  for (int dim=0; dim < num_dims; ++dim) {
    for (auto c: mesh.cells()) {
      for ( auto m : velocity_handle.entries(c) ) {
        velocity[c + dim*mesh.num_cells()] = velocity_handle(c,m)[dim];
      }
    }
  }

  // Create a vector of strings that correspond to the names of the variables
  //   that will be remapped
  std::vector<std::string> var_names;

  // Create the state wrapper objects
  flecsi_state_t<mesh_t> state_wrapper_a(mesh);
  flecsi_state_t<mesh_t> state_wrapper_b(mesh);

  // Add the fields that need to be remapped to the state wrappers
  // Special care is taken to handle each dimension of the velocity field
  var_names.push_back(std::string{"density"});
  state_wrapper_a.add_field( "density", &density_handle(0), "CELL" );
  state_wrapper_b.add_field( "density", remap_density.data(), "CELL" );

  static char coordinate[] = {'x', 'y', 'z'};    
  for ( int dim=0; dim<num_dims; ++dim ) {
    std::string name = "vel_";
    name += coordinate[dim];
    var_names.emplace_back(name);
    state_wrapper_a.add_field( name, velocity.data() + mesh.num_cells()*dim, "CELL" );
    state_wrapper_b.add_field( name, remap_velocity.data() + mesh.num_cells()*dim, "CELLL" );
  }

  
  std::cout << " Before make_remapper" << std::endl;
  for ( auto v: var_names) {
    state_wrapper_a.check_map(v);
    state_wrapper_b.check_map(v);
  }
  std::cout << " Before make_remapper" << std::endl;

  auto remapper = make_remapper(
             mesh_wrapper_a,
             state_wrapper_a,
             mesh_wrapper_b,
             state_wrapper_b,
             var_names);
  

  // Assign the remap varaible names for the portage driver
  remapper.set_remap_var_names(var_names);
  //remapper.set_limiter( Portage::Limiter_type::BARTH_JESPERSEN );
  remapper.set_limiter( Portage::Limiter_type::NOLIMITER );

  std::cout << " After make_remapper" << std::endl;
  for ( auto v: var_names) {
    state_wrapper_a.check_map(v);
    state_wrapper_b.check_map(v);
  }
  std::cout << " After make_remapper" << std::endl;

  // Do the remap 
  // Argument is boolean: distributed (true or false)
  auto mpi_comm = MPI_COMM_WORLD;
  Wonton::MPIExecutor_type mpiexecutor(mpi_comm);
  remapper.run(&mpiexecutor);

  //---------------------------------------------------------------------------
  // Swap old for new data

  // Checking conservation for remap test
  real_t total_density{0};
  vector_t total_velocity{0};

  for (auto c: mesh.cells(flecsi::owned)) {
    total_density +=  c->volume() *  density_handle(c);
    density_handle(c) = remap_density[c];
    for ( auto m : velocity_handle.entries(c) ) {
      auto & vel = velocity_handle(c,m);
      for (int dim=0; dim < num_dims; ++dim) {
        total_velocity[dim] += c->volume() * vel[dim];
        vel[dim] = remap_velocity[c + dim * mesh.num_cells()];
      }
    }
  }

  // Apply the new coordinates to the mesh and update its geometry
  for ( auto v : mesh.vertices() ) v->coordinates() = new_vertex_coords(v);
  mesh.update_geometry();
  
  //---------------------------------------------------------------------------
  // Post process

  // Check conservation for remap test
  real_t total_remap_density{0};
  vector_t total_remap_velocity{0};
  // Calculate the L1 and L2 norms
  real_t total_vol = 0.0;
  real_t L1 = 0.0;
  real_t L2 = 0.0;

  for (auto c : mesh.cells(flecsi::owned) ) {

    // remapped sums
    auto actual = density_handle(c);
    total_remap_density += c->volume() * density_handle(c);
    for ( auto m : velocity_handle.entries(c) ) {
      const auto & vel = velocity_handle(c,m);
      for (int dim=0; dim<num_dims; ++dim) {
        total_remap_velocity[dim] += c->volume() * vel[dim];
      }

      // expected answer
      real_t sum{0}, expected;
#ifdef CUBIC
      for ( int i=0; i < num_dims; ++i)
        sum += ((c->centroid())[i]+0.5);
      expected = pow(sum,3.0);
#elif LINEAR
      for ( int i=0; i < num_dims; ++i)
        sum += ((c->centroid())[i]+0.5);
      expected = sum * 2.0;
#elif COSINE
      real_t arg = 10.0 * std::sqrt( pow((c->centroid())[0], 2) + pow((c->centroid())[1], 2) );
      expected = 1.0 + (std::cos(arg) / 3.0);
#else
      expected = 1.0;
      EXPECT_NEAR( expected, actual, epsilon );
#endif

      // L1,L2 errors
      auto each_error = std::abs(actual - expected);
      L1 +=  c->volume() * each_error; // L1
      L2 +=  c->volume() * std::pow(each_error, 2); // L2
      total_vol += c->volume();
    }
  }

  // output L1, L2 errors for debug
  real_t total_vol_sum{0}, L1_sum{0}, L2_sum{0};
  int ret;
  ret = MPI_Allreduce( &total_vol, &total_vol_sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  ret = MPI_Allreduce( &L1, &L1_sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  ret = MPI_Allreduce( &L2, &L2_sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  L1_sum = L1_sum/total_vol_sum; // L1
  L2_sum = pow(L2_sum/total_vol_sum, 0.5); // L2

  if ( comm_rank == 0 ) {
    printf("L1 Norm: %.20f, Total Volume: %f \n", L1_sum, total_vol_sum);
    printf("L2 Norm: %.20f, Total Volume: %f \n", L2_sum, total_vol_sum);
  }

  // Verify conservation
  real_t total_density_sum{0}, total_remap_density_sum{0};
  vector_t total_velocity_sum{0}, total_remap_velocity_sum{0};
  ret = MPI_Allreduce( &total_density, &total_density_sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  ret = MPI_Allreduce( &total_remap_density, &total_remap_density_sum, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
  ret = MPI_Allreduce( total_velocity.data(), total_velocity_sum.data(), total_velocity.size(), MPI_DOUBLE, MPI_SUM, mpi_comm);
  ret = MPI_Allreduce( total_remap_velocity.data(), total_remap_velocity_sum.data(), total_remap_velocity.size(), MPI_DOUBLE, MPI_SUM, mpi_comm);

  EXPECT_NEAR( total_density_sum, total_remap_density_sum, epsilon);
  for ( int i=0; i<num_dims; ++i)
    EXPECT_NEAR( total_velocity_sum[i], total_remap_velocity_sum[i], epsilon);

}


////////////////////////////////////////////////////////////////////////////////
/// \brief Test the remap capabilities of portage
/// \param [in] mesh       the mesh object
/// \param [in] mat_state  a densely populated set of data
////////////////////////////////////////////////////////////////////////////////
void initialize_flat(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_w__<real_t> density,
  utils::sparse_mutator__<vector_t> velocity
) {

  int m = 0;
  for (auto c: mesh.cells()){
    density(c) = 1.0;
    velocity(c,m) = 1.0;
  }
}

void initialize_linear(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_w__<real_t> density,
  utils::sparse_mutator__<vector_t> velocity
) {

  real_t val;
  int m = 0;
  for (auto c: mesh.cells()){
    real_t sum = 0.0;
    for ( int i=0; i < mesh_t::num_dimensions; i++)
      sum += ((c->centroid())[i]+0.5);
    val = sum * 2.0;
    density(c) = val;
    velocity(c,m) = val;
  }
}

void initialize_cubic(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_w__<real_t> density,
  utils::sparse_mutator__<vector_t> velocity
) {

  int m = 0;
  for (auto c: mesh.cells()){
    real_t sum = 0.0;
    for ( int i=0; i < mesh_t::num_dimensions; i++)
      sum += pow((c->centroid())[i]+0.5, 3.0);
    auto val = sum;
    density(c) = val;
    velocity(c,m) = val;
  }
}

void initialize_cosine(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_w__<real_t> density,
  utils::sparse_mutator__<vector_t> velocity
) {

  int m = 0;
  for (auto c: mesh.cells()){
    real_t sum = 0.0;
    for ( int i=0; i < mesh_t::num_dimensions; i++)
      sum += pow((c->centroid())[i], 2);
    auto arg = 10.0 * std::sqrt( sum );
    auto val = 1.0 + (std::cos(arg) / 3.0);
    density(c) = val;
    velocity(c,m) = val;
  }
}

void initialize_step(
  utils::client_handle_r__<mesh_t> mesh,
  utils::dense_handle_w__<real_t> density,
  utils::sparse_mutator__<vector_t> velocity
) {

  int m = 0;
  for (auto c: mesh.cells()){
    real_t val;
    if ( c < mesh.num_cells()/2 ){
      val = 0;
    } else {
      val = 1;
    }
    density(c) = val;
    velocity(c,m) = val;
  }
}

////////////////////////////////////////////////////////////////////////////////
//! \brief modify the coordinates & solution
//!
//! \param [in] mesh the mesh object
//! \param [out] coord0  storage for the mesh coordinates
////////////////////////////////////////////////////////////////////////////////
void modify( 
  utils::client_handle_r__<mesh_t>  mesh,
  utils::dense_handle_w__<vector_t> coord
)
{
  constexpr auto num_dims = mesh_t::num_dimensions;

  // Loop over vertices
  const auto & vs = mesh.vertices();
  auto num_verts = vs.size();

  std::mt19937 mt_rand(0);
  auto real_rand = std::bind(
    std::uniform_real_distribution<double>(0, 0.2),
    mt_rand);

  //Randomly modifies each of the vertices such that 
  // no vertex moves by more than 40% of the spacing
  // in any direction. These vertices are then used
  // for the new mesh.

  double spacing = 1/std::pow(num_verts, 0.5);
  for (auto v : vs ) {
    auto vertex = v->coordinates();
    
    for ( int j=0; j < num_dims; j++) {
      if ( vertex[j] < 0.5 && vertex[j] > 0) {
        auto perturbation = spacing * real_rand();
        vertex[j] -= perturbation;
      } else if (vertex[j] <= 0 && vertex[j]> -0.5) {
        auto perturbation = spacing * real_rand();
        vertex[j] += perturbation;
      }
    }
    
    coord(v) = vertex;
  }
} 

////////////////////////////////////////////////////////////////////////////////
//! \brief restore the coordinates
//!
//! \param [in,out] mesh  the mesh object
//! \param [in] coord0  the mesh coordinates to restore
////////////////////////////////////////////////////////////////////////////////
void restore( 
  utils::client_handle_r__<mesh_t>  mesh,
  utils::dense_handle_r__<vector_t> coord
)
{

  // Loop over vertices
  auto vs = mesh.vertices();
  auto num_verts = vs.size();

  for ( mesh_t::counter_t i=0; i<num_verts; i++ ) {
    auto vt = vs[i];
    vt->coordinates() = coord(vt);
  }

}// TEST_F

flecsi_register_mpi_task(remap_test, flecsi_sp::burton::test);
flecsi_register_task(output, flecsi_sp::burton::test, loc,
  single|flecsi::leaf);

// Different Initialization Tasks
flecsi_register_task(initialize_flat, flecsi_sp::burton::test, loc,
  single|flecsi::leaf);
flecsi_register_task(initialize_linear, flecsi_sp::burton::test, loc,
  single|flecsi::leaf);
flecsi_register_task(initialize_cubic, flecsi_sp::burton::test, loc,
  single|flecsi::leaf);
flecsi_register_task(initialize_cosine, flecsi_sp::burton::test, loc,
  single|flecsi::leaf);
flecsi_register_task(initialize_step, flecsi_sp::burton::test, loc,
  single|flecsi::leaf);

flecsi_register_task(restore, flecsi_sp::burton::test, loc,
         single|flecsi::leaf);
flecsi_register_task(modify, flecsi_sp::burton::test, loc,
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



  size_t time_cnt{0};

  auto xn = flecsi_get_handle(mesh_handle, hydro, node_coordinates, vector_t, dense, 0);
  auto density_handle = flecsi_get_handle(mesh_handle, hydro, density, real_t, dense, 0);
  auto velocity_handle = flecsi_get_handle(mesh_handle, hydro, velocity, vector_t, sparse, 0);
  auto velocity_mutator = flecsi_get_mutator(mesh_handle, hydro, velocity, vector_t, sparse, 0, 1);

#ifdef CUBIC
  flecsi_execute_task(
          initialize_cubic,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          density_handle,
          velocity_mutator);
#elif LINEAR
  flecsi_execute_task(
          initialize_linear,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          density_handle,
          velocity_mutator);
#elif STEP
  flecsi_execute_task(
          initialize_step,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          density_handle,
          velocity_mutator);
#elif COSINE
  flecsi_execute_task(
          initialize_cosine,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          density_handle,
          velocity_mutator);
#else
  flecsi_execute_task(
          initialize_flat,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          density_handle,
          velocity_mutator);
#endif
  
  time_cnt++;
  flecsi_execute_task(
            output,
            flecsi_sp::burton::test,
            single,
            mesh_handle, time_cnt, density_handle);

  flecsi_execute_task(
          modify,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          xn);

#if 1
  flecsi_execute_mpi_task(
          remap_test, 
          flecsi_sp::burton::test, 
          mesh_handle,
          density_handle,
          velocity_handle,
          xn);
#else
  flecsi_execute_task(
          restore,
          flecsi_sp::burton::test,
          single,
          mesh_handle,
          xn);
#endif

  time_cnt++;
  flecsi_execute_task(
            output,
            flecsi_sp::burton::test,
            single,
            mesh_handle, time_cnt, density_handle);

} // driver

} // namespace execution
} // namespace flecsi

////////////////////////////////////////////////////////////////////////////////
//! \brief Only here so test runs?
////////////////////////////////////////////////////////////////////////////////
TEST(burton, portage) {}

