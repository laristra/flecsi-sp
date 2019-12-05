#include <boost/program_options.hpp>
#include <exodusII.h>
#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace po = boost::program_options;

using integer_t = int;
using unsigned_integer_t = std::size_t;
using real_t = double;
    
template <typename T, typename Container>
auto write(
  int exo_id,
  ex_entity_id blk_id,
  const Container & cell_vertices,
  unsigned_integer_t num_verts_per_cell
) {
  
  std::string desc;
  
  switch(num_verts_per_cell) {
    case 4:
      desc = "quad4";
      break;
    case 8:
      desc = "hex8";
      break;
    default:
      std::cerr
        << "Unknown element type with " << num_verts_per_cell << " verts"
        << std::endl;
      return -1;
  }
  
  auto num_elem = cell_vertices.size() / num_verts_per_cell;
  auto status = ex_put_block( exo_id, EX_ELEM_BLOCK, blk_id, desc.c_str(),
    num_elem, num_verts_per_cell, 0, 0, 0 );
  if (status) return status;

  status = ex_put_name(exo_id, EX_ELEM_BLOCK, blk_id, "element block");
  if (status) return status;

  if ( std::is_same<typename Container::value_type, T>::value )
    status = ex_put_elem_conn( exo_id, blk_id, cell_vertices.data() );
  else {
    std::vector<T> vec(cell_vertices.begin(), cell_vertices.end());
    status = ex_put_elem_conn( exo_id, blk_id, vec.data() );
  }
  
  return status;
}

template <typename T, typename Container>
auto write_mapping(
  int exo_id,
  const Container & cell_mapping,
  const Container & vert_mapping
) {
  int status;
  if ( std::is_same<typename Container::value_type, T>::value ) {
    status = ex_put_elem_num_map( exo_id, cell_mapping.data() );
    if (status) return status;
    status = ex_put_node_num_map( exo_id, vert_mapping.data() );
    if (status) return status;
  }
  else {
    std::vector<T> vec(cell_mapping.begin(), cell_mapping.end());
    status = ex_put_elem_num_map( exo_id, vec.data() );
    if (status) return status;
    vec.assign(vert_mapping.begin(), vert_mapping.end());
    status = ex_put_node_num_map( exo_id, vec.data() );
    if (status) return status;
  }
  return status;
}

void subdivide(
    unsigned_integer_t nelem,
    unsigned_integer_t npart,
    std::vector<unsigned_integer_t> & dist )
{

  size_t quot = nelem / npart;
  size_t rem = nelem % npart;

  dist.clear();
  dist.reserve(npart + 1);
  dist.push_back(0);

  // Set the distributions for each rank. This happens on all ranks.
  // Each rank gets the average number of indices, with higher ranks
  // getting an additional index for non-zero remainders.
  for(size_t r(0); r < npart; ++r) {
    const size_t indices = quot + ((r >= (npart - rem)) ? 1 : 0);
    dist.push_back(dist[r] + indices);
  } // for
};

int main( int argc, char* argv[] )
{
  
  bool is_verbose = false;
  bool force_int64 = false;

  MPI_Init(&argc, &argv);

  int comm_rank, comm_size;
  MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("dimensions,d", po::value<std::vector<integer_t>>()->multitoken(), "Box dimensions.")
      ("partitions,p", po::value<std::vector<integer_t>>()->multitoken(), "Number partitions in each direction.")
      ("lower,l", po::value<std::vector<real_t>>()->multitoken(), "Box lower coordinates.")
      ("upper,u", po::value<std::vector<real_t>>()->multitoken(), "Box upper coordinates.")
      ("output-file,o", po::value<std::string>(), "output file")
      ("verbose,v", "Print extra debug info.")
      ("large-integer,i", "Use large integer support.")
  ;
  
  po::positional_options_description p;
  p.add("output-file,o", -1);

  po::variables_map vm;
  auto extras = po::command_line_style::unix_style ^ po::command_line_style::allow_short;
  po::store( po::command_line_parser(argc, argv).style(extras).options(desc).positional(p).run(), vm);
  
  if (!vm.count("dimensions")) {
    if (comm_rank == 0) {
      std::cout << "No dimensions provided!" << std::endl;
      std::cout << desc << "\n";
    }
    MPI_Finalize();
    return 1;
  }
  
  if (!vm.count("output-file")) {
    if (comm_rank == 0) {
      std::cout << "No output file provided!" << std::endl;
      std::cout << desc << "\n";
    }
    MPI_Finalize();
    return 1;
  }

  if (vm.count("help")) {
    if (comm_rank == 0) {
      std::cout << desc << "\n";
    }
    MPI_Finalize();
    return 1;
  }

  po::notify(vm);
 
  bool print_help = false;

  auto name = vm["output-file"].as<std::string>();

  const auto & args_dims = vm["dimensions"].as<std::vector<integer_t>>();
  std::vector<unsigned_integer_t> global_dims( args_dims.begin(), args_dims.end() );
  auto num_dims = global_dims.size();

  std::vector<real_t> global_lower(num_dims, -0.5);
  std::vector<real_t> global_upper(num_dims,  0.5);
  
  if (num_dims > 3 || num_dims <= 1) {
    if (comm_rank == 0) {
      std::cout << "Only two to three dimensions supported!" << std::endl;
    }
    print_help = true;
  }

  if (vm.count("lower")) {
    const auto & args = vm["lower"].as<std::vector<real_t>>();
    auto n = args.size();
    if ( n != num_dims ) {
      if (comm_rank == 0) {
        std::cout
          << "Number of lower box coordinates must match the number of "
          << "dimensions, " << n << " provided, " << num_dims << " expected."
          << std::endl;
      }
      print_help = true;
    }
    global_lower = args;
  }

  if (vm.count("upper")) {
    const auto & args = vm["upper"].as<std::vector<real_t>>();
    auto n = args.size();
    if ( n != num_dims ) {
      if (comm_rank == 0) {
        std::cout
          << "Number of upper box coordinates must match the number of "
          << "dimensions, " << n << " provided, " << num_dims << " expected."
          << std::endl;
      }
      print_help = true;
    }
    global_upper = args;
  }

  if (print_help) {
    if (comm_rank == 0) {
      std::cout << desc << "\n";
    }
    MPI_Finalize();
    return -1;
  }

  std::vector<unsigned_integer_t> block_sizes(num_dims, 1);

  if (vm.count("partitions")) {
    if (comm_rank==0) std::cout << "Using provided partitioning" << std::endl;
    const auto & args = vm["partitions"].as<std::vector<int>>();
    auto n = args.size();
    if ( n != num_dims ) {
      if ( comm_rank == 0 ) {
        std::cout
          << "Number of partition must match the number of "
          << "dimensions, " << n << " privided, " << num_dims << " expected."
          << std::endl;
      }
      MPI_Finalize();
      return -1;
    }
    unsigned_integer_t total_partitions = 1;
    for ( unsigned_integer_t i=0; i<num_dims; ++i ) {
      block_sizes[i] = args[i];
      total_partitions *= args[i];
    }
    if ( total_partitions != comm_size ) {
      if ( comm_rank == 0 ) {
        std::cout
          << "Total number of partitions does not match the total number of "
          << " mpi ranks" << std::endl;
      }
      MPI_Finalize();
      return -1;
    }
  }
  else {

    auto is_even = (comm_size % 2 == 0);

    if (is_even) {
      unsigned_integer_t avg_block_size = std::pow<real_t>( comm_size, 1./num_dims );
      unsigned_integer_t tot_parts = 1;
      for ( int i=0; i<num_dims-1; ++i ) {
        block_sizes[i] = avg_block_size;
        tot_parts *= avg_block_size;
      }
      block_sizes[num_dims-1] = comm_size / tot_parts;
    }

    else {
      block_sizes[0] = comm_size;
    }
      
    if (comm_rank==0) {
      std::cout << "Automatically partitioned block sizes: ";
      for ( auto i : block_sizes ) std::cout << i << " ";
      std::cout << std::endl;
    }

  } // partitioning
  
  if (vm.count("large-integer")) {
    force_int64 = true;
    if (comm_rank == 0) {
      std::cout << "Forcing 64 bit integers" << std::endl;
    }
  }

  if (vm.count("verbose")) {
    is_verbose = true;
    if (comm_rank == 0) {
      std::cout << "Adding extra exodus debug." << std::endl;
    }
    ex_opts(EX_ABORT | EX_VERBOSE);
  }

  if (comm_rank == 0) {
    std::cout << num_dims << "-dimensional mesh selected." << std::endl;
  }

  // build a mesh
  std::vector<real_t> delta(num_dims);
  for ( unsigned_integer_t i=0; i<num_dims; ++i ) {
    if ( global_upper[i] < global_lower[i] ) {
      if (comm_rank == 0) {
        std::cout
          << "Bounding box is invalid: for i=" << i << ", " << global_lower[i]
          << " should be less than " << global_upper[i]
          << std::endl;
      }
      MPI_Finalize();
      return -1;
    }
    delta[i] =  (global_upper[i] - global_lower[i]) / global_dims[i];
  }

  std::vector<unsigned_integer_t> block_ijk(num_dims);
  auto temp = comm_rank;
  for ( int i=0; i<num_dims-1; ++i ) {
    block_ijk[i] = temp % block_sizes[i];
    block_ijk[i+1] = (temp - block_ijk[i]) / block_sizes[i];
    temp -= block_ijk[i];
    temp /= block_sizes[i];
  }

  if (is_verbose) {
    std::cout << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    for ( unsigned_integer_t i=0; i<comm_size; ++i )
    {
      if ( comm_rank == i ) {
        std::cout << "Rank " << comm_rank << " has block dimensions: ";
        for ( unsigned_integer_t j=0; j<num_dims; ++j ) std::cout << block_ijk[j] << " ";
        std::cout << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  
  std::vector< std::vector<unsigned_integer_t> > block_displ(num_dims);
  for ( int i=0; i<num_dims; ++i )
    subdivide( global_dims[i], block_sizes[i], block_displ[i] );

  unsigned_integer_t num_verts = 1;
  unsigned_integer_t num_cells = 1;
  unsigned_integer_t global_verts = 1;
  unsigned_integer_t global_cells = 1;
  unsigned_integer_t num_verts_per_cell = 1;

  std::vector<unsigned_integer_t> local_dims(num_dims);
  std::vector<real_t> local_lower(num_dims);
  std::vector<real_t> local_upper(num_dims);

  for ( unsigned_integer_t i=0; i<num_dims; ++i ) {
    global_verts *= global_dims[i]+1;
    global_cells *= global_dims[i];
    num_verts_per_cell *= 2;
    auto jblock = block_ijk[i];
    const auto & displ = block_displ[i];
    local_dims[i] = displ[jblock+1] - displ[jblock];
    num_cells *= local_dims[i];
    num_verts *= local_dims[i]+1;
    local_lower[i] = global_lower[i] + (displ[jblock])*delta[i] ;
    local_upper[i] = global_lower[i] + (displ[jblock+1])*delta[i];
  }


  MPI_Datatype mpi_type_t;
  if ( sizeof(unsigned_integer_t) == sizeof(unsigned int) )
    mpi_type_t = MPI_UNSIGNED;
  else if ( sizeof(unsigned_integer_t) == sizeof(unsigned long long) )
    mpi_type_t = MPI_UNSIGNED_LONG_LONG;
  else {
    if (comm_rank == 0) std::cout << "UNKNOWN MPI DATA TYPE." << std::endl;
    MPI_Finalize();
    return -1;
  }

  unsigned_integer_t temp_cells;
  MPI_Allreduce( &num_cells, &temp_cells, 1, mpi_type_t, MPI_SUM, MPI_COMM_WORLD );
  if ( temp_cells != global_cells ) {
    if ( comm_rank == 0) {
      std::cout << "Summed cell count should be " << global_cells
        << " got " << temp_cells << std::endl;
    }
    MPI_Finalize();
    return -1;
  }

  std::vector<real_t> coordinates(num_dims * num_verts);
  std::vector<unsigned_integer_t> cell_vertices(num_cells*num_verts_per_cell);
  std::vector<unsigned_integer_t> global_cell_id(num_cells);
  std::vector<unsigned_integer_t> global_vertex_id(num_verts);

  if ( num_dims == 2 ) {

    auto vertex_local_id = [&](auto i, auto j) {
      return i + (local_dims[0]+1)*j;
    };
    auto vertex_global_id = [&](auto i, auto j) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      return iglobal + (global_dims[0]+1)*jglobal;
    };

    for ( unsigned_integer_t j=0; j<local_dims[1]+1; ++j ) {
      for ( unsigned_integer_t i=0; i<local_dims[0]+1; ++i ) {
        auto id = vertex_local_id(i, j);
        coordinates[id              ] = local_lower[0] + i*delta[0];
        coordinates[id +   num_verts] = local_lower[1] + j*delta[1];
        global_vertex_id[id] = vertex_global_id(i,j) + 1;
      }
    }
    
    auto cell_local_id = [&](auto i, auto j) {
      return i + local_dims[0]*j;
    };
    auto cell_global_id = [&](auto i, auto j) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      return iglobal + global_dims[0]*jglobal;
    };

    for ( unsigned_integer_t j=0, id=0; j<local_dims[1]; ++j ) {
      for ( unsigned_integer_t i=0; i<local_dims[0]; ++i ) {
        auto local_id = cell_local_id(i, j);
        global_cell_id[local_id] = cell_global_id(i, j) + 1;
        cell_vertices[id++] = vertex_local_id(i  ,   j) + 1;
        cell_vertices[id++] = vertex_local_id(i+1,   j) + 1;
        cell_vertices[id++] = vertex_local_id(i+1, j+1) + 1;
        cell_vertices[id++] = vertex_local_id(i  , j+1) + 1;
      }
    }

  } // two-dimensional

  else if ( num_dims == 3 ) {

    auto vertex_local_id = [&](auto i, auto j, auto k) {
      return i + (local_dims[0]+1)*(j + (local_dims[1]+1)*k);
    };
    auto vertex_global_id = [&](auto i, auto j, auto k) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      
      auto kblk = block_ijk[2];
      auto kstart = block_displ[2][ kblk ];
      auto kglobal = kstart + k;
      return iglobal + (global_dims[0]+1)*(jglobal + (global_dims[1]+1)*kglobal);
    };
    
    auto cell_local_id = [&](auto i, auto j, auto k) {
      return i + local_dims[0]*(j + local_dims[1]*k);
    };
    auto cell_global_id = [&](auto i, auto j, auto k) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      
      auto kblk = block_ijk[2];
      auto kstart = block_displ[2][ kblk ];
      auto kglobal = kstart + k;
      return iglobal + global_dims[0]*(jglobal + global_dims[1]*kglobal);
    };

    for ( unsigned_integer_t k=0; k<local_dims[2]+1; ++k ) {
      for ( unsigned_integer_t j=0; j<local_dims[1]+1; ++j ) {
        for ( unsigned_integer_t i=0; i<local_dims[0]+1; ++i ) {
          auto id = vertex_local_id(i, j, k);
          coordinates[id              ] = local_lower[0] + i*delta[0];
          coordinates[id +   num_verts] = local_lower[1] + j*delta[1];
          coordinates[id + 2*num_verts] = local_lower[2] + k*delta[2];
          global_vertex_id[id] = vertex_global_id(i,j,k) + 1;
        }
      }
    }

    for ( unsigned_integer_t k=0, id=0; k<local_dims[2]; ++k ) {
      for ( unsigned_integer_t j=0; j<local_dims[1]; ++j ) {
        for ( unsigned_integer_t i=0; i<local_dims[0]; ++i ) {
          auto local_id = cell_local_id(i, j, k);
          global_cell_id[local_id] = cell_global_id(i, j, k) + 1;
          cell_vertices[id++] = vertex_local_id(i  ,   j, k) + 1;
          cell_vertices[id++] = vertex_local_id(i+1,   j, k) + 1;
          cell_vertices[id++] = vertex_local_id(i+1, j+1, k) + 1;
          cell_vertices[id++] = vertex_local_id(i  , j+1, k) + 1;

          cell_vertices[id++] = vertex_local_id(i  ,   j, k+1) + 1;
          cell_vertices[id++] = vertex_local_id(i+1,   j, k+1) + 1;
          cell_vertices[id++] = vertex_local_id(i+1, j+1, k+1) + 1;
          cell_vertices[id++] = vertex_local_id(i  , j+1, k+1) + 1;
        }
      }
    }
  
  } // three-dimensional

  // check if we need 64 bit integers
  constexpr auto max_int = std::numeric_limits<int>::max();
  auto is_int64 = 
    force_int64 ||
    (cell_vertices.size() >= max_int/2) ||
    global_verts >= max_int ||
    global_cells >= max_int;

  // size of floating point variables used in app.
  int app_word_size = sizeof(real_t);

  // size of floating point to be stored in file.
  // change to float to save space
  int exo_word_size = sizeof(real_t);

  // determine the file creation mode
  auto cmode = EX_CLOBBER;
  if (is_int64) {
    cmode |= EX_ALL_INT64_DB | EX_ALL_INT64_API;
    if (comm_rank==0) std::cout << "Using 64-bit integers." << std::endl;
  }
  else
    if (comm_rank==0) std::cout << "Using 32-bit integers." << std::endl;

  // figure out mesh file name
  std::string filename;
  if (comm_size == 1 )
    filename = name;
  else {
    auto number = comm_size;
    unsigned int num_digits = 1;
    while (number) {
      number /= 10;
      num_digits++;
    }
    std::stringstream ss;
    ss
      << name
      << "."
      << std::setfill('0') << std::setw(num_digits) << comm_size
      << "."
      << std::setfill('0') << std::setw(num_digits) << comm_rank;
    filename = ss.str();
  }

  // create file
  auto exo_id =
      ex_create(filename.c_str(), cmode, &app_word_size, &exo_word_size);
  if (exo_id < 0) {
    if ( comm_rank == 0 ) {
      std::cerr << "Problem writing exodus file, ex_create() returned " << exo_id << std::endl;
    }
    MPI_Finalize();
    return exo_id;
  }
  else {
    if (comm_rank==0) std::cout << "Opened file for writing: " << name << std::endl;
  }

  
  ex_init_params exopar;
  std::strcpy(exopar.title, "Exodus II output from make_mesh.");
  exopar.num_dim = num_dims;
  exopar.num_nodes = num_verts;
  exopar.num_edge = 0;
  exopar.num_edge_blk = 0;
  exopar.num_face = 0;
  exopar.num_face_blk = 0;
  exopar.num_elem = num_cells;
  exopar.num_elem_blk = 1;
  exopar.num_node_sets = 0;
  exopar.num_edge_sets = 0;
  exopar.num_face_sets = 0;
  exopar.num_side_sets = 0;
  exopar.num_elem_sets = 0;
  exopar.num_node_maps = 0;
  exopar.num_edge_maps = 0;
  exopar.num_face_maps = 0;
  exopar.num_elem_maps = 0;

  auto status = ex_put_init_ext(exo_id, &exopar);
  if (status) {
    if (comm_rank==0)
      std::cerr << "Problem putting exodus file parameters, ex_put_init_ext() returned "
          << status << std::endl;
    MPI_Finalize();
    return status;
  }

  // exodus is kind enough to fetch the data in the real type we ask for
  status = ex_put_coord(
      exo_id, coordinates.data(), coordinates.data() + num_verts,
      coordinates.data() + 2 * num_verts);

  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting vertex coordinates to exodus file, "
          << " ex_put_coord() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0)
      std::cout << "Wrote " << num_verts << " vertex coordinates" << std::endl;
  }

  const char *coord_names[3] = {"x", "y", "z"};
  status = ex_put_coord_names (exo_id, const_cast<char**>(coord_names));
  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting coordinate names to exodus file, "
          << " ex_put_coord_names() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }

  ex_entity_id elem_blk_id = 1;
  if (is_int64) 
    status = write<long long>(exo_id, elem_blk_id, cell_vertices, num_verts_per_cell);
  else
    status = write<int>(exo_id, elem_blk_id, cell_vertices, num_verts_per_cell);
  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting elements to exodus file, "
          << " ex_put_block() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0) std::cout << "Wrote " << num_cells << " cells" << std::endl;
  }

  if (is_int64) 
    status = write_mapping<long long>(exo_id, global_cell_id, global_vertex_id);
  else
    status = write_mapping<int>(exo_id, global_cell_id, global_vertex_id);
  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting mapping to exodus file, "
          << " ex_put_num_map() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0) std::cout << "Wrote mapping." << std::endl;
  }
    
  status = ex_close(exo_id);
  if (status) {
    if (comm_rank==0)
      std::cerr
        << "Problem closing exodus file, ex_close() returned " << exo_id
        << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0) std::cout << "Closed exodus file!" << std::endl;
  }

  return MPI_Finalize();

}
