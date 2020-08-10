#ifndef FLECSI_SP_BURTON_PORTAGE_HELPERS
#define FLECSI_SP_BURTON_PORTAGE_HELPERS

#include <flecsi-sp/burton/portage_mesh_wrapper.h>
#include <flecsi-sp/burton/portage_state_wrapper.h>
#include <flecsi-sp/burton/portage_mm_state_wrapper.h>

#include <portage/driver/uberdriver.h>
#include <tangram/reconstruct/VOF.h>
#include <tangram/intersect/split_r3d.h>
#include <tangram/intersect/split_r2d.h>

namespace flecsi_sp {
namespace burton {

template<
	int DIM,
  typename mesh_wrapper_a_t,
  typename state_wrapper_a_t,
  typename mesh_wrapper_b_t,
  typename state_wrapper_b_t,
	std::enable_if_t<(DIM==2)>* = nullptr
>
auto make_remapper(
    mesh_wrapper_a_t & mesh_wrapper_a,
    state_wrapper_a_t & state_wrapper_a,
    mesh_wrapper_b_t & mesh_wrapper_b,
    state_wrapper_b_t & state_wrapper_b)
{

    Portage::UberDriver<
      mesh_t::num_dimensions,
      mesh_wrapper_a_t,
      state_wrapper_a_t,
      mesh_wrapper_b_t,
      state_wrapper_b_t,
      Tangram::VOF,
      Tangram::SplitR2D,
      Tangram::ClipR2D > remapper(
          mesh_wrapper_a,
          state_wrapper_a,
          mesh_wrapper_b,
          state_wrapper_b);
    return std::move(remapper);
}

template<
	int DIM,
  typename mesh_wrapper_a_t,
  typename state_wrapper_a_t,
  typename mesh_wrapper_b_t,
  typename state_wrapper_b_t,
	typename = typename std::enable_if_t<(DIM==3)>
>
auto make_remapper(
    mesh_wrapper_a_t & mesh_wrapper_a,
    state_wrapper_a_t & state_wrapper_a,
    mesh_wrapper_b_t & mesh_wrapper_b,
    state_wrapper_b_t & state_wrapper_b)
{
    Portage::UberDriver<
      mesh_t::num_dimensions,
      mesh_wrapper_a_t,
      state_wrapper_a_t,
      mesh_wrapper_b_t,
      state_wrapper_b_t,
      Tangram::VOF,
      Tangram::SplitR3D,
      Tangram::ClipR3D > remapper(
          mesh_wrapper_a,
          state_wrapper_a,
          mesh_wrapper_b,
          state_wrapper_b);
    return std::move(remapper);
}

template< typename T, typename U, typename V, typename W>
auto distrubute_mesh(
 	const T & sourceMeshWrapper,
  const U & sourceStateWrapper,
 	const V & targetMeshWrapper,
  const W & targetStateWrapper,
	const std::vector<std::string> & varNames)
{ 
  // create the flat mesh
  auto source_mesh_flat = std::make_unique<Wonton::Flat_Mesh_Wrapper<>>();
  source_mesh_flat->initialize(sourceMeshWrapper);

  // create the flat state
  auto source_state_flat =
    std::make_unique<Wonton::Flat_State_Wrapper<Wonton::Flat_Mesh_Wrapper<>>>(*source_mesh_flat);

  // mat_volfracs and mat_centroids are always imported from the state wrapper
  source_state_flat->initialize(sourceStateWrapper, varNames);

  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);

  // Use a bounding box distributor to send the source cells to the target
  // partitions where they are needed
  Portage::MPI_Bounding_Boxes distributor(&mpiexecutor);
  distributor.distribute(*source_mesh_flat, *source_state_flat, targetMeshWrapper, 
		targetStateWrapper);
	return std::move(std::make_pair(std::move(source_mesh_flat), std::move(source_state_flat)));
}


template<int D, typename T>
std::enable_if_t<(D==2)> compute_weights(T & remapper)
{
	remapper.template compute_interpolation_weights<Portage::SearchKDTree, Portage::IntersectR2D>();
}

template<int D, typename T, std::enable_if_t<(D==3)>* = nullptr >
void compute_weights(T & remapper)
{
	remapper.template compute_interpolation_weights<Portage::SearchKDTree, Portage::IntersectR3D>();
}


} // namespace
} // namespace


#endif // FLECSI_SP_BURTON_PORTAGE_HELPERS
