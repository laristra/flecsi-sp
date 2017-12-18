/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Defines the mesh traits for the burton_mesh_t class.
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsi-sp-config.h>

#include <flecsi/data/data.h>
#include <ristra/geometry/point.h>
#include <ristra/geometry/shapes/geometric_shapes.h>
#include <ristra/math/vector.h>
#include <ristra/utils/fixed_vector.h>

// system includes
#include<bitset>

namespace flecsi_sp {
namespace burton {

////////////////////////////////////////////////////////////////////////////////
/// \brief The mesh traits
/// \tparam N  The number of mesh dimensions.
////////////////////////////////////////////////////////////////////////////////
template< std::size_t N >
struct burton_config_t {

  //! the size type
  using size_t = std::size_t;

  //! The dimension of the burton mesh.
  static constexpr size_t num_dimensions = N;

  //! The number of mesh domains in the burton mesh.
  static constexpr size_t num_domains = 2;

  //! the constant string type
  using const_string_t = flecsi::utils::const_string_t;

  //! the bitfield type
  using bitfield_t = std::bitset<8>;

  //! A type used for loop indexing.
  using counter_t = flecsi_sp::config::counter_t;

  //! The type for floating-point values.
  using real_t = flecsi_sp::config::real_t;

  //! The type for integer values.
  using integer_t = flecsi_sp::config::integer_t;

  //! A point type with real_t data and mesh dimension.
  using point_t = ristra::geometry::point<real_t, num_dimensions>;

  //! A space ("physics") vector type with real_t data and mesh dimension.
  using vector_t = ristra::math::vector<real_t, num_dimensions>;

  //! \brief The locations of different bits that we set as flags
  enum bits : uint8_t
  {
    boundary = 0 //!< the boundary bit
  };

  // the tags type
  using tag_t = uint8_t;
  using tag_list_t = ristra::utils::fixed_vector< tag_t, N*N >;

  //! \brief the shape type
  using shape_t = ristra::geometry::shapes::geometric_shapes_t;

  //! the flecsi id type
  using id_t = flecsi::utils::id_t;

  //! the flecsi mesh topology storage type
  using mesh_storage_t = 
    flecsi::topology::mesh_storage__<num_dimensions, num_domains>;

  //! the flecsi mesh topology type
  using mesh_topology_base_t = 
    flecsi::topology::mesh_topology_base__< mesh_storage_t >;

  //! the base type for the entities
  using mesh_entity_base_t = flecsi::topology::mesh_entity_base__<num_domains>;
  
  //! The flecsi domain connectivity type.
  using connectivity_t = flecsi::topology::domain_connectivity__<num_dimensions>;

};

} // namespace burton
} // namespace flecsale

