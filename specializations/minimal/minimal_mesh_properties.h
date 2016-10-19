/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_specializations_minimal_mesh_properties_h
#define flecsi_specializations_minimal_mesh_properties_h

///
// \file minimal_mesh_properties.h
// \authors bergen bird
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace specializations {

#ifndef FLECSI_MESH_DIMENSION
#define FLECSI_MESH_DIMENSION 2
#endif // FLECSI_MESH_DIMENSION

struct minimal_mesh_properties_t {

  /// The dimension of the burton mesh.
  static constexpr size_t num_dimensions = FLECSI_MESH_DIMENSION;

  /// The number of mesh domains.
  static constexpr size_t num_domains = 1;

}; // minimal_mesh_properties_t

} // namespace specializations
} // namespace flecsi

#endif // flecsi_specializations_minimal_mesh_properties_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
