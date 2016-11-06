/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_sp_minimal_config_h
#define flecsi_sp_minimal_config_h

#include "flecsi-sp/geometry/point.h"

///
// \file minimal_basic_types.h
// \authors bergen bird
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace sp {

#ifndef FLECSI_MESH_DIMENSION
#define FLECSI_MESH_DIMENSION 2
#endif // FLECSI_MESH_DIMENSION

///
// \class minimal_config_t minimal_types.h
// \brief minimal_config_t provides...
///
struct minimal_config_t
{
  //--------------------------------------------------------------------------//
  // Define local properties to satisfy mesh_topology requirements.
  //--------------------------------------------------------------------------//

  /// The dimension of the mesh
  static constexpr size_t num_dimensions = FLECSI_MESH_DIMENSION;

  /// The number of domains
  static constexpr size_t num_domains = 1;

  /// Floating-point type
  using real_t = double;

  /// Point type
  using point_t = point<real_t, num_dimensions>;

}; // class minimal_config_t

} // namespace sp
} // namespace flecsi

#endif // flecsi_sp_minimal_config_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
