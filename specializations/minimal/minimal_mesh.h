/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_specializations_minimal_mesh_h
#define flecsi_specializations_minimal_mesh_h

#include "flecsi-specializations/minimal/minimal_types.h"

///
// \file minimal_mesh.h
// \authors bergen
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace specializations {

///
// \class minimal_mesh_t minimal_mesh.h
// \brief minimal_mesh_t provides...
///
class minimal_mesh_t
  : topology::mesh_topology_t<minimal_types_t>
{
public:

  /// Default constructor
  minimal_mesh_t() {}

  /// Copy constructor (disabled)
  minimal_mesh_t(const minimal_mesh_t &) = delete;

  /// Assignment operator (disabled)
  minimal_mesh_t & operator = (const minimal_mesh_t &) = delete;

  /// Destructor
   ~minimal_mesh_t() {}

private:

}; // class minimal_mesh_t

} // namespace specializations
} // namespace flecsi

#endif // flecsi_specializations_minimal_mesh_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
