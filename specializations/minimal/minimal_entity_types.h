/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_specializations_minimal_entity_types_h
#define flecsi_specializations_minimal_entity_types_h

#include <flecsi/topology/mesh_types.h>

#include "flecsi-specializations/minimal/minimal_mesh_properties.h"

///
// \file minimal_entity_types.h
// \authors bergen bird
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace specializations {

///
// \struct minimal_vertex_t
// \breif FIXME
///
struct minimal_vertex_t
  : topology::mesh_entity_t<0, minimal_mesh_properties_t::num_domains>
{
  ///
  // Constructor.
  //
  // \param mesh FIXME
  ///
  minimal_vertex_t(topology::mesh_topology_base_t & mesh)
    : mesh_(mesh) {}

private:

  topology::mesh_topology_base_t & mesh_;

}; // class minimal_vertex_t

///
// \struct minimal_edge_t
// \breif FIXME
///
struct minimal_edge_t
  : topology::mesh_entity_t<0, minimal_mesh_properties_t::num_domains>
{
  ///
  // Constructor.
  //
  // \param mesh FIXME
  ///
  minimal_edge_t(topology::mesh_topology_base_t & mesh)
    : mesh_(mesh) {}

private:

  topology::mesh_topology_base_t & mesh_;

}; // class minimal_edge_t

///
// \struct minimal_face_t
// \breif FIXME
///
struct minimal_face_t
  : topology::mesh_entity_t<0, minimal_mesh_properties_t::num_domains>
{
  ///
  // Constructor.
  //
  // \param mesh FIXME
  ///
  minimal_face_t(topology::mesh_topology_base_t & mesh)
    : mesh_(mesh) {}

private:

  topology::mesh_topology_base_t & mesh_;

}; // class minimal_face_t

///
// \struct minimal_cell_t
// \breif FIXME
///
struct minimal_cell_t
  : topology::mesh_entity_t<0, minimal_mesh_properties_t::num_domains>
{
  ///
  // Constructor.
  //
  // \param mesh FIXME
  ///
  minimal_cell_t(topology::mesh_topology_base_t & mesh)
    : mesh_(mesh) {}

private:

  topology::mesh_topology_base_t & mesh_;

}; // class minimal_cell_t

} // namespace specializations
} // namespace flecsi

#endif // flecsi_specializations_minimal_entity_types_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
