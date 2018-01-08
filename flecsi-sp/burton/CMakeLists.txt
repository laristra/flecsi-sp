#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(burton_HEADERS
  burton_config.h
  burton_mesh.h
  burton_mesh_topology.h
  burton_types.h
  burton_corner.h
  burton_element.h
  burton_vertex.h
  burton_wedge.h

  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

set(burton_SOURCES
  burton_types.cc
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

# get a relative path
cinch_add_unit(burton_2d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/2d/square/square_32x32.g
  SOURCES 
    test/burton_2d.cc
    burton_specialization_default_driver.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_RUNTIME_DRIVER}
  DEFINES
    -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
    -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
  POLICY ${FLECSI_SP_UNIT_POLICY}
  LIBRARIES FleCSI-SP
)
