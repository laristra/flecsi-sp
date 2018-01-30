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

set(FLECSI_SP_BURTON_SPECIALIZATION_INIT
  ${FLECSI_SP_SHARE_DIR}/burton_specialization_init.cc PARENT_SCOPE)
install(FILES burton_specialization_init.cc DESTINATION ${FLECSI_SP_SHARE_DIR})

cinch_add_unit(burton_2d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/2d/square/square_32x32.g
  SOURCES 
    test/burton.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_RUNTIME_DRIVER}
  INPUTS
   test/burton_2d_dump_rank0.out.std
   test/burton_2d_dump_rank1.out.std
   test/burton_2d_connectivity_rank0.out.std
   test/burton_2d_connectivity_rank1.out.std
   test/burton_2d_geometry_rank0.out.std
   test/burton_2d_geometry_rank1.out.std
  DEFINES
    -DFLECSI_SP_BURTON_MESH_DIMENSION=2
    -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
    -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
  POLICY ${FLECSI_SP_UNIT_POLICY}
  THREADS 2
  LIBRARIES FleCSI-SP
)

cinch_add_unit(burton_3d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.g
  SOURCES 
    test/burton.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_RUNTIME_DRIVER}
  INPUTS
   test/burton_3d_dump_rank0.out.std
   test/burton_3d_dump_rank1.out.std
   test/burton_3d_connectivity_rank0.out.std
   test/burton_3d_connectivity_rank1.out.std
   test/burton_3d_geometry_rank0.out.std
   test/burton_3d_geometry_rank1.out.std
  DEFINES
    -DFLECSI_SP_BURTON_MESH_DIMENSION=3
    -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
    -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
  POLICY ${FLECSI_SP_UNIT_POLICY}
  THREADS 2
  LIBRARIES FleCSI-SP
)
