#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(burton_HEADERS
  burton_config.h
  burton_corner.h
  burton_element.h
  burton_mesh.h
  burton_mesh_topology.h
  burton_specialization_init.h
  burton_specialization_arguments.h
  burton_types.h
  burton_vertex.h
  burton_wedge.h

  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

set(burton_SOURCES
  burton_types.cc
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

# install common files
install(FILES burton_specialization_init.cc DESTINATION ${FLECSI_SP_SHARE_DIR})
    
cinch_add_unit(flecsi_sp_burton_2d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/2d/square/square_32x32.g
  SOURCES 
    test/burton.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
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

cinch_add_unit(flecsi_sp_burton_3d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.g
  SOURCES 
    test/burton.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
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

cinch_add_unit(flecsi_sp_larsen_thin_2d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/2d/larsen-thin.g
  SOURCES 
    test/larsen-thin.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
  DEFINES
    -DFLECSI_SP_BURTON_MESH_DIMENSION=2
    -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
    -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
    -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
  POLICY ${FLECSI_SP_UNIT_POLICY}
  THREADS 4
  LIBRARIES FleCSI-SP
)

