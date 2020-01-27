#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set_property(SOURCE burton_specialization_init.cc APPEND PROPERTY COMPILE_FLAGS "-xcuda --cuda-gpu-arch=sm_70" )
# burton requires parmetis
if (NOT FLECSI_SP_ENABLE_PARMETIS)
  return()
endif()

set(burton_HEADERS
  burton_config.h
  burton_element.h
  burton_extras.h
  burton_mesh.h
  burton_mesh_topology.h
  burton_specialization_init.h
  burton_specialization_arguments.h
  burton_types.h
  burton_vertex.h

  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

set(burton_SOURCES
  burton_types.cc
  PARENT_SCOPE # THIS NEEDS TO BE HERE
)

# install common files
install(FILES burton_specialization_init.cc DESTINATION ${FLECSI_SP_SHARE_DIR})

# cinch_add_unit(flecsi_sp_burton_1d
#   ARGUMENTS
#   -m ${FLECSI_SP_DATA_DIR}/meshes/1d/uniform/uniform_32.exo
#   SOURCES
#     test/burton.cc
#     burton_specialization_init.cc
#     ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
#     ${FLECSI_SP_RUNTIME_DRIVER}
#   INPUTS
#    test/burton_1d_dump_rank0.out.std
#    test/burton_1d_dump_rank1.out.std
#    test/burton_1d_connectivity_rank0.out.std
#    test/burton_1d_connectivity_rank1.out.std
#    test/burton_1d_geometry_rank0.out.std
#    test/burton_1d_geometry_rank1.out.std
#   DEFINES
#     -DFLECSI_SP_BURTON_MESH_DIMENSION=1
#     -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
#     -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
#     -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
#   POLICY ${FLECSI_SP_UNIT_POLICY}
#   THREADS 2
#   LIBRARIES FleCSI-SP
# )
# 
# cinch_add_unit(flecsi_sp_burton_extras_1d
#   ARGUMENTS
#   -m ${FLECSI_SP_DATA_DIR}/meshes/1d/uniform/uniform_32.exo
#   SOURCES
#     test/burton.cc
#     burton_specialization_init.cc
#     ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
#     ${FLECSI_SP_RUNTIME_DRIVER}
#   INPUTS
#    test/burton_extras_1d_dump_rank0.out.std
#    test/burton_extras_1d_dump_rank1.out.std
#    test/burton_extras_1d_connectivity_rank0.out.std
#    test/burton_extras_1d_connectivity_rank1.out.std
#    test/burton_extras_1d_geometry_rank0.out.std
#    test/burton_extras_1d_geometry_rank1.out.std
#   DEFINES
#     -DFLECSI_SP_BURTON_MESH_DIMENSION=1
#     -DFLECSI_SP_BURTON_MESH_EXTRAS
#     -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
#     -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
#     -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
#   POLICY ${FLECSI_SP_UNIT_POLICY}
#   THREADS 2
#   LIBRARIES FleCSI-SP
# )


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

cinch_add_unit(flecsi_sp_burton_extras_2d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/2d/mixed/mixed.g
  SOURCES 
    test/burton.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
  INPUTS
   test/burton_extras_2d_dump_rank0.out.std
   test/burton_extras_2d_dump_rank1.out.std
   test/burton_extras_2d_connectivity_rank0.out.std
   test/burton_extras_2d_connectivity_rank1.out.std
   test/burton_extras_2d_geometry_rank0.out.std
   test/burton_extras_2d_geometry_rank1.out.std
  DEFINES
    -DFLECSI_SP_BURTON_MESH_DIMENSION=2
    -DFLECSI_SP_BURTON_MESH_EXTRAS
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

cinch_add_unit(flecsi_sp_burton_extras_3d
  ARGUMENTS
  -m ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.g
  SOURCES 
    test/burton.cc
    burton_specialization_init.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
  INPUTS
   test/burton_extras_3d_dump_rank0.out.std
   test/burton_extras_3d_dump_rank1.out.std
   test/burton_extras_3d_connectivity_rank0.out.std
   test/burton_extras_3d_connectivity_rank1.out.std
   test/burton_extras_3d_geometry_rank0.out.std
   test/burton_extras_3d_geometry_rank1.out.std
  DEFINES
    -DFLECSI_SP_BURTON_MESH_DIMENSION=3
    -DFLECSI_SP_BURTON_MESH_EXTRAS
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

if (FLECSI_SP_ENABLE_PORTAGE)
  cinch_add_unit(flecsi_sp_burton_portage_2d
    ARGUMENTS
    -m ${FLECSI_SP_DATA_DIR}/meshes/2d/square/square_32x32.g
    SOURCES
      test/portage.cc
      burton_specialization_init.cc
      ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
      ${FLECSI_SP_RUNTIME_DRIVER}
    DEFINES
      -DFLECSI_SP_BURTON_MESH_DIMENSION=2
      -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
      -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
      -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
    POLICY ${FLECSI_SP_UNIT_POLICY}
    THREADS 2
    LIBRARIES FleCSI-SP
  )

  cinch_add_unit(flecsi_sp_burton_portage_tangram_2d
    ARGUMENTS
    -m ${FLECSI_SP_DATA_DIR}/meshes/2d/square/square_32x32.g
    SOURCES
      test/portage_tangram.cc
      burton_specialization_init.cc
      ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
      ${FLECSI_SP_RUNTIME_DRIVER}
    DEFINES
      -DFLECSI_SP_BURTON_MESH_DIMENSION=2
      -DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT
      -DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
      -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
    POLICY ${FLECSI_SP_UNIT_POLICY}
    THREADS 1
    LIBRARIES FleCSI-SP
  )
endif()
