#------------------------------------------------------------------------------#
#  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
# /@@/////  /@@          @@////@@ @@////// /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@
# /@@       /@@/@@//// //@@    @@       /@@/@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
# //       ///  //////   //////  ////////  //
#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#------------------------------------------------------------------------------#

set(io_HEADERS
  exodus_definition.h
  mpas_definition.h
  io_utils.h
  PARENT_SCOPE
  )

if ( FLECSI_SP_ENABLE_MPAS )
  cinch_add_unit(flecsi_sp_mpas_definition
    SOURCES
    test/mpas_definition.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
    INPUTS
    ${FLECSI_SP_DATA_DIR}/meshes/mpas/init.h5
    DEFINES
    -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
    POLICY ${FLECSI_SP_UNIT_POLICY}
    LIBRARIES
    ${FLECSI_SP_LIBRARIES}
  )
endif()
if ( FLECSI_SP_ENABLE_EXODUS )
  cinch_add_unit(flecsi_sp_exodus_definition
    SOURCES test/exodus_definition.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
    INPUTS
      ${FLECSI_SP_DATA_DIR}/meshes/1d/uniform/uniform_32.exo
      ${FLECSI_SP_DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${FLECSI_SP_DATA_DIR}/meshes/2d/voro/voro.g
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.exo
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-tet.exo
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-nfaced.exo
    DEFINES
    -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
    POLICY ${FLECSI_SP_UNIT_POLICY}
    LIBRARIES
    ${FLECSI_SP_LIBRARIES}
    )


  cinch_add_unit(flecsi_sp_color_exodus
    SOURCES test/color_exodus.cc
    ${FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER}
    ${FLECSI_SP_RUNTIME_DRIVER}
    INPUTS
      ${FLECSI_SP_DATA_DIR}/meshes/1d/uniform/uniform_32.exo
      ${FLECSI_SP_DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.exo
    LIBRARIES
    ${FLECSI_SP_LIBRARIES}
    DEFINES
	  -DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER
    POLICY ${FLECSI_SP_UNIT_POLICY}
    # POLICY MPI
    THREADS 5
    )
endif()
