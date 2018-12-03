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
  io_exodus.h
  exodus_definition.h
  PARENT_SCOPE
)

if ( FLECSI_SP_ENABLE_EXODUS )
  cinch_add_unit(flecsi_sp_exodus_definition
    SOURCES test/exodus_definition.cc
    INPUTS 
      ${FLECSI_SP_DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${FLECSI_SP_DATA_DIR}/meshes/2d/voro/voro.g
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.exo
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-tet.exo
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-nfaced.exo
    LIBRARIES
      ${FLECSI_SP_LIBRARIES}
  )


  cinch_add_unit(flecsi_sp_color_exodus
    SOURCES test/color_exodus.cc
    INPUTS
      ${FLECSI_SP_DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${FLECSI_SP_DATA_DIR}/meshes/3d/box/box-hex.exo
    LIBRARIES
      ${FLECSI_SP_LIBRARIES}
    POLICY MPI
    THREADS 5
  )

endif()

