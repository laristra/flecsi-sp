#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

file(GLOB headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/detail/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/common/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/geom/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/geom/shapes/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/math/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/math/detail/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/std/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/detail/*.h
)

file(GLOB sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/detail/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/common/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/geom/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/geom/shapes/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/math/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/math/detail/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/std/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/detail/*.cc
)

set(burton_HEADERS
  ${headers}
  PARENT_SCOPE
)


set(burton_SOURCES
  ${sources} 
  PARENT_SCOPE
)


# get a relative path
file(RELATIVE_PATH DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ${FleCSI_SP_DATA_DIR})

cinch_add_unit(burton
    SOURCES 
      test/burton_test_base.h
      test/burton_create_test.h
      test/burton_2d_test.h
      test/burton_3d_test.h
      test/burton_io_test.h
      test/burton_voro_test.h

      test/burton_create.cc
      test/burton_2d.cc
      test/burton_3d.cc
      test/burton_io.cc
      test/burton_voro.cc
      test/burton_grad.cc
    INPUTS 
      test/burton_2d.dump.blessed
      test/burton_2d.mesh.blessed
      test/burton_2d.geometry.blessed
      test/burton_2d.copy.blessed
      test/burton_3d.dump.blessed
      test/burton_3d.mesh.blessed
      test/burton_3d.geometry.blessed
      test/burton_create.faces.blessed
      test/burton_create.minimal.blessed
      test/burton_create.points.blessed
      ${DATA_DIR}/meshes/2d/mixed/mixed.g
      ${DATA_DIR}/meshes/2d/mixed/mixed.exo
      ${DATA_DIR}/meshes/2d/mixed/mixed.vtu
      ${DATA_DIR}/meshes/2d/mixed/mixed.vtm
      ${DATA_DIR}/meshes/2d/mixed/mixed_0_0.vtu
      ${DATA_DIR}/meshes/2d/mixed/mixed_1_0.vtu
      ${DATA_DIR}/meshes/2d/mixed/mixed.vtk
      ${DATA_DIR}/meshes/2d/voro/voro.g
      ${DATA_DIR}/meshes/3d/box/box-hex.exo
      ${DATA_DIR}/meshes/3d/box/box-hex.g
      ${DATA_DIR}/meshes/3d/box/box-nfaced.exo
      ${DATA_DIR}/meshes/3d/box/box-tet.exo
      ${DATA_DIR}/meshes/3d/box/box-tet.g
      ${DATA_DIR}/meshes/3d/voro/sedov-0_05cm.exo
      LIBRARIES ${FleCSI_SP_LIBRARIES}
)
