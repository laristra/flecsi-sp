#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(burton_HEADERS
  burton.h

  burton_io_exodus.h
  burton_io_tecplot.h
  burton_io_vtk.h
  burton_io_vtm.h
  burton_io_vtu.h

  burton_config.h
  burton_mesh.h
  burton_mesh_topology.h
  burton_types.h

  burton_corner.h
  burton_element.h
  burton_vertex.h
  burton_wedge.h

  burton_triangle.h
  burton_quadrilateral.h
  burton_polygon.h

  burton_tetrahedron.h
  burton_hexahedron.h
  burton_polyhedron.h

  factory.h
  voronoi.h  detail/voronoi_impl.h
  vtk_utils.h

  PARENT_SCOPE
)


set(burton_SOURCES
  burton_element.cc
  burton_hexahedron.cc
  burton_polygon.cc
  burton_polyhedron.cc
  burton_quadrilateral.cc
  burton_tetrahedron.cc
  burton_triangle.cc
  burton_wedge.cc
  burton_vertex.cc
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
