#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(burton_HEADERS
  burton.h
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
  SOURCES 
    test/burton_2d.cc
    ${FLECSI_RUNTIME_DRIVER}
  LIBRARIES FleCSI-SP
)
