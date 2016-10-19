#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(minimal_HEADERS
  minimal_types.h
  minimal_mesh_properties.h
  minimal_entity_types.h
  minimal_mesh.h
  PARENT_SCOPE
)

set(minimal_SOURCES
  PARENT_SCOPE
)

cinch_add_unit(properties
 SOURCES test/properties.cc
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
