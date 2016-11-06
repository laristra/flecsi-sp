#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(minimal_HEADERS
  minimal_config.h
  minimal_entity_types.h
  minimal_mesh.h
  minimal_types.h
  PARENT_SCOPE
)

set(minimal_SOURCES
  PARENT_SCOPE
)

cinch_add_unit(minimal
  SOURCES test/minimal.cc
  LIBRARIES ${FleCSI_LIBRARY}
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
