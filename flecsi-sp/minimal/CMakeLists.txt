#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(minimal_HEADERS
  config.h
  entity_types.h
  mesh.h
  types.h
  PARENT_SCOPE
)

set(minimal_SOURCES
  PARENT_SCOPE
)

cinch_add_unit(minimal
  SOURCES test/minimal.cc
  LIBRARIES ${FleCSI_LIBRARIES}
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
