#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(example_HEADERS
  utils.h
  PARENT_SCOPE
)

set(example_SOURCES
  utils.cc
  PARENT_SCOPE
)

cinch_add_unit(sanity
  SOURCES test/unit.cc
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#