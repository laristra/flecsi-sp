#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Configure header with version information
#------------------------------------------------------------------------------#

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/header.tex.in
	${CMAKE_BINARY_DIR}/doc/header.tex)

#------------------------------------------------------------------------------#
# Pandoc options for user guide
#------------------------------------------------------------------------------#

set(pandoc_options
    "--toc"
    "--include-in-header=${CMAKE_SOURCE_DIR}/cinch/tex/addtolength.tex"
    "--include-in-header=${CMAKE_BINARY_DIR}/doc/header.tex"
)

#------------------------------------------------------------------------------#
# Add user guide target
#------------------------------------------------------------------------------#

cinch_add_doc(user-guide ugconfig.py src
    user-guide-${${PROJECT_NAME}_VERSION}.pdf
    PANDOC_OPTIONS ${pandoc_options} IMAGE_GLOB "*.pdf")

#------------------------------------------------------------------------------#
# Add developer guide target
#------------------------------------------------------------------------------#

cinch_add_doc(developer-guide dgconfig.py src
    developer-guide-${${PROJECT_NAME}_VERSION}.pdf
    PANDOC_OPTIONS ${pandoc_options} IMAGE_GLOB "*.pdf")

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#