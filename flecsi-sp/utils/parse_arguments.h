/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utilities for parsing arguments.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include <ristra/assertions/errors.h>
#include <ristra/utils/string_utils.h>

// system libraries
#include <cstring>
#include <getopt.h>
#include <map>
#include <string>

namespace flecsi_sp {
namespace utils {


///////////////////////////////////////////////////////////////////////////////
//! \brief Parse the argument list
///////////////////////////////////////////////////////////////////////////////
static auto parse_arguments(
  int argc, 
  char** argv, 
  const option * long_options, 
  const char * short_options
) {

  std::map<std::string, std::string> key_value_pair;

  // reset getopts global variable
  optind = 0;

    
  // getopt_long stores the option index here.
  int option_index = 0;

  // if no arguments, set c to -1 to skip while lop
  int c = ( argc > 1 ) ? 0 : -1;

  // make a copy to avoid getopt permuting invalid options
  std::vector<char*> argvcopy(argv, argv + argc);

  while (c != -1) {
    c = getopt_long(argc, argvcopy.data(), short_options, long_options, &option_index);
    auto c_char = static_cast<char>(c);
    auto c_str = ristra::utils::to_string( c_char );
    // finished with arguments
    if ( c == -1 )
      break;
    // long options that set a flag 
    else if (c == 0)
      key_value_pair[long_options[option_index].name] =
        optarg ? optarg : "";
    // standard short/long option
    else
      key_value_pair[c_str] = optarg ? optarg : "";
    // make sure we have not gone past argc
    if (optind > argc)
      THROW_RUNTIME_ERROR( "Expected argument after options" );
  }

  return key_value_pair;
}

} // namespace
} // namespace
