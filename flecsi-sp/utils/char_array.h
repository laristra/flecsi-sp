/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some commonly used utilities.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include <flecsi-sp-config.h>
#include <ristra/utils/trivial_string.h>

namespace flecsi_sp {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! \brief A constexpr, trivially copyable string.
////////////////////////////////////////////////////////////////////////////////
using char_array_t = ristra::utils::trivial_string_u< config::max_char_length >;

///////////////////////////////////////////////////////////////////////////////
//! \brief Convert a std::string to a character array
//! \tparam N  The maximum length of the string.
////////////////////////////////////////////////////////////////////////////////
inline auto to_char_array( const std::string & str )
{
  char_array_t tmp;
  strcpy( tmp.data(), str.c_str() );
  return tmp;
}

} // namespace
} // namespace
