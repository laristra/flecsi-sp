/*~--------------------------------------------------------------------------~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_sp_point_h
#define flecsi_sp_point_h

#include "flecsi-sp/math/vector.h"

#include <flecsi/utils/common.h>

#include <cmath>

///
// \file point.h
// \authors bergen
// \date Initial file creation: Sep 23, 2015
///

namespace flecsi {
namespace sp {
namespace geometry {

///
// \class point point.h
// \brief point defines an interface for storing and manipulating
// coordinate data associated with a geometric domain.
//
// The point type is implemented using \ref dimensioned_array.  Look there
// for more information on the point interface.
///
template <typename T, std::size_t D> 
using point = math::vector<T,D>;

///
// \function distance
///
template <typename T, size_t D>
T distance(const point<T, D> & a, const point<T, D> & b)
{
  T sum(0);
  for (size_t d(0); d < D; ++d) {
    sum += flecsi::utils::square(a[d] - b[d]);
  } // for

  return std::sqrt(sum);
} // distance

///
// \function midpoint
///
template <typename T, size_t D>
point<T, D> midpoint(const point<T, D> & a, const point<T, D> & b)
{
  return point<T, D>((a + b) / 2.0);
} // distance

///
// Compute the centroid of a list of points.
//
// \param[in] cell The cell to return the centroid for.
// \return a point that is the centroid.
///
template <template <typename...> class LIST, typename T, size_t D>
auto centroid(const LIST<point<T, D>> & vert_list)
{
  point<T, D> tmp(0.0);
  for (auto v : vert_list)
    tmp += v;
  tmp /= vert_list.size();
  return tmp;
} // centroid

template <typename T, size_t D>
auto centroid(std::initializer_list<point<T, D>> vert_list)
{
  point<T, D> tmp(0.0);
  for (auto v : vert_list)
    tmp += v;
  tmp /= vert_list.size();
  return tmp;
} // centroid

} // namespace sp
} // namespace flecsi
} // namespace flecsi

#endif // flecsi_sp_point_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
