
#pragma once

namespace flecsi_sp {
namespace utils {

// the handle type
template<typename T>
using dense_handle_w__ =
  flecsi::dense_accessor<T, flecsi::wo, flecsi::wo, flecsi::ro>;

template<typename T>
using dense_handle_rw__ =
  flecsi::dense_accessor<T, flecsi::rw, flecsi::rw, flecsi::ro>;

template<typename T>
using dense_handle_r__ =
  flecsi::dense_accessor<T, flecsi::ro, flecsi::ro, flecsi::ro>;

template<typename DC>
using client_handle_w__ = flecsi::data_client_handle__<DC, flecsi::wo>;

template<typename DC>
using client_handle_r__ = flecsi::data_client_handle__<DC, flecsi::ro>;

} // namespace
} // namespace
