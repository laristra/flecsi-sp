
#pragma once

namespace flecsi_sp {
namespace utils {

// the handle type
#if FLECSI_RUNTIME_MODEL == FLECSI_RUNTIME_MODEL_legion

template<typename T>
using dense_handle_w__ =
  flecsi::data::legion::dense_handle_t<T, flecsi::wo, flecsi::wo, flecsi::ro>;

template<typename T>
using dense_handle_rw__ =
  flecsi::data::legion::dense_handle_t<T, flecsi::rw, flecsi::rw, flecsi::ro>;

template<typename T>
using dense_handle_r__ =
  flecsi::data::legion::dense_handle_t<T, flecsi::ro, flecsi::ro, flecsi::ro>;

#elif FLECSI_RUNTIME_MODEL == FLECSI_RUNTIME_MODEL_mpi

template<typename T>
using dense_handle_w__ =
  flecsi::data::mpi::dense_handle_t<T, flecsi::wo, flecsi::wo, flecsi::ro>;

template<typename T>
using dense_handle_rw__ =
  flecsi::data::mpi::dense_handle_t<T, flecsi::rw, flecsi::rw, flecsi::ro>;

template<typename T>
using dense_handle_r__ =
  flecsi::data::mpi::dense_handle_t<T, flecsi::ro, flecsi::ro, flecsi::ro>;

#else

#error Unknown runtime model.

#endif

template<typename DC>
using client_handle_w__ = flecsi::data_client_handle__<DC, flecsi::wo>;

template<typename DC>
using client_handle_r__ = flecsi::data_client_handle__<DC, flecsi::ro>;

} // namespace
} // namespace
