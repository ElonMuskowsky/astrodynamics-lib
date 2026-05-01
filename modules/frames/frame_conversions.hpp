#ifndef FRAME_CONVERSIONS_HPP
#define FRAME_CONVERSIONS_HPP

#include "core/state.hpp"
#include "core/epoch.hpp"

namespace astrodynamics_lib {

CartState icrf_to_itrf(CartState state, UtcJd epoch);
CartState itrf_to_icrf(CartState state, UtcJd epoch);

} // namespace astrodynamics_lib

#endif // FRAME_CONVERSIONS_HPP
