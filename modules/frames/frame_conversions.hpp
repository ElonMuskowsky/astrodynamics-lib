#ifndef FRAME_CONVERSIONS_HPP
#define FRAME_CONVERSIONS_HPP

#include "core/state.hpp"
#include "core/epoch.hpp"

namespace astrodynamics_lib {

// Uses SOFA iauC2t00b (IAU 2000B, ~1 mas). EOP pole offsets and dUT1 zeroed.
CartState icrf_to_itrf(CartState state, UtcJd epoch);

// Uses SOFA iauC2t00b (IAU 2000B, ~1 mas). EOP pole offsets and dUT1 zeroed.
CartState itrf_to_icrf(CartState state, UtcJd epoch);

} // namespace astrodynamics_lib

#endif // FRAME_CONVERSIONS_HPP
