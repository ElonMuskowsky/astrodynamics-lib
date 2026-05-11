#ifndef FRAME_CONVERSIONS_HPP
#define FRAME_CONVERSIONS_HPP

#include "core/state.hpp"
#include "core/epoch.hpp"

namespace astrodynamics_lib {

// Uses SOFA iauC2t00b (IAU 2000B, ~1 mas). EOP pole offsets and dUT1 zeroed.
CartState icrf_to_itrf(const CartState& state, UtcJd epoch);

// Uses SOFA iauC2t00b (IAU 2000B, ~1 mas). EOP pole offsets and dUT1 zeroed.
CartState itrf_to_icrf(const CartState& state, UtcJd epoch);

// Uses SOFA iauGst00b + iauPom00 (IAU 2000B). EOP pole offsets and dUT1 zeroed.
CartState tod_to_itrf(const CartState& state, UtcJd epoch);

// Uses SOFA iauGst00b + iauPom00 (IAU 2000B). EOP pole offsets and dUT1 zeroed.
CartState itrf_to_tod(const CartState& state, UtcJd epoch);

// Uses SOFA iauPn00b (IAU 2000B). Pure rotation — no Earth-rate correction.
CartState icrf_to_tod(const CartState& state, UtcJd epoch);

// Uses SOFA iauPn00b (IAU 2000B). Pure rotation — no Earth-rate correction.
CartState tod_to_icrf(const CartState& state, UtcJd epoch);

} // namespace astrodynamics_lib

#endif // FRAME_CONVERSIONS_HPP
