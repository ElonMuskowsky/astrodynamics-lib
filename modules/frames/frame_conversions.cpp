#include "frames/frame_conversions.hpp"

#include <cmath>
#include "core/constants.hpp"
#include "sofa.h"

namespace astrodynamics_lib {

static void build_rc2t(UtcJd epoch, double rc2t[3][3]) {
    double utc1 = std::floor(epoch.jd);
    double utc2 = epoch.jd - utc1;

    double tai1, tai2;
    iauUtctai(utc1, utc2, &tai1, &tai2);

    double tt1, tt2;
    iauTaitt(tai1, tai2, &tt1, &tt2);

    int iy, im, id;
    double fd;
    iauJd2cal(epoch.jd, 0.0, &iy, &im, &id, &fd);
    double nls;
    iauDat(iy, im, id, fd, &nls);

    // EOP zeroed — no table lookup
    double dut1 = 0.0;
    double ut11, ut12;
    iauTaiut1(tai1, tai2, dut1 - nls, &ut11, &ut12);

    iauC2t00b(tt1, tt2, ut11, ut12, 0.0, 0.0, rc2t);
}

CartState icrf_to_itrf(CartState state, UtcJd epoch) {
    double rc2t[3][3];
    build_rc2t(epoch, rc2t);

    double r_eci[3] = {state.pos.x, state.pos.y, state.pos.z};
    double v_eci[3] = {state.vel.x, state.vel.y, state.vel.z};

    double r_ecef[3], v_rot[3];
    iauRxp(rc2t, r_eci, r_ecef);
    iauRxp(rc2t, v_eci, v_rot);

    return CartState{
        {r_ecef[0], r_ecef[1], r_ecef[2]},
        {v_rot[0] + omega_earth * r_ecef[1],
         v_rot[1] - omega_earth * r_ecef[0],
         v_rot[2]}
    };
}

CartState itrf_to_icrf(CartState state, UtcJd epoch) {
    double rc2t[3][3];
    build_rc2t(epoch, rc2t);

    double r_ecef[3] = {state.pos.x, state.pos.y, state.pos.z};
    double v_ecef[3] = {state.vel.x, state.vel.y, state.vel.z};

    double v_corr[3] = {
        v_ecef[0] - omega_earth * r_ecef[1],
        v_ecef[1] + omega_earth * r_ecef[0],
        v_ecef[2]
    };

    double r_eci[3], v_eci[3];
    iauTrxp(rc2t, r_ecef, r_eci);
    iauTrxp(rc2t, v_corr, v_eci);

    return CartState{
        {r_eci[0], r_eci[1], r_eci[2]},
        {v_eci[0], v_eci[1], v_eci[2]}
    };
}

static void build_rt2t(UtcJd epoch, double rt2t[3][3]) {
    double utc1 = std::floor(epoch.jd);
    double utc2 = epoch.jd - utc1;

    double tai1, tai2;
    iauUtctai(utc1, utc2, &tai1, &tai2);

    double tt1, tt2;
    iauTaitt(tai1, tai2, &tt1, &tt2);

    int iy, im, id;
    double fd;
    iauJd2cal(epoch.jd, 0.0, &iy, &im, &id, &fd);
    double nls;
    iauDat(iy, im, id, fd, &nls);

    // EOP zeroed — no table lookup
    double ut11, ut12;
    iauTaiut1(tai1, tai2, -nls, &ut11, &ut12);

    double gast = iauGst00b(ut11, ut12);

    double rpom[3][3];
    iauPom00(0.0, 0.0, iauSp00(tt1, tt2), rpom);

    // TOD -> ITRF: R = RPOM * R3(GAST), identity BPN since we start in TOD
    double rbpn_id[3][3];
    iauIr(rbpn_id);
    iauC2teqx(rbpn_id, gast, rpom, rt2t);
}

CartState tod_to_itrf(CartState state, UtcJd epoch) {
    double rt2t[3][3];
    build_rt2t(epoch, rt2t);

    double r_tod[3] = {state.pos.x, state.pos.y, state.pos.z};
    double v_tod[3] = {state.vel.x, state.vel.y, state.vel.z};

    double r_itrf[3], v_rot[3];
    iauRxp(rt2t, r_tod, r_itrf);
    iauRxp(rt2t, v_tod, v_rot);

    return CartState{
        {r_itrf[0], r_itrf[1], r_itrf[2]},
        {v_rot[0] + omega_earth * r_itrf[1],
         v_rot[1] - omega_earth * r_itrf[0],
         v_rot[2]}
    };
}

CartState itrf_to_tod(CartState state, UtcJd epoch) {
    double rt2t[3][3];
    build_rt2t(epoch, rt2t);

    double r_itrf[3] = {state.pos.x, state.pos.y, state.pos.z};
    double v_itrf[3] = {state.vel.x, state.vel.y, state.vel.z};

    double v_corr[3] = {
        v_itrf[0] - omega_earth * r_itrf[1],
        v_itrf[1] + omega_earth * r_itrf[0],
        v_itrf[2]
    };

    double r_tod[3], v_tod[3];
    iauTrxp(rt2t, r_itrf, r_tod);
    iauTrxp(rt2t, v_corr, v_tod);

    return CartState{
        {r_tod[0], r_tod[1], r_tod[2]},
        {v_tod[0], v_tod[1], v_tod[2]}
    };
}

} // namespace astrodynamics_lib
