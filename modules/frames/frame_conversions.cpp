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

void TODtoITRF(double epoch, Vector3& posTOD, Vector3& velTOD,
               Vector3& posITRF, Vector3& velITRF) {
    /* Вход:
     *    double     epoch,    UTC as JD
     *    Vector3    posTOD,  [km]
     *    Vector3    velTOD,  [km/s]
     * Выход:
     *    Vector3    posITRF, [km]
     *    Vector3    velITRF, [km/s]
     */

    // EOP (как и в ICRFtoITRF — нулевые)
    double xPole = 0.0, yPole = 0.0, dUT1 = 0.0;

    // UTC -> TAI
    double UTC1 = std::floor(epoch),
           UTC2 = epoch - UTC1,
           TAI1, TAI2;
    iauUtctai(UTC1, UTC2, &TAI1, &TAI2);

    // TAI -> TT
    double TT1, TT2;
    iauTaitt(TAI1, TAI2, &TT1, &TT2);

    // TAI -> UT1
    int iy, im, id;
    double fd;
    iauJd2cal(epoch, 0.0, &iy, &im, &id, &fd);
    double NLS;
    iauDat(iy, im, id, fd, &NLS);
    double DTA = -NLS + dUT1;
    double UT11, UT12;
    iauTaiut1(TAI1, TAI2, DTA, &UT11, &UT12);

    // GAST согласованный с IAU 2000B (как и iauPn00b в (TOD<->ICRF))
    double GAST = iauGst00b(UT11, UT12);

    // Матрица движения полюса (с TIO locator s' от iauSp00)
    double RPOM[3][3];
    iauPom00(xPole, yPole, iauSp00(TT1, TT2), RPOM);

    // Матрица TOD -> ITRF: RT2T = RPOM * R3(GAST) * I
    double RBPN_id[3][3];
    iauIr(RBPN_id);                       // единичная — мы уже в TOD
    double RT2T[3][3];
    iauC2teqx(RBPN_id, GAST, RPOM, RT2T); // собираем RPOM * R3(GAST)

    // Позиция
    double rTOD[3] = { posTOD.x * 1000.0, posTOD.y * 1000.0, posTOD.z * 1000.0 };
    double rTer[3];
    iauRxp(RT2T, rTOD, rTer);
    posITRF.x = rTer[0] / 1000.0;
    posITRF.y = rTer[1] / 1000.0;
    posITRF.z = rTer[2] / 1000.0;

    // Скорость: v_ITRF = R · v_TOD - omega x r_ITRF
    double vTOD[3] = { velTOD.x * 1000.0, velTOD.y * 1000.0, velTOD.z * 1000.0 };
    double vTer[3];
    iauRxp(RT2T, vTOD, vTer);
    velITRF.x = (vTer[0] + Constants::omegaEearth * rTer[1]) / 1000.0;
    velITRF.y = (vTer[1] - Constants::omegaEearth * rTer[0]) / 1000.0;
    velITRF.z =  vTer[2] / 1000.0;
}

void ITRFtoTOD(double epoch, Vector3& posITRF, Vector3& velITRF,
               Vector3& posTOD,  Vector3& velTOD) {
    /* Вход:
     *    double     epoch,    UTC as JD
     *    Vector3    posITRF, [km]
     *    Vector3    velITRF, [km/s]
     * Выход:
     *    Vector3    posTOD,  [km]
     *    Vector3    velTOD,  [km/s]
     */

    double xPole = 0.0, yPole = 0.0, dUT1 = 0.0;

    // UTC -> TAI
    double UTC1 = std::floor(epoch),
           UTC2 = epoch - UTC1,
           TAI1, TAI2;
    iauUtctai(UTC1, UTC2, &TAI1, &TAI2);

    // TAI -> TT
    double TT1, TT2;
    iauTaitt(TAI1, TAI2, &TT1, &TT2);

    // TAI -> UT1
    int iy, im, id;
    double fd;
    iauJd2cal(epoch, 0.0, &iy, &im, &id, &fd);
    double NLS;
    iauDat(iy, im, id, fd, &NLS);
    double DTA = -NLS + dUT1;
    double UT11, UT12;
    iauTaiut1(TAI1, TAI2, DTA, &UT11, &UT12);

    // GAST + полюс
    double GAST = iauGst00b(UT11, UT12);
    double RPOM[3][3];
    iauPom00(xPole, yPole, iauSp00(TT1, TT2), RPOM);

    // Та же матрица TOD -> ITRF, обратное преобразование — её транспонированная
    double RBPN_id[3][3];
    iauIr(RBPN_id);
    double RT2T[3][3];
    iauC2teqx(RBPN_id, GAST, RPOM, RT2T);

    // Позиция: r_TOD = R^T · r_ITRF
    double rTer[3] = { posITRF.x * 1000.0, posITRF.y * 1000.0, posITRF.z * 1000.0 };
    double rTODm[3];
    iauTrxp(RT2T, rTer, rTODm);
    posTOD.x = rTODm[0] / 1000.0;
    posTOD.y = rTODm[1] / 1000.0;
    posTOD.z = rTODm[2] / 1000.0;

    // Скорость: v_TOD = R^T · (v_ITRF + omega x r_ITRF)
    double vITRFm[3] = { velITRF.x * 1000.0, velITRF.y * 1000.0, velITRF.z * 1000.0 };
    double vTer[3];
    vTer[0] = (vITRFm[0] - Constants::omegaEearth * rTer[1]) / 1000.0;
    vTer[1] = (vITRFm[1] + Constants::omegaEearth * rTer[0]) / 1000.0;
    vTer[2] =  vITRFm[2] / 1000.0;

    double vTODm[3];
    iauTrxp(RT2T, vTer, vTODm);
    velTOD.x = vTODm[0];
    velTOD.y = vTODm[1];
    velTOD.z = vTODm[2];
}

} // namespace astrodynamics_lib
