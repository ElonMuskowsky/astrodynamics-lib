#include <gtest/gtest.h>
#include <cmath>

#include "frames/frame_conversions.hpp"

using namespace astrodynamics_lib;

// ICRF <-> ITRF: IAU 2000B matches STK at ~1 mm level
static constexpr double POS_TOL    = 1e-3;   // km  (~1 m)    — reference tests
static constexpr double VEL_TOL    = 1e-6;   // km/s (~1 mm/s) — reference tests

// TOD <-> ITRF / TOD <-> ICRF: ~40 mas model gap between IAU 2000B and
// STK's IAU 1976/1980 produces ~2 m position and ~2 µm/s velocity offsets
static constexpr double POS_TOL_TOD = 2e-3;  // km  (~2 m)
static constexpr double VEL_TOL_TOD = 5e-6;  // km/s (~5 mm/s)

static constexpr double POS_TOL_RT  = 1e-9;  // km  (~1 µm)   — round-trips
static constexpr double VEL_TOL_RT  = 1e-12; // km/s           — round-trips

// ============================================================
// icrf_to_itrf
// ============================================================

TEST(FrameConversions, icrf_to_itrf_reference_leo)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState eci{{5459.2396713670341342, -3527.2148870837145296, 2214.7311273388722839},
                   {-2.6136887550416064, 0.4567270970933447, 7.1456293866401950}};
    UtcJd epoch{2460676.5};

    CartState ecef = icrf_to_itrf(eci, epoch);

    EXPECT_NEAR(ecef.pos.x, -4468.6325575372111416, POS_TOL);
    EXPECT_NEAR(ecef.pos.y, -4713.5451335432780979, POS_TOL);
    EXPECT_NEAR(ecef.pos.z, 2227.8671314208245349, POS_TOL);
    EXPECT_NEAR(ecef.vel.x, 0.5880513773723662, VEL_TOL);
    EXPECT_NEAR(ecef.vel.y, 2.8283696437551167, VEL_TOL);
    EXPECT_NEAR(ecef.vel.z, 7.1392737686105887, VEL_TOL);
}

TEST(FrameConversions, icrf_to_itrf_position_magnitude_preserved)
{
    // Pure rotation — position magnitude must be invariant.
    CartState eci{{5459.2396713670341342, -3527.2148870837145296, 2214.7311273388722839},
                   {-2.6136887550416064, 0.4567270970933447, 7.1456293866401950}};
    UtcJd epoch{2460676.5};

    CartState ecef = icrf_to_itrf(eci, epoch);

    EXPECT_NEAR(norm(ecef.pos), norm(eci.pos), 1e-9);
}

// ============================================================
// itrf_to_icrf
// ============================================================

TEST(FrameConversions, itrf_to_icrf_reference_leo)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState ecef{{-5344.8007498102224417, -974.8552535676452635, -3706.3868933148032738},
                    {-1.3992290458120542, -6.3277729610962936, 3.6984098055584855}};
    UtcJd epoch{2460676.5};

    CartState eci = itrf_to_icrf(ecef, epoch);

    EXPECT_NEAR(eci.pos.x, 1930.5528299901454830, POS_TOL);
    EXPECT_NEAR(eci.pos.y, -5075.0985731923992716, POS_TOL);
    EXPECT_NEAR(eci.pos.z, -3710.9138856235445019, POS_TOL);
    EXPECT_NEAR(eci.vel.x, 6.8561429466375774, VEL_TOL);
    EXPECT_NEAR(eci.vel.y, -0.0721353352807623, VEL_TOL);
    EXPECT_NEAR(eci.vel.z, 3.6817654595132039, VEL_TOL);
}

// ============================================================
// Round-trips: ICRF <-> ITRF
// ============================================================

TEST(FrameConversions, roundtrip_icrf_itrf_icrf)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState orig{{5459.2396713670341342, -3527.2148870837145296, 2214.7311273388722839},
                   {-2.6136887550416064, 0.4567270970933447, 7.1456293866401950}};
    UtcJd epoch{2460676.5};

    CartState back = itrf_to_icrf(icrf_to_itrf(orig, epoch), epoch);

    EXPECT_NEAR(back.pos.x, orig.pos.x, POS_TOL_RT);
    EXPECT_NEAR(back.pos.y, orig.pos.y, POS_TOL_RT);
    EXPECT_NEAR(back.pos.z, orig.pos.z, POS_TOL_RT);
    EXPECT_NEAR(back.vel.x, orig.vel.x, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.y, orig.vel.y, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.z, orig.vel.z, VEL_TOL_RT);
}

TEST(FrameConversions, roundtrip_itrf_icrf_itrf)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState orig{{-5344.8007498102224417, -974.8552535676452635, -3706.3868933148032738},
                    {-1.3992290458120542, -6.3277729610962936, 3.6984098055584855}};
    UtcJd epoch{2451545.0};  // J2000

    CartState back = icrf_to_itrf(itrf_to_icrf(orig, epoch), epoch);

    EXPECT_NEAR(back.pos.x, orig.pos.x, POS_TOL_RT);
    EXPECT_NEAR(back.pos.y, orig.pos.y, POS_TOL_RT);
    EXPECT_NEAR(back.pos.z, orig.pos.z, POS_TOL_RT);
    EXPECT_NEAR(back.vel.x, orig.vel.x, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.y, orig.vel.y, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.z, orig.vel.z, VEL_TOL_RT);
}

// ============================================================
// tod_to_itrf
// ============================================================

TEST(FrameConversions, tod_to_itrf_reference_leo)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState tod{{5473.4814187764450253, -3496.7393123455508430, 2227.8671314291636918},
                   {-2.6335535085359751, 0.4417622156983962, 7.1392737686073797}};
    UtcJd epoch{2460676.5};

    CartState ecef = tod_to_itrf(tod, epoch);

    EXPECT_NEAR(ecef.pos.x, -4468.6325575361288429, POS_TOL_TOD);
    EXPECT_NEAR(ecef.pos.y, -4713.5451335403531630, POS_TOL_TOD);
    EXPECT_NEAR(ecef.pos.z, 2227.8671314291636918, POS_TOL_TOD);
    EXPECT_NEAR(ecef.vel.x, 0.5880513773790149, VEL_TOL_TOD);
    EXPECT_NEAR(ecef.vel.y, 2.8283696437618220, VEL_TOL_TOD);
    EXPECT_NEAR(ecef.vel.z, 7.1392737686073797, VEL_TOL_TOD);
}

TEST(FrameConversions, tod_to_itrf_position_magnitude_preserved)
{
    // Pure rotation — position magnitude must be invariant.
    CartState tod{{5200.0, -3800.0, 2900.0}, {4.1, 5.6, -3.8}};
    UtcJd epoch{2460676.5};

    CartState ecef = tod_to_itrf(tod, epoch);

    EXPECT_NEAR(norm(ecef.pos), norm(tod.pos), 1e-9);
}

// ============================================================
// itrf_to_tod
// ============================================================

TEST(FrameConversions, itrf_to_tod_reference_leo)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState ecef{{-4468.6325575361288429 , -4713.5451335403558915, 2227.8671314291568706}, {0.5880513773790095, 2.8283696437618140, 7.1392737686073833}};
    UtcJd epoch{2460676.5};

    CartState tod = itrf_to_tod(ecef, epoch);

    EXPECT_NEAR(tod.pos.x, 5473.4814187764477538, POS_TOL_TOD);
    EXPECT_NEAR(tod.pos.y, -3496.7393123455503883, POS_TOL_TOD);
    EXPECT_NEAR(tod.pos.z, 2227.8671314291568706, POS_TOL_TOD);
    EXPECT_NEAR(tod.vel.x, -2.6335535085359663, VEL_TOL_TOD);
    EXPECT_NEAR(tod.vel.y, 0.4417622156983927, VEL_TOL_TOD);
    EXPECT_NEAR(tod.vel.z, 7.1392737686073833, VEL_TOL_TOD);
}

// ============================================================
// Round-trips: TOD <-> ITRF
// ============================================================

TEST(FrameConversions, roundtrip_tod_itrf_tod)
{
    CartState orig{{5200.0, -3800.0, 2900.0}, {4.1, 5.6, -3.8}};
    UtcJd epoch{2460676.5};

    CartState back = itrf_to_tod(tod_to_itrf(orig, epoch), epoch);

    EXPECT_NEAR(back.pos.x, orig.pos.x, POS_TOL_RT);
    EXPECT_NEAR(back.pos.y, orig.pos.y, POS_TOL_RT);
    EXPECT_NEAR(back.pos.z, orig.pos.z, POS_TOL_RT);
    EXPECT_NEAR(back.vel.x, orig.vel.x, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.y, orig.vel.y, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.z, orig.vel.z, VEL_TOL_RT);
}

TEST(FrameConversions, roundtrip_itrf_tod_itrf)
{
    CartState orig{{4200.0, 5100.0, -2700.0}, {-5.3, 3.8, -3.1}};
    UtcJd epoch{2451545.0};  // J2000

    CartState back = tod_to_itrf(itrf_to_tod(orig, epoch), epoch);

    EXPECT_NEAR(back.pos.x, orig.pos.x, POS_TOL_RT);
    EXPECT_NEAR(back.pos.y, orig.pos.y, POS_TOL_RT);
    EXPECT_NEAR(back.pos.z, orig.pos.z, POS_TOL_RT);
    EXPECT_NEAR(back.vel.x, orig.vel.x, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.y, orig.vel.y, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.z, orig.vel.z, VEL_TOL_RT);
}

// ============================================================
// icrf_to_tod
// ============================================================

TEST(FrameConversions, icrf_to_tod_reference_leo)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState eci{{5459.2396713670341342, -3527.2148870837145296, 2214.7311273388722839},
                  {-2.6136887550416064, 0.4567270970933447, 7.1456293866401950}};
    UtcJd epoch{2460676.5};

    CartState tod = icrf_to_tod(eci, epoch);

    EXPECT_NEAR(tod.pos.x, 5473.4814187764477538, POS_TOL_TOD);
    EXPECT_NEAR(tod.pos.y, -3496.7393123455503883, POS_TOL_TOD);
    EXPECT_NEAR(tod.pos.z, 2227.8671314291568706, POS_TOL_TOD);
    EXPECT_NEAR(tod.vel.x, -2.6335535085359663, VEL_TOL_TOD);
    EXPECT_NEAR(tod.vel.y, 0.4417622156983927, VEL_TOL_TOD);
    EXPECT_NEAR(tod.vel.z, 7.1392737686073833, VEL_TOL_TOD);
}

TEST(FrameConversions, icrf_to_tod_position_magnitude_preserved)
{
    // Pure rotation — position magnitude must be invariant.
    CartState eci{{5459.2396713670341342, -3527.2148870837145296, 2214.7311273388722839},
                  {-2.6136887550416064, 0.4567270970933447, 7.1456293866401950}};
    UtcJd epoch{2460676.5};

    CartState tod = icrf_to_tod(eci, epoch);

    EXPECT_NEAR(norm(tod.pos), norm(eci.pos), 1e-9);
}

// ============================================================
// tod_to_icrf
// ============================================================

TEST(FrameConversions, tod_to_icrf_reference_leo)
{
    // Epoch: 2025-01-01 00:00:00 UTC (JD 2460676.5)
    CartState tod{{5473.4814187764477538, -3496.7393123455503883, 2227.8671314291568706},
                  {-2.6335535085359663, 0.4417622156983927, 7.1392737686073833}};
    UtcJd epoch{2460676.5};

    CartState eci = tod_to_icrf(tod, epoch);

    EXPECT_NEAR(eci.pos.x, 5459.2396713639827794, POS_TOL_TOD);
    EXPECT_NEAR(eci.pos.y, -3527.2148870831865679, POS_TOL_TOD);
    EXPECT_NEAR(eci.pos.z, 2214.7311273472118955, POS_TOL_TOD);
    EXPECT_NEAR(eci.vel.x, -2.6136887550494552, VEL_TOL_TOD);
    EXPECT_NEAR(eci.vel.y, 0.4567270970984216, VEL_TOL_TOD);
    EXPECT_NEAR(eci.vel.z, 7.1456293866370082, VEL_TOL_TOD);
}

// ============================================================
// Round-trips: ICRF <-> TOD
// ============================================================

TEST(FrameConversions, roundtrip_icrf_tod_icrf)
{
    CartState orig{{5459.2396713670341342, -3527.2148870837145296, 2214.7311273388722839},
                   {-2.6136887550416064, 0.4567270970933447, 7.1456293866401950}};
    UtcJd epoch{2460676.5};

    CartState back = tod_to_icrf(icrf_to_tod(orig, epoch), epoch);

    EXPECT_NEAR(back.pos.x, orig.pos.x, POS_TOL_RT);
    EXPECT_NEAR(back.pos.y, orig.pos.y, POS_TOL_RT);
    EXPECT_NEAR(back.pos.z, orig.pos.z, POS_TOL_RT);
    EXPECT_NEAR(back.vel.x, orig.vel.x, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.y, orig.vel.y, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.z, orig.vel.z, VEL_TOL_RT);
}

TEST(FrameConversions, roundtrip_tod_icrf_tod)
{
    CartState orig{{5473.4814187764450253, -3496.7393123455508430, 2227.8671314291636918},
                   {-2.6335535085359751, 0.4417622156983962, 7.1392737686073797}};
    UtcJd epoch{2451545.0};  // J2000

    CartState back = icrf_to_tod(tod_to_icrf(orig, epoch), epoch);

    EXPECT_NEAR(back.pos.x, orig.pos.x, POS_TOL_RT);
    EXPECT_NEAR(back.pos.y, orig.pos.y, POS_TOL_RT);
    EXPECT_NEAR(back.pos.z, orig.pos.z, POS_TOL_RT);
    EXPECT_NEAR(back.vel.x, orig.vel.x, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.y, orig.vel.y, VEL_TOL_RT);
    EXPECT_NEAR(back.vel.z, orig.vel.z, VEL_TOL_RT);
}
