/*
  test_kernels.cpp
  -----------------
  Quick sanity check for the Voyager 1 / DSN CSPICE kernel pool.

  What it checks, layer by layer:

    1. DSS-43 (Canberra, 399043) geodetic position
       -> exercises: leapseconds, pck00011.tpc (Earth radii),
                      earthstns_itrf93_201023.bsp, earth_assoc_itrf93.tf,
                      and whichever Earth-orientation .bpc is loaded.
       Known truth: DSS-43 is at roughly lat -35.4 deg, lon 148.9 deg E,
       altitude a few hundred meters.

    2. Voyager 1 -> DSS-43 range/range-rate at the Jupiter encounter
       (1979-03-05) -> this ONLY works if a *historical* Earth
       orientation kernel is loaded (earth_620120_*.bpc or the
       combined kernel). With only earth_latest_high_prec.bpc
       (coverage starts 2000-01-01) this call will fail with a
       SPICE(NOFRAMECONNECT)-style error.
       Sanity check: Voyager 1 was about 5-6 AU from Earth at Jupiter.

    3. Voyager 1 -> DSS-43 range/range-rate "today" -> exercises the
       recent end of the merged trajectory file and the
       high-precision/predict EOP coverage.
       Sanity check: Voyager 1 is currently ~165+ AU out, light time
       roughly 22-24 hours.

  Compile (adjust paths to your CSPICE install):
    g++ -std=c++17 test_kernels.cpp -I/path/to/cspice/include \
        -L/path/to/cspice/lib -lcspice -o test_kernels

  Run:
    ./test_kernels

  Note: CSPICE's default error action is to print a diagnostic and
  abort, which is exactly what you want for a first "does this load"
  check -- if a kernel is missing or coverage is insufficient, you'll
  get a clear SPICE error message telling you which call failed.
*/

#include <iostream>
#include <iomanip>
#include "SpiceUsr.h"

// Point this at your meta-kernel (the \begindata/\begintext file
// listing naif0012.tls, pck00011.tpc, de442.bsp, the Voyager SPKs,
// gm_de440.tpc, earthstns_itrf93_201023.bsp, the Earth orientation
// .bpc, and earth_assoc_itrf93.tf).
static const char* METAKR = "../Kernels.tm";

static const double KM_PER_AU = 1.495978707e8;
static const char*  VOYAGER1  = "-31";       // NAIF ID, Voyager 1
static const char*  DSS43     = "399043";    // Canberra 70m

void checkStationGeometry()
{
    SpiceDouble et;
    str2et_c("2026-06-21T00:00:00", &et);

    SpiceDouble pos[3], lt;
    spkpos_c(DSS43, et, "ITRF93", "NONE", "EARTH", pos, &lt);

    SpiceInt dim;
    SpiceDouble radii[3];
    bodvrd_c("EARTH", "RADII", 3, &dim, radii);
    SpiceDouble re = radii[0];
    SpiceDouble rp = radii[2];
    SpiceDouble f  = (re - rp) / re;

    SpiceDouble lon, lat, alt;
    recgeo_c(pos, re, f, &lon, &lat, &alt);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[1] DSS-43 geodetic position (ITRF93 chain)\n"
              << "    lat = " << lat * dpr_c() << " deg   (expect ~ -35.4)\n"
              << "    lon = " << lon * dpr_c() << " deg E (expect ~ 148.9)\n"
              << "    alt = " << alt            << " km    (expect ~ 0.7)\n\n";
}

void checkVoyagerRange(const char* label, const char* utc)
{
    SpiceDouble et;
    str2et_c(utc, &et);

    SpiceDouble state[6], lt;
    spkezr_c(VOYAGER1, et, "J2000", "NONE", DSS43, state, &lt);

    SpiceDouble range     = vnorm_c(state);
    SpiceDouble rangeRate = vdot_c(state, state + 3) / range;

    std::cout << "[" << label << "] epoch " << utc << "\n"
              << "    range      = " << range / KM_PER_AU << " AU\n"
              << "    range-rate = " << rangeRate          << " km/s\n"
              << "    light time = " << lt / 3600.0         << " hours\n\n";
}

int main()
{
    furnsh_c(METAKR);

    checkStationGeometry();
    checkVoyagerRange("2] Jupiter encounter", "1979-03-05T00:00:00");
    checkVoyagerRange("3] Present day       ", "2026-06-21T00:00:00");

    std::cout << "All calls completed without a SPICE error -> "
                 "kernel pool covers both epochs correctly.\n";

    return 0;
}