#pragma once

#include "problemsize.h"
#include "mod_debug.h"
#include "mod_thrmdyn.h"
#include "mod_vadv1d.h"

using namespace PROBLEM_SIZE;
using namespace DEBUG;
using namespace THRMDYN;
using namespace VADV1D;

namespace PRECIP_TRANSPORT {

void precip_transport_new(  double rhog               [kdim][ijdim],
                            double rhogvx             [kdim][ijdim],
                            double rhogvy             [kdim][ijdim],
                            double rhogvz             [kdim][ijdim],
                            double rhogw              [kdim][ijdim],
                            double rhoge              [kdim][ijdim],
                            double rhogq              [nqmax][kdim][ijdim],
                            double rho                [kdim][ijdim],
                            double tem                [kdim][ijdim],
                            double pre                [kdim][ijdim],
                            double vx                 [kdim][ijdim],
                            double vy                 [kdim][ijdim],
                            double vz                 [kdim][ijdim],
                            double w                  [kdim][ijdim],
                            double q                  [nqmax][kdim][ijdim],
                            double qd                 [kdim][ijdim],
                            double z                  [kdim][ijdim],
                            double Vterm              [nqmax][kdim][ijdim],
                            bool   precipitating_flag [nqmax],
                            double precip             [2][ijdim],
                            double precip_rhoe        [ijdim],
                            double precip_lh_heat     [ijdim],
                            double precip_rhophi      [ijdim],
                            double precip_rhokin      [ijdim],
                            double frain              [kdim][ijdim],
                            double gsgam2             [kdim][ijdim],
                            double gsgam2h            [kdim][ijdim],
                            double rgs                [kdim][ijdim],
                            double rgsh               [kdim][ijdim],
                            double ix                 [ijdim],
                            double iy                 [ijdim],
                            double iz                 [ijdim],
                            double jx                 [ijdim],
                            double jy                 [ijdim],
                            double jz                 [ijdim],
                            double dt,
                            double **precip_trc      // precip[nqmax][ijdim]
                            );

}
