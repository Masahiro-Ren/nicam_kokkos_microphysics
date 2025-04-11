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
/**
 * Kokkos ver.
 */
void precip_transport_new(  const View<double**>&  rhog               ,
                            const View<double**>&  rhogvx             ,
                            const View<double**>&  rhogvy             ,
                            const View<double**>&  rhogvz             ,
                            const View<double**>&  rhogw              ,
                            const View<double**>&  rhoge              ,
                            const View<double***>& rhogq              ,
                            const View<double**>&  rho                ,
                            const View<double**>&  tem                ,
                            const View<double**>&  pre                ,
                            const View<double**>&  vx                 ,
                            const View<double**>&  vy                 ,
                            const View<double**>&  vz                 ,
                            const View<double**>&  w                  ,
                            const View<double***>& q                  ,
                            const View<double**>&  qd                 ,
                            const View<double**>&  z                  ,
                            const View<double***>& Vterm              ,
                            bool             precipitating_flag [nqmax],
                            const View<double**>&  precip             ,
                            const View<double*>&   precip_rhoe        ,
                            const View<double*>&   precip_lh_heat     ,
                            const View<double*>&   precip_rhophi      ,
                            const View<double*>&   precip_rhokin      ,
                            const View<double**>&  frain              ,
                            const View<double**>&  gsgam2             ,
                            const View<double**>&  gsgam2h            ,
                            const View<double**>&  rgs                ,
                            const View<double**>&  rgsh               ,
                            const View<double*>&   ix                 ,
                            const View<double*>&   iy                 ,
                            const View<double*>&   iz                 ,
                            const View<double*>&   jx                 ,
                            const View<double*>&   jy                 ,
                            const View<double*>&   jz                 ,
                            double dt,
                            double **precip_trc      // precip[nqmax][ijdim]
                            );
}
