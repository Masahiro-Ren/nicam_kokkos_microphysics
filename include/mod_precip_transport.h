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
void precip_transport_new(  const View2D<double, DEFAULT_MEM>&  rhog               ,
                            const View2D<double, DEFAULT_MEM>&  rhogvx             ,
                            const View2D<double, DEFAULT_MEM>&  rhogvy             ,
                            const View2D<double, DEFAULT_MEM>&  rhogvz             ,
                            const View2D<double, DEFAULT_MEM>&  rhogw              ,
                            const View2D<double, DEFAULT_MEM>&  rhoge              ,
                            const View3D<double, DEFAULT_MEM>&  rhogq              ,
                            const View2D<double, DEFAULT_MEM>&  rho                ,
                            const View2D<double, DEFAULT_MEM>&  tem                ,
                            const View2D<double, DEFAULT_MEM>&  pre                ,
                            const View2D<double, DEFAULT_MEM>&  vx                 ,
                            const View2D<double, DEFAULT_MEM>&  vy                 ,
                            const View2D<double, DEFAULT_MEM>&  vz                 ,
                            const View2D<double, DEFAULT_MEM>&  w                  ,
                            const View3D<double, DEFAULT_MEM>&  q                  ,
                            const View2D<double, DEFAULT_MEM>&  qd                 ,
                            const View2D<double, DEFAULT_MEM>&  z                  ,
                            const View3D<double, DEFAULT_MEM>&  Vterm              ,
                            bool             precipitating_flag [nqmax],
                            const View2D<double, DEFAULT_MEM>&  precip             ,
                            const View1D<double, DEFAULT_MEM>&  precip_rhoe        ,
                            const View1D<double, DEFAULT_MEM>&  precip_lh_heat     ,
                            const View1D<double, DEFAULT_MEM>&  precip_rhophi      ,
                            const View1D<double, DEFAULT_MEM>&  precip_rhokin      ,
                            const View2D<double, DEFAULT_MEM>&  frain              ,
                            const View2D<double, DEFAULT_MEM>&  gsgam2             ,
                            const View2D<double, DEFAULT_MEM>&  gsgam2h            ,
                            const View2D<double, DEFAULT_MEM>&  rgs                ,
                            const View2D<double, DEFAULT_MEM>&  rgsh               ,
                            const View1D<double, DEFAULT_MEM>&  ix                 ,
                            const View1D<double, DEFAULT_MEM>&  iy                 ,
                            const View1D<double, DEFAULT_MEM>&  iz                 ,
                            const View1D<double, DEFAULT_MEM>&  jx                 ,
                            const View1D<double, DEFAULT_MEM>&  jy                 ,
                            const View1D<double, DEFAULT_MEM>&  jz                 ,
                            double h_dt,
                            double **precip_trc      // precip[nqmax][ijdim]
                            );
}
