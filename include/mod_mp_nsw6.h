#pragma once

#include "problemsize.h"
#include "mod_debug.h"
#include "mod_satadjust.h"


using namespace PROBLEM_SIZE;
using namespace DEBUG;
using namespace SATADJUST;

namespace MP_NSW6{

void mp_nsw6_init();

template <size_t ijdim, size_t kdim>
void mp_nsw6(
            int l_region,
            double rhog          [kdim][ijdim],
            double rhogvx        [kdim][ijdim],
            double rhogvy        [kdim][ijdim],
            double rhogvz        [kdim][ijdim],
            double rhogw         [kdim][ijdim],
            double rhoge         [kdim][ijdim],
            double rhogq         [nqmax][kdim][ijdim],
            double vx            [kdim][ijdim],
            double vy            [kdim][ijdim],
            double vz            [kdim][ijdim],
            double w             [kdim][ijdim],
            double UNCCN         [kdim][ijdim],
            double rho           [kdim][ijdim],
            double tem           [kdim][ijdim],
            double pre           [kdim][ijdim],
            double q             [nqmax][kdim][ijdim],
            double qd            [kdim][ijdim],
            double precip        [2][ijdim],
            double precip_rhoe   [ijdim],
            double precip_lh_heat[ijdim],
            double precip_rhophi [ijdim],
            double precip_rhokin [ijdim],
            double gprec         [kdim][ijdim],
            double rceff         [kdim][ijdim],
            double rctop         [1][ijdim],
            double rwtop         [1][ijdim],
            double tctop         [1][ijdim],
            double rceff_cld     [kdim][ijdim],
            double rctop_cld     [1][ijdim],
            double rwtop_cld     [1][ijdim],
            double tctop_cld     [1][ijdim],
            double gsgam2        [kdim][ijdim],
            double gsgam2h       [kdim][ijdim],
            double gam2          [kdim][ijdim],
            double gam2h         [kdim][ijdim],
            double ix			 [ijdim],
            double iy			 [ijdim],
            double iz			 [ijdim],
            double jx			 [ijdim],
            double jy			 [ijdim],
            double jz			 [ijdim],
            double z			 [kdim][ijdim],
            double dt
            );
}
