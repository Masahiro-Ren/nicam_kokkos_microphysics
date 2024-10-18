#pragma once

#include "problemsize.h"
#include "mod_debug.h"

using namespace PROBLEM_SIZE;
using namespace DEBUG;

namespace DRIVER{

void mp_init(const std::string& MP_TYPE_in);

template<size_t ijdim, size_t kdim>
void mp_driver( size_t l_region,
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
				double uccn          [kdim][ijdim],
				double rho           [kdim][ijdim],
				double tem           [kdim][ijdim],
				double pre           [kdim][ijdim],
				double q             [nqmax][kdim][ijdim],
				double qd            [kdim][ijdim],
				double precip        [2][ijdim],
				double ISO1_precip   [2][ijdim],
				double ISO2_precip   [2][ijdim],
				double precip_rhoe   [ijdim],
				double precip_lh_heat[ijdim],
				double precip_rhophi [ijdim],
				double precip_rhokin [ijdim],
				double re_liquid	 [kdim][ijdim],
				double re_solid		 [kdim][ijdim],
				double re_cld		 [kdim][ijdim],
				double rctop 		 [1][ijdim],
				double rwtop 		 [1][ijdim],
				double tctop 		 [1][ijdim],
				double frhoge_af     [kdim][ijdim],
				double frhogqv_af    [kdim][ijdim],
				double frhoge_rad    [kdim][ijdim],
				double rhogqke       [kdim][ijdim],
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
				double zh			 [kdim][ijdim],
				double dt,
				double ct,
				double GDCLW		 [kdim][ijdim],
				double GDCFRC		 [kdim][ijdim],
				double GPREC		 [kdim][ijdim],
				double CBMFX		 [kdim][ijdim]
				);

}