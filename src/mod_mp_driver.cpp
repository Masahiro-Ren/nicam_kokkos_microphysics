#include "mod_mp_driver.h"

bool opt_radius_explicit = false;
bool opt_volume_explicit = false;

double TSICE = 273.15;
double TWICE = 258.15;

namespace MP_DRIVER{

void mp_init(const std::string& MP_TYPE_in)
{
	if(MP_TYPE_in == "NSW6")
	{
		std::cout << "*** microphysics type: NSW6 *** \n";
		// std::cout << "Call mp_nsw6_init \n";
		mp_nsw6_init();
	}
	else
	{
		std::cerr << __PRETTY_FUNCTION__ << " NOT appropriate type. Type:  " << MP_TYPE << std::endl;
		ADM_Proc_stop();
	}
}

void mp_driver( int l_region,
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
				double unccn         [kdim][ijdim],
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
				)
{
	std::cout << __PRETTY_FUNCTION__ << std::endl;

	double precip_trc        [nqmax][ijdim];

	double precip_sum		 [2][ijdim];
	double ISO1_precip_sum   [2][ijdim];
	double ISO2_precip_sum   [2][ijdim];
	double precip_rhoe_sum   [ijdim];
	double precip_lh_heat_sum[ijdim];
	double precip_rhophi_sum [ijdim];
	double precip_rhokin_sum [ijdim];
	double precip_trc_sum	 [nqmax][ijdim];
	double GDCLW_sum		 [kdim][ijdim];
	double GDCFRC_sum		 [kdim][ijdim];
	double GPREC_sum		 [kdim][ijdim];

	double re_rain           [kdim][ijdim]; // Effective Radius
	double re_ice            [kdim][ijdim]; // Effective Radius
	double re_snow           [kdim][ijdim]; // Effective Radius
	double re_graupel        [kdim][ijdim]; // Effective Radius
	double rctop_cld         [1][ijdim];    // Effective Radius of Cloud Top
	double rwtop_cld         [1][ijdim];    // Effective Radius of Warm-Cloud Top
	double tctop_cld         [1][ijdim];    // Cloud Top Temperature

	double qke				 [kdim][ijdim];

	double fraction_mp = 1.0 / double(MP_DIV_NUM); // 1 / MP_DIV_NUM
	double dt_mp = dt * fraction_mp;

	for(int k = 0; k < 2; k++)
	{
		for(int ij = 0; ij < ijdim; ij++)
		{
			precip_sum[k][ij] = 0.0;
			ISO1_precip_sum[k][ij] = 0.0;
			ISO2_precip_sum[k][ij] = 0.0;
		}
	}
	for(int nq = 0; nq < nqmax; nq++)
	{
		for(int ij = 0; ij < ijdim; ij++)
		{
			precip_trc_sum[nq][ij] = 0.0;
		}
	}
	for(int k = 0; k < kdim; k++)
	{
		for(int ij = 0; ij < ijdim; ij++)
		{
			GDCLW_sum[k][ij] = 0.0;
			GDCFRC_sum[k][ij] = 0.0;
			GPREC_sum[k][ij] = 0.0;
		}
	}
	for(int ij = 0; ij < ijdim; ij++)
	{
		precip_rhoe_sum[ij] = 0.0;
		precip_lh_heat_sum[ij] = 0.0;
		precip_rhophi_sum [ij] = 0.0;
		precip_rhokin_sum [ij] = 0.0;
	}

	// memset(precip_sum,         0, sizeof(precip_sum));
	// memset(ISO1_precip_sum,    0, sizeof(ISO1_precip_sum));
	// memset(ISO2_precip_sum,    0, sizeof(ISO2_precip_sum));
	// memset(precip_rhoe_sum,    0, sizeof(precip_rhoe_sum));
	// memset(precip_lh_heat_sum, 0, sizeof(precip_lh_heat_sum));
	// memset(precip_rhophi_sum,  0, sizeof(precip_rhophi_sum));
	// memset(precip_rhokin_sum,  0, sizeof(precip_rhokin_sum));
	// memset(precip_trc_sum,     0, sizeof(precip_trc_sum));
	// memset(GDCLW_sum,          0, sizeof(GDCLW_sum));
	// memset(GDCFRC_sum,         0, sizeof(GDCFRC_sum));
	// memset(GPREC_sum,          0, sizeof(GPREC_sum));

	if(MP_TYPE != "NDW6")
	{
		for(int k = 0; k < kdim; k++)
		{
			for(int ij = 0; ij < ijdim; ij++)
			{
				re_rain[k][ij]    = CONST_UNDEF;
				re_ice[k][ij]     = CONST_UNDEF;
				re_snow[k][ij]    = CONST_UNDEF;
				re_graupel[k][ij] = CONST_UNDEF;
			}
		}
	}


	for(int m = 1; m < MP_DIV_NUM; m++)
	{
		for(int k = 0; k < 2; k++)
		{
			for(int ij = 0; ij < ijdim; ij++)
			{
				precip_sum[k][ij] = 0.0;
				ISO1_precip_sum[k][ij] = 0.0;
				ISO2_precip_sum[k][ij] = 0.0;
			}
		}
		for(int nq = 0; nq < nqmax; nq++)
		{
			for(int ij = 0; ij < ijdim; ij++)
			{
				precip_trc_sum[nq][ij] = 0.0;
			}
		}
		for(int k = 0; k < kdim; k++)
		{
			for(int ij = 0; ij < ijdim; ij++)
			{
				GDCLW_sum[k][ij] = 0.0;
				GDCFRC_sum[k][ij] = 0.0;
				GPREC_sum[k][ij] = 0.0;
			}
		}
		for(int ij = 0; ij < ijdim; ij++)
		{
			precip_rhoe_sum[ij] = 0.0;
			precip_lh_heat_sum[ij] = 0.0;
			precip_rhophi_sum [ij] = 0.0;
			precip_rhokin_sum [ij] = 0.0;
		}

		if(MP_TYPE == "NONE")
		{
			for(int k = 0; k < kdim; k++)
			{
				for(int ij = 0; ij < ijdim; ij++)
				{
					re_liquid[k][ij] = 10.0E-6;
					re_solid[k][ij] = 20.0E-6;
				}
			}
		}
		else if(MP_TYPE == "NSW6")
		{
			mp_nsw6(l_region,
					rhog,
					rhogvx,
					rhogvy,
					rhogvz,
					rhogw,
					rhoge,
					rhogq,
					vx,
					vy,
					vz,
					w,
					unccn,
					rho,
					tem,
					pre,
					q,
					qd,
					precip,
					precip_rhoe,
					precip_lh_heat,
					precip_rhophi,
					precip_rhokin,
					GPREC,
					re_liquid,
					rctop,
					rwtop,
					tctop,
					re_cld,
					rctop_cld,
					rwtop_cld,
					tctop_cld,
					gsgam2,
					gsgam2h,
					gam2,
					gam2h,
					ix,
					iy,
					iz,
					jx,
					jy,
					jz,
					z,
					dt_mp);

			for(int k = 0; k < kdim; k++)
				for(int ij = 0; ij < ijdim; ij++)
					re_solid[k][ij] = 20.0E-6;
		}

		for(int k = 0; k < 2; k++)
		{
			for(int ij = 0; ij < ijdim; ij++)
			{
				precip_sum[k][ij]      += precip[k][ij];
				ISO1_precip_sum[k][ij] += ISO1_precip[k][ij];
				ISO2_precip_sum[k][ij] += ISO2_precip[k][ij];
			}
		}
		for(int ij = 0; ij < ijdim; ij++)
		{
			precip_rhoe_sum[ij]    += precip_rhoe[ij];
			precip_lh_heat_sum[ij] += precip_lh_heat[ij];
			precip_rhophi_sum [ij] += precip_rhophi[ij];
			precip_rhokin_sum [ij] += precip_rhokin[ij];
		}

		for(int nq = 0; nq < nqmax; nq++)
			for(int ij = 0; ij < ijdim; ij++)
				precip_trc_sum[nq][ij] += precip_trc[nq][ij];

		for(int k = 0; k < kdim; k++)
		{
			for(int ij = 0; ij < ijdim; ij++)
			{
				GDCLW_sum [k][ij] += GDCLW[k][ij];
				GDCFRC_sum[k][ij] += GDCFRC[k][ij];
				GPREC_sum [k][ij] += GPREC[k][ij];
			}
		}
	}

	for(int k = 0; k < 2; k++)
	{
		for(int ij = 0; ij < ijdim; ij++)
		{
			precip[k][ij]      = precip_sum[k][ij] * fraction_mp;
			ISO1_precip[k][ij] = ISO1_precip_sum[k][ij] * fraction_mp;
			ISO2_precip[k][ij] = ISO2_precip_sum[k][ij] * fraction_mp;
		}
	}
	for(int ij = 0; ij < ijdim; ij++)
	{
		precip_rhoe[ij]    = precip_rhoe_sum[ij] * fraction_mp;
		precip_lh_heat[ij] = precip_lh_heat_sum[ij] * fraction_mp;
		precip_rhophi [ij] = precip_rhophi_sum[ij] * fraction_mp;
		precip_rhokin [ij] = precip_rhokin_sum[ij] * fraction_mp;
	}
	for(int nq = 0; nq < nqmax; nq++)
		for(int ij = 0; ij < ijdim; ij++)
			precip_trc[nq][ij] += precip_trc_sum[nq][ij] * fraction_mp;

	for(int k = 0; k < kdim; k++)
	{
		for(int ij = 0; ij < ijdim; ij++)
		{
			GDCLW [k][ij] = GDCLW[k][ij] * fraction_mp;
			GDCFRC[k][ij] = GDCFRC[k][ij] * fraction_mp;
			GPREC [k][ij] = GPREC[k][ij] * fraction_mp;
		}
	}

}

}