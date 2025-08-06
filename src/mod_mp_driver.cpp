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

// void mp_driver( int l_region,
// 				double rhog          [kdim][ijdim],
// 				double rhogvx        [kdim][ijdim],
// 				double rhogvy        [kdim][ijdim],
// 				double rhogvz        [kdim][ijdim],
// 				double rhogw         [kdim][ijdim],
// 				double rhoge         [kdim][ijdim],
// 				double rhogq         [nqmax][kdim][ijdim],
// 				double vx            [kdim][ijdim],
// 				double vy            [kdim][ijdim],
// 				double vz            [kdim][ijdim],
// 				double w             [kdim][ijdim],
// 				double unccn         [kdim][ijdim],
// 				double rho           [kdim][ijdim],
// 				double tem           [kdim][ijdim],
// 				double pre           [kdim][ijdim],
// 				double q             [nqmax][kdim][ijdim],
// 				double qd            [kdim][ijdim],
// 				double precip        [2][ijdim],
// 				double ISO1_precip   [2][ijdim],
// 				double ISO2_precip   [2][ijdim],
// 				double precip_rhoe   [ijdim],
// 				double precip_lh_heat[ijdim],
// 				double precip_rhophi [ijdim],
// 				double precip_rhokin [ijdim],
// 				double re_liquid	 [kdim][ijdim],
// 				double re_solid		 [kdim][ijdim],
// 				double re_cld		 [kdim][ijdim],
// 				double rctop 		 [1][ijdim],
// 				double rwtop 		 [1][ijdim],
// 				double tctop 		 [1][ijdim],
// 				double frhoge_af     [kdim][ijdim],
// 				double frhogqv_af    [kdim][ijdim],
// 				double frhoge_rad    [kdim][ijdim],
// 				double rhogqke       [kdim][ijdim],
// 				double gsgam2        [kdim][ijdim],
// 				double gsgam2h       [kdim][ijdim],
// 				double gam2          [kdim][ijdim],
// 				double gam2h         [kdim][ijdim],
// 				double ix			 [ijdim],
// 				double iy			 [ijdim],
// 				double iz			 [ijdim],
// 				double jx			 [ijdim],
// 				double jy			 [ijdim],
// 				double jz			 [ijdim],
// 				double z			 [kdim][ijdim],
// 				double zh			 [kdim][ijdim],
// 				double dt,
// 				double ct,
// 				double GDCLW		 [kdim][ijdim],
// 				double GDCFRC		 [kdim][ijdim],
// 				double GPREC		 [kdim][ijdim],
// 				double CBMFX		 [kdim][ijdim]
// 				)
// {
// #ifdef ENABLE_DEBUG
// 	std::cout << __PRETTY_FUNCTION__ << std::endl;
// #endif

// 	double precip_trc        [nqmax][ijdim];

// 	double precip_sum		 [2][ijdim];
// 	double ISO1_precip_sum   [2][ijdim];
// 	double ISO2_precip_sum   [2][ijdim];
// 	double precip_rhoe_sum   [ijdim];
// 	double precip_lh_heat_sum[ijdim];
// 	double precip_rhophi_sum [ijdim];
// 	double precip_rhokin_sum [ijdim];
// 	double precip_trc_sum	 [nqmax][ijdim];
// 	double GDCLW_sum		 [kdim][ijdim];
// 	double GDCFRC_sum		 [kdim][ijdim];
// 	double GPREC_sum		 [kdim][ijdim];

// 	double re_rain           [kdim][ijdim]; // Effective Radius
// 	double re_ice            [kdim][ijdim]; // Effective Radius
// 	double re_snow           [kdim][ijdim]; // Effective Radius
// 	double re_graupel        [kdim][ijdim]; // Effective Radius
// 	double rctop_cld         [1][ijdim];    // Effective Radius of Cloud Top
// 	double rwtop_cld         [1][ijdim];    // Effective Radius of Warm-Cloud Top
// 	double tctop_cld         [1][ijdim];    // Cloud Top Temperature

// 	double qke				 [kdim][ijdim];

// 	double fraction_mp = 1.0 / double(MP_DIV_NUM); // 1 / MP_DIV_NUM
// 	double dt_mp = dt * fraction_mp;

// 	for(int k = 0; k < 2; k++)
// 	{
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			precip_sum[k][ij] = 0.0;
// 			ISO1_precip_sum[k][ij] = 0.0;
// 			ISO2_precip_sum[k][ij] = 0.0;
// 		}
// 	}
// 	for(int nq = 0; nq < nqmax; nq++)
// 	{
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			precip_trc_sum[nq][ij] = 0.0;
// 		}
// 	}
// 	for(int k = 0; k < kdim; k++)
// 	{
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			GDCLW_sum[k][ij] = 0.0;
// 			GDCFRC_sum[k][ij] = 0.0;
// 			GPREC_sum[k][ij] = 0.0;
// 		}
// 	}
// 	for(int ij = 0; ij < ijdim; ij++)
// 	{
// 		precip_rhoe_sum[ij] = 0.0;
// 		precip_lh_heat_sum[ij] = 0.0;
// 		precip_rhophi_sum [ij] = 0.0;
// 		precip_rhokin_sum [ij] = 0.0;
// 	}

// 	// memset(precip_sum,         0, sizeof(precip_sum));
// 	// memset(ISO1_precip_sum,    0, sizeof(ISO1_precip_sum));
// 	// memset(ISO2_precip_sum,    0, sizeof(ISO2_precip_sum));
// 	// memset(precip_rhoe_sum,    0, sizeof(precip_rhoe_sum));
// 	// memset(precip_lh_heat_sum, 0, sizeof(precip_lh_heat_sum));
// 	// memset(precip_rhophi_sum,  0, sizeof(precip_rhophi_sum));
// 	// memset(precip_rhokin_sum,  0, sizeof(precip_rhokin_sum));
// 	// memset(precip_trc_sum,     0, sizeof(precip_trc_sum));
// 	// memset(GDCLW_sum,          0, sizeof(GDCLW_sum));
// 	// memset(GDCFRC_sum,         0, sizeof(GDCFRC_sum));
// 	// memset(GPREC_sum,          0, sizeof(GPREC_sum));

// 	if(MP_TYPE != "NDW6")
// 	{
// 		for(int k = 0; k < kdim; k++)
// 		{
// 			for(int ij = 0; ij < ijdim; ij++)
// 			{
// 				re_rain[k][ij]    = CONST_UNDEF;
// 				re_ice[k][ij]     = CONST_UNDEF;
// 				re_snow[k][ij]    = CONST_UNDEF;
// 				re_graupel[k][ij] = CONST_UNDEF;
// 			}
// 		}
// 	}


// 	for(int m = 0; m < MP_DIV_NUM; m++)
// 	{
// 		for(int k = 0; k < 2; k++)
// 		{
// 			for(int ij = 0; ij < ijdim; ij++)
// 			{
// 				precip     [k][ij] = 0.0;
// 				ISO1_precip[k][ij] = 0.0;
// 				ISO2_precip[k][ij] = 0.0;
// 			}
// 		}
// 		for(int nq = 0; nq < nqmax; nq++)
// 		{
// 			for(int ij = 0; ij < ijdim; ij++)
// 			{
// 				precip_trc[nq][ij] = 0.0;
// 			}
// 		}
// 		for(int k = 0; k < kdim; k++)
// 		{
// 			for(int ij = 0; ij < ijdim; ij++)
// 			{
// 				GDCLW [k][ij] = 0.0;
// 				GDCFRC[k][ij] = 0.0;
// 				GPREC [k][ij] = 0.0;
// 			}
// 		}
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			precip_rhoe   [ij] = 0.0;
// 			precip_lh_heat[ij] = 0.0;
// 			precip_rhophi [ij] = 0.0;
// 			precip_rhokin [ij] = 0.0;
// 		}

// 		if(MP_TYPE == "NONE")
// 		{
// 			for(int k = 0; k < kdim; k++)
// 			{
// 				for(int ij = 0; ij < ijdim; ij++)
// 				{
// 					re_liquid[k][ij] = 10.0E-6;
// 					re_solid[k][ij] = 20.0E-6;
// 				}
// 			}
// 		}
// 		else if(MP_TYPE == "NSW6")
// 		{
// 			mp_nsw6(l_region,
// 					rhog,
// 					rhogvx,
// 					rhogvy,
// 					rhogvz,
// 					rhogw,
// 					rhoge,
// 					rhogq,
// 					vx,
// 					vy,
// 					vz,
// 					w,
// 					unccn,
// 					rho,
// 					tem,
// 					pre,
// 					q,
// 					qd,
// 					precip,
// 					precip_rhoe,
// 					precip_lh_heat,
// 					precip_rhophi,
// 					precip_rhokin,
// 					GPREC,
// 					re_liquid,
// 					rctop,
// 					rwtop,
// 					tctop,
// 					re_cld,
// 					rctop_cld,
// 					rwtop_cld,
// 					tctop_cld,
// 					gsgam2,
// 					gsgam2h,
// 					gam2,
// 					gam2h,
// 					ix,
// 					iy,
// 					iz,
// 					jx,
// 					jy,
// 					jz,
// 					z,
// 					dt_mp);

// 			for(int k = 0; k < kdim; k++)
// 				for(int ij = 0; ij < ijdim; ij++)
// 					re_solid[k][ij] = 20.0E-6;
// 		}

// 		for(int k = 0; k < 2; k++)
// 		{
// 			for(int ij = 0; ij < ijdim; ij++)
// 			{
// 				precip_sum[k][ij]      += precip[k][ij];
// 				ISO1_precip_sum[k][ij] += ISO1_precip[k][ij];
// 				ISO2_precip_sum[k][ij] += ISO2_precip[k][ij];
// 			}
// 		}
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			precip_rhoe_sum[ij]    += precip_rhoe[ij];
// 			precip_lh_heat_sum[ij] += precip_lh_heat[ij];
// 			precip_rhophi_sum [ij] += precip_rhophi[ij];
// 			precip_rhokin_sum [ij] += precip_rhokin[ij];
// 		}

// 		for(int nq = 0; nq < nqmax; nq++)
// 			for(int ij = 0; ij < ijdim; ij++)
// 				precip_trc_sum[nq][ij] += precip_trc[nq][ij];

// 		for(int k = 0; k < kdim; k++)
// 		{
// 			for(int ij = 0; ij < ijdim; ij++)
// 			{
// 				GDCLW_sum [k][ij] += GDCLW[k][ij];
// 				GDCFRC_sum[k][ij] += GDCFRC[k][ij];
// 				GPREC_sum [k][ij] += GPREC[k][ij];
// 			}
// 		}
// 	}

// 	for(int k = 0; k < 2; k++)
// 	{
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			precip[k][ij]      = precip_sum[k][ij] * fraction_mp;
// 			ISO1_precip[k][ij] = ISO1_precip_sum[k][ij] * fraction_mp;
// 			ISO2_precip[k][ij] = ISO2_precip_sum[k][ij] * fraction_mp;
// 		}
// 	}
// 	for(int ij = 0; ij < ijdim; ij++)
// 	{
// 		precip_rhoe[ij]    = precip_rhoe_sum[ij] * fraction_mp;
// 		precip_lh_heat[ij] = precip_lh_heat_sum[ij] * fraction_mp;
// 		precip_rhophi [ij] = precip_rhophi_sum[ij] * fraction_mp;
// 		precip_rhokin [ij] = precip_rhokin_sum[ij] * fraction_mp;
// 	}
// 	for(int nq = 0; nq < nqmax; nq++)
// 		for(int ij = 0; ij < ijdim; ij++)
// 			precip_trc[nq][ij] += precip_trc_sum[nq][ij] * fraction_mp;

// 	for(int k = 0; k < kdim; k++)
// 	{
// 		for(int ij = 0; ij < ijdim; ij++)
// 		{
// 			GDCLW [k][ij] = GDCLW_sum[k][ij] * fraction_mp;
// 			GDCFRC[k][ij] = GDCFRC_sum[k][ij] * fraction_mp;
// 			GPREC [k][ij] = GPREC_sum[k][ij] * fraction_mp;
// 		}
// 	}

// }

// Kokkos ver.
void mp_driver( int l_region,
				const View2D<double, DEFAULT_MEM>&  rhog          ,
				const View2D<double, DEFAULT_MEM>&  rhogvx        ,
				const View2D<double, DEFAULT_MEM>&  rhogvy        ,
				const View2D<double, DEFAULT_MEM>&  rhogvz        ,
				const View2D<double, DEFAULT_MEM>&  rhogw         ,
				const View2D<double, DEFAULT_MEM>&  rhoge         ,
				const View3D<double, DEFAULT_MEM>&  rhogq         ,
				const View2D<double, DEFAULT_MEM>&  vx            ,
				const View2D<double, DEFAULT_MEM>&  vy            ,
				const View2D<double, DEFAULT_MEM>&  vz            ,
				const View2D<double, DEFAULT_MEM>&  w             ,
				const View2D<double, DEFAULT_MEM>&  unccn         ,
				const View2D<double, DEFAULT_MEM>&  rho           ,
				const View2D<double, DEFAULT_MEM>&  tem           ,
				const View2D<double, DEFAULT_MEM>&  pre           ,
				const View3D<double, DEFAULT_MEM>&  q             ,
				const View2D<double, DEFAULT_MEM>&  qd            ,
				const View2D<double, DEFAULT_MEM>&  precip        ,
				const View2D<double, DEFAULT_MEM>&  ISO1_precip   ,
				const View2D<double, DEFAULT_MEM>&  ISO2_precip   ,
				const View1D<double, DEFAULT_MEM>&  precip_rhoe   ,
				const View1D<double, DEFAULT_MEM>&  precip_lh_heat,
				const View1D<double, DEFAULT_MEM>&  precip_rhophi ,
				const View1D<double, DEFAULT_MEM>&  precip_rhokin ,
				const View2D<double, DEFAULT_MEM>&  re_liquid	 ,
				const View2D<double, DEFAULT_MEM>&  re_solid		 ,
				const View2D<double, DEFAULT_MEM>&  re_cld		 ,
				const View2D<double, DEFAULT_MEM>&  rctop 		 ,
				const View2D<double, DEFAULT_MEM>&  rwtop 		 ,
				const View2D<double, DEFAULT_MEM>&  tctop 		 ,
				const View2D<double, DEFAULT_MEM>&  frhoge_af     ,
				const View2D<double, DEFAULT_MEM>&  frhogqv_af    ,
				const View2D<double, DEFAULT_MEM>&  frhoge_rad    ,
				const View2D<double, DEFAULT_MEM>&  rhogqke       ,
				const View2D<double, DEFAULT_MEM>&  gsgam2        ,
				const View2D<double, DEFAULT_MEM>&  gsgam2h       ,
				const View2D<double, DEFAULT_MEM>&  gam2          ,
				const View2D<double, DEFAULT_MEM>&  gam2h         ,
				const View1D<double, DEFAULT_MEM>&  ix			 ,
				const View1D<double, DEFAULT_MEM>&  iy			 ,
				const View1D<double, DEFAULT_MEM>&  iz			 ,
				const View1D<double, DEFAULT_MEM>&  jx			 ,
				const View1D<double, DEFAULT_MEM>&  jy			 ,
				const View1D<double, DEFAULT_MEM>&  jz			 ,
				const View2D<double, DEFAULT_MEM>&  z			 ,
				const View2D<double, DEFAULT_MEM>&  zh			 ,
				double 				                dt			 ,
				double 				                ct			 ,
				const View2D<double, DEFAULT_MEM>&  GDCLW		 ,
				const View2D<double, DEFAULT_MEM>&  GDCFRC		 ,
				const View2D<double, DEFAULT_MEM>&  GPREC		 ,
				const View2D<double, DEFAULT_MEM>&  CBMFX		 
				)
{
#ifdef ENABLE_DEBUG
	std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

	View2D<double, DEFAULT_MEM> precip_trc("precip_trc", nqmax, ijdim);

	View2D<double, DEFAULT_MEM> precip_sum		 ("precip_sum", 2, ijdim);
	View2D<double, DEFAULT_MEM> ISO1_precip_sum   ("ISO1_precip_sum", 2, ijdim);
	View2D<double, DEFAULT_MEM> ISO2_precip_sum   ("ISO2_precip_sum", 2, ijdim);
	View1D<double, DEFAULT_MEM> precip_rhoe_sum   ("precip_rhoe_sum   ", ijdim);
	View1D<double, DEFAULT_MEM> precip_lh_heat_sum("precip_lh_heat_sum", ijdim);
	View1D<double, DEFAULT_MEM> precip_rhophi_sum ("precip_rhophi_sum ", ijdim);
	View1D<double, DEFAULT_MEM> precip_rhokin_sum ("precip_rhokin_sum ", ijdim);
	View2D<double, DEFAULT_MEM> precip_trc_sum	 ("precip_trc_sum", nqmax, ijdim);
	View2D<double, DEFAULT_MEM> GDCLW_sum		 ("GDCLW_sum", kdim, ijdim);
	View2D<double, DEFAULT_MEM> GDCFRC_sum		 ("GDCFRC_sum", kdim, ijdim);
	View2D<double, DEFAULT_MEM> GPREC_sum		 ("GPREC_sum", kdim, ijdim);

	View2D<double, DEFAULT_MEM> re_rain           ("re_rain   ", kdim, ijdim); // Effective Radius
	View2D<double, DEFAULT_MEM> re_ice            ("re_ice    ", kdim, ijdim); // Effective Radius
	View2D<double, DEFAULT_MEM> re_snow           ("re_snow   ", kdim, ijdim); // Effective Radius
	View2D<double, DEFAULT_MEM> re_graupel        ("re_graupel", kdim, ijdim); // Effective Radius
	View2D<double, DEFAULT_MEM> rctop_cld         ("rctop_cld", 1, ijdim);    // Effective Radius of Cloud Top
	View2D<double, DEFAULT_MEM> rwtop_cld         ("rwtop_cld", 1, ijdim);    // Effective Radius of Warm-Cloud Top
	View2D<double, DEFAULT_MEM> tctop_cld         ("tctop_cld", 1, ijdim);    // Cloud Top Temperature

	View2D<double, DEFAULT_MEM> qke				 ("qke", kdim, ijdim);

	View<double> fraction_mp("fraction_mp");

	double h_fraction_mp = 1.0 / double(MP_DIV_NUM); // 1 / MP_DIV_NUM
	double dt_mp = dt * h_fraction_mp;

	Kokkos::deep_copy(fraction_mp, h_fraction_mp);

	Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{2,ijdim}), 
	KOKKOS_LAMBDA(const size_t k, const size_t ij){
		precip_sum     (k,ij) = 0.0;
		ISO1_precip_sum(k,ij) = 0.0;
		ISO2_precip_sum(k,ij) = 0.0;	
	});

	Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), 
	KOKKOS_LAMBDA(const size_t nq, const size_t ij){
		precip_trc_sum(nq,ij) = 0.0;
	});

	Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
	KOKKOS_LAMBDA(const size_t k, const size_t ij){
		GDCLW_sum (k,ij) = 0.0;
		GDCFRC_sum(k,ij) = 0.0;
		GPREC_sum (k,ij) = 0.0;
	});

	Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
		precip_rhoe_sum   (ij) = 0.0;
		precip_lh_heat_sum(ij) = 0.0;
		precip_rhophi_sum (ij) = 0.0;
		precip_rhokin_sum (ij) = 0.0;
	});

	if(MP_TYPE != "NDW6")
	{
		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
		KOKKOS_LAMBDA(const size_t k, const size_t ij){
			re_rain(k,ij)    = CONST_UNDEF;
			re_ice(k,ij)     = CONST_UNDEF;
			re_snow(k,ij)    = CONST_UNDEF;
			re_graupel(k,ij) = CONST_UNDEF;	
		});
	}

	for(int m = 0; m < MP_DIV_NUM; m++)
	{
		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{2,ijdim}), 
		KOKKOS_LAMBDA(const size_t k, const size_t ij){
			precip     (k,ij) = 0.0;
			ISO1_precip(k,ij) = 0.0;
			ISO2_precip(k,ij) = 0.0;
		});

		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), 
		KOKKOS_LAMBDA(const size_t nq, const size_t ij){
			precip_trc(nq,ij) = 0.0;
		});

		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
		KOKKOS_LAMBDA(const size_t k, const size_t ij){
			GDCLW (k,ij) = 0.0;
			GDCFRC(k,ij) = 0.0;
			GPREC (k,ij) = 0.0;
		});

		Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
			precip_rhoe   (ij) = 0.0;
			precip_lh_heat(ij) = 0.0;
			precip_rhophi (ij) = 0.0;
			precip_rhokin (ij) = 0.0;
		});

		if(MP_TYPE == "NONE")
		{
			Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
			KOKKOS_LAMBDA(const size_t k, const size_t ij){
				re_liquid(k,ij) = 10.0E-6;
				re_solid (k,ij) = 20.0E-6;
			});
		}
		else if(MP_TYPE == "NSW6")
		{
			Kokkos::fence();
			Kokkos::Timer timer;
			mp_nsw6_new mp_nsw6(l_region,
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
					dt_mp,
				    para);
				
			mp_nsw6.compute();

			Kokkos::fence();
			double elapsed = timer.seconds();
			std::cout << "NSW6 Elapsed Time: " << elapsed << std::endl;
			
			
			Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
			KOKKOS_LAMBDA(const size_t k, const size_t ij){
				re_solid (k,ij) = 20.0E-6;
			});
		}

		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{2,ijdim}), 
		KOKKOS_LAMBDA(const size_t k, const size_t ij){
			precip_sum(k,ij)      += precip(k,ij);
			ISO1_precip_sum(k,ij) += ISO1_precip(k,ij);
			ISO2_precip_sum(k,ij) += ISO2_precip(k,ij);
		});

		Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
			precip_rhoe_sum(ij)    += precip_rhoe(ij);
			precip_lh_heat_sum(ij) += precip_lh_heat(ij);
			precip_rhophi_sum (ij) += precip_rhophi(ij);
			precip_rhokin_sum (ij) += precip_rhokin(ij);
		});

		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), 
		KOKKOS_LAMBDA(const size_t nq, const size_t ij){
			precip_trc_sum(nq,ij) += precip_trc(nq,ij);
		});

		Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
		KOKKOS_LAMBDA(const size_t k, const size_t ij){
			GDCLW_sum (k,ij) += GDCLW (k,ij);
			GDCFRC_sum(k,ij) += GDCFRC(k,ij);
			GPREC_sum (k,ij) += GPREC (k,ij);	
		});
	}

	Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{2,ijdim}), 
	KOKKOS_LAMBDA(const size_t k, const size_t ij){
		precip(k,ij)      = precip_sum(k,ij) * fraction_mp();
		ISO1_precip(k,ij) = ISO1_precip_sum(k,ij) * fraction_mp();
		ISO2_precip(k,ij) = ISO2_precip_sum(k,ij) * fraction_mp();
	});

	Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
		precip_rhoe(ij)    = precip_rhoe_sum(ij) * fraction_mp();
		precip_lh_heat(ij) = precip_lh_heat_sum(ij) * fraction_mp();
		precip_rhophi (ij) = precip_rhophi_sum (ij) * fraction_mp();
		precip_rhokin (ij) = precip_rhokin_sum (ij) * fraction_mp();
	});

	Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), 
	KOKKOS_LAMBDA(const size_t nq, const size_t ij){
		precip_trc(nq,ij) += precip_trc_sum(nq,ij) * fraction_mp();
	});

	Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
	KOKKOS_LAMBDA(const size_t k, const size_t ij){
		GDCLW (k,ij) = GDCLW_sum (k,ij) * fraction_mp();
		GDCFRC(k,ij) = GDCFRC_sum(k,ij) * fraction_mp();
		GPREC (k,ij) = GPREC_sum (k,ij) * fraction_mp();
	});
}

}
