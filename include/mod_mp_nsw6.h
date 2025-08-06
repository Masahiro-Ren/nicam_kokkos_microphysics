#pragma once

#include "problemsize.h"
#include "mod_debug.h"
#include "mod_satadjust.h"
#include "mod_thrmdyn.h"
#include "mod_precip_transport.h"


using namespace PROBLEM_SIZE;
using namespace DEBUG;
using namespace SATADJUST;
using namespace THRMDYN;
using namespace PRECIP_TRANSPORT;

struct GlobalParams {
double N0r    = 8.E+6;  //< intercept parameter for rain    [1/m4]
double N0s    = 3.E+6;  //< intercept parameter for snow    [1/m4]
double N0g    = 4.E+6;  //< intercept parameter for graupel [1/m4]

double rho_s  = 100.0;  //< density of snow    [kg/m3]
double rho_g  = 400.0;  //< density of graupel [kg/m3]
                        //   graupel : 400
                        //   hail    : 917

double C_d    =   0.6;  //< drag coefficient for graupel
double Cr     = 130.0;
double Cs     =  4.84;

// Empirical parameter
double Ar, As, Ag;
double Br, Bs, Bg;
double Cg        ;
double Dr, Ds, Dg;

// GAMMA function
double GAM, GAM_2, GAM_3;

double GAM_1br, GAM_2br, GAM_3br;
double GAM_3dr;
double GAM_6dr;
double GAM_1brdr;
double GAM_5dr_h;

double GAM_1bs, GAM_2bs, GAM_3bs;
double GAM_3ds;
double GAM_1bsds;
double GAM_5ds_h;

double GAM_1bg, GAM_3dg;
double GAM_1bgdg;
double GAM_5dg_h;

//---< Khairoutdinov and Kogan (2000) >---
double sw_kk2000  = 0.0; //< switch for k-k scheme

//---< Roh and Satoh (2014) >---
double sw_roh2014 = 0.0; //< switch for Roh scheme
double ln10;             //< log(10)
double coef_a[10] = {5.065339, -0.062659, -3.032362, 0.029469, -0.000285, 0.31255, 0.000204, 0.003199, 0.0, -0.015952};
double coef_b[10] = {0.476221, -0.015896,  0.165977, 0.007468, -0.000141, 0.060366, 0.000079, 0.000594, 0.0, -0.003577};

// Accretion parameter
double Eiw        = 1.0;      //< collection efficiency of cloud ice for cloud water
double Erw        = 1.0;      //< collection efficiency of rain    for cloud water
double Esw        = 1.0;      //< collection efficiency of snow    for cloud water
double Egw        = 1.0;      //< collection efficiency of graupel for cloud water
double Eri        = 1.0;      //< collection efficiency of rain    for cloud ice
double Esi        = 1.0;      //< collection efficiency of snow    for cloud ice
double Egi        = 0.1;      //< collection efficiency of graupel for cloud ice
double Esr        = 1.0;      //< collection efficiency of snow    for rain
double Egr        = 1.0;      //< collection efficiency of graupel for rain
double Egs        = 1.0;      //< collection efficiency of graupel for snow
double gamma_sacr = 25.E-3;   //< effect of low temperature for Esi
double gamma_gacs = 90.E-3;   //< effect of low temperature for Egs
double mi         = 4.19E-13; //< mass of one cloud ice crystal [kg]

// Auto-conversion parameter
// static const double Nc_lnd     = 2000.0;   //< number concentration of cloud water (land)  [1/cc]
// static const double Nc_ocn     =   50.0;   //< number concentration of cloud water (ocean) [1/cc]
double Nc_def;                             //< number concentration of cloud water         [1/cc]

double beta_saut  =  1.E-3;   //< auto-conversion factor beta  for ice
double gamma_saut = 25.E-3;   //< auto-conversion factor gamma for ice
double beta_gaut  =  1.E-3;   //< auto-conversion factor beta  for snow
double gamma_gaut = 90.E-3;   //< auto-conversion factor gamma for snow
double qicrt_saut =  0.0;     //< mixing ratio threshold for Psaut [kg/kg]
double qscrt_gaut =  6.E-4;   //< mixing ratio threshold for Pgaut [kg/kg]

// // Evaporation, Sublimation parameter
// static const double Ka0        = 2.428E-2; //< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
// static const double dKa_dT     =  7.47E-5; //< Coefficient of Ka depending on temperature      [J/m/s/K/K]
// static const double Kd0        = 2.222E-5; //< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
// static const double dKd_dT     =  1.37E-7; //< Coefficient of Dw depending on temperature                 [m2/s/K]
// static const double nu0        = 1.718E-5; //< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
// static const double dnu_dT     =  5.28E-8; //< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

double f1r        = 0.78;     //< ventilation factor 1 for rain
double f2r        = 0.27;     //< ventilation factor 2 for rain
double f1s        = 0.65;     //< ventilation factor 1 for snow
double f2s        = 0.39;     //< ventilation factor 2 for snow
double f1g        = 0.78;     //< ventilation factor 1 for graupel
double f2g        = 0.27;     //< ventilation factor 2 for graupel

// Freezing parameter
double A_frz      = 0.66;     //< freezing factor [/K]
double B_frz      = 100.0;    //< freezing factor [/m3/s]

// Bergeron process parameter
double mi40       = 2.46E-10; //< mass              of a 40 micron ice crystal [kg]
double mi50       = 4.80E-10; //< mass              of a 50 micron ice crystal [kg]
double vti50      = 1.0;      //< terminal velocity of a 50 micron ice crystal [m/s]
double Ri50       = 5.E-5;    //< radius            of a 50 micron ice crystal [m]

// Explicit ice generation
double sw_expice  = 0.0;      //< switch for explicit ice generation
double Nc_ihtr    = 300.0;    //< cloud number concentration for heterogeneous ice nucleation [1/cc]
double Di_max     = 500.E-6;
double Di_a       = 11.9;

double sw_constVti = 0.0;
double CONST_Vti;           // force constant terminal velocity for ice
};

extern GlobalParams para;

namespace MP_NSW6{

void mp_nsw6_init();

// void mp_nsw6(
//             int l_region,
//             double rhog          [kdim][ijdim],
//             double rhogvx        [kdim][ijdim],
//             double rhogvy        [kdim][ijdim],
//             double rhogvz        [kdim][ijdim],
//             double rhogw         [kdim][ijdim],
//             double rhoge         [kdim][ijdim],
//             double rhogq         [nqmax][kdim][ijdim],
//             double vx            [kdim][ijdim],
//             double vy            [kdim][ijdim],
//             double vz            [kdim][ijdim],
//             double w             [kdim][ijdim],
//             double UNCCN         [kdim][ijdim],
//             double rho           [kdim][ijdim],
//             double tem           [kdim][ijdim],
//             double pre           [kdim][ijdim],
//             double q             [nqmax][kdim][ijdim],
//             double qd            [kdim][ijdim],
//             double precip        [2][ijdim],
//             double precip_rhoe   [ijdim],
//             double precip_lh_heat[ijdim],
//             double precip_rhophi [ijdim],
//             double precip_rhokin [ijdim],
//             double gprec         [kdim][ijdim],
//             double rceff         [kdim][ijdim],
//             double rctop         [1][ijdim],
//             double rwtop         [1][ijdim],
//             double tctop         [1][ijdim],
//             double rceff_cld     [kdim][ijdim],
//             double rctop_cld     [1][ijdim],
//             double rwtop_cld     [1][ijdim],
//             double tctop_cld     [1][ijdim],
//             double gsgam2        [kdim][ijdim],
//             double gsgam2h       [kdim][ijdim],
//             double gam2          [kdim][ijdim],
//             double gam2h         [kdim][ijdim],
//             double ix			 [ijdim],
//             double iy			 [ijdim],
//             double iz			 [ijdim],
//             double jx			 [ijdim],
//             double jy			 [ijdim],
//             double jz			 [ijdim],
//             double z			 [kdim][ijdim],
//             double dt
//             );

/**
 * Real Kokkos ver.
 */
struct mp_nsw6_new{
private:
    //-----< Gloabal parameters >------
    GlobalParams paras;
    // Views from mp_nsw6
    int l_region                 ;
    View2D<double, DEFAULT_MEM> rhog          ;
    View2D<double, DEFAULT_MEM> rhogvx        ;
    View2D<double, DEFAULT_MEM> rhogvy        ;
    View2D<double, DEFAULT_MEM> rhogvz        ;
    View2D<double, DEFAULT_MEM> rhogw         ;
    View2D<double, DEFAULT_MEM> rhoge         ;
    View3D<double, DEFAULT_MEM> rhogq         ;
    View2D<double, DEFAULT_MEM> vx            ;
    View2D<double, DEFAULT_MEM> vy            ;
    View2D<double, DEFAULT_MEM> vz            ;
    View2D<double, DEFAULT_MEM> w             ;
    View2D<double, DEFAULT_MEM> UNCCN         ;
    View2D<double, DEFAULT_MEM> rho           ;
    View2D<double, DEFAULT_MEM> tem           ;
    View2D<double, DEFAULT_MEM> pre           ;
    View3D<double, DEFAULT_MEM> q             ;
    View2D<double, DEFAULT_MEM> qd            ;
    View2D<double, DEFAULT_MEM> precip        ;
    View1D<double, DEFAULT_MEM> precip_rhoe   ;
    View1D<double, DEFAULT_MEM> precip_lh_heat;
    View1D<double, DEFAULT_MEM> precip_rhophi ;
    View1D<double, DEFAULT_MEM> precip_rhokin ;
    View2D<double, DEFAULT_MEM> gprec         ;
    View2D<double, DEFAULT_MEM> rceff         ;
    View2D<double, DEFAULT_MEM> rctop         ;
    View2D<double, DEFAULT_MEM> rwtop         ;
    View2D<double, DEFAULT_MEM> tctop         ;
    View2D<double, DEFAULT_MEM> rceff_cld     ;
    View2D<double, DEFAULT_MEM> rctop_cld     ;
    View2D<double, DEFAULT_MEM> rwtop_cld     ;
    View2D<double, DEFAULT_MEM> tctop_cld     ;
    View2D<double, DEFAULT_MEM> gsgam2        ;
    View2D<double, DEFAULT_MEM> gsgam2h       ;
    View2D<double, DEFAULT_MEM> gam2          ;
    View2D<double, DEFAULT_MEM> gam2h         ;
    View1D<double, DEFAULT_MEM>  ix			 ;
    View1D<double, DEFAULT_MEM>  iy			 ;
    View1D<double, DEFAULT_MEM>  iz			 ;
    View1D<double, DEFAULT_MEM>  jx			 ;
    View1D<double, DEFAULT_MEM>  jy			 ;
    View1D<double, DEFAULT_MEM>  jz			 ;
    View2D<double, DEFAULT_MEM> z			 ;
    double dt;
    // Views for working
    View2D<double, DEFAULT_MEM> drhogqv;
    View2D<double, DEFAULT_MEM> drhogqc;
    View2D<double, DEFAULT_MEM> drhogqi;
    View2D<double, DEFAULT_MEM> drhogqr;
    View2D<double, DEFAULT_MEM> drhogqs;
    View2D<double, DEFAULT_MEM> drhogqg;

    View2D<double, DEFAULT_MEM> psatl;
    View2D<double, DEFAULT_MEM> qsatl;                //< saturated water vapor for liquid water [kg/kg]
    View2D<double, DEFAULT_MEM> psati;
    View2D<double, DEFAULT_MEM> qsati;                //< saturated water vapor for ice water    [kg/kg]
    View2D<double, DEFAULT_MEM> Nc   ;                //< Number concentration of cloud water [1/cc]


    //---< Bergeron process >---
    View2D<double, DEFAULT_MEM> a1      ;           //<
    View2D<double, DEFAULT_MEM> a2      ;           //<
    View2D<double, DEFAULT_MEM> ma2     ;           //< 1-a2


    //---< Effective radius >---
    static constexpr double r2_min = 1.0E-10;
    static constexpr double q_min  = 1.0E-5;
    double dgamma_a, GAM_dgam23, GAM_dgam;
    double coef_dgam, coef_xf;

    //---< Precipitation >---
    bool preciptation_flag[nqmax] = {false};

    View3D<double, DEFAULT_MEM> Vt    ;
    View2D<double, DEFAULT_MEM> cva   ;
    View2D<double, DEFAULT_MEM> rgs   ;
    View2D<double, DEFAULT_MEM> rgsh  ;

    View2D<double, DEFAULT_MEM> ml_Pconv;
    View2D<double, DEFAULT_MEM> ml_Pconw;
    View2D<double, DEFAULT_MEM> ml_Pconi;

    double UNDEF, EPS, PI, Rvap, LHV0, LHS0, LHF0, PRE00;
public:
    // constructor
    mp_nsw6_new( int l_region_            ,
            View2D<double, DEFAULT_MEM> rhog_          ,
            View2D<double, DEFAULT_MEM> rhogvx_        ,
            View2D<double, DEFAULT_MEM> rhogvy_        ,
            View2D<double, DEFAULT_MEM> rhogvz_        ,
            View2D<double, DEFAULT_MEM> rhogw_         ,
            View2D<double, DEFAULT_MEM> rhoge_         ,
            View3D<double, DEFAULT_MEM> rhogq_         ,
            View2D<double, DEFAULT_MEM> vx_            ,
            View2D<double, DEFAULT_MEM> vy_            ,
            View2D<double, DEFAULT_MEM> vz_            ,
            View2D<double, DEFAULT_MEM> w_             ,
            View2D<double, DEFAULT_MEM> UNCCN_         ,
            View2D<double, DEFAULT_MEM> rho_           ,
            View2D<double, DEFAULT_MEM> tem_           ,
            View2D<double, DEFAULT_MEM> pre_           ,
            View3D<double, DEFAULT_MEM> q_             ,
            View2D<double, DEFAULT_MEM> qd_            ,
            View2D<double, DEFAULT_MEM> precip_        ,
            View1D<double, DEFAULT_MEM> precip_rhoe_   ,
            View1D<double, DEFAULT_MEM> precip_lh_heat_,
            View1D<double, DEFAULT_MEM> precip_rhophi_ ,
            View1D<double, DEFAULT_MEM> precip_rhokin_ ,
            View2D<double, DEFAULT_MEM> gprec_         ,
            View2D<double, DEFAULT_MEM> rceff_         ,
            View2D<double, DEFAULT_MEM> rctop_         ,
            View2D<double, DEFAULT_MEM> rwtop_         ,
            View2D<double, DEFAULT_MEM> tctop_         ,
            View2D<double, DEFAULT_MEM> rceff_cld_     ,
            View2D<double, DEFAULT_MEM> rctop_cld_     ,
            View2D<double, DEFAULT_MEM> rwtop_cld_     ,
            View2D<double, DEFAULT_MEM> tctop_cld_     ,
            View2D<double, DEFAULT_MEM> gsgam2_        ,
            View2D<double, DEFAULT_MEM> gsgam2h_       ,
            View2D<double, DEFAULT_MEM> gam2_          ,
            View2D<double, DEFAULT_MEM> gam2h_         ,
            View1D<double, DEFAULT_MEM> ix_			,
            View1D<double, DEFAULT_MEM> iy_			,
            View1D<double, DEFAULT_MEM> iz_			,
            View1D<double, DEFAULT_MEM> jx_			,
            View1D<double, DEFAULT_MEM> jy_			,
            View1D<double, DEFAULT_MEM> jz_			,
            View2D<double, DEFAULT_MEM> z_			    , 
            double dt_,
            GlobalParams& paras_) : l_region(l_region_),
                          rhog  (rhog_),
                          rhogvx(rhogvx_),
                          rhogvy(rhogvy_),
                          rhogvz(rhogvz_),
                          rhogw (rhogw_),
                          rhoge (rhoge_),
                          rhogq (rhogq_),
                          vx    (vx_),
                          vy    (vy_),
                          vz    (vz_),
                          w     (w_ ),
                          UNCCN (UNCCN_),
                          rho   (rho_),
                          tem   (tem_),
                          pre   (pre_),
                          q     (q_  ),
                          qd    (qd_ ),
                          precip        (precip_        ),
                          precip_rhoe   (precip_rhoe_   ),
                          precip_lh_heat(precip_lh_heat_),
                          precip_rhophi (precip_rhophi_ ),
                          precip_rhokin (precip_rhokin_ ),
                          gprec(gprec_)         ,
                          rceff(rceff_)         ,
                          rctop(rctop_)         ,
                          rwtop(rwtop_)         ,
                          tctop(tctop_)         ,
                          rceff_cld(rceff_cld_)     ,
                          rctop_cld(rctop_cld_)     ,
                          rwtop_cld(rwtop_cld_)     ,
                          tctop_cld(tctop_cld_)     ,
                          gsgam2(gsgam2_)        ,
                          gsgam2h(gsgam2h_)       ,
                          gam2(gam2_)          ,
                          gam2h(gam2h_)         ,
                          ix(ix_)			,
                          iy(iy_)			,
                          iz(iz_)			,
                          jx(jx_)			,
                          jy(jy_)			,
                          jz(jz_)			,
                          z(z_)			    , 
                          dt(dt_),
                          paras(paras_)           
    { init(); }
    void init();

public:
    // ----- < Loop Tags > ------
    struct TagAerosolT{};
    struct TagAerosolF{};
    struct TagSATWater{};
    struct TagBigLoop{};
    struct TagUpdateRHOGQ{};
    struct TagUpdateRHOGE{};
    struct TagUpdateCLD{};
    struct TagPrecip{};
    struct TagUpdateMassConc{};
    struct TagUpdateTEM{};
    struct TagUpdatePCON{};
    struct TagUpdatePCON2{};

    void compute();

    KOKKOS_FUNCTION
    void operator()(TagAerosolT, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagAerosolF, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagSATWater, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagBigLoop, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdateRHOGQ, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdateRHOGE, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdateCLD, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagPrecip, const size_t nq, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdateMassConc, const size_t nq, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdateTEM, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdatePCON, const size_t k, const size_t ij) const;

    KOKKOS_FUNCTION
    void operator()(TagUpdatePCON2, const size_t k, const size_t ij) const;
};

/**
 * Naive Kokkos Ver.
 */
// void mp_nsw6(
//             int l_region,
//             const View<double**>&  rhog          ,
//             const View<double**>&  rhogvx        ,
//             const View<double**>&  rhogvy        ,
//             const View<double**>&  rhogvz        ,
//             const View<double**>&  rhogw         ,
//             const View<double**>&  rhoge         ,
//             const View<double***>& rhogq         ,
//             const View<double**>&  vx            ,
//             const View<double**>&  vy            ,
//             const View<double**>&  vz            ,
//             const View<double**>&  w             ,
//             const View<double**>&  UNCCN         ,
//             const View<double**>&  rho           ,
//             const View<double**>&  tem           ,
//             const View<double**>&  pre           ,
//             const View<double***>& q             ,
//             const View<double**>&  qd            ,
//             const View<double**>&  precip        ,
//             const View<double*>&   precip_rhoe   ,
//             const View<double*>&   precip_lh_heat,
//             const View<double*>&   precip_rhophi ,
//             const View<double*>&   precip_rhokin ,
//             const View<double**>&  gprec         ,
//             const View<double**>&  rceff         ,
//             const View<double**>&  rctop         ,
//             const View<double**>&  rwtop         ,
//             const View<double**>&  tctop         ,
//             const View<double**>&  rceff_cld     ,
//             const View<double**>&  rctop_cld     ,
//             const View<double**>&  rwtop_cld     ,
//             const View<double**>&  tctop_cld     ,
//             const View<double**>&  gsgam2        ,
//             const View<double**>&  gsgam2h       ,
//             const View<double**>&  gam2          ,
//             const View<double**>&  gam2h         ,
//             const View<double*>&   ix			   ,
//             const View<double*>&   iy			   ,
//             const View<double*>&   iz			   ,
//             const View<double*>&   jx			   ,
//             const View<double*>&   jy			   ,
//             const View<double*>&   jz			   ,
//             const View<double**>&  z			   ,
//             double dt
//             );
}

/**
 * Deprecated variables
 */

// double dens    ;                         //< density
// double temp    ;                         //< T [K]
// double qv      ;                         //< mixing ratio of water vapor  [kg/kg]
// double qc      ;                         //< mixing ratio of liquid water [kg/kg]
// double qr      ;                         //< mixing ratio of rain         [kg/kg]
// double qi      ;                         //< mixing ratio of ice water    [kg/kg]
// double qs      ;                         //< mixing ratio of snow         [kg/kg]
// double qg      ;                         //< mixing ratio of graupel      [kg/kg]
// double qv_t    ;                         //< tendency     of water vapor  [kg/kg/s]
// double qc_t    ;                         //< tendency     of liquid water [kg/kg/s]
// double qr_t    ;                         //< tendency     of rain         [kg/kg/s]
// double qi_t    ;                         //< tendency     of ice water    [kg/kg/s]
// double qs_t    ;                         //< tendency     of snow         [kg/kg/s]
// double qg_t    ;                         //< tendency     of graupel      [kg/kg/s]
// double Sliq    ;                         //< saturated ratio S for liquid water [0-1]
// double Sice    ;                         //< saturated ratio S for ice water    [0-1]
// double Rdens   ;                         //< 1 / density
// double rho_fact;                         //< density factor
// double temc    ;                         //< T - T0 [K]

// double RLMDr, RLMDr_2, RLMDr_3        ;
// double RLMDs, RLMDs_2, RLMDs_3        ;
// double RLMDg, RLMDg_2, RLMDg_3        ;
// double RLMDr_1br, RLMDr_2br, RLMDr_3br;
// double RLMDs_1bs, RLMDs_2bs, RLMDs_3bs;
// double RLMDr_dr, RLMDr_3dr, RLMDr_5dr ;
// double RLMDs_ds, RLMDs_3ds, RLMDs_5ds ;
// double RLMDg_dg, RLMDg_3dg, RLMDg_5dg ;
// double RLMDr_7                        ;
// double RLMDr_6dr                      ;

//---< Roh and Satoh (2014) >---
// double tems, Xs2                     ;
// double MOMs_0, MOMs_1, MOMs_2        ;
// double MOMs_0bs, MOMs_1bs, MOMs_2bs  ;
// double MOMs_2ds, MOMs_5ds_h, RMOMs_Vt;
// double coef_at[4]                    ;
// double coef_bt[4]                    ;
// double loga_, b_, nm                 ;

// double Vti, Vtr, Vts, Vtg   ;           //< terminal velocity
// double Esi_mod, Egs_mod     ;           //< modified accretion efficiency
// double rhoqc                ;           //< rho * qc
// double Pracw_orig,  Pracw_kk;           //< accretion       term by orig  & k-k scheme
// double Praut_berry, Praut_kk;           //< auto-conversion term by berry & k-k scheme
// double Dc                   ;           //< relative variance
// double betai, betas         ;           //< sticky parameter for auto-conversion
// double Ka                   ;           //< thermal diffusion coefficient of air
// double Kd                   ;           //< diffusion coefficient of water vapor in air
// double Nu                   ;           //< kinematic viscosity of air
// double Glv, Giv, Gil        ;           //< thermodynamic function
// double ventr, vents, ventg  ;           //< ventilation factor
// double net, fac, fac_sw     ;
// double zerosw, tmp          ;

//---< Bergeron process >---
// double sw_bergeron     ;                 //< if 0C<T<30C, sw=1
// double dt1             ;                 //< time during which the an ice particle of 40um grows to 50um
// double Ni50            ;                 //< number concentration of ice particle of 50um

//---< Explicit ice generation >---
// double sw, rhoqi, XNi, XMi, Di, Ni0, Qi0;

//---< Effective radius >---
// double xf_qc;                 //< mean mass   of qc [kg]
// double rf_qc;                 //< mean radius of qc [m]
// double r2_qc;                 //< r^2 moment  of qc
// double r2_qr;                 //< r^2 moment  of qr
// double r3_qc;                 //< r^3 moment  of qc
// double r3_qr;                 //< r^3 moment  of qr

// double wk[wk_nmax];
// double ml_wk   [wk_nmax][kdim][ijdim]; tentatively remove due to the stack overflow