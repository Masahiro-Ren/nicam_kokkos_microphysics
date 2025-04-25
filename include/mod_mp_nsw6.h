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

namespace MP_NSW6{

void mp_nsw6_init();

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

/**
 * Real Kokkos ver.
 */
struct mp_nsw6_new{
private:
    struct BigLoopTag{};
    // Views from mp_nsw6
    int l_region                 ;
    View<double**> rhog          ;
    View<double**> rhogvx        ;
    View<double**> rhogvy        ;
    View<double**> rhogvz        ;
    View<double**> rhogw         ;
    View<double**> rhoge         ;
    View<double***>rhogq         ;
    View<double**> vx            ;
    View<double**> vy            ;
    View<double**> vz            ;
    View<double**> w             ;
    View<double**> UNCCN         ;
    View<double**> rho           ;
    View<double**> tem           ;
    View<double**> pre           ;
    View<double***>q             ;
    View<double**> qd            ;
    View<double**> precip        ;
    View<double*>  precip_rhoe   ;
    View<double*>  precip_lh_heat;
    View<double*>  precip_rhophi ;
    View<double*>  precip_rhokin ;
    View<double**> gprec         ;
    View<double**> rceff         ;
    View<double**> rctop         ;
    View<double**> rwtop         ;
    View<double**> tctop         ;
    View<double**> rceff_cld     ;
    View<double**> rctop_cld     ;
    View<double**> rwtop_cld     ;
    View<double**> tctop_cld     ;
    View<double**> gsgam2        ;
    View<double**> gsgam2h       ;
    View<double**> gam2          ;
    View<double**> gam2h         ;
    View<double*>  ix			 ;
    View<double*>  iy			 ;
    View<double*>  iz			 ;
    View<double*>  jx			 ;
    View<double*>  jy			 ;
    View<double*>  jz			 ;
    View<double**> z			 ;
    double dt;
    // Views for working
    View<double**> drhogqv;
    View<double**> drhogqc;
    View<double**> drhogqi;
    View<double**> drhogqr;
    View<double**> drhogqs;
    View<double**> drhogqg;

    View<double**> psatl;
    View<double**> qsatl;                //< saturated water vapor for liquid water [kg/kg]
    View<double**> psati;
    View<double**> qsati;                //< saturated water vapor for ice water    [kg/kg]
    View<double**> Nc   ;                //< Number concentration of cloud water [1/cc]


    //---< Bergeron process >---
    View<double**> a1      ;           //<
    View<double**> a2      ;           //<
    View<double**> ma2     ;           //< 1-a2


    //---< Effective radius >---
    const double r2_min = 1.0E-10;
    const double q_min  = 1.0E-5;
    double dgamma_a, GAM_dgam23, GAM_dgam;
    double coef_dgam, coef_xf;

    //---< Precipitation >---
    bool preciptation_flag[nqmax] = {false};

    View<double***> Vt    ;
    View<double**>  cva   ;
    View<double**>  rgs   ;
    View<double**>  rgsh  ;

    View<double**> ml_Pconv;
    View<double**> ml_Pconw;
    View<double**> ml_Pconi;

    double UNDEF, EPS, PI, Rvap, LHV0, LHS0, LHF0, PRE00;
public:
    // constructor
    mp_nsw6_new( int l_region_            ,
            View<double**> rhog_          ,
            View<double**> rhogvx_        ,
            View<double**> rhogvy_        ,
            View<double**> rhogvz_        ,
            View<double**> rhogw_         ,
            View<double**> rhoge_         ,
            View<double***>rhogq_         ,
            View<double**> vx_            ,
            View<double**> vy_            ,
            View<double**> vz_            ,
            View<double**> w_             ,
            View<double**> UNCCN_         ,
            View<double**> rho_           ,
            View<double**> tem_           ,
            View<double**> pre_           ,
            View<double***>q_             ,
            View<double**> qd_            ,
            View<double**> precip_        ,
            View<double*>  precip_rhoe_   ,
            View<double*>  precip_lh_heat_,
            View<double*>  precip_rhophi_ ,
            View<double*>  precip_rhokin_ ,
            View<double**> gprec_         ,
            View<double**> rceff_         ,
            View<double**> rctop_         ,
            View<double**> rwtop_         ,
            View<double**> tctop_         ,
            View<double**> rceff_cld_     ,
            View<double**> rctop_cld_     ,
            View<double**> rwtop_cld_     ,
            View<double**> tctop_cld_     ,
            View<double**> gsgam2_        ,
            View<double**> gsgam2h_       ,
            View<double**> gam2_          ,
            View<double**> gam2h_         ,
            View<double*>  ix_			,
            View<double*>  iy_			,
            View<double*>  iz_			,
            View<double*>  jx_			,
            View<double*>  jy_			,
            View<double*>  jz_			,
            View<double**> z_			    , 
            double dt_) : l_region(l_region_),
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
                          dt(dt_)           
    { init(); }
    void init();
public:
    void compute();

    KOKKOS_FUNCTION
    void operator()(BigLoopTag, const size_t k, const size_t ij) const;
};

/**
 * Naive Kokkos Ver.
 */
void mp_nsw6(
            int l_region,
            const View<double**>&  rhog          ,
            const View<double**>&  rhogvx        ,
            const View<double**>&  rhogvy        ,
            const View<double**>&  rhogvz        ,
            const View<double**>&  rhogw         ,
            const View<double**>&  rhoge         ,
            const View<double***>& rhogq         ,
            const View<double**>&  vx            ,
            const View<double**>&  vy            ,
            const View<double**>&  vz            ,
            const View<double**>&  w             ,
            const View<double**>&  UNCCN         ,
            const View<double**>&  rho           ,
            const View<double**>&  tem           ,
            const View<double**>&  pre           ,
            const View<double***>& q             ,
            const View<double**>&  qd            ,
            const View<double**>&  precip        ,
            const View<double*>&   precip_rhoe   ,
            const View<double*>&   precip_lh_heat,
            const View<double*>&   precip_rhophi ,
            const View<double*>&   precip_rhokin ,
            const View<double**>&  gprec         ,
            const View<double**>&  rceff         ,
            const View<double**>&  rctop         ,
            const View<double**>&  rwtop         ,
            const View<double**>&  tctop         ,
            const View<double**>&  rceff_cld     ,
            const View<double**>&  rctop_cld     ,
            const View<double**>&  rwtop_cld     ,
            const View<double**>&  tctop_cld     ,
            const View<double**>&  gsgam2        ,
            const View<double**>&  gsgam2h       ,
            const View<double**>&  gam2          ,
            const View<double**>&  gam2h         ,
            const View<double*>&   ix			   ,
            const View<double*>&   iy			   ,
            const View<double*>&   iz			   ,
            const View<double*>&   jx			   ,
            const View<double*>&   jy			   ,
            const View<double*>&   jz			   ,
            const View<double**>&  z			   ,
            double dt
            );
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