#include "mod_mp_nsw6.h"

/**
 * Private parameters
 */
static const double dens00 = 1.28;
static const int wk_nmax = 49;

enum{
I_dqv_dt  , // 
I_dqc_dt  , // 
I_dqr_dt  , // 
I_dqi_dt  , // 
I_dqs_dt  , // 
I_dqg_dt  , // 
I_delta1  , // separation switch for r->s,g
I_delta2  , // separation switch for s->g
I_spsati  , // separation switch for ice sublimation
I_iceflg  , // separation switch for T > 0
I_RLMDr   ,
I_RLMDs   ,
I_RLMDg   ,
I_Piacr   , // r->s,g
I_Psacr   , // r->s,g
I_Praci   , // i->s,g
I_Pigen   , // v->i
I_Pidep   , // v->i
I_Psdep   , // v->s
I_Pgdep   , // v->g
I_Praut   , // c->r
I_Pracw   , // c->r
I_Pihom   , // c->i
I_Pihtr   , // c->i
I_Psacw   , // c->s
I_Psfw    , // c->s
I_Pgacw   , // c->g
I_Prevp   , // r->v
I_Piacr_s , // r->s
I_Psacr_s , // r->s
I_Piacr_g , // r->g
I_Psacr_g , // r->g
I_Pgacr   , // r->g
I_Pgfrz   , // r->g
I_Pisub   , // i->v
I_Pimlt   , // i->c
I_Psaut   , // i->s
I_Praci_s , // i->s
I_Psaci   , // i->s
I_Psfi    , // i->s
I_Praci_g , // i->g
I_Pgaci   , // i->g
I_Pssub   , // s->v
I_Psmlt   , // s->r
I_Pgaut   , // s->g
I_Pracs   , // s->g
I_Pgacs   , // s->g
I_Pgsub   , // g->v
I_Pgmlt     // g->r
};

std::string w_name[wk_nmax] = {
"dqv_dt",
"dqc_dt",
"dqr_dt",
"dqi_dt",
"dqs_dt",
"dqg_dt",
"delta1",
"delta2",
"spsati",
"iceflg",
"RLMDr",
"RLMDs",
"RLMDg",
"Piacr",
"Psacr",
"Praci",
"Pigen",
"Pidep",
"Psdep",
"Pgdep",
"Praut",
"Pracw",
"Pihom",
"Pihtr",
"Psacw",
"Psfw ",
"Pgacw",
"Prevp",
"Piacr_s",
"Psacr_s",
"Piacr_g",
"Psacr_g",
"Pgacr",
"Pgfrz",
"Pisub",
"Pimlt",
"Psaut",
"Praci_s",
"Psaci",
"Psfi",
"Praci_g",
"Pgaci",
"Pssub",
"Psmlt",
"Pgaut",
"Pracs",
"Pgacs",
"Pgsub",
"Pgmlt"
};

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
static const double Nc_lnd     = 2000.0;   //< number concentration of cloud water (land)  [1/cc]
static const double Nc_ocn     =   50.0;   //< number concentration of cloud water (ocean) [1/cc]
double Nc_def;                             //< number concentration of cloud water         [1/cc]

double beta_saut  =  1.E-3;   //< auto-conversion factor beta  for ice
double gamma_saut = 25.E-3;   //< auto-conversion factor gamma for ice
double beta_gaut  =  1.E-3;   //< auto-conversion factor beta  for snow
double gamma_gaut = 90.E-3;   //< auto-conversion factor gamma for snow
double qicrt_saut =  0.0;     //< mixing ratio threshold for Psaut [kg/kg]
double qscrt_gaut =  6.E-4;   //< mixing ratio threshold for Pgaut [kg/kg]

// Evaporation, Sublimation parameter
static const double Ka0        = 2.428E-2; //< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
static const double dKa_dT     =  7.47E-5; //< Coefficient of Ka depending on temperature      [J/m/s/K/K]
static const double Kd0        = 2.222E-5; //< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
static const double dKd_dT     =  1.37E-7; //< Coefficient of Dw depending on temperature                 [m2/s/K]
static const double nu0        = 1.718E-5; //< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
static const double dnu_dT     =  5.28E-8; //< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

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


std::string precip_transport_type = "3WATER"; 
std::string precip_scheme_type = "Default";

bool OPT_EXPLICIT_ICEGEN   = false;   // enable explicit ice generation?
bool OPT_INDIR             = false;   // enable aerosol indirect effect?
bool Roh_flag              = false;   // enable setting by Roh and Satoh (2014)?

double sw_constVti = 0.0;
double CONST_Vti;           // force constant terminal velocity for ice

namespace MP_NSW6{

};