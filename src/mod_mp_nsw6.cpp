#include "mod_mp_nsw6.h"

GlobalParams para;

/**
 * Private parameters
 */
constexpr double dens00 = 1.28;
// static const int wk_nmax = 49;

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

// ---------- < Parameters no need to copy > ----------------------
// Auto-conversion parameter
static constexpr double Nc_lnd     = 2000.0;   //< number concentration of cloud water (land)  [1/cc]
static constexpr double Nc_ocn     =   50.0;   //< number concentration of cloud water (ocean) [1/cc]
// Evaporation, Sublimation parameter
static constexpr double Ka0        = 2.428E-2; //< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
static constexpr double dKa_dT     =  7.47E-5; //< Coefficient of Ka depending on temperature      [J/m/s/K/K]
static constexpr double Kd0        = 2.222E-5; //< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
static constexpr double dKd_dT     =  1.37E-7; //< Coefficient of Dw depending on temperature                 [m2/s/K]
static constexpr double nu0        = 1.718E-5; //< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
static constexpr double dnu_dT     =  5.28E-8; //< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

std::string precip_transport_type = "3WATER"; 
std::string precip_scheme_type = "Default";

bool OPT_EXPLICIT_ICEGEN   = false;   // enable explicit ice generation?
bool OPT_INDIR             = false;   // enable aerosol indirect effect?
bool Roh_flag              = false;   // enable setting by Roh and Satoh (2014)?
// ---------- < END Parameters no need to copy > -------------------

// double N0r    = 8.E+6;  //< intercept parameter for rain    [1/m4]
// double N0s    = 3.E+6;  //< intercept parameter for snow    [1/m4]
// double N0g    = 4.E+6;  //< intercept parameter for graupel [1/m4]

// double rho_s  = 100.0;  //< density of snow    [kg/m3]
// double rho_g  = 400.0;  //< density of graupel [kg/m3]
//                         //   graupel : 400
//                         //   hail    : 917

// double C_d    =   0.6;  //< drag coefficient for graupel
// double Cr     = 130.0;
// double Cs     =  4.84;

// // Empirical parameter
// double Ar, As, Ag;
// double Br, Bs, Bg;
// double Cg        ;
// double Dr, Ds, Dg;

// // GAMMA function
// double GAM, GAM_2, GAM_3;

// double GAM_1br, GAM_2br, GAM_3br;
// double GAM_3dr;
// double GAM_6dr;
// double GAM_1brdr;
// double GAM_5dr_h;

// double GAM_1bs, GAM_2bs, GAM_3bs;
// double GAM_3ds;
// double GAM_1bsds;
// double GAM_5ds_h;

// double GAM_1bg, GAM_3dg;
// double GAM_1bgdg;
// double GAM_5dg_h;

// //---< Khairoutdinov and Kogan (2000) >---
// double sw_kk2000  = 0.0; //< switch for k-k scheme

// //---< Roh and Satoh (2014) >---
// double sw_roh2014 = 0.0; //< switch for Roh scheme
// double ln10;             //< log(10)
// double coef_a[10] = {5.065339, -0.062659, -3.032362, 0.029469, -0.000285, 0.31255, 0.000204, 0.003199, 0.0, -0.015952};
// double coef_b[10] = {0.476221, -0.015896,  0.165977, 0.007468, -0.000141, 0.060366, 0.000079, 0.000594, 0.0, -0.003577};

// // Accretion parameter
// double Eiw        = 1.0;      //< collection efficiency of cloud ice for cloud water
// double Erw        = 1.0;      //< collection efficiency of rain    for cloud water
// double Esw        = 1.0;      //< collection efficiency of snow    for cloud water
// double Egw        = 1.0;      //< collection efficiency of graupel for cloud water
// double Eri        = 1.0;      //< collection efficiency of rain    for cloud ice
// double Esi        = 1.0;      //< collection efficiency of snow    for cloud ice
// double Egi        = 0.1;      //< collection efficiency of graupel for cloud ice
// double Esr        = 1.0;      //< collection efficiency of snow    for rain
// double Egr        = 1.0;      //< collection efficiency of graupel for rain
// double Egs        = 1.0;      //< collection efficiency of graupel for snow
// double gamma_sacr = 25.E-3;   //< effect of low temperature for Esi
// double gamma_gacs = 90.E-3;   //< effect of low temperature for Egs
// double mi         = 4.19E-13; //< mass of one cloud ice crystal [kg]

// // Auto-conversion parameter
// // static const double Nc_lnd     = 2000.0;   //< number concentration of cloud water (land)  [1/cc]
// // static const double Nc_ocn     =   50.0;   //< number concentration of cloud water (ocean) [1/cc]
// double Nc_def;                             //< number concentration of cloud water         [1/cc]

// double beta_saut  =  1.E-3;   //< auto-conversion factor beta  for ice
// double gamma_saut = 25.E-3;   //< auto-conversion factor gamma for ice
// double beta_gaut  =  1.E-3;   //< auto-conversion factor beta  for snow
// double gamma_gaut = 90.E-3;   //< auto-conversion factor gamma for snow
// double qicrt_saut =  0.0;     //< mixing ratio threshold for Psaut [kg/kg]
// double qscrt_gaut =  6.E-4;   //< mixing ratio threshold for Pgaut [kg/kg]

// // // Evaporation, Sublimation parameter
// // static const double Ka0        = 2.428E-2; //< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
// // static const double dKa_dT     =  7.47E-5; //< Coefficient of Ka depending on temperature      [J/m/s/K/K]
// // static const double Kd0        = 2.222E-5; //< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
// // static const double dKd_dT     =  1.37E-7; //< Coefficient of Dw depending on temperature                 [m2/s/K]
// // static const double nu0        = 1.718E-5; //< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
// // static const double dnu_dT     =  5.28E-8; //< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

// double f1r        = 0.78;     //< ventilation factor 1 for rain
// double f2r        = 0.27;     //< ventilation factor 2 for rain
// double f1s        = 0.65;     //< ventilation factor 1 for snow
// double f2s        = 0.39;     //< ventilation factor 2 for snow
// double f1g        = 0.78;     //< ventilation factor 1 for graupel
// double f2g        = 0.27;     //< ventilation factor 2 for graupel

// // Freezing parameter
// double A_frz      = 0.66;     //< freezing factor [/K]
// double B_frz      = 100.0;    //< freezing factor [/m3/s]

// // Bergeron process parameter
// double mi40       = 2.46E-10; //< mass              of a 40 micron ice crystal [kg]
// double mi50       = 4.80E-10; //< mass              of a 50 micron ice crystal [kg]
// double vti50      = 1.0;      //< terminal velocity of a 50 micron ice crystal [m/s]
// double Ri50       = 5.E-5;    //< radius            of a 50 micron ice crystal [m]

// // Explicit ice generation
// double sw_expice  = 0.0;      //< switch for explicit ice generation
// double Nc_ihtr    = 300.0;    //< cloud number concentration for heterogeneous ice nucleation [1/cc]
// double Di_max     = 500.E-6;
// double Di_a       = 11.9;

// double sw_constVti = 0.0;
// double CONST_Vti;           // force constant terminal velocity for ice

void negative_filter( double rhog  [kdim][ijdim],
                      double rhoge [kdim][ijdim],
                      double rhogq [nqmax][kdim][ijdim],
                      double rho   [kdim][ijdim],
                      double tem   [kdim][ijdim],
                      double pre   [kdim][ijdim],
                      double q     [nqmax][kdim][ijdim],
                      double gsgam2[kdim][ijdim]  );

void Bergeron_param( double tem[kdim][ijdim],
                     double a1 [kdim][ijdim],
                     double a2 [kdim][ijdim],
                     double ma2[kdim][ijdim]);

/**
 * Kokkos ver.
 */
void negative_filter( const View2D<double, DEFAULT_MEM>&  rhog,
                      const View2D<double, DEFAULT_MEM>&  rhoge, 
                      const View3D<double, DEFAULT_MEM>&  rhogq,
                      const View2D<double, DEFAULT_MEM>&  rho,   
                      const View2D<double, DEFAULT_MEM>&  tem,   
                      const View2D<double, DEFAULT_MEM>&  pre,   
                      const View3D<double, DEFAULT_MEM>&  q,     
                      const View2D<double, DEFAULT_MEM>&  gsgam2);

void Bergeron_param( const View2D<double, DEFAULT_MEM> tem,
                     const View2D<double, DEFAULT_MEM> a1 ,
                     const View2D<double, DEFAULT_MEM> a2 ,
                     const View2D<double, DEFAULT_MEM> ma2);


namespace MP_NSW6{

void mp_nsw6_init()
{
    bool INC_PGAUT = false;
    double PSAUT_BETA0;

    std::string qr_aut_acc_type = "Default";

    para.Nc_def = Nc_ocn;
    para.CONST_Vti = PROBLEM_SIZE::UNDEF;
    PSAUT_BETA0 = para.beta_saut;

    para.gamma_sacr = 6.0E-2;
    PSAUT_BETA0 = 5.0E-3;
    para.gamma_saut = 6.0E-2;
    precip_scheme_type = "Flux-Semilag_new";
    Roh_flag = true;
    qr_aut_acc_type = "KK2000";

    para.beta_saut = PSAUT_BETA0;

    std::cout << std::endl;
    std::cout << "*** Calculation flag of sedimentation: \n";
    std::cout << "*** QV => NO \n";
    std::cout << "*** QC => NO \n";
    std::cout << "*** QR => YES \n";

    if(precip_transport_type == "3WATER")
    {
        std::cout << "*** QI => NO \n";
    }
    else if(precip_transport_type == "4WATER")
    {
        std::cout << "*** QI => YES \n";
    }

    std::cout << "*** QS => YES \n";
    std::cout << "*** QG => YES \n\n";

    std::cout << "*** Precipitation(sedimentation) scheme: \n";
    if(precip_scheme_type == "Upwind-Euler")
    {
        std::cout << "*** => Upwind-Euler\n";
    }
    else if (precip_scheme_type == "Flux-Semilag_new")
    {
        std::cout << "*** => Flux-Semilag_new\n";
    }
    else
    {
        std::cout << "*** => Default(Flux-Semilag) \n";
    }

    //--- empirical coefficients A, B, C, D
    para.Ar = PI * rho_w / 6.0;
    para.As = PI * para.rho_s / 6.0;
    para.Ag = PI * para.rho_g / 6.0;

    para.Br = 3.0;
    para.Bs = 3.0;
    para.Bg = 3.0;
    
    para.Cg = std::sqrt( ( 4.0 * para.rho_g * GRAV ) / ( 3.0 * dens00 * para.C_d ) );

    para.Dr = 0.50;
    para.Ds = 0.25;
    para.Dg = 0.50;

    std::cout << std::endl;
    std::cout << "*** Use setting of Roh and Satoh(2014)?: \n";
    if(Roh_flag)
    {
        std::cout << "*** => Yes\n";
        OPT_EXPLICIT_ICEGEN = true;
        para.sw_roh2014 = 1.0;
        para.N0g        = 4.0E+8;
        para.As         = 0.069;
        para.Bs         = 2.0;
        para.Esi        = 0.25;
        para.Egi        = 0.0;
        para.Egs        = 0.0;
    }
    else
    {
        std::cout << "*** => No: default. \n";
    }
    
    std::cout << std::endl;
    std::cout << "*** Use explicit ice generation scheme?:\n";
    if(OPT_EXPLICIT_ICEGEN)
    {
        std::cout << "*** => Yes \n";
        para.sw_expice = 1.0;
    }
    else
    {
        std::cout << "*** => No: default.\n";
        para.sw_expice = 0.0;
    }

    std::cout << std::endl;
    std::cout << "*** Autoconversion & Accretion scheme for QC->Qr:\n";
    if(qr_aut_acc_type == "Default")
    {
        std::cout << "*** => Berry(1968) : default \n";
        para.sw_kk2000 = 0.0;
    }
    else if (qr_aut_acc_type == "KK2000")
    {
        std::cout << "*** => Khairoutdinov and Kogan(2000) \n";
        para.sw_kk2000 = 1.0;
    }
    else
    {
        std::cerr << "Not appropriate qr_aut_acc_type. Type: " << qr_aut_acc_type << std::endl;
        std::cerr << "STOP! \n";
        ADM_Proc_stop();
    }

    if(!INC_PGAUT)
    {
        para.beta_gaut = 0.0;
    }

    if(para.CONST_Vti != UNDEF)
    {
        para.CONST_Vti = std::abs(para.CONST_Vti);
        para.sw_constVti = 1.0;
    }
    else
    {
        para.sw_constVti = 0.0;
    }

    para.GAM       = 1.0; // =0!
    para.GAM_2     = 1.0; // =1!
    para.GAM_3     = 2.0; // =2!

    para.GAM_1br   = MISC_gammafunc( 1.0 + para.Br ); // = 4!
    para.GAM_2br   = MISC_gammafunc( 2.0 + para.Br ); // = 5!
    para.GAM_3br   = MISC_gammafunc( 3.0 + para.Br ); // = 6!
    para.GAM_3dr   = MISC_gammafunc( 3.0 + para.Dr );
    para.GAM_6dr   = MISC_gammafunc( 6.0 + para.Dr );
    para.GAM_1brdr = MISC_gammafunc( 1.0 + para.Br + para.Dr );
    para.GAM_5dr_h = MISC_gammafunc( 0.5 * (5.0 + para.Dr) );

    para.GAM_1bs   = MISC_gammafunc( 1.0 + para.Bs ); // = 4!
    para.GAM_2bs   = MISC_gammafunc( 2.0 + para.Bs ); // = 5!
    para.GAM_3bs   = MISC_gammafunc( 3.0 + para.Bs ); // = 6!
    para.GAM_3ds   = MISC_gammafunc( 3.0 + para.Ds );
    para.GAM_1bsds = MISC_gammafunc( 1.0 + para.Bs + para.Ds );
    para.GAM_5ds_h = MISC_gammafunc( 0.5 * (5.0 + para.Ds) );

    para.GAM_1bg   = MISC_gammafunc( 1.0 + para.Bg ); // = 4!
    para.GAM_3dg   = MISC_gammafunc( 3.0 + para.Dg );
    para.GAM_1bgdg = MISC_gammafunc( 1.0 + para.Bg + para.Dg);
    para.GAM_5dg_h = MISC_gammafunc( 0.5 * (5.0 + para.Dg) );

    para.ln10 = std::log(10.0);
}


// ================================ mp_nsw6_new Kokkos Implementation =========================

void mp_nsw6_new::init()
{
    drhogqv  = View2D<double, DEFAULT_MEM>("drhogqv", kdim, ijdim);
    drhogqc  = View2D<double, DEFAULT_MEM>("drhogqc", kdim, ijdim);
    drhogqi  = View2D<double, DEFAULT_MEM>("drhogqi", kdim, ijdim);
    drhogqr  = View2D<double, DEFAULT_MEM>("drhogqr", kdim, ijdim);
    drhogqs  = View2D<double, DEFAULT_MEM>("drhogqs", kdim, ijdim);
    drhogqg  = View2D<double, DEFAULT_MEM>("drhogqg", kdim, ijdim);
    psatl    = View2D<double, DEFAULT_MEM>("psatl", kdim, ijdim);
    qsatl    = View2D<double, DEFAULT_MEM>("qsatl", kdim, ijdim);                
    psati    = View2D<double, DEFAULT_MEM>("psati", kdim, ijdim);
    qsati    = View2D<double, DEFAULT_MEM>("qsati", kdim, ijdim);                
    Nc       = View2D<double, DEFAULT_MEM>("Nc", kdim, ijdim);                 
    a1       = View2D<double, DEFAULT_MEM>("a1", kdim, ijdim);                 //<
    a2       = View2D<double, DEFAULT_MEM>("a2", kdim, ijdim);                 //<
    ma2      = View2D<double, DEFAULT_MEM>("ma2", kdim, ijdim);                //< 1-a2
    Vt       = View3D<double, DEFAULT_MEM>("Vt  ", nqmax, kdim, ijdim);
    cva      = View2D<double, DEFAULT_MEM>("cva ", kdim, ijdim);
    rgs      = View2D<double, DEFAULT_MEM>("rgs ", kdim, ijdim);
    rgsh     = View2D<double, DEFAULT_MEM>("rgsh", kdim, ijdim);
    ml_Pconv = View2D<double, DEFAULT_MEM>("ml_Pconv", kdim, ijdim);
    ml_Pconw = View2D<double, DEFAULT_MEM>("ml_Pconw", kdim, ijdim);
    ml_Pconi = View2D<double, DEFAULT_MEM>("ml_Pconi", kdim, ijdim);

    //---------------------------------------------------------------------------
    UNDEF = CONST_UNDEF;
    EPS   = CONST_EPS  ;
    PI    = CONST_PI   ;
    Rvap  = CONST_Rvap ;
    LHV0  = CONST_LHV0 ;
    LHS0  = CONST_LHS0 ;
    LHF0  = CONST_LHF0 ;
    PRE00 = CONST_Pstd ;

    negative_filter(rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2);

    if(OPT_INDIR) // aerosol indirect effect
    {
        // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        // KOKKOS_LAMBDA(const size_t k, const size_t ij){
        //     Nc(k,ij) = max( UNCCN(k,ij) * 1.0E-6, Nc_def );
        // });
        Kokkos::parallel_for(MDRangePolicy<TagAerosolT, Kokkos::Rank<2>>({0,0},{kdim,ijdim}), *this);
    }
    else
    {
        // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        // KOKKOS_LAMBDA(const size_t k, const size_t ij){
        //     Nc(k,ij) = Nc_def;
        // });
        Kokkos::parallel_for(MDRangePolicy<TagAerosolF, Kokkos::Rank<2>>({0,0},{kdim,ijdim}), *this);
    }

    // saturation water contensts
    SATURATION_psat_liq(tem, psatl);
    SATURATION_psat_ice(tem, psati);

    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    // KOKKOS_LAMBDA(const size_t k, const size_t ij){
    //     qsatl(k,ij) = psatl(k,ij) / ( rho(k,ij) * Rvap * tem(k,ij) );
    //     qsati(k,ij) = psati(k,ij) / ( rho(k,ij) * Rvap * tem(k,ij) );
    // });
    Kokkos::parallel_for(MDRangePolicy<TagSATWater, Kokkos::Rank<2>>({0,0},{kdim,ijdim}), *this);

    // Bergeron process parameters
    Bergeron_param(tem, a1, a2, ma2);

    // work for Effective Radius of Liquid Water
    // dgamma_a is parameter of Gamma-Dist.
    // "Berry and Reinhardt(1974-a) eq.(16-a,b)"
    // 0.38 * {var1/2x} =~ {var1/2 r}
    //         dgamma_a = 1.0_RP /{var x}
    double Dc = 0.146 - 5.964E-2 * std::log(para.Nc_def / 2000.0);
    dgamma_a = 0.1444 / (Dc * Dc); // Dc**2 in fortran
    GAM_dgam   = MISC_gammafunc(dgamma_a);
    GAM_dgam23 = MISC_gammafunc(dgamma_a + 2.0/3.0);
    coef_dgam  = GAM_dgam23 / GAM_dgam * std::pow(dgamma_a, (-2.0/3.0));
    coef_xf    = 3.0 / 4.0 / PI / rho_w;

}

void mp_nsw6_new::compute()
{
    // nsw6 big loop
    Kokkos::parallel_for(MDRangePolicy<TagBigLoop, Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), *this);

    // update rhogq
    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), 
    // KOKKOS_LAMBDA(const size_t k, const size_t ij){
    //     rhogq(I_QV,k,ij) += drhogqv(k,ij);
    //     rhogq(I_QC,k,ij) += drhogqc(k,ij);
    //     rhogq(I_QR,k,ij) += drhogqr(k,ij);
    //     rhogq(I_QI,k,ij) += drhogqi(k,ij);
    //     rhogq(I_QS,k,ij) += drhogqs(k,ij);
    //     rhogq(I_QG,k,ij) += drhogqg(k,ij);
    // });
    Kokkos::parallel_for(MDRangePolicy<TagUpdateRHOGQ, Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), *this);

    // update rhoge
    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), 
    // KOKKOS_LAMBDA(const size_t k, const size_t ij){
    //     rhoge(k,ij) = rhoge(k,ij) - LHV * drhogqv(k,ij)
    //                             + LHF * drhogqi(k,ij)
    //                             + LHF * drhogqs(k,ij)
    //                             + LHF * drhogqg(k,ij);
    // });
    Kokkos::parallel_for(MDRangePolicy<TagUpdateRHOGE, Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), *this);

    // Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
    //     rceff    (kmin-1,ij) = 0.0;
    //     rceff_cld(kmin-1,ij) = 0.0;
    //     rceff    (kmax+1,ij) = 0.0;
    //     rceff_cld(kmax+1,ij) = 0.0;

    //     double sw = ( 0.5 + std::copysign(0.5, tctop(0,ij) - TEM00 ) ); // if T > 0C, sw=1

    //     rwtop(0,ij) =  (       sw ) * rctop(0,ij)
    //                 + ( 1.0 - sw ) * UNDEF;

    //     sw = ( 0.5 + std::copysign(0.5, rctop_cld(0,ij) - TEM00 ) ); // if T > 0C, sw=1

    //     rwtop_cld(0,ij) =  (       sw ) * rctop_cld(0,ij)
    //                     + ( 1.0 - sw ) * UNDEF;
    // });
    Kokkos::parallel_for(RangePolicy<TagUpdateCLD>(0,ijdim), *this);
    

    //---< preciptation（降水量） >---
    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), 
    // KOKKOS_LAMBDA(const size_t nq, const size_t ij){
    //     Vt(nq,kmin-1,ij) = 0.0;
    //     Vt(nq,kmax+1,ij) = 0.0;
    // });
    Kokkos::parallel_for(MDRangePolicy<TagPrecip, Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), *this);
    //--- update mass concentration（質量濃度）
    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<3>>({NQW_STR,kmin,IDX_ZERO},{NQW_END+1,kmax+1,ijdim}), 
    // KOKKOS_LAMBDA(const size_t nq, const size_t k, const size_t ij){
    //     q(nq,k,ij) = rhogq(nq,k,ij) / rhog(k,ij);
    // });
    Kokkos::parallel_for(MDRangePolicy<TagUpdateMassConc, Kokkos::Rank<3>>({NQW_STR,kmin,IDX_ZERO},{NQW_END+1,kmax+1,ijdim}), *this);

    THRMDYN_qd(q, qd);
    THRMDYN_cv(qd, q, cva);

    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), 
    // KOKKOS_LAMBDA(const size_t k, const size_t ij){
    //         tem (k,ij) = rhoge(k,ij) / ( rhog(k,ij) * cva(k,ij) );
    //         rgs (k,ij) = gam2 (k,ij) / gsgam2 (k,ij);
    //         rgsh(k,ij) = gam2h(k,ij) / gsgam2h(k,ij);
    // });
    Kokkos::parallel_for(MDRangePolicy<TagUpdateTEM, Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), *this);

    // for(int nq = 0; nq < nqmax; nq++)
    //     preciptation_flag[nq] = false;
    preciptation_flag[I_QR] = true;
    preciptation_flag[I_QI] = true;
    preciptation_flag[I_QS] = true;
    preciptation_flag[I_QG] = true;
    if(precip_transport_type == "3WATER")
    {
        preciptation_flag[I_QI] = false;
    }

    precip_transport_new( rhog,
                        rhogvx,
                        rhogvy,
                        rhogvz,
                        rhogw,
                        rhoge,
                        rhogq,
                        rho,
                        tem,
                        pre,
                        vx,
                        vy,
                        vz,
                        w,
                        q,
                        qd,
                        z,
                        Vt,
                        preciptation_flag,
                        precip,precip_rhoe,
                        precip_lh_heat,
                        precip_rhophi,
                        precip_rhokin,
                        gprec,
                        gsgam2,
                        gsgam2h,
                        rgs,
                        rgsh,
                        ix,
                        iy,
                        iz,
                        jx,
                        jy,
                        jz,
                        dt,
                        nullptr); // precip_trc**

    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    // KOKKOS_LAMBDA(const size_t k, const size_t ij){
    //         ml_Pconv(k,ij) = q(I_QV,k,ij);
    //         ml_Pconw(k,ij) = q(I_QC,k,ij);
    //         ml_Pconi(k,ij) = q(I_QI,k,ij);
    // });
    Kokkos::parallel_for(MDRangePolicy<TagUpdatePCON, Kokkos::Rank<2>>({0,0},{kdim,ijdim}), *this);

    bool ice_adjust;
    if(OPT_EXPLICIT_ICEGEN)
    {
        ice_adjust = false;
        SATURATION_Setrange(100.0, 90.0);

        SATURATION_adjustment(rhog, rhoge, rhogq, tem, q, qd, gsgam2, ice_adjust);
    }
    else
    {
        ice_adjust = true;
        SATURATION_Setrange(273.16, 233.16);

        SATURATION_adjustment(rhog, rhoge, rhogq, tem , q, qd, gsgam2, ice_adjust);
    }

    // Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    // KOKKOS_LAMBDA(const size_t k, const size_t ij){
    //     ml_Pconv(k,ij) = q(I_QV,k,ij) - ml_Pconv(k,ij);
    //     ml_Pconw(k,ij) = q(I_QC,k,ij) - ml_Pconw(k,ij);
    //     ml_Pconi(k,ij) = q(I_QI,k,ij) - ml_Pconi(k,ij);
    // });
    Kokkos::parallel_for(MDRangePolicy<TagUpdatePCON2, Kokkos::Rank<2>>({0,0},{kdim,ijdim}), *this);

    negative_filter(rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2);

    // End of nsw6
}

void mp_nsw6_new::operator()(TagAerosolT, const size_t k, const size_t ij) const
{
    Nc(k,ij) = Kokkos::max( UNCCN(k,ij) * 1.0E-6, paras.Nc_def );
}

void mp_nsw6_new::operator()(TagAerosolF, const size_t k, const size_t ij) const
{
    Nc(k,ij) = paras.Nc_def;
}

void mp_nsw6_new::operator()(TagSATWater, const size_t k, const size_t ij) const
{
    qsatl(k,ij) = psatl(k,ij) / ( rho(k,ij) * Rvap * tem(k,ij) );
    qsati(k,ij) = psati(k,ij) / ( rho(k,ij) * Rvap * tem(k,ij) );
}

void mp_nsw6_new::operator()(TagUpdateRHOGQ, const size_t k, const size_t ij) const
{
    rhogq(I_QV,k,ij) += drhogqv(k,ij);
    rhogq(I_QC,k,ij) += drhogqc(k,ij);
    rhogq(I_QR,k,ij) += drhogqr(k,ij);
    rhogq(I_QI,k,ij) += drhogqi(k,ij);
    rhogq(I_QS,k,ij) += drhogqs(k,ij);
    rhogq(I_QG,k,ij) += drhogqg(k,ij);
}

void mp_nsw6_new::operator()(TagUpdateRHOGE, const size_t k, const size_t ij) const
{
    rhoge(k,ij) = rhoge(k,ij) - LHV * drhogqv(k,ij)
                            + LHF * drhogqi(k,ij)
                            + LHF * drhogqs(k,ij)
                            + LHF * drhogqg(k,ij);
}

void mp_nsw6_new::operator()(TagUpdateCLD, const size_t ij) const
{
    rceff    (kmin-1,ij) = 0.0;
    rceff_cld(kmin-1,ij) = 0.0;
    rceff    (kmax+1,ij) = 0.0;
    rceff_cld(kmax+1,ij) = 0.0;

    double sw = ( 0.5 + Kokkos::copysign(0.5, tctop(0,ij) - TEM00 ) ); // if T > 0C, sw=1

    rwtop(0,ij) =  (       sw ) * rctop(0,ij)
                + ( 1.0 - sw ) * UNDEF;

    sw = ( 0.5 + Kokkos::copysign(0.5, rctop_cld(0,ij) - TEM00 ) ); // if T > 0C, sw=1

    rwtop_cld(0,ij) =  (       sw ) * rctop_cld(0,ij)
                    + ( 1.0 - sw ) * UNDEF;
}

void mp_nsw6_new::operator()(TagPrecip, const size_t nq, const size_t ij) const
{
    Vt(nq,kmin-1,ij) = 0.0;
    Vt(nq,kmax+1,ij) = 0.0;
}

void mp_nsw6_new::operator()(TagUpdateMassConc, const size_t nq, const size_t k, const size_t ij) const
{
    q(nq,k,ij) = rhogq(nq,k,ij) / rhog(k,ij);
}

void mp_nsw6_new::operator()(TagUpdateTEM, const size_t k, const size_t ij) const
{
    tem (k,ij) = rhoge(k,ij) / ( rhog(k,ij) * cva(k,ij) );
    rgs (k,ij) = gam2 (k,ij) / gsgam2 (k,ij);
    rgsh(k,ij) = gam2h(k,ij) / gsgam2h(k,ij);
}

void mp_nsw6_new::operator()(TagUpdatePCON, const size_t k, const size_t ij) const
{
    ml_Pconv(k,ij) = q(I_QV,k,ij);
    ml_Pconw(k,ij) = q(I_QC,k,ij);
    ml_Pconi(k,ij) = q(I_QI,k,ij);
}

void mp_nsw6_new::operator()(TagUpdatePCON2, const size_t k, const size_t ij) const
{
    ml_Pconv(k,ij) = q(I_QV,k,ij) - ml_Pconv(k,ij);
    ml_Pconw(k,ij) = q(I_QC,k,ij) - ml_Pconw(k,ij);
    ml_Pconi(k,ij) = q(I_QI,k,ij) - ml_Pconi(k,ij);
}

void mp_nsw6_new::operator()(TagBigLoop, const size_t k, const size_t ij) const
{
    // define private variables
    double dens    ;                         //< density
    double temp    ;                         //< T [K]
    double qv      ;                         //< mixing ratio of water vapor  [kg/kg]
    double qc      ;                         //< mixing ratio of liquid water [kg/kg]
    double qr      ;                         //< mixing ratio of rain         [kg/kg]
    double qi      ;                         //< mixing ratio of ice water    [kg/kg]
    double qs      ;                         //< mixing ratio of snow         [kg/kg]
    double qg      ;                         //< mixing ratio of graupel      [kg/kg]
    double qv_t    ;                         //< tendency     of water vapor  [kg/kg/s]
    double qc_t    ;                         //< tendency     of liquid water [kg/kg/s]
    double qr_t    ;                         //< tendency     of rain         [kg/kg/s]
    double qi_t    ;                         //< tendency     of ice water    [kg/kg/s]
    double qs_t    ;                         //< tendency     of snow         [kg/kg/s]
    double qg_t    ;                         //< tendency     of graupel      [kg/kg/s]
    double Sliq    ;                         //< saturated ratio S for liquid water [0-1]
    double Sice    ;                         //< saturated ratio S for ice water    [0-1]
    double Rdens   ;                         //< 1 / density
    double rho_fact;                         //< density factor
    double temc    ;                         //< T - T0 [K]

    double RLMDr, RLMDr_2, RLMDr_3        ;
    double RLMDs, RLMDs_2, RLMDs_3        ;
    double RLMDg, RLMDg_2, RLMDg_3        ;
    double RLMDr_1br, RLMDr_2br, RLMDr_3br;
    double RLMDs_1bs, RLMDs_2bs, RLMDs_3bs;
    double RLMDr_dr, RLMDr_3dr, RLMDr_5dr ;
    double RLMDs_ds, RLMDs_3ds, RLMDs_5ds ;
    double RLMDg_dg, RLMDg_3dg, RLMDg_5dg ;
    double RLMDr_7                        ;
    double RLMDr_6dr                      ;

    // ---< Roh and Satoh (2014) >---
    double tems, Xs2                     ;
    double MOMs_0, MOMs_1, MOMs_2        ;
    double MOMs_0bs, MOMs_1bs, MOMs_2bs  ;
    double MOMs_2ds, MOMs_5ds_h, RMOMs_Vt;
    double coef_at[4]                    ;
    double coef_bt[4]                    ;
    double loga_, b_, nm                 ;

    double Vti, Vtr, Vts, Vtg   ;           //< terminal velocity
    double Esi_mod, Egs_mod     ;           //< modified accretion efficiency
    double rhoqc                ;           //< rho * qc
    double Pracw_orig,  Pracw_kk;           //< accretion       term by orig  & k-k scheme
    double Praut_berry, Praut_kk;           //< auto-conversion term by berry & k-k scheme
    double Dc                   ;           //< relative variance
    double betai, betas         ;           //< sticky parameter for auto-conversion
    double Ka                   ;           //< thermal diffusion coefficient of air
    double Kd                   ;           //< diffusion coefficient of water vapor in air
    double Nu                   ;           //< kinematic viscosity of air
    double Glv, Giv, Gil        ;           //< thermodynamic function
    double ventr, vents, ventg  ;           //< ventilation factor
    double net, fac, fac_sw     ;
    double zerosw, tmp          ;
    double wk[wk_nmax]          ;
    //---< Bergeron process >---
    double sw_bergeron     ;                 //< if 0C<T<30C, sw=1
    double dt1             ;                 //< time during which the an ice particle of 40um grows to 50um
    double Ni50            ;                 //< number concentration of ice particle of 50um
    //---< Explicit ice generation >---
    double sw, rhoqi, XNi, XMi, Di, Ni0, Qi0;
    //---< Effective radius >---
    double xf_qc;                 //< mean mass   of qc [kg]
    double rf_qc;                 //< mean radius of qc [m]
    double r2_qc;                 //< r^2 moment  of qc
    double r2_qr;                 //< r^2 moment  of qr
    double r3_qc;                 //< r^3 moment  of qc
    double r3_qr;                 //< r^3 moment  of qr

    // ---- < create alias from global params > ----
    const double Ar = paras.Ar;
    const double As = paras.As;
    const double Ag = paras.Ag;
    const double Cr = paras.Cr;
    const double Cs = paras.Cs;
    const double Cg = paras.Cg;
    const double N0r = paras.N0r;
    const double N0s = paras.N0s;
    const double N0g = paras.N0g;

    const double GAM     = paras.GAM  ;
    const double GAM_2   = paras.GAM_2;
    const double GAM_3   = paras.GAM_3;

    const double GAM_1br = paras.GAM_1br;
    const double GAM_2br = paras.GAM_2br;
    const double GAM_3br = paras.GAM_3br;
    const double GAM_3dr = paras.GAM_3dr;
    const double GAM_6dr = paras.GAM_6dr;
    const double GAM_1brdr = paras.GAM_1brdr;
    const double GAM_5dr_h = paras.GAM_5dr_h;

    const double GAM_1bs = paras.GAM_1bs;
    const double GAM_2bs = paras.GAM_2bs;
    const double GAM_3bs = paras.GAM_3bs;
    const double GAM_3ds = paras.GAM_3ds;
    const double GAM_1bsds = paras.GAM_1bsds;
    const double GAM_5ds_h = paras.GAM_5ds_h;

    const double GAM_1bg = paras.GAM_1bg;
    const double GAM_3dg   = paras.GAM_3dg;
    const double GAM_1bgdg = paras.GAM_1bgdg;
    const double GAM_5dg_h = paras.GAM_5dg_h;

    const double Eiw        = paras.Eiw       ; 
    const double Erw        = paras.Erw       ; 
    const double Esw        = paras.Esw       ; 
    const double Egw        = paras.Egw       ; 
    const double Eri        = paras.Eri       ; 
    const double Esi        = paras.Esi       ; 
    const double Egi        = paras.Egi       ; 
    const double Esr        = paras.Esr       ; 
    const double Egr        = paras.Egr       ; 
    const double Egs        = paras.Egs       ; 
    const double gamma_sacr = paras.gamma_sacr; 
    const double gamma_gacs = paras.gamma_gacs; 
    const double mi         = paras.mi        ; 

    const double beta_saut  = paras.beta_saut ;  
    const double gamma_saut = paras.gamma_saut;  
    const double beta_gaut  = paras.beta_gaut ;  
    const double gamma_gaut = paras.gamma_gaut;  
    const double qicrt_saut = paras.qicrt_saut;  
    const double qscrt_gaut = paras.qscrt_gaut;  

    const double f1r        = paras.f1r;     
    const double f2r        = paras.f2r;     
    const double f1s        = paras.f1s;     
    const double f2s        = paras.f2s;     
    const double f1g        = paras.f1g;     
    const double f2g        = paras.f2g;     

    const double A_frz      = paras.A_frz;
    const double B_frz      = paras.B_frz;

    const double mi40       = paras.mi40 ; 
    const double mi50       = paras.mi50 ; 
    const double vti50      = paras.vti50; 
    const double Ri50       = paras.Ri50 ; 

    const double sw_expice  = paras.sw_expice;
    const double Nc_ihtr    = paras.Nc_ihtr  ;
    const double Di_max     = paras.Di_max;
    const double Di_a       = paras.Di_a  ;

    const double sw_kk2000 = paras.sw_kk2000;
    const double sw_constVti = paras.sw_constVti;
    const double CONST_Vti = paras.CONST_Vti;

    const double sw_roh2014 = paras.sw_roh2014;
    const double ln10 = paras.ln10;

    const double dens00 = ::dens00;
    // ---- < END create alias from global params > ----

    dens = rho(k,ij);
    temp = tem(k,ij);
    qv   = Kokkos::max( q(I_QV, k, ij), 0.0 );
    qc   = Kokkos::max( q(I_QC, k, ij), 0.0 );
    qr   = Kokkos::max( q(I_QR, k, ij), 0.0 );
    qi   = Kokkos::max( q(I_QI, k, ij), 0.0 );
    qs   = Kokkos::max( q(I_QS, k, ij), 0.0 );
    qg   = Kokkos::max( q(I_QG, k, ij), 0.0 );

    // saturation ration S
    Sliq = qv / (Kokkos::max(qsatl(k,ij), EPS));
    Sice = qv / (Kokkos::max(qsati(k,ij), EPS));

    Rdens = 1.0 / dens;
    rho_fact = Kokkos::sqrt(dens00 * Rdens);
    temc = temp - TEM00;

    wk[I_delta1] = ( 0.5 + Kokkos::copysign(0.5, qr - 1.0E-4) );

    wk[I_delta2] = ( 0.5 + Kokkos::copysign(0.5, 1.0E-4 - qr) )
                    * ( 0.5 + Kokkos::copysign(0.5, 1.0E-4 - qs) );
    
    wk[I_spsati] = 0.5 + Kokkos::copysign(0.5, Sice - 1.0);

    wk[I_iceflg] = 0.5 - Kokkos::copysign(0.5, temc); // 0: warm, 1: ice

    wk[I_dqv_dt] = qv / dt;
    wk[I_dqc_dt] = qc / dt;
    wk[I_dqr_dt] = qr / dt;
    wk[I_dqi_dt] = qi / dt;
    wk[I_dqs_dt] = qs / dt;
    wk[I_dqg_dt] = qg / dt;

    sw_bergeron = ( 0.5 + Kokkos::copysign(0.5, temc + 30.0) )
                    * ( 0.5 + Kokkos::copysign(0.5, 0.0 - temc) )
                    * ( 1.0 - sw_expice );
    
    // slope parameter lambda (Rain)
    zerosw = 0.5 - Kokkos::copysign(0.5, qr - 1.0E-12);
    RLMDr = Kokkos::sqrt( Kokkos::sqrt( dens * qr / (Ar * N0r * GAM_1br) + zerosw ) ) * (1.0 - zerosw);

    RLMDr_dr  = Kokkos::sqrt(RLMDr);   // **Dr
    RLMDr_2   = Kokkos::pow(RLMDr, 2);
    RLMDr_3   = Kokkos::pow(RLMDr, 3);
    RLMDr_7   = Kokkos::pow(RLMDr, 7);
    RLMDr_1br = Kokkos::pow(RLMDr, 4); // (1+Br)
    RLMDr_2br = Kokkos::pow(RLMDr, 5); // (2+Br)
    RLMDr_3br = Kokkos::pow(RLMDr, 6); // (3+Br)
    RLMDr_3dr = Kokkos::pow(RLMDr, 3) * RLMDr_dr;
    RLMDr_5dr = Kokkos::pow(RLMDr, 5) * RLMDr_dr;
    RLMDr_6dr = Kokkos::pow(RLMDr, 6) * RLMDr_dr;

    // slope parameter lambda (Snow)
    zerosw = 0.5 - Kokkos::copysign(0.5, qs - 1.0E-12);
    RLMDs  = Kokkos::sqrt( Kokkos::sqrt( dens * qs / ( As * N0s * GAM_1bs ) + zerosw ) ) * ( 1.0 - zerosw );

    RLMDs_ds  = Kokkos::sqrt( Kokkos::sqrt(RLMDs) ); // **Ds
    RLMDs_2   = Kokkos::pow(RLMDs, 2);
    RLMDs_3   = Kokkos::pow(RLMDs, 3);
    RLMDs_1bs = Kokkos::pow(RLMDs, 4); // (1+Bs)
    RLMDs_2bs = Kokkos::pow(RLMDs, 5); // (2+Bs)
    RLMDs_3bs = Kokkos::pow(RLMDs, 6); // (3+Bs)
    RLMDs_3ds = Kokkos::pow(RLMDs, 3) * RLMDs_ds;
    RLMDs_5ds = Kokkos::pow(RLMDs, 5) * RLMDs_ds;

    MOMs_0     = N0s * GAM       * RLMDs    ;       // Ns * 0th moment
    MOMs_1     = N0s * GAM_2     * RLMDs_2  ;       // Ns * 1st moment
    MOMs_2     = N0s * GAM_3     * RLMDs_3  ;       // Ns * 2nd moment
    MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs;       // Ns * 0+bs
    MOMs_1bs   = N0s * GAM_2bs   * RLMDs_2bs;       // Ns * 1+bs
    MOMs_2bs   = N0s * GAM_3bs   * RLMDs_3bs;       // Ns * 2+bs
    MOMs_2ds   = N0s * GAM_3ds   * RLMDs_3ds;       // Ns * 2+ds
    MOMs_5ds_h = N0s * GAM_5ds_h * Kokkos::sqrt(RLMDs_5ds); // Ns * (5+ds)/2
    RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds;

    //---< modification by Roh and Satoh (2014) >---

    // bimodal size distribution of snow
    Xs2 = dens * qs / As;
    zerosw = 0.5 - Kokkos::copysign(0.5, Xs2 - 1.0E-12);

    tems = Kokkos::min(-0.1, temc);
    coef_at[0] = paras.coef_a[0] + tems * ( paras.coef_a[1] + tems * ( paras.coef_a[4] + tems * paras.coef_a[8] ) );
    coef_at[1] = paras.coef_a[2] + tems * ( paras.coef_a[3] + tems *   paras.coef_a[6] );
    coef_at[2] = paras.coef_a[5] + tems *   paras.coef_a[7];
    coef_at[3] = paras.coef_a[9];
    coef_bt[0] = paras.coef_b[0] + tems * ( paras.coef_b[1] + tems * ( paras.coef_b[4] + tems * paras.coef_b[8] ) );
    coef_bt[1] = paras.coef_b[2] + tems * ( paras.coef_b[3] + tems *   paras.coef_b[6] );
    coef_bt[2] = paras.coef_b[5] + tems *   paras.coef_b[7];
    coef_bt[3] = paras.coef_b[9];

    // 0th moment
    loga_ = coef_at[0];
    b_    = coef_bt[0];
    MOMs_0 = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
                + ( 1.0 - sw_roh2014 ) * MOMs_0;
    // 1st moment
    nm = 1.0;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    MOMs_1 = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw ) * b_) * ( 1.0 - zerosw )
                + ( 1.0 - sw_roh2014 ) * MOMs_1;
    // 2nd moment
    MOMs_2 = sw_roh2014 * Xs2 + (1.0 - sw_roh2014) * MOMs_2;
    // 0 + Bs(=2) moment
    nm = 2.0;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    MOMs_0bs = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
                + ( 1.0 - sw_roh2014 ) * MOMs_0bs;
    // 1 + Bs(=2) moment
    nm = 3.0;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    MOMs_1bs = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
                + ( 1.0 - sw_roh2014 ) * MOMs_1bs;
    // 2 + Bs(=2) moment
    nm = 4.0;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    MOMs_2bs = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
                + ( 1.0 - sw_roh2014 ) * MOMs_2bs;
    // 2 + Ds(=0.25) moment
    nm = 2.25;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    MOMs_2ds = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
                + ( 1.0 - sw_roh2014 ) * MOMs_2ds;
    // ( 3 + Ds(=0.25) ) / 2  moment
    nm = 1.625;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    MOMs_5ds_h = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
                    + ( 1.0 - sw_roh2014 ) * MOMs_5ds_h;
    // Bs(=2) + Ds(=0.25) moment
    nm = 2.25;
    loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
    b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
    RMOMs_Vt = sw_roh2014 * Kokkos::exp(ln10 * loga_) * Kokkos::exp(Kokkos::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw ) / (MOMs_0bs + zerosw)
                + ( 1.0 - sw_roh2014 ) * RMOMs_Vt;

    // slope parameter lambda (Graupel)
    zerosw = 0.5 - Kokkos::copysign(0.5, qg - 1.0E-12);
    RLMDg  = Kokkos::sqrt( Kokkos::sqrt( dens * qg / ( Ag * N0g * GAM_1bg ) + zerosw ) ) * ( 1.0 - zerosw );
    RLMDg_dg  = Kokkos::sqrt( RLMDg );       // **Dg
    RLMDg_2   = Kokkos::pow(RLMDg, 2);
    RLMDg_3   = Kokkos::pow(RLMDg, 3);
    RLMDg_3dg = Kokkos::pow(RLMDg, 3) * RLMDg_dg;
    RLMDg_5dg = Kokkos::pow(RLMDg, 5) * RLMDg_dg;

    wk[I_RLMDr] = RLMDr;
    wk[I_RLMDs] = RLMDs;
    wk[I_RLMDg] = RLMDg;

    //---< terminal velocity >---
    zerosw = 0.5 - Kokkos::copysign(0.5, qi - 1.0E-8 );
    Vti = ( 1.0 - sw_constVti ) * (-3.29) * Kokkos::exp( Kokkos::log( dens * qi + zerosw ) * 0.16 ) * ( 1.0 - zerosw ) 
            + ( sw_constVti ) * (-CONST_Vti);
    Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr;
    Vts = -Cs * rho_fact * RMOMs_Vt;
    Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg;

    //---< Nucleation >---
    // [Pigen] ice nucleation
    Ni0 = Kokkos::max( Kokkos::exp(-0.1 * temc), 1.0 ) * 1000.0;
    Qi0 = 4.92E-11 * Kokkos::exp( Kokkos::log(Ni0) * 1.33 ) * Rdens;

    wk[I_Pigen] = Kokkos::max( Kokkos::min( Qi0 - qi, qv - qsati(k,ij) ), 0.0 ) / dt;

    //---< Accretion >---
    Esi_mod = Kokkos::min( Esi, Esi * Kokkos::exp( gamma_sacr * temc ) );
    Egs_mod = Kokkos::min( Egs, Egs * Kokkos::exp( gamma_gacs * temc ) );

    // [Pracw] accretion rate of cloud water by rain
    Pracw_orig = qc * 0.25 * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact;

    zerosw     = 0.5 - Kokkos::copysign(0.5, qc * qr - 1.0E-12 );
    Pracw_kk   = 67.0 * Kokkos::exp( Kokkos::log( qc * qr + zerosw ) * 1.15 ) * ( 1.0 - zerosw ); // eq.(33) in KK(2000)

    // switch orig / k-k scheme
    wk[I_Pracw] = ( 1.0 - sw_kk2000 ) * Pracw_orig
                + (       sw_kk2000 ) * Pracw_kk;

    // [Psacw] accretion rate of cloud water by snow
    wk[I_Psacw] = qc * 0.25 * PI * Esw * Cs * MOMs_2ds * rho_fact;

    // [Pgacw] accretion rate of cloud water by graupel
    wk[I_Pgacw] = qc * 0.25 * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact;

    // [Praci] accretion rate of cloud ice by rain
    wk[I_Praci] = qi * 0.25 * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact;

    // [Psaci] accretion rate of cloud ice by snow
    wk[I_Psaci] = qi * 0.25 * PI * Esi_mod * Cs * MOMs_2ds * rho_fact;

    // [Pgaci] accretion rate of cloud ice by grupel
    wk[I_Pgaci] = qi * 0.25 * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact;

    // [Piacr] accretion rate of rain by cloud ice
    wk[I_Piacr] = qi * Ar / mi * 0.25 * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact;

    // [Psacr] accretion rate of rain by snow
    wk[I_Psacr] = Ar * 0.25 * PI * Rdens * Esr * N0r * Kokkos::abs(Vtr - Vts)
                        * ( GAM_1br * RLMDr_1br * MOMs_2          
                            + 2.0 * GAM_2br * RLMDr_2br * MOMs_1          
                            + GAM_3br * RLMDr_3br * MOMs_0 );

    // [Pgacr] accretion rate of rain by graupel
    wk[I_Pgacr] = Ar * 0.25 * PI * Rdens * Egr * N0g * N0r * Kokkos::abs(Vtg - Vtr)
                        * ( GAM_1br * RLMDr_1br * GAM_3 * RLMDg_3
                            + 2.0 * GAM_2br * RLMDr_2br * GAM_2 * RLMDg_2
                            + GAM_3br * RLMDr_3br * GAM * RLMDg );

    // [Pracs] accretion rate of snow by rain
    wk[I_Pracs] = As * 0.25 * PI * Rdens * Esr *  N0r * Kokkos::abs(Vtr - Vts)
                        * ( MOMs_0bs * GAM_3 * RLMDr_3
                            + 2.0 * MOMs_1bs * GAM_2 * RLMDr_2
                            + MOMs_2bs * GAM * RLMDr );

    // [Pgacs] accretion rate of snow by graupel
    wk[I_Pgacs] = As * 0.25 * PI * Rdens * Egs_mod * N0g * Kokkos::abs(Vtg - Vts)
                        * ( MOMs_0bs * GAM_3 * RLMDg_3
                            + 2.0 * MOMs_1bs * GAM_2 * RLMDg_2
                            + MOMs_2bs * GAM * RLMDg );

    //---< Autoconversion >---            
    // [Praut] auto-conversion rate from cloud water to rain
    rhoqc = dens * qc * 1000.0; // [g/m3]
    Dc    = 0.146 - 5.964E-2 * Kokkos::log( Nc(k,ij) / 2000.0 );
    Praut_berry = Rdens * 1.67E-5 * rhoqc * rhoqc / ( 5.0 + 3.66E-2 * Nc(k,ij) / ( Dc * rhoqc + EPS ) );

    zerosw   = 0.5 - Kokkos::copysign(0.5, qc - 1.0E-12 );
    Praut_kk = 1350.0                                           
                * Kokkos::exp( Kokkos::log( qc + zerosw ) * 2.47 ) * ( 1.0 - zerosw ) 
                * Kokkos::exp( Kokkos::log( Nc(k,ij) ) * (-1.79) );                     // eq.(29) in KK(2000)

    Praut_kk = 1350.0 * Kokkos::pow(qc, 2.47) * Kokkos::pow(Nc(k,ij), -1.79);


    // switch berry / k-k scheme
    wk[I_Praut] = ( 1.0 - sw_kk2000 ) * Praut_berry
                    + sw_kk2000 * Praut_kk;

    // [Psaut] auto-conversion rate from cloud ice to snow
    betai = Kokkos::min( beta_saut, beta_saut * Kokkos::exp( gamma_saut * temc ) );
    wk[I_Psaut] = Kokkos::max( betai * (qi - qicrt_saut), 0.0 );

    // [Pgaut] auto-conversion rate from snow to graupel
    betas = Kokkos::min( beta_gaut, beta_gaut * Kokkos::exp( gamma_gaut * temc ) );
    wk[I_Pgaut] = Kokkos::max( betas * (qs - qscrt_gaut), 0.0 );

    //---< Evaporation, Sublimation, Melting, and Freezing >---
    Ka  = ( Ka0 + dKa_dT * temc );
    Kd  = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(k,ij);
    Nu  = ( nu0 + dnu_dT * temc ) * Rdens;

    Glv = 1.0 / ( LHV0/(Ka * temp) * ( LHV0/(Rvap * temp) - 1.0 ) + 1.0 / (Kd * dens * qsatl(k,ij)) );
    Giv = 1.0 / ( LHS0/(Ka * temp) * ( LHS0/(Rvap * temp) - 1.0 ) + 1.0 / (Kd * dens * qsati(k,ij)) );
    Gil = 1.0 / ( LHF0/(Ka * temc) );

    // [Prevp] evaporation rate of rain
    ventr = f1r * GAM_2 * RLMDr_2 + f2r * Kokkos::sqrt( Cr * rho_fact / Nu * RLMDr_5dr ) * GAM_5dr_h;

    wk[I_Prevp] = 2.0 * PI * Rdens * N0r * ( 1.0 - Kokkos::min(Sliq, 1.0) ) * Glv * ventr;

    // [Pidep,Pisub] deposition/sublimation rate for ice
    rhoqi = Kokkos::max(dens * qi, EPS);
    XNi   = Kokkos::min( Kokkos::max( 5.38E+7 * Kokkos::exp( Kokkos::log(rhoqi) * 0.75 ), 1.0E+3 ), 1.0E+6 );
    XMi   = rhoqi / XNi;
    Di    = Kokkos::min( Di_a * Kokkos::sqrt(XMi), Di_max );

    tmp = 4.0 * Di * XNi * Rdens * ( Sice - 1.0 ) * Giv;

    wk[I_Pidep] = (       wk[I_spsati] ) * ( tmp); // Sice > 2
    wk[I_Pisub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1

    // [Pihom] homogenious freezing at T < -40C
    sw = ( 0.5 - Kokkos::copysign(0.5, temc + 40.0) ); // if T < -40C, sw=1

    wk[I_Pihom] = sw * qc / dt;

    // [Pihtr] heteroginous freezing at -40C < T < 0C
    sw =   ( 0.5 + Kokkos::copysign(0.5, temc + 40.0) )
         * ( 0.5 - Kokkos::copysign(0.5, temc       ) ); // if -40C < T < 0C, sw=1

    wk[I_Pihtr] = sw * ( dens / rho_w * (qc*qc) / ( Nc_ihtr * 1.0E+6 ) )
                     * B_frz * ( Kokkos::exp(-A_frz * temc) - 1.0 );

    // [Pimlt] ice melting at T > 0C
    sw = ( 0.5 + Kokkos::copysign(0.5, temc) ); // if T > 0C, sw=1

    wk[I_Pimlt] = sw * qi / dt;

    // [Psdep,Pssub] deposition/sublimation rate for snow
    vents = f1s * MOMs_1 + f2s * Kokkos::sqrt( Cs * rho_fact / Nu ) * MOMs_5ds_h;

    tmp = 2.0 * PI * Rdens * ( Sice - 1.0 ) * Giv * vents;

    wk[I_Psdep] = (       wk[I_spsati] ) * ( tmp); // Sice > 1
    wk[I_Pssub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1            

    // [Psmlt] melting rate of snow
    wk[I_Psmlt] = 2.0 * PI * Rdens * Gil * vents
                    + CL * temc / LHF0 * ( wk[I_Psacw] + wk[I_Psacr] );
    wk[I_Psmlt] = Kokkos::max( wk[I_Psmlt], 0.0 );

    // [Pgdep/pgsub] deposition/sublimation rate for graupel
    ventg = f1g * GAM_2 * RLMDg_2 + f2g * Kokkos::sqrt( Cg * rho_fact / Nu * RLMDg_5dg ) * GAM_5dg_h;

    tmp = 2.0 * PI * Rdens * N0g * ( Sice - 1.0 ) * Giv * ventg;

    wk[I_Pgdep] = (       wk[I_spsati] ) * ( tmp); // Sice > 1
    wk[I_Pgsub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1

    // [Pgmlt] melting rate of graupel
    wk[I_Pgmlt] = 2.0 * PI * Rdens * N0g * Gil * ventg
                    + CL * temc / LHF0 * ( wk[I_Pgacw] + wk[I_Pgacr] );
    wk[I_Pgmlt] = Kokkos::max( wk[I_Pgmlt], 0.0 );

    // [Pgfrz] freezing rate of graupel
    wk[I_Pgfrz] = 2.0 * PI * Rdens * N0r * 60.0 * B_frz * Ar * ( Kokkos::exp(-A_frz * temc) - 1.0 ) * RLMDr_7;

    // [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
    dt1  = ( Kokkos::exp( Kokkos::log(mi50) * ma2(k,ij) )
           - Kokkos::exp( Kokkos::log(mi40) * ma2(k,ij) ) ) / ( a1(k,ij) * ma2(k,ij) );
    Ni50 = qi * dt / ( mi50 * dt1 );

    wk[I_Psfw] = Ni50 * ( a1(k,ij) * Kokkos::pow(mi50, a2(k,ij))
                            + PI * Eiw * dens * qc * Ri50 * Ri50 * vti50 );
    wk[I_Psfi] = qi / dt1;

    //---< limiter >---
    wk[I_Pigen] = Kokkos::min( wk[I_Pigen], wk[I_dqv_dt] ) * ( wk[I_iceflg] ) * sw_expice;
    wk[I_Pidep] = Kokkos::min( wk[I_Pidep], wk[I_dqv_dt] ) * ( wk[I_iceflg] ) * sw_expice;
    wk[I_Psdep] = Kokkos::min( wk[I_Psdep], wk[I_dqv_dt] ) * ( wk[I_iceflg] )            ;
    wk[I_Pgdep] = Kokkos::min( wk[I_Pgdep], wk[I_dqv_dt] ) * ( wk[I_iceflg] )            ;

    wk[I_Pracw] = wk[I_Pracw]                           
                + wk[I_Psacw] * ( 1.0 - wk[I_iceflg] )   // c->r by s
                + wk[I_Pgacw] * ( 1.0 - wk[I_iceflg] );  // c->r by g

    wk[I_Praut] = Kokkos::min( wk[I_Praut], wk[I_dqc_dt] ); 
    wk[I_Pracw] = Kokkos::min( wk[I_Pracw], wk[I_dqc_dt] ); 
    wk[I_Pihom] = Kokkos::min( wk[I_Pihom], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_expice;
    wk[I_Pihtr] = Kokkos::min( wk[I_Pihtr], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_expice;
    wk[I_Psacw] = Kokkos::min( wk[I_Psacw], wk[I_dqc_dt] ) * (        wk[I_iceflg] )            ;
    wk[I_Psfw ] = Kokkos::min( wk[I_Psfw ], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_bergeron;
    wk[I_Pgacw] = Kokkos::min( wk[I_Pgacw], wk[I_dqc_dt] ) * (        wk[I_iceflg] )            ;

    wk[I_Prevp] = Kokkos::min( wk[I_Prevp], wk[I_dqr_dt] );
    wk[I_Piacr] = Kokkos::min( wk[I_Piacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
    wk[I_Psacr] = Kokkos::min( wk[I_Psacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
    wk[I_Pgacr] = Kokkos::min( wk[I_Pgacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
    wk[I_Pgfrz] = Kokkos::min( wk[I_Pgfrz], wk[I_dqr_dt] ) * (        wk[I_iceflg] );

    wk[I_Pisub] = Kokkos::min( wk[I_Pisub], wk[I_dqi_dt] ) * (       wk[I_iceflg] ) * sw_expice;
    wk[I_Pimlt] = Kokkos::min( wk[I_Pimlt], wk[I_dqi_dt] ) * ( 1.0 - wk[I_iceflg] ) * sw_expice;
    wk[I_Psaut] = Kokkos::min( wk[I_Psaut], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
    wk[I_Praci] = Kokkos::min( wk[I_Praci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
    wk[I_Psaci] = Kokkos::min( wk[I_Psaci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
    wk[I_Psfi ] = Kokkos::min( wk[I_Psfi ], wk[I_dqi_dt] ) * (       wk[I_iceflg] ) * sw_bergeron;
    wk[I_Pgaci] = Kokkos::min( wk[I_Pgaci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;

    wk[I_Pssub] = Kokkos::min( wk[I_Pssub], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
    wk[I_Psmlt] = Kokkos::min( wk[I_Psmlt], wk[I_dqs_dt] ) * ( 1.0 - wk[I_iceflg] );
    wk[I_Pgaut] = Kokkos::min( wk[I_Pgaut], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
    wk[I_Pracs] = Kokkos::min( wk[I_Pracs], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
    wk[I_Pgacs] = Kokkos::min( wk[I_Pgacs], wk[I_dqs_dt] );

    wk[I_Pgsub] = Kokkos::min( wk[I_Pgsub], wk[I_dqg_dt] ) * (       wk[I_iceflg] );
    wk[I_Pgmlt] = Kokkos::min( wk[I_Pgmlt], wk[I_dqg_dt] ) * ( 1.0 - wk[I_iceflg] );

    wk[I_Piacr_s] = ( 1.0 - wk[I_delta1] ) * wk[I_Piacr];
    wk[I_Piacr_g] = (       wk[I_delta1] ) * wk[I_Piacr];
    wk[I_Praci_s] = ( 1.0 - wk[I_delta1] ) * wk[I_Praci];
    wk[I_Praci_g] = (       wk[I_delta1] ) * wk[I_Praci];
    wk[I_Psacr_s] = (       wk[I_delta2] ) * wk[I_Psacr];
    wk[I_Psacr_g] = ( 1.0 - wk[I_delta2] ) * wk[I_Psacr];
    wk[I_Pracs  ] = ( 1.0 - wk[I_delta2] ) * wk[I_Pracs];

    // [QC]
    net = 
        + wk[I_Pimlt]  // [prod] i->c
        - wk[I_Praut]  // [loss] c->r
        - wk[I_Pracw]  // [loss] c->r
        - wk[I_Pihom]  // [loss] c->i
        - wk[I_Pihtr]  // [loss] c->i
        - wk[I_Psacw]  // [loss] c->s
        - wk[I_Psfw ]  // [loss] c->s
        - wk[I_Pgacw]; // [loss] c->g

    fac_sw = 0.5 + Kokkos::copysign( 0.5, net + EPS ); // if production > loss , fac_sw=1
    fac    = fac_sw + ( 1.0 - fac_sw ) 
                * Kokkos::min( -wk[I_dqc_dt]/(net - fac_sw), 1.0 ); // loss limiter

    wk[I_Pimlt] = wk[I_Pimlt] * fac;
    wk[I_Praut] = wk[I_Praut] * fac;
    wk[I_Pracw] = wk[I_Pracw] * fac;
    wk[I_Pihom] = wk[I_Pihom] * fac;
    wk[I_Pihtr] = wk[I_Pihtr] * fac;
    wk[I_Psacw] = wk[I_Psacw] * fac;
    wk[I_Psfw ] = wk[I_Psfw ] * fac;
    wk[I_Pgacw] = wk[I_Pgacw] * fac;

    // [QI]
    net = 
        + wk[I_Pigen  ]  // [prod] v->i
        + wk[I_Pidep  ]  // [prod] v->i
        + wk[I_Pihom  ]  // [prod] c->i
        + wk[I_Pihtr  ]  // [prod] c->i
        - wk[I_Pisub  ]  // [loss] i->v
        - wk[I_Pimlt  ]  // [loss] i->c
        - wk[I_Psaut  ]  // [loss] i->s
        - wk[I_Praci_s]  // [loss] i->s
        - wk[I_Psaci  ]  // [loss] i->s
        - wk[I_Psfi   ]  // [loss] i->s
        - wk[I_Praci_g]  // [loss] i->g
        - wk[I_Pgaci  ]; // [loss] i->g

    fac_sw = 0.5 + Kokkos::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
    fac = fac_sw + ( 1.0 - fac_sw ) 
            * Kokkos::min( -wk[I_dqi_dt]/(net - fac_sw), 1.0 ); // loss limiter

    wk[I_Pigen  ] = wk[I_Pigen  ] * fac;
    wk[I_Pidep  ] = wk[I_Pidep  ] * fac;
    wk[I_Pihom  ] = wk[I_Pihom  ] * fac;
    wk[I_Pihtr  ] = wk[I_Pihtr  ] * fac;
    wk[I_Pisub  ] = wk[I_Pisub  ] * fac;
    wk[I_Pimlt  ] = wk[I_Pimlt  ] * fac;
    wk[I_Psaut  ] = wk[I_Psaut  ] * fac;
    wk[I_Praci_s] = wk[I_Praci_s] * fac;
    wk[I_Psaci  ] = wk[I_Psaci  ] * fac;
    wk[I_Psfi   ] = wk[I_Psfi   ] * fac;
    wk[I_Praci_g] = wk[I_Praci_g] * fac;
    wk[I_Pgaci  ] = wk[I_Pgaci  ] * fac;


    // [QR]
    net = 
        + wk[I_Praut  ]  // [prod] c->r
        + wk[I_Pracw  ]  // [prod] c->r
        + wk[I_Psmlt  ]  // [prod] s->r
        + wk[I_Pgmlt  ]  // [prod] g->r
        - wk[I_Prevp  ]  // [loss] r->v
        - wk[I_Piacr_s]  // [loss] r->s
        - wk[I_Psacr_s]  // [loss] r->s
        - wk[I_Piacr_g]  // [loss] r->g
        - wk[I_Psacr_g]  // [loss] r->g
        - wk[I_Pgacr  ]  // [loss] r->g
        - wk[I_Pgfrz  ]; // [loss] r->g

    fac_sw = 0.5 + Kokkos::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
    fac    = fac_sw + ( 1.0 - fac_sw ) 
                * Kokkos::min( -wk[I_dqr_dt]/(net - fac_sw), 1.0 ); // loss limiter

    wk[I_Praut  ] = wk[I_Praut  ] * fac;
    wk[I_Pracw  ] = wk[I_Pracw  ] * fac;
    wk[I_Psmlt  ] = wk[I_Psmlt  ] * fac;
    wk[I_Pgmlt  ] = wk[I_Pgmlt  ] * fac;
    wk[I_Prevp  ] = wk[I_Prevp  ] * fac;
    wk[I_Piacr_s] = wk[I_Piacr_s] * fac;
    wk[I_Psacr_s] = wk[I_Psacr_s] * fac;
    wk[I_Piacr_g] = wk[I_Piacr_g] * fac;
    wk[I_Psacr_g] = wk[I_Psacr_g] * fac;
    wk[I_Pgacr  ] = wk[I_Pgacr  ] * fac;
    wk[I_Pgfrz  ] = wk[I_Pgfrz  ] * fac;

    // [QV]
    net =
        + wk[I_Prevp]  // [prod] r->v
        + wk[I_Pisub]  // [prod] i->v
        + wk[I_Pssub]  // [prod] s->v
        + wk[I_Pgsub]  // [prod] g->v
        - wk[I_Pigen]  // [loss] v->i
        - wk[I_Pidep]  // [loss] v->i
        - wk[I_Psdep]  // [loss] v->s
        - wk[I_Pgdep]; // [loss] v->g

    fac_sw = 0.5 + Kokkos::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
    fac = fac_sw + ( 1.0 - fac_sw ) 
            * Kokkos::min( -wk[I_dqv_dt]/(net - fac_sw), 1.0 ); // loss limiter

    wk[I_Prevp] = wk[I_Prevp] * fac;
    wk[I_Pisub] = wk[I_Pisub] * fac;
    wk[I_Pssub] = wk[I_Pssub] * fac;
    wk[I_Pgsub] = wk[I_Pgsub] * fac;
    wk[I_Pigen] = wk[I_Pigen] * fac;
    wk[I_Pidep] = wk[I_Pidep] * fac;
    wk[I_Psdep] = wk[I_Psdep] * fac;
    wk[I_Pgdep] = wk[I_Pgdep] * fac;

    // [QS]
    net =
        + wk[I_Psdep  ]  // [prod] v->s
        + wk[I_Psacw  ]  // [prod] c->s
        + wk[I_Psfw   ]  // [prod] c->s
        + wk[I_Piacr_s]  // [prod] r->s
        + wk[I_Psacr_s]  // [prod] r->s
        + wk[I_Psaut  ]  // [prod] i->s
        + wk[I_Praci_s]  // [prod] i->s
        + wk[I_Psaci  ]  // [prod] i->s
        + wk[I_Psfi   ]  // [prod] i->s
        - wk[I_Pssub  ]  // [loss] s->v
        - wk[I_Psmlt  ]  // [loss] s->r
        - wk[I_Pgaut  ]  // [loss] s->g
        - wk[I_Pracs  ]  // [loss] s->g
        - wk[I_Pgacs  ]; // [loss] s->g

    fac_sw = 0.5 + Kokkos::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
    fac = fac_sw + ( 1.0 - fac_sw ) 
            * Kokkos::min( -wk[I_dqs_dt]/(net - fac_sw), 1.0 ); // loss limiter

    wk[I_Psdep  ] = wk[I_Psdep  ] * fac;
    wk[I_Psacw  ] = wk[I_Psacw  ] * fac;
    wk[I_Psfw   ] = wk[I_Psfw   ] * fac;
    wk[I_Piacr_s] = wk[I_Piacr_s] * fac;
    wk[I_Psacr_s] = wk[I_Psacr_s] * fac;
    wk[I_Psaut  ] = wk[I_Psaut  ] * fac;
    wk[I_Praci_s] = wk[I_Praci_s] * fac;
    wk[I_Psaci  ] = wk[I_Psaci  ] * fac;
    wk[I_Psfi   ] = wk[I_Psfi   ] * fac;
    wk[I_Pssub  ] = wk[I_Pssub  ] * fac;
    wk[I_Psmlt  ] = wk[I_Psmlt  ] * fac;
    wk[I_Pgaut  ] = wk[I_Pgaut  ] * fac;
    wk[I_Pracs  ] = wk[I_Pracs  ] * fac;
    wk[I_Pgacs  ] = wk[I_Pgacs  ] * fac;

    // [QG]
    net =
        + wk[I_Pgdep  ]  // [prod] v->g
        + wk[I_Pgacw  ]  // [prod] c->g
        + wk[I_Piacr_g]  // [prod] r->g
        + wk[I_Psacr_g]  // [prod] r->g
        + wk[I_Pgacr  ]  // [prod] r->g
        + wk[I_Pgfrz  ]  // [prod] r->g
        + wk[I_Praci_g]  // [prod] i->g
        + wk[I_Pgaci  ]  // [prod] i->g
        + wk[I_Pgaut  ]  // [prod] s->g
        + wk[I_Pracs  ]  // [prod] s->g
        + wk[I_Pgacs  ]  // [prod] s->g
        - wk[I_Pgsub  ]  // [loss] g->v
        - wk[I_Pgmlt  ]; // [loss] g->r

    fac_sw = 0.5 + Kokkos::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
    fac = fac_sw + ( 1.0 - fac_sw ) 
            * Kokkos::min( -wk[I_dqg_dt]/(net - fac_sw), 1.0 ); // loss limiter

    wk[I_Pgdep  ] = wk[I_Pgdep  ] * fac;
    wk[I_Pgacw  ] = wk[I_Pgacw  ] * fac;
    wk[I_Piacr_g] = wk[I_Piacr_g] * fac;
    wk[I_Psacr_g] = wk[I_Psacr_g] * fac;
    wk[I_Pgacr  ] = wk[I_Pgacr  ] * fac;
    wk[I_Pgfrz  ] = wk[I_Pgfrz  ] * fac;
    wk[I_Praci_g] = wk[I_Praci_g] * fac;
    wk[I_Pgaci  ] = wk[I_Pgaci  ] * fac;
    wk[I_Pgaut  ] = wk[I_Pgaut  ] * fac;
    wk[I_Pracs  ] = wk[I_Pracs  ] * fac;
    wk[I_Pgacs  ] = wk[I_Pgacs  ] * fac;
    wk[I_Pgsub  ] = wk[I_Pgsub  ] * fac;
    wk[I_Pgmlt  ] = wk[I_Pgmlt  ] * fac;

    qc_t = + wk[I_Pimlt]  // [prod] i->c
            - wk[I_Praut]  // [loss] c->r
            - wk[I_Pracw]  // [loss] c->r
            - wk[I_Pihom]  // [loss] c->i
            - wk[I_Pihtr]  // [loss] c->i
            - wk[I_Psacw]  // [loss] c->s
            - wk[I_Psfw ]  // [loss] c->s
            - wk[I_Pgacw]; // [loss] c->g

    qr_t = + wk[I_Praut  ]  // [prod] c->r
            + wk[I_Pracw  ]  // [prod] c->r
            + wk[I_Psmlt  ]  // [prod] s->r
            + wk[I_Pgmlt  ]  // [prod] g->r
            - wk[I_Prevp  ]  // [loss] r->v
            - wk[I_Piacr_s]  // [loss] r->s
            - wk[I_Psacr_s]  // [loss] r->s
            - wk[I_Piacr_g]  // [loss] r->g
            - wk[I_Psacr_g]  // [loss] r->g
            - wk[I_Pgacr  ]  // [loss] r->g
            - wk[I_Pgfrz  ]; // [loss] r->g

    qi_t = + wk[I_Pigen  ]  // [prod] v->i
            + wk[I_Pidep  ]  // [prod] v->i
            + wk[I_Pihom  ]  // [prod] c->i
            + wk[I_Pihtr  ]  // [prod] c->i
            - wk[I_Pisub  ]  // [loss] i->v
            - wk[I_Pimlt  ]  // [loss] i->c
            - wk[I_Psaut  ]  // [loss] i->s
            - wk[I_Praci_s]  // [loss] i->s
            - wk[I_Psaci  ]  // [loss] i->s
            - wk[I_Psfi   ]  // [loss] i->s
            - wk[I_Praci_g]  // [loss] i->g
            - wk[I_Pgaci  ]; // [loss] i->g

    qs_t = + wk[I_Psdep  ]  // [prod] v->s
            + wk[I_Psacw  ]  // [prod] c->s
            + wk[I_Psfw   ]  // [prod] c->s
            + wk[I_Piacr_s]  // [prod] r->s
            + wk[I_Psacr_s]  // [prod] r->s
            + wk[I_Psaut  ]  // [prod] i->s
            + wk[I_Praci_s]  // [prod] i->s
            + wk[I_Psaci  ]  // [prod] i->s
            + wk[I_Psfi   ]  // [prod] i->s
            - wk[I_Pssub  ]  // [loss] s->v
            - wk[I_Psmlt  ]  // [loss] s->r
            - wk[I_Pgaut  ]  // [loss] s->g
            - wk[I_Pracs  ]  // [loss] s->g
            - wk[I_Pgacs  ]; // [loss] s->g

    qg_t = + wk[I_Pgdep  ]  // [prod] v->g
            + wk[I_Pgacw  ]  // [prod] c->g
            + wk[I_Piacr_g]  // [prod] r->g
            + wk[I_Psacr_g]  // [prod] r->g
            + wk[I_Pgacr  ]  // [prod] r->g
            + wk[I_Pgfrz  ]  // [prod] r->g
            + wk[I_Praci_g]  // [prod] i->g
            + wk[I_Pgaci  ]  // [prod] i->g
            + wk[I_Pgaut  ]  // [prod] s->g
            + wk[I_Pracs  ]  // [prod] s->g
            + wk[I_Pgacs  ]  // [prod] s->g
            - wk[I_Pgsub  ]  // [loss] g->v
            - wk[I_Pgmlt  ]; // [loss] g->r

    qc_t = Kokkos::max( qc_t, -wk[I_dqc_dt] );
    qr_t = Kokkos::max( qr_t, -wk[I_dqr_dt] );
    qi_t = Kokkos::max( qi_t, -wk[I_dqi_dt] );
    qs_t = Kokkos::max( qs_t, -wk[I_dqs_dt] );
    qg_t = Kokkos::max( qg_t, -wk[I_dqg_dt] );

    qv_t = - ( qc_t 
                + qr_t 
                + qi_t 
                + qs_t 
                + qg_t );

    drhogqv(k,ij) = rhog(k,ij) * qv_t * dt;
    drhogqc(k,ij) = rhog(k,ij) * qc_t * dt;
    drhogqr(k,ij) = rhog(k,ij) * qr_t * dt;
    drhogqi(k,ij) = rhog(k,ij) * qi_t * dt;
    drhogqs(k,ij) = rhog(k,ij) * qs_t * dt;
    drhogqg(k,ij) = rhog(k,ij) * qg_t * dt;

    Vt(I_QR,k,ij) = Vtr;
    Vt(I_QI,k,ij) = Vti;
    Vt(I_QS,k,ij) = Vts;
    Vt(I_QG,k,ij) = Vtg;

    //--- Effective Radius of Liquid Water
    xf_qc = rhoqc / Nc(k,ij) * 1.0E-9;                            // mean mass   of qc [kg]
    rf_qc = Kokkos::pow(( coef_xf * xf_qc + EPS ), 0.33333333);       // mean radius of qc [m]
    r2_qc = coef_dgam * (rf_qc * rf_qc) * ( Nc(k,ij) * 1.0E+6 );  // r^2 moment  of qc
    r2_qr = 0.25 * N0r * GAM_3 * RLMDr_3;                          // r^2 moment  of qr
    r3_qc = coef_xf * dens * qc;                                   // r^3 moment  of qc
    r3_qr = coef_xf * dens * qr;                                   // r^3 moment  of qr

    zerosw = 0.5 - Kokkos::copysign(0.5, (r2_qc + r2_qr) - r2_min );
    rceff(k,ij) = ( r3_qc + r3_qr ) / ( r2_qc + r2_qr + zerosw ) * ( 1.0 - zerosw );

    sw = ( 0.5 + Kokkos::copysign(0.5, (qc + qr) - q_min ) ) * zerosw; // if qc+qr > 0.01[g/kg], sw=1

    rctop(0,ij) = (       sw ) * rceff(k,ij)
                    + ( 1.0 - sw ) * UNDEF;
    tctop(0,ij) = (       sw ) * temp
                    + ( 1.0 - sw ) * UNDEF;

    zerosw = 0.5 - Kokkos::copysign( 0.5, r2_qc - r2_min );
    rceff_cld(k,ij) = r3_qc / ( r2_qc + zerosw ) * ( 1.0 - zerosw );

    sw = ( 0.5 + Kokkos::copysign( 0.5, qc - q_min ) ) * zerosw; // if qc > 0.01[g/kg], sw=1

    rctop_cld(0,ij) = (       sw ) * rceff_cld(k,ij)
                        + ( 1.0 - sw ) * UNDEF;
    tctop_cld(0,ij) = (       sw ) * temp
                        + ( 1.0 - sw ) * UNDEF;

}
// ================================ mp_nsw6_new Kokkos Implementation End =========================





}

/**
 * ======================= assistant function =======================
 */
void negative_filter( double rhog  [kdim][ijdim],
                      double rhoge [kdim][ijdim],
                      double rhogq [nqmax][kdim][ijdim],
                      double rho   [kdim][ijdim],
                      double tem   [kdim][ijdim],
                      double pre   [kdim][ijdim],
                      double q     [nqmax][kdim][ijdim],
                      double gsgam2[kdim][ijdim]  )
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double qd[kdim][ijdim];
    double cva[kdim][ijdim];
    double Rdry = CONST_Rdry;
    double Rvap = CONST_Rvap;

    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            double diffq = 0.0;
            for(int nq = NQW_STR + 1; nq <= NQW_END; nq++)
            {
                // total hydrometeor (before correction)
                diffq += rhogq[nq][k][ij];
                // remove negative value of hydrometeors (mass)
                rhogq[nq][k][ij] = std::max(rhogq[nq][k][ij], 0.0);
            }

            for(int nq = NQW_STR + 1; nq <= NQW_END; nq++)
            {
                // difference between before and after correction
                diffq -= rhogq[nq][k][ij];
            }

            // Compensate for the lack of hydrometeors by the water vapor
            rhogq[I_QV][k][ij] += diffq;
        }
    }


    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            double diffq = rhogq[I_QV][k][ij];
            // remove negative value of water vapor (mass)
            rhogq[I_QV][k][ij] = std::max(rhogq[I_QV][k][ij], 0.0);

            diffq -= rhogq[I_QV][k][ij];

            // Apply correction to total density
            rhog[k][ij] = rhog[k][ij] * (1.0 - diffq);
            rho[k][ij] = rhog[k][ij] / gsgam2[k][ij];
        }
    }

    for(int nq = NQW_STR; nq <= NQW_END; nq++)
    {
        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                //--- update mass concentration
                q[nq][k][ij] = rhogq[nq][k][ij] / rhog[k][ij];
            }
        }
    }

    /**
     * q  [IN]
     * qd [OUT]
     */
    THRMDYN_qd(q, qd);
    /**
     * qd  [IN]
     * q   [IN]
     * cva [OUT]
     */
    THRMDYN_cv(qd, q, cva);

    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rhoge[k][ij] = tem[k][ij] * rhog[k][ij] * cva[k][ij];
            pre[k][ij] = rho[k][ij] * ( qd[k][ij] * Rdry + q[I_QV][k][ij] * Rvap ) * tem[k][ij];
        }
    }

}

void Bergeron_param( double tem[kdim][ijdim],
                     double a1 [kdim][ijdim],
                     double a2 [kdim][ijdim],
                     double ma2[kdim][ijdim])
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    int itemc;
    double temc;
    double fact;
    double a1_tab[32] = {0.0001E-7, 0.7939E-7, 0.7841E-6, 0.3369E-5, 0.4336E-5,
                         0.5285E-5, 0.3728E-5, 0.1852E-5, 0.2991E-6, 0.4248E-6,
                         0.7434E-6, 0.1812E-5, 0.4394E-5, 0.9145E-5, 0.1725E-4,
                         0.3348E-4, 0.1725E-4, 0.9175E-5, 0.4412E-5, 0.2252E-5,
                         0.9115E-6, 0.4876E-6, 0.3473E-6, 0.4758E-6, 0.6306E-6,
                         0.8573E-6, 0.7868E-6, 0.7192E-6, 0.6513E-6, 0.5956E-6,
                         0.5333E-6, 0.4834E-6};

    double a2_tab[32] = { 0.0100, 0.4006, 0.4831, 0.5320, 0.5307,
                          0.5319, 0.5249, 0.4888, 0.3849, 0.4047,
                          0.4318, 0.4771, 0.5183, 0.5463, 0.5651,
                          0.5813, 0.5655, 0.5478, 0.5203, 0.4906,
                          0.4447, 0.4126, 0.3960, 0.4149, 0.4320,
                          0.4506, 0.4483, 0.4460, 0.4433, 0.4413,
                          0.4382, 0.4361 };

    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            temc = std::min( std::max( tem[k][ij] - TEM00, -30.99 ), 0.0 );
            // itemc = int(-temc) + 1; in fortran 
            itemc = int(-temc);
            fact = -(temc + double(itemc - 1));
            a1[k][ij] = (1.0 - fact) * a1_tab[itemc] +
                        (fact) * a1_tab[itemc + 1];
            a2[k][ij] = (1.0 - fact) * a2_tab[itemc] +
                        (fact) * a2_tab[itemc + 1];
            ma2[k][ij] = 1.0 - a2[k][ij];

            a1[k][ij] = a1[k][ij] * std::pow(1.0E-3, ma2[k][ij]);
        }
    } 
}

/**
 * ================================= Kokkos ver. ==============================================
 */
void negative_filter( const View2D<double, DEFAULT_MEM>&  rhog,
                      const View2D<double, DEFAULT_MEM>&  rhoge, 
                      const View3D<double, DEFAULT_MEM>&  rhogq,
                      const View2D<double, DEFAULT_MEM>&  rho,   
                      const View2D<double, DEFAULT_MEM>&  tem,   
                      const View2D<double, DEFAULT_MEM>&  pre,   
                      const View3D<double, DEFAULT_MEM>&  q,     
                      const View2D<double, DEFAULT_MEM>&  gsgam2)
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    View2D<double, DEFAULT_MEM> qd ("qd", kdim, ijdim);
    View2D<double, DEFAULT_MEM> cva("cva", kdim, ijdim);
    // View<double> Rdry;
    // View<double> Rvap; 
    // Kokkos::deep_copy(Rdry, CONST_Rdry);
    // Kokkos::deep_copy(Rvap, CONST_Rvap); 

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        double diffq = 0.0;
        for(int nq = NQW_STR + 1; nq <= NQW_END; nq++)
        {
            // total hydrometeor (before correction)
            diffq += rhogq(nq,k,ij);
            // remove negative value of hydrometeors (mass)
            rhogq(nq,k,ij) = Kokkos::max(rhogq(nq,k,ij), 0.0);
        }

        for(int nq = NQW_STR + 1; nq <= NQW_END; nq++)
        {
            // difference between before and after correction
            diffq -= rhogq(nq,k,ij);
        }

        // Compensate for the lack of hydrometeors by the water vapor
        rhogq(I_QV,k,ij) += diffq; 
    });
    
    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
       double diffq = rhogq(I_QV,k,ij);
       // remove negative value of water vapor (mass)
       rhogq(I_QV,k,ij) = Kokkos::max(rhogq(I_QV,k,ij), 0.0);

       diffq -= rhogq(I_QV,k,ij);

       // Apply correction to total density
       rhog(k,ij) = rhog(k,ij) * (1.0 - diffq);
       rho (k,ij) = rhog(k,ij) / gsgam2(k,ij);
    });

    for(int nq = NQW_STR; nq <= NQW_END; nq++)
    {
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            //--- update mass concentration
            q(nq,k,ij) = rhogq(nq,k,ij) / rhog(k,ij);
        });
    }

    /**
     * q  [IN]
     * qd [OUT]
     */
    THRMDYN_qd(q, qd);
    /**
     * qd  [IN]
     * q   [IN]
     * cva [OUT]
     */
    THRMDYN_cv(qd, q, cva);

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        const double Rdry = CONST_Rdry;
        const double Rvap = CONST_Rvap;
        rhoge(k,ij) = tem(k,ij) * rhog(k,ij) * cva(k,ij);
        pre  (k,ij) = rho(k,ij) * ( qd(k,ij) * Rdry + q(I_QV,k,ij) * Rvap ) * tem(k,ij);
    });
}

void Bergeron_param( const View2D<double, DEFAULT_MEM> tem,
                     const View2D<double, DEFAULT_MEM> a1 ,
                     const View2D<double, DEFAULT_MEM> a2 ,
                     const View2D<double, DEFAULT_MEM> ma2)
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    // int itemc;
    // double temc;
    // double fact;
    double arr_a1_tab[32] = {0.0001E-7, 0.7939E-7, 0.7841E-6, 0.3369E-5, 0.4336E-5,
                             0.5285E-5, 0.3728E-5, 0.1852E-5, 0.2991E-6, 0.4248E-6,
                             0.7434E-6, 0.1812E-5, 0.4394E-5, 0.9145E-5, 0.1725E-4,
                             0.3348E-4, 0.1725E-4, 0.9175E-5, 0.4412E-5, 0.2252E-5,
                             0.9115E-6, 0.4876E-6, 0.3473E-6, 0.4758E-6, 0.6306E-6,
                             0.8573E-6, 0.7868E-6, 0.7192E-6, 0.6513E-6, 0.5956E-6,
                             0.5333E-6, 0.4834E-6};

    double arr_a2_tab[32] = { 0.0100, 0.4006, 0.4831, 0.5320, 0.5307,
                              0.5319, 0.5249, 0.4888, 0.3849, 0.4047,
                              0.4318, 0.4771, 0.5183, 0.5463, 0.5651,
                              0.5813, 0.5655, 0.5478, 0.5203, 0.4906,
                              0.4447, 0.4126, 0.3960, 0.4149, 0.4320,
                              0.4506, 0.4483, 0.4460, 0.4433, 0.4413,
                              0.4382, 0.4361 };
    
    View1D<double, DEFAULT_MEM> a1_tab("a1_tab", 32);
    View1D<double, DEFAULT_MEM> a2_tab("a2_tab", 32);

    auto h_a1_tab = Kokkos::create_mirror_view(a1_tab);
    auto h_a2_tab = Kokkos::create_mirror_view(a2_tab);
    for(int i = 0; i < 32; i++)
    {
        h_a1_tab(i) = arr_a1_tab[i];
        h_a2_tab(i) = arr_a2_tab[i];
    }
    Kokkos::deep_copy(a1_tab, h_a1_tab);
    Kokkos::deep_copy(a2_tab, h_a2_tab);

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        const double tem00 = TEM00;
        double temc = Kokkos::min( Kokkos::max( tem(k,ij) - tem00, -30.99 ), 0.0 );
        // itemc = int(-temc) + 1; in fortran 
        int itemc = int(-temc);
        double fact = -(temc + double(itemc - 1));
        a1(k,ij) = (1.0 - fact) * a1_tab[itemc] +
                    (fact) * a1_tab[itemc + 1];
        a2(k,ij) = (1.0 - fact) * a2_tab[itemc] +
                    (fact) * a2_tab[itemc + 1];
        ma2(k,ij) = 1.0 - a2(k,ij);
        
        a1(k,ij) = a1(k,ij) * Kokkos::pow(1.0E-3, ma2(k,ij));
    });
}


// /**
//  * Naive Kokkos Ver.
//  */
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
//             )
// {
// #ifdef ENABLE_DEBUG
//     std::cout << __PRETTY_FUNCTION__ << std::endl;
// #endif
//     // working
//     View<double**> drhogqv("drhogqv", kdim, ijdim);
//     View<double**> drhogqc("drhogqc", kdim, ijdim);
//     View<double**> drhogqi("drhogqi", kdim, ijdim);
//     View<double**> drhogqr("drhogqr", kdim, ijdim);
//     View<double**> drhogqs("drhogqs", kdim, ijdim);
//     View<double**> drhogqg("drhogqg", kdim, ijdim);

//     View<double**> psatl("psatl", kdim, ijdim);
//     View<double**> qsatl("qsatl", kdim, ijdim);                //< saturated water vapor for liquid water [kg/kg]
//     View<double**> psati("psati", kdim, ijdim);
//     View<double**> qsati("qsati", kdim, ijdim);                //< saturated water vapor for ice water    [kg/kg]
//     View<double**> Nc   ("Nc   ", kdim, ijdim);                //< Number concentration of cloud water [1/cc]

//     double dens    ;                         //< density
//     double temp    ;                         //< T [K]
//     double qv      ;                         //< mixing ratio of water vapor  [kg/kg]
//     double qc      ;                         //< mixing ratio of liquid water [kg/kg]
//     double qr      ;                         //< mixing ratio of rain         [kg/kg]
//     double qi      ;                         //< mixing ratio of ice water    [kg/kg]
//     double qs      ;                         //< mixing ratio of snow         [kg/kg]
//     double qg      ;                         //< mixing ratio of graupel      [kg/kg]
//     double qv_t    ;                         //< tendency     of water vapor  [kg/kg/s]
//     double qc_t    ;                         //< tendency     of liquid water [kg/kg/s]
//     double qr_t    ;                         //< tendency     of rain         [kg/kg/s]
//     double qi_t    ;                         //< tendency     of ice water    [kg/kg/s]
//     double qs_t    ;                         //< tendency     of snow         [kg/kg/s]
//     double qg_t    ;                         //< tendency     of graupel      [kg/kg/s]
//     double Sliq    ;                         //< saturated ratio S for liquid water [0-1]
//     double Sice    ;                         //< saturated ratio S for ice water    [0-1]
//     double Rdens   ;                         //< 1 / density
//     double rho_fact;                         //< density factor
//     double temc    ;                         //< T - T0 [K]

//     double RLMDr, RLMDr_2, RLMDr_3        ;
//     double RLMDs, RLMDs_2, RLMDs_3        ;
//     double RLMDg, RLMDg_2, RLMDg_3        ;
//     double RLMDr_1br, RLMDr_2br, RLMDr_3br;
//     double RLMDs_1bs, RLMDs_2bs, RLMDs_3bs;
//     double RLMDr_dr, RLMDr_3dr, RLMDr_5dr ;
//     double RLMDs_ds, RLMDs_3ds, RLMDs_5ds ;
//     double RLMDg_dg, RLMDg_3dg, RLMDg_5dg ;
//     double RLMDr_7                        ;
//     double RLMDr_6dr                      ;

//     //---< Roh and Satoh (2014) >---
//     double tems, Xs2                     ;
//     double MOMs_0, MOMs_1, MOMs_2        ;
//     double MOMs_0bs, MOMs_1bs, MOMs_2bs  ;
//     double MOMs_2ds, MOMs_5ds_h, RMOMs_Vt;
//     double coef_at[4]                    ;
//     double coef_bt[4]                    ;
//     double loga_, b_, nm                 ;

//     double Vti, Vtr, Vts, Vtg   ;           //< terminal velocity
//     double Esi_mod, Egs_mod     ;           //< modified accretion efficiency
//     double rhoqc                ;           //< rho * qc
//     double Pracw_orig,  Pracw_kk;           //< accretion       term by orig  & k-k scheme
//     double Praut_berry, Praut_kk;           //< auto-conversion term by berry & k-k scheme
//     double Dc                   ;           //< relative variance
//     double betai, betas         ;           //< sticky parameter for auto-conversion
//     double Ka                   ;           //< thermal diffusion coefficient of air
//     double Kd                   ;           //< diffusion coefficient of water vapor in air
//     double Nu                   ;           //< kinematic viscosity of air
//     double Glv, Giv, Gil        ;           //< thermodynamic function
//     double ventr, vents, ventg  ;           //< ventilation factor
//     double net, fac, fac_sw     ;
//     double zerosw, tmp          ;

//     //---< Bergeron process >---
//     double sw_bergeron     ;                 //< if 0C<T<30C, sw=1
//     View<double**> a1 ("a1", kdim, ijdim);                 //<
//     View<double**> a2 ("a2", kdim, ijdim);                 //<
//     View<double**> ma2("ma2", kdim, ijdim);                //< 1-a2
//     double dt1             ;                 //< time during which the an ice particle of 40um grows to 50um
//     double Ni50            ;                 //< number concentration of ice particle of 50um

//     //---< Explicit ice generation >---
//     double sw, rhoqi, XNi, XMi, Di, Ni0, Qi0;

//     //---< Effective radius >---
//     constexpr double r2_min = 1.0E-10;
//     constexpr double q_min  = 1.0E-5;
//     double xf_qc;                 //< mean mass   of qc [kg]
//     double rf_qc;                 //< mean radius of qc [m]
//     double r2_qc;                 //< r^2 moment  of qc
//     double r2_qr;                 //< r^2 moment  of qr
//     double r3_qc;                 //< r^3 moment  of qc
//     double r3_qr;                 //< r^3 moment  of qr
//     double dgamma_a, GAM_dgam23, GAM_dgam;
//     double coef_dgam, coef_xf;

//     //---< Precipitation >---
//     bool preciptation_flag[nqmax] = {false};

//     View<double***> Vt    ("Vt  ", nqmax, kdim, ijdim);
//     View<double**>  cva   ("cva ", kdim, ijdim);
//     View<double**>  rgs   ("rgs ", kdim, ijdim);
//     View<double**>  rgsh  ("rgsh", kdim, ijdim);

//     double wk    [wk_nmax];
//     //double ml_wk   [wk_nmax][kdim][ijdim]; tentatively remove due to the stack overflow
//     View<double**> ml_Pconv("ml_Pconv", kdim, ijdim);
//     View<double**> ml_Pconw("ml_Pconw", kdim, ijdim);
//     View<double**> ml_Pconi("ml_Pconi", kdim, ijdim);

//     double UNDEF, EPS, PI, Rvap, LHV0, LHS0, LHF0, PRE00;

//     // constexpr int simdlen = 8;
//     // int blk, vec, veclen;

//     //---------------------------------------------------------------------------
//     UNDEF = CONST_UNDEF;
//     EPS   = CONST_EPS  ;
//     PI    = CONST_PI   ;
//     Rvap  = CONST_Rvap ;
//     LHV0  = CONST_LHV0 ;
//     LHS0  = CONST_LHS0 ;
//     LHF0  = CONST_LHF0 ;
//     PRE00 = CONST_Pstd ;

//     negative_filter(rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2);

//     if(OPT_INDIR) // aerosol indirect effect
//     {
//         Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
//         KOKKOS_LAMBDA(const size_t k, const size_t ij){
//             Nc(k,ij) = std::max( UNCCN(k,ij) * 1.0E-6, Nc_def );
//         });
//     }
//     else
//     {
//         Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
//         KOKKOS_LAMBDA(const size_t k, const size_t ij){
//             Nc(k,ij) = Nc_def;
//         });
//     }

//     // saturation water contensts
//     SATURATION_psat_liq(tem, psatl);
//     SATURATION_psat_ice(tem, psati);

//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
//     KOKKOS_LAMBDA(const size_t k, const size_t ij){
//         qsatl(k,ij) = psatl(k,ij) / ( rho(k,ij) * Rvap * tem(k,ij) );
//         qsati(k,ij) = psati(k,ij) / ( rho(k,ij) * Rvap * tem(k,ij) );
//     });

//     // Bergeron process parameters
//     Bergeron_param(tem, a1, a2, ma2);

//     // work for Effective Radius of Liquid Water
//     // dgamma_a is parameter of Gamma-Dist.
//     // "Berry and Reinhardt(1974-a) eq.(16-a,b)"
//     // 0.38 * {var1/2x} =~ {var1/2 r}
//     //         dgamma_a = 1.0_RP /{var x}
//     Dc = 0.146 - 5.964E-2 * std::log(Nc_def / 2000.0);
//     dgamma_a = 0.1444 / (Dc * Dc); // Dc**2 in fortran
//     GAM_dgam   = MISC_gammafunc(dgamma_a);
//     GAM_dgam23 = MISC_gammafunc(dgamma_a + 2.0/3.0);
//     coef_dgam  = GAM_dgam23 / GAM_dgam * std::pow(dgamma_a, (-2.0/3.0));
//     coef_xf    = 3.0 / 4.0 / PI / rho_w;

//     // !! Big loop start here !!
//     for(int k = kmin; k <= kmax; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             dens = rho(k,ij);
//             temp = tem(k,ij);
//             qv   = std::max( q(I_QV, k, ij), 0.0 );
//             qc   = std::max( q(I_QC, k, ij), 0.0 );
//             qr   = std::max( q(I_QR, k, ij), 0.0 );
//             qi   = std::max( q(I_QI, k, ij), 0.0 );
//             qs   = std::max( q(I_QS, k, ij), 0.0 );
//             qg   = std::max( q(I_QG, k, ij), 0.0 );

//             // saturation ration S
//             Sliq = qv / (std::max(qsatl(k,ij), EPS));
//             Sice = qv / (std::max(qsati(k,ij), EPS));

//             Rdens = 1.0 / dens;
//             rho_fact = std::sqrt(dens00 * Rdens);
//             temc = temp - TEM00;

//             wk[I_delta1] = ( 0.5 + std::copysign(0.5, qr - 1.0E-4) );

//             wk[I_delta2] = ( 0.5 + std::copysign(0.5, 1.0E-4 - qr) )
//                             * ( 0.5 + std::copysign(0.5, 1.0E-4 - qs) );
            
//             wk[I_spsati] = 0.5 + std::copysign(0.5, Sice - 1.0);

//             wk[I_iceflg] = 0.5 - std::copysign(0.5, temc); // 0: warm, 1: ice

//             wk[I_dqv_dt] = qv / dt;
//             wk[I_dqc_dt] = qc / dt;
//             wk[I_dqr_dt] = qr / dt;
//             wk[I_dqi_dt] = qi / dt;
//             wk[I_dqs_dt] = qs / dt;
//             wk[I_dqg_dt] = qg / dt;

//             sw_bergeron = ( 0.5 + std::copysign(0.5, temc + 30.0) )
//                             * ( 0.5 + std::copysign(0.5, 0.0 - temc) )
//                             * ( 1.0 - sw_expice );
            
//             // slope parameter lambda (Rain)
//             zerosw = 0.5 - std::copysign(0.5, qr - 1.0E-12);
//             RLMDr = std::sqrt( std::sqrt( dens * qr / (Ar * N0r * GAM_1br) + zerosw ) ) * (1.0 - zerosw);

//             RLMDr_dr  = std::sqrt(RLMDr);   // **Dr
//             RLMDr_2   = std::pow(RLMDr, 2);
//             RLMDr_3   = std::pow(RLMDr, 3);
//             RLMDr_7   = std::pow(RLMDr, 7);
//             RLMDr_1br = std::pow(RLMDr, 4); // (1+Br)
//             RLMDr_2br = std::pow(RLMDr, 5); // (2+Br)
//             RLMDr_3br = std::pow(RLMDr, 6); // (3+Br)
//             RLMDr_3dr = std::pow(RLMDr, 3) * RLMDr_dr;
//             RLMDr_5dr = std::pow(RLMDr, 5) * RLMDr_dr;
//             RLMDr_6dr = std::pow(RLMDr, 6) * RLMDr_dr;

//             // slope parameter lambda (Snow)
//             zerosw = 0.5 - std::copysign(0.5, qs - 1.0E-12);
//             RLMDs  = std::sqrt( std::sqrt( dens * qs / ( As * N0s * GAM_1bs ) + zerosw ) ) * ( 1.0 - zerosw );

//             RLMDs_ds  = std::sqrt( std::sqrt(RLMDs) ); // **Ds
//             RLMDs_2   = std::pow(RLMDs, 2);
//             RLMDs_3   = std::pow(RLMDs, 3);
//             RLMDs_1bs = std::pow(RLMDs, 4); // (1+Bs)
//             RLMDs_2bs = std::pow(RLMDs, 5); // (2+Bs)
//             RLMDs_3bs = std::pow(RLMDs, 6); // (3+Bs)
//             RLMDs_3ds = std::pow(RLMDs, 3) * RLMDs_ds;
//             RLMDs_5ds = std::pow(RLMDs, 5) * RLMDs_ds;

//             MOMs_0     = N0s * GAM       * RLMDs    ;       // Ns * 0th moment
//             MOMs_1     = N0s * GAM_2     * RLMDs_2  ;       // Ns * 1st moment
//             MOMs_2     = N0s * GAM_3     * RLMDs_3  ;       // Ns * 2nd moment
//             MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs;       // Ns * 0+bs
//             MOMs_1bs   = N0s * GAM_2bs   * RLMDs_2bs;       // Ns * 1+bs
//             MOMs_2bs   = N0s * GAM_3bs   * RLMDs_3bs;       // Ns * 2+bs
//             MOMs_2ds   = N0s * GAM_3ds   * RLMDs_3ds;       // Ns * 2+ds
//             MOMs_5ds_h = N0s * GAM_5ds_h * std::sqrt(RLMDs_5ds); // Ns * (5+ds)/2
//             RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds;

//             //---< modification by Roh and Satoh (2014) >---

//             // bimodal size distribution of snow
//             Xs2 = dens * qs / As;
//             zerosw = 0.5 - std::copysign(0.5, Xs2 - 1.0E-12);

//             tems = std::min(-0.1, temc);
//             coef_at[0] = coef_a[0] + tems * ( coef_a[1] + tems * ( coef_a[4] + tems * coef_a[8] ) );
//             coef_at[1] = coef_a[2] + tems * ( coef_a[3] + tems *   coef_a[6] );
//             coef_at[2] = coef_a[5] + tems *   coef_a[7];
//             coef_at[3] = coef_a[9];
//             coef_bt[0] = coef_b[0] + tems * ( coef_b[1] + tems * ( coef_b[4] + tems * coef_b[8] ) );
//             coef_bt[1] = coef_b[2] + tems * ( coef_b[3] + tems *   coef_b[6] );
//             coef_bt[2] = coef_b[5] + tems *   coef_b[7];
//             coef_bt[3] = coef_b[9];

//             // 0th moment
//             loga_ = coef_at[0];
//             b_    = coef_bt[0];
//             MOMs_0 = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                      + ( 1.0 - sw_roh2014 ) * MOMs_0;
//             // 1st moment
//             nm = 1.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_1 = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw ) * b_) * ( 1.0 - zerosw )
//                      + ( 1.0 - sw_roh2014 ) * MOMs_1;
//             // 2nd moment
//             MOMs_2 = sw_roh2014 * Xs2 + (1.0 - sw_roh2014) * MOMs_2;
//             // 0 + Bs(=2) moment
//             nm = 2.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_0bs = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_0bs;
//             // 1 + Bs(=2) moment
//             nm = 3.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_1bs = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_1bs;
//             // 2 + Bs(=2) moment
//             nm = 4.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_2bs = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_2bs;
//             // 2 + Ds(=0.25) moment
//             nm = 2.25;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_2ds = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_2ds;
//             // ( 3 + Ds(=0.25) ) / 2  moment
//             nm = 1.625;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_5ds_h = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                          + ( 1.0 - sw_roh2014 ) * MOMs_5ds_h;
//             // Bs(=2) + Ds(=0.25) moment
//             nm = 2.25;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             RMOMs_Vt = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw ) / (MOMs_0bs + zerosw)
//                        + ( 1.0 - sw_roh2014 ) * RMOMs_Vt;

//             // slope parameter lambda (Graupel)
//             zerosw = 0.5 - std::copysign(0.5, qg - 1.0E-12);
//             RLMDg  = std::sqrt( std::sqrt( dens * qg / ( Ag * N0g * GAM_1bg ) + zerosw ) ) * ( 1.0 - zerosw );
//             RLMDg_dg  = std::sqrt( RLMDg );       // **Dg
//             RLMDg_2   = std::pow(RLMDg, 2);
//             RLMDg_3   = std::pow(RLMDg, 3);
//             RLMDg_3dg = std::pow(RLMDg, 3) * RLMDg_dg;
//             RLMDg_5dg = std::pow(RLMDg, 5) * RLMDg_dg;

//             wk[I_RLMDr] = RLMDr;
//             wk[I_RLMDs] = RLMDs;
//             wk[I_RLMDg] = RLMDg;

//             //---< terminal velocity >---
//             zerosw = 0.5 - std::copysign(0.5, qi - 1.0E-8 );
//             Vti = ( 1.0 - sw_constVti ) * (-3.29) * std::exp( std::log( dens * qi + zerosw ) * 0.16 ) * ( 1.0 - zerosw ) 
//                   + ( sw_constVti ) * (-CONST_Vti);
//             Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr;
//             Vts = -Cs * rho_fact * RMOMs_Vt;
//             Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg;

//             //---< Nucleation >---
//             // [Pigen] ice nucleation
//             Ni0 = std::max( std::exp(-0.1 * temc), 1.0 ) * 1000.0;
//             Qi0 = 4.92E-11 * std::exp( std::log(Ni0) * 1.33 ) * Rdens;

//             wk[I_Pigen] = std::max( std::min( Qi0 - qi, qv - qsati(k,ij) ), 0.0 ) / dt;

//             //---< Accretion >---
//             Esi_mod = std::min( Esi, Esi * std::exp( gamma_sacr * temc ) );
//             Egs_mod = std::min( Egs, Egs * std::exp( gamma_gacs * temc ) );

//             // [Pracw] accretion rate of cloud water by rain
//             Pracw_orig = qc * 0.25 * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact;

//             zerosw     = 0.5 - std::copysign(0.5, qc * qr - 1.0E-12 );
//             Pracw_kk   = 67.0 * std::exp( std::log( qc * qr + zerosw ) * 1.15 ) * ( 1.0 - zerosw ); // eq.(33) in KK(2000)

//             // switch orig / k-k scheme
//             wk[I_Pracw] = ( 1.0 - sw_kk2000 ) * Pracw_orig
//                         + (       sw_kk2000 ) * Pracw_kk;

//             // [Psacw] accretion rate of cloud water by snow
//             wk[I_Psacw] = qc * 0.25 * PI * Esw * Cs * MOMs_2ds * rho_fact;

//             // [Pgacw] accretion rate of cloud water by graupel
//             wk[I_Pgacw] = qc * 0.25 * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact;

//             // [Praci] accretion rate of cloud ice by rain
//             wk[I_Praci] = qi * 0.25 * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact;

//             // [Psaci] accretion rate of cloud ice by snow
//             wk[I_Psaci] = qi * 0.25 * PI * Esi_mod * Cs * MOMs_2ds * rho_fact;

//             // [Pgaci] accretion rate of cloud ice by grupel
//             wk[I_Pgaci] = qi * 0.25 * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact;

//             // [Piacr] accretion rate of rain by cloud ice
//             wk[I_Piacr] = qi * Ar / mi * 0.25 * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact;

//             // [Psacr] accretion rate of rain by snow
//             wk[I_Psacr] = Ar * 0.25 * PI * Rdens * Esr * N0r * std::abs(Vtr - Vts)
//                              * ( GAM_1br * RLMDr_1br * MOMs_2          
//                                  + 2.0 * GAM_2br * RLMDr_2br * MOMs_1          
//                                  + GAM_3br * RLMDr_3br * MOMs_0 );

//             // [Pgacr] accretion rate of rain by graupel
//             wk[I_Pgacr] = Ar * 0.25 * PI * Rdens * Egr * N0g * N0r * std::abs(Vtg - Vtr)
//                              * ( GAM_1br * RLMDr_1br * GAM_3 * RLMDg_3
//                                  + 2.0 * GAM_2br * RLMDr_2br * GAM_2 * RLMDg_2
//                                  + GAM_3br * RLMDr_3br * GAM * RLMDg );

//             // [Pracs] accretion rate of snow by rain
//             wk[I_Pracs] = As * 0.25 * PI * Rdens * Esr *  N0r * std::abs(Vtr - Vts)
//                              * ( MOMs_0bs * GAM_3 * RLMDr_3
//                                  + 2.0 * MOMs_1bs * GAM_2 * RLMDr_2
//                                  + MOMs_2bs * GAM * RLMDr );

//             // [Pgacs] accretion rate of snow by graupel
//             wk[I_Pgacs] = As * 0.25 * PI * Rdens * Egs_mod * N0g * std::abs(Vtg - Vts)
//                              * ( MOMs_0bs * GAM_3 * RLMDg_3
//                                  + 2.0 * MOMs_1bs * GAM_2 * RLMDg_2
//                                  + MOMs_2bs * GAM * RLMDg );

//             //---< Autoconversion >---            
//             // [Praut] auto-conversion rate from cloud water to rain
//             rhoqc = dens * qc * 1000.0; // [g/m3]
//             Dc    = 0.146 - 5.964E-2 * std::log( Nc(k,ij) / 2000.0 );
//             Praut_berry = Rdens * 1.67E-5 * rhoqc * rhoqc / ( 5.0 + 3.66E-2 * Nc(k,ij) / ( Dc * rhoqc + EPS ) );

//             zerosw   = 0.5 - std::copysign(0.5, qc - 1.0E-12 );
//             Praut_kk = 1350.0                                           
//                        * std::exp( std::log( qc + zerosw ) * 2.47 ) * ( 1.0 - zerosw ) 
//                        * std::exp( std::log( Nc(k,ij) ) * (-1.79) );                     // eq.(29) in KK(2000)

//             Praut_kk = 1350.0 * std::pow(qc, 2.47) * std::pow(Nc(k,ij), -1.79);


//             // switch berry / k-k scheme
//             wk[I_Praut] = ( 1.0 - sw_kk2000 ) * Praut_berry
//                           + sw_kk2000 * Praut_kk;

//             // [Psaut] auto-conversion rate from cloud ice to snow
//             betai = std::min( beta_saut, beta_saut * std::exp( gamma_saut * temc ) );
//             wk[I_Psaut] = std::max( betai * (qi - qicrt_saut), 0.0 );

//             // [Pgaut] auto-conversion rate from snow to graupel
//             betas = std::min( beta_gaut, beta_gaut * std::exp( gamma_gaut * temc ) );
//             wk[I_Pgaut] = std::max( betas * (qs - qscrt_gaut), 0.0 );

//             //---< Evaporation, Sublimation, Melting, and Freezing >---
//             Ka  = ( Ka0 + dKa_dT * temc );
//             Kd  = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(k,ij);
//             Nu  = ( nu0 + dnu_dT * temc ) * Rdens;

//             Glv = 1.0 / ( LHV0/(Ka * temp) * ( LHV0/(Rvap * temp) - 1.0 ) + 1.0 / (Kd * dens * qsatl(k,ij)) );
//             Giv = 1.0 / ( LHS0/(Ka * temp) * ( LHS0/(Rvap * temp) - 1.0 ) + 1.0 / (Kd * dens * qsati(k,ij)) );
//             Gil = 1.0 / ( LHF0/(Ka * temc) );

//             // [Prevp] evaporation rate of rain
//             ventr = f1r * GAM_2 * RLMDr_2 + f2r * std::sqrt( Cr * rho_fact / Nu * RLMDr_5dr ) * GAM_5dr_h;

//             wk[I_Prevp] = 2.0 * PI * Rdens * N0r * ( 1.0 - std::min(Sliq, 1.0) ) * Glv * ventr;

//             // [Pidep,Pisub] deposition/sublimation rate for ice
//             rhoqi = std::max(dens * qi, EPS);
//             XNi   = std::min( std::max( 5.38E+7 * std::exp( std::log(rhoqi) * 0.75 ), 1.0E+3 ), 1.0E+6 );
//             XMi   = rhoqi / XNi;
//             Di    = std::min( Di_a * std::sqrt(XMi), Di_max );

//             tmp = 4.0 * Di * XNi * Rdens * ( Sice - 1.0 ) * Giv;

//             wk[I_Pidep] = (       wk[I_spsati] ) * ( tmp); // Sice > 2
//             wk[I_Pisub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1

//             // [Pihom] homogenious freezing at T < -40C
//             sw = ( 0.5 - std::copysign(0.5, temc + 40.0) ); // if T < -40C, sw=1

//             wk[I_Pihom] = sw * qc / dt;

//             // [Pihtr] heteroginous freezing at -40C < T < 0C
//             sw =   ( 0.5 + std::copysign(0.5, temc + 40.0) )
//                  * ( 0.5 - std::copysign(0.5, temc       ) ); // if -40C < T < 0C, sw=1

//             wk[I_Pihtr] = sw * ( dens / rho_w * (qc*qc) / ( Nc_ihtr * 1.0E+6 ) )
//                              * B_frz * ( std::exp(-A_frz * temc) - 1.0 );

//             // [Pimlt] ice melting at T > 0C
//             sw = ( 0.5 + std::copysign(0.5, temc) ); // if T > 0C, sw=1

//             wk[I_Pimlt] = sw * qi / dt;

//             // [Psdep,Pssub] deposition/sublimation rate for snow
//             vents = f1s * MOMs_1 + f2s * std::sqrt( Cs * rho_fact / Nu ) * MOMs_5ds_h;

//             tmp = 2.0 * PI * Rdens * ( Sice - 1.0 ) * Giv * vents;

//             wk[I_Psdep] = (       wk[I_spsati] ) * ( tmp); // Sice > 1
//             wk[I_Pssub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1            

//             // [Psmlt] melting rate of snow
//             wk[I_Psmlt] = 2.0 * PI * Rdens * Gil * vents
//                           + CL * temc / LHF0 * ( wk[I_Psacw] + wk[I_Psacr] );
//             wk[I_Psmlt] = std::max( wk[I_Psmlt], 0.0 );

//             // [Pgdep/pgsub] deposition/sublimation rate for graupel
//             ventg = f1g * GAM_2 * RLMDg_2 + f2g * std::sqrt( Cg * rho_fact / Nu * RLMDg_5dg ) * GAM_5dg_h;

//             tmp = 2.0 * PI * Rdens * N0g * ( Sice - 1.0 ) * Giv * ventg;

//             wk[I_Pgdep] = (       wk[I_spsati] ) * ( tmp); // Sice > 1
//             wk[I_Pgsub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1

//             // [Pgmlt] melting rate of graupel
//             wk[I_Pgmlt] = 2.0 * PI * Rdens * N0g * Gil * ventg
//                           + CL * temc / LHF0 * ( wk[I_Pgacw] + wk[I_Pgacr] );
//             wk[I_Pgmlt] = std::max( wk[I_Pgmlt], 0.0 );

//             // [Pgfrz] freezing rate of graupel
//             wk[I_Pgfrz] = 2.0 * PI * Rdens * N0r * 60.0 * B_frz * Ar * ( std::exp(-A_frz * temc) - 1.0 ) * RLMDr_7;

//             // [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
//             dt1  = ( std::exp( std::log(mi50) * ma2(k,ij) )
//                    - std::exp( std::log(mi40) * ma2(k,ij) ) ) / ( a1(k,ij) * ma2(k,ij) );
//             Ni50 = qi * dt / ( mi50 * dt1 );

//             wk[I_Psfw] = Ni50 * ( a1(k,ij) * std::pow(mi50, a2(k,ij))
//                                   + PI * Eiw * dens * qc * Ri50 * Ri50 * vti50 );
//             wk[I_Psfi] = qi / dt1;

//             //---< limiter >---
//             wk[I_Pigen] = std::min( wk[I_Pigen], wk[I_dqv_dt] ) * ( wk[I_iceflg] ) * sw_expice;
//             wk[I_Pidep] = std::min( wk[I_Pidep], wk[I_dqv_dt] ) * ( wk[I_iceflg] ) * sw_expice;
//             wk[I_Psdep] = std::min( wk[I_Psdep], wk[I_dqv_dt] ) * ( wk[I_iceflg] )            ;
//             wk[I_Pgdep] = std::min( wk[I_Pgdep], wk[I_dqv_dt] ) * ( wk[I_iceflg] )            ;

//             wk[I_Pracw] = wk[I_Pracw]                           
//                         + wk[I_Psacw] * ( 1.0 - wk[I_iceflg] )   // c->r by s
//                         + wk[I_Pgacw] * ( 1.0 - wk[I_iceflg] );  // c->r by g

//             wk[I_Praut] = std::min( wk[I_Praut], wk[I_dqc_dt] ); 
//             wk[I_Pracw] = std::min( wk[I_Pracw], wk[I_dqc_dt] ); 
//             wk[I_Pihom] = std::min( wk[I_Pihom], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_expice;
//             wk[I_Pihtr] = std::min( wk[I_Pihtr], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_expice;
//             wk[I_Psacw] = std::min( wk[I_Psacw], wk[I_dqc_dt] ) * (        wk[I_iceflg] )            ;
//             wk[I_Psfw ] = std::min( wk[I_Psfw ], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_bergeron;
//             wk[I_Pgacw] = std::min( wk[I_Pgacw], wk[I_dqc_dt] ) * (        wk[I_iceflg] )            ;

//             wk[I_Prevp] = std::min( wk[I_Prevp], wk[I_dqr_dt] );
//             wk[I_Piacr] = std::min( wk[I_Piacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
//             wk[I_Psacr] = std::min( wk[I_Psacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
//             wk[I_Pgacr] = std::min( wk[I_Pgacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
//             wk[I_Pgfrz] = std::min( wk[I_Pgfrz], wk[I_dqr_dt] ) * (        wk[I_iceflg] );

//             wk[I_Pisub] = std::min( wk[I_Pisub], wk[I_dqi_dt] ) * (       wk[I_iceflg] ) * sw_expice;
//             wk[I_Pimlt] = std::min( wk[I_Pimlt], wk[I_dqi_dt] ) * ( 1.0 - wk[I_iceflg] ) * sw_expice;
//             wk[I_Psaut] = std::min( wk[I_Psaut], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
//             wk[I_Praci] = std::min( wk[I_Praci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
//             wk[I_Psaci] = std::min( wk[I_Psaci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
//             wk[I_Psfi ] = std::min( wk[I_Psfi ], wk[I_dqi_dt] ) * (       wk[I_iceflg] ) * sw_bergeron;
//             wk[I_Pgaci] = std::min( wk[I_Pgaci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;

//             wk[I_Pssub] = std::min( wk[I_Pssub], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
//             wk[I_Psmlt] = std::min( wk[I_Psmlt], wk[I_dqs_dt] ) * ( 1.0 - wk[I_iceflg] );
//             wk[I_Pgaut] = std::min( wk[I_Pgaut], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
//             wk[I_Pracs] = std::min( wk[I_Pracs], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
//             wk[I_Pgacs] = std::min( wk[I_Pgacs], wk[I_dqs_dt] );

//             wk[I_Pgsub] = std::min( wk[I_Pgsub], wk[I_dqg_dt] ) * (       wk[I_iceflg] );
//             wk[I_Pgmlt] = std::min( wk[I_Pgmlt], wk[I_dqg_dt] ) * ( 1.0 - wk[I_iceflg] );

//             wk[I_Piacr_s] = ( 1.0 - wk[I_delta1] ) * wk[I_Piacr];
//             wk[I_Piacr_g] = (       wk[I_delta1] ) * wk[I_Piacr];
//             wk[I_Praci_s] = ( 1.0 - wk[I_delta1] ) * wk[I_Praci];
//             wk[I_Praci_g] = (       wk[I_delta1] ) * wk[I_Praci];
//             wk[I_Psacr_s] = (       wk[I_delta2] ) * wk[I_Psacr];
//             wk[I_Psacr_g] = ( 1.0 - wk[I_delta2] ) * wk[I_Psacr];
//             wk[I_Pracs  ] = ( 1.0 - wk[I_delta2] ) * wk[I_Pracs];

//             // [QC]
//             net = 
//                + wk[I_Pimlt]  // [prod] i->c
//                - wk[I_Praut]  // [loss] c->r
//                - wk[I_Pracw]  // [loss] c->r
//                - wk[I_Pihom]  // [loss] c->i
//                - wk[I_Pihtr]  // [loss] c->i
//                - wk[I_Psacw]  // [loss] c->s
//                - wk[I_Psfw ]  // [loss] c->s
//                - wk[I_Pgacw]; // [loss] c->g

//             fac_sw = 0.5 + std::copysign( 0.5, net + EPS ); // if production > loss , fac_sw=1
//             fac    = fac_sw + ( 1.0 - fac_sw ) 
//                      * std::min( -wk[I_dqc_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Pimlt] = wk[I_Pimlt] * fac;
//             wk[I_Praut] = wk[I_Praut] * fac;
//             wk[I_Pracw] = wk[I_Pracw] * fac;
//             wk[I_Pihom] = wk[I_Pihom] * fac;
//             wk[I_Pihtr] = wk[I_Pihtr] * fac;
//             wk[I_Psacw] = wk[I_Psacw] * fac;
//             wk[I_Psfw ] = wk[I_Psfw ] * fac;
//             wk[I_Pgacw] = wk[I_Pgacw] * fac;

//             // [QI]
//             net = 
//                + wk[I_Pigen  ]  // [prod] v->i
//                + wk[I_Pidep  ]  // [prod] v->i
//                + wk[I_Pihom  ]  // [prod] c->i
//                + wk[I_Pihtr  ]  // [prod] c->i
//                - wk[I_Pisub  ]  // [loss] i->v
//                - wk[I_Pimlt  ]  // [loss] i->c
//                - wk[I_Psaut  ]  // [loss] i->s
//                - wk[I_Praci_s]  // [loss] i->s
//                - wk[I_Psaci  ]  // [loss] i->s
//                - wk[I_Psfi   ]  // [loss] i->s
//                - wk[I_Praci_g]  // [loss] i->g
//                - wk[I_Pgaci  ]; // [loss] i->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqi_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Pigen  ] = wk[I_Pigen  ] * fac;
//             wk[I_Pidep  ] = wk[I_Pidep  ] * fac;
//             wk[I_Pihom  ] = wk[I_Pihom  ] * fac;
//             wk[I_Pihtr  ] = wk[I_Pihtr  ] * fac;
//             wk[I_Pisub  ] = wk[I_Pisub  ] * fac;
//             wk[I_Pimlt  ] = wk[I_Pimlt  ] * fac;
//             wk[I_Psaut  ] = wk[I_Psaut  ] * fac;
//             wk[I_Praci_s] = wk[I_Praci_s] * fac;
//             wk[I_Psaci  ] = wk[I_Psaci  ] * fac;
//             wk[I_Psfi   ] = wk[I_Psfi   ] * fac;
//             wk[I_Praci_g] = wk[I_Praci_g] * fac;
//             wk[I_Pgaci  ] = wk[I_Pgaci  ] * fac;


//             // [QR]
//             net = 
//                + wk[I_Praut  ]  // [prod] c->r
//                + wk[I_Pracw  ]  // [prod] c->r
//                + wk[I_Psmlt  ]  // [prod] s->r
//                + wk[I_Pgmlt  ]  // [prod] g->r
//                - wk[I_Prevp  ]  // [loss] r->v
//                - wk[I_Piacr_s]  // [loss] r->s
//                - wk[I_Psacr_s]  // [loss] r->s
//                - wk[I_Piacr_g]  // [loss] r->g
//                - wk[I_Psacr_g]  // [loss] r->g
//                - wk[I_Pgacr  ]  // [loss] r->g
//                - wk[I_Pgfrz  ]; // [loss] r->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac    = fac_sw + ( 1.0 - fac_sw ) 
//                      * std::min( -wk[I_dqr_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Praut  ] = wk[I_Praut  ] * fac;
//             wk[I_Pracw  ] = wk[I_Pracw  ] * fac;
//             wk[I_Psmlt  ] = wk[I_Psmlt  ] * fac;
//             wk[I_Pgmlt  ] = wk[I_Pgmlt  ] * fac;
//             wk[I_Prevp  ] = wk[I_Prevp  ] * fac;
//             wk[I_Piacr_s] = wk[I_Piacr_s] * fac;
//             wk[I_Psacr_s] = wk[I_Psacr_s] * fac;
//             wk[I_Piacr_g] = wk[I_Piacr_g] * fac;
//             wk[I_Psacr_g] = wk[I_Psacr_g] * fac;
//             wk[I_Pgacr  ] = wk[I_Pgacr  ] * fac;
//             wk[I_Pgfrz  ] = wk[I_Pgfrz  ] * fac;

//             // [QV]
//             net =
//                + wk[I_Prevp]  // [prod] r->v
//                + wk[I_Pisub]  // [prod] i->v
//                + wk[I_Pssub]  // [prod] s->v
//                + wk[I_Pgsub]  // [prod] g->v
//                - wk[I_Pigen]  // [loss] v->i
//                - wk[I_Pidep]  // [loss] v->i
//                - wk[I_Psdep]  // [loss] v->s
//                - wk[I_Pgdep]; // [loss] v->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqv_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Prevp] = wk[I_Prevp] * fac;
//             wk[I_Pisub] = wk[I_Pisub] * fac;
//             wk[I_Pssub] = wk[I_Pssub] * fac;
//             wk[I_Pgsub] = wk[I_Pgsub] * fac;
//             wk[I_Pigen] = wk[I_Pigen] * fac;
//             wk[I_Pidep] = wk[I_Pidep] * fac;
//             wk[I_Psdep] = wk[I_Psdep] * fac;
//             wk[I_Pgdep] = wk[I_Pgdep] * fac;

//             // [QS]
//             net =
//                + wk[I_Psdep  ]  // [prod] v->s
//                + wk[I_Psacw  ]  // [prod] c->s
//                + wk[I_Psfw   ]  // [prod] c->s
//                + wk[I_Piacr_s]  // [prod] r->s
//                + wk[I_Psacr_s]  // [prod] r->s
//                + wk[I_Psaut  ]  // [prod] i->s
//                + wk[I_Praci_s]  // [prod] i->s
//                + wk[I_Psaci  ]  // [prod] i->s
//                + wk[I_Psfi   ]  // [prod] i->s
//                - wk[I_Pssub  ]  // [loss] s->v
//                - wk[I_Psmlt  ]  // [loss] s->r
//                - wk[I_Pgaut  ]  // [loss] s->g
//                - wk[I_Pracs  ]  // [loss] s->g
//                - wk[I_Pgacs  ]; // [loss] s->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqs_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Psdep  ] = wk[I_Psdep  ] * fac;
//             wk[I_Psacw  ] = wk[I_Psacw  ] * fac;
//             wk[I_Psfw   ] = wk[I_Psfw   ] * fac;
//             wk[I_Piacr_s] = wk[I_Piacr_s] * fac;
//             wk[I_Psacr_s] = wk[I_Psacr_s] * fac;
//             wk[I_Psaut  ] = wk[I_Psaut  ] * fac;
//             wk[I_Praci_s] = wk[I_Praci_s] * fac;
//             wk[I_Psaci  ] = wk[I_Psaci  ] * fac;
//             wk[I_Psfi   ] = wk[I_Psfi   ] * fac;
//             wk[I_Pssub  ] = wk[I_Pssub  ] * fac;
//             wk[I_Psmlt  ] = wk[I_Psmlt  ] * fac;
//             wk[I_Pgaut  ] = wk[I_Pgaut  ] * fac;
//             wk[I_Pracs  ] = wk[I_Pracs  ] * fac;
//             wk[I_Pgacs  ] = wk[I_Pgacs  ] * fac;

//             // [QG]
//             net =
//                + wk[I_Pgdep  ]  // [prod] v->g
//                + wk[I_Pgacw  ]  // [prod] c->g
//                + wk[I_Piacr_g]  // [prod] r->g
//                + wk[I_Psacr_g]  // [prod] r->g
//                + wk[I_Pgacr  ]  // [prod] r->g
//                + wk[I_Pgfrz  ]  // [prod] r->g
//                + wk[I_Praci_g]  // [prod] i->g
//                + wk[I_Pgaci  ]  // [prod] i->g
//                + wk[I_Pgaut  ]  // [prod] s->g
//                + wk[I_Pracs  ]  // [prod] s->g
//                + wk[I_Pgacs  ]  // [prod] s->g
//                - wk[I_Pgsub  ]  // [loss] g->v
//                - wk[I_Pgmlt  ]; // [loss] g->r

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqg_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Pgdep  ] = wk[I_Pgdep  ] * fac;
//             wk[I_Pgacw  ] = wk[I_Pgacw  ] * fac;
//             wk[I_Piacr_g] = wk[I_Piacr_g] * fac;
//             wk[I_Psacr_g] = wk[I_Psacr_g] * fac;
//             wk[I_Pgacr  ] = wk[I_Pgacr  ] * fac;
//             wk[I_Pgfrz  ] = wk[I_Pgfrz  ] * fac;
//             wk[I_Praci_g] = wk[I_Praci_g] * fac;
//             wk[I_Pgaci  ] = wk[I_Pgaci  ] * fac;
//             wk[I_Pgaut  ] = wk[I_Pgaut  ] * fac;
//             wk[I_Pracs  ] = wk[I_Pracs  ] * fac;
//             wk[I_Pgacs  ] = wk[I_Pgacs  ] * fac;
//             wk[I_Pgsub  ] = wk[I_Pgsub  ] * fac;
//             wk[I_Pgmlt  ] = wk[I_Pgmlt  ] * fac;

//             qc_t = + wk[I_Pimlt]  // [prod] i->c
//                    - wk[I_Praut]  // [loss] c->r
//                    - wk[I_Pracw]  // [loss] c->r
//                    - wk[I_Pihom]  // [loss] c->i
//                    - wk[I_Pihtr]  // [loss] c->i
//                    - wk[I_Psacw]  // [loss] c->s
//                    - wk[I_Psfw ]  // [loss] c->s
//                    - wk[I_Pgacw]; // [loss] c->g
 
//             qr_t = + wk[I_Praut  ]  // [prod] c->r
//                    + wk[I_Pracw  ]  // [prod] c->r
//                    + wk[I_Psmlt  ]  // [prod] s->r
//                    + wk[I_Pgmlt  ]  // [prod] g->r
//                    - wk[I_Prevp  ]  // [loss] r->v
//                    - wk[I_Piacr_s]  // [loss] r->s
//                    - wk[I_Psacr_s]  // [loss] r->s
//                    - wk[I_Piacr_g]  // [loss] r->g
//                    - wk[I_Psacr_g]  // [loss] r->g
//                    - wk[I_Pgacr  ]  // [loss] r->g
//                    - wk[I_Pgfrz  ]; // [loss] r->g

//             qi_t = + wk[I_Pigen  ]  // [prod] v->i
//                    + wk[I_Pidep  ]  // [prod] v->i
//                    + wk[I_Pihom  ]  // [prod] c->i
//                    + wk[I_Pihtr  ]  // [prod] c->i
//                    - wk[I_Pisub  ]  // [loss] i->v
//                    - wk[I_Pimlt  ]  // [loss] i->c
//                    - wk[I_Psaut  ]  // [loss] i->s
//                    - wk[I_Praci_s]  // [loss] i->s
//                    - wk[I_Psaci  ]  // [loss] i->s
//                    - wk[I_Psfi   ]  // [loss] i->s
//                    - wk[I_Praci_g]  // [loss] i->g
//                    - wk[I_Pgaci  ]; // [loss] i->g

//             qs_t = + wk[I_Psdep  ]  // [prod] v->s
//                    + wk[I_Psacw  ]  // [prod] c->s
//                    + wk[I_Psfw   ]  // [prod] c->s
//                    + wk[I_Piacr_s]  // [prod] r->s
//                    + wk[I_Psacr_s]  // [prod] r->s
//                    + wk[I_Psaut  ]  // [prod] i->s
//                    + wk[I_Praci_s]  // [prod] i->s
//                    + wk[I_Psaci  ]  // [prod] i->s
//                    + wk[I_Psfi   ]  // [prod] i->s
//                    - wk[I_Pssub  ]  // [loss] s->v
//                    - wk[I_Psmlt  ]  // [loss] s->r
//                    - wk[I_Pgaut  ]  // [loss] s->g
//                    - wk[I_Pracs  ]  // [loss] s->g
//                    - wk[I_Pgacs  ]; // [loss] s->g

//             qg_t = + wk[I_Pgdep  ]  // [prod] v->g
//                    + wk[I_Pgacw  ]  // [prod] c->g
//                    + wk[I_Piacr_g]  // [prod] r->g
//                    + wk[I_Psacr_g]  // [prod] r->g
//                    + wk[I_Pgacr  ]  // [prod] r->g
//                    + wk[I_Pgfrz  ]  // [prod] r->g
//                    + wk[I_Praci_g]  // [prod] i->g
//                    + wk[I_Pgaci  ]  // [prod] i->g
//                    + wk[I_Pgaut  ]  // [prod] s->g
//                    + wk[I_Pracs  ]  // [prod] s->g
//                    + wk[I_Pgacs  ]  // [prod] s->g
//                    - wk[I_Pgsub  ]  // [loss] g->v
//                    - wk[I_Pgmlt  ]; // [loss] g->r

//             qc_t = std::max( qc_t, -wk[I_dqc_dt] );
//             qr_t = std::max( qr_t, -wk[I_dqr_dt] );
//             qi_t = std::max( qi_t, -wk[I_dqi_dt] );
//             qs_t = std::max( qs_t, -wk[I_dqs_dt] );
//             qg_t = std::max( qg_t, -wk[I_dqg_dt] );

//             qv_t = - ( qc_t 
//                      + qr_t 
//                      + qi_t 
//                      + qs_t 
//                      + qg_t );

//             drhogqv(k,ij) = rhog(k,ij) * qv_t * dt;
//             drhogqc(k,ij) = rhog(k,ij) * qc_t * dt;
//             drhogqr(k,ij) = rhog(k,ij) * qr_t * dt;
//             drhogqi(k,ij) = rhog(k,ij) * qi_t * dt;
//             drhogqs(k,ij) = rhog(k,ij) * qs_t * dt;
//             drhogqg(k,ij) = rhog(k,ij) * qg_t * dt;

//             Vt(I_QR,k,ij) = Vtr;
//             Vt(I_QI,k,ij) = Vti;
//             Vt(I_QS,k,ij) = Vts;
//             Vt(I_QG,k,ij) = Vtg;

//             //--- Effective Radius of Liquid Water
//             xf_qc = rhoqc / Nc(k,ij) * 1.0E-9;                            // mean mass   of qc [kg]
//             rf_qc = std::pow(( coef_xf * xf_qc + EPS ), 0.33333333);       // mean radius of qc [m]
//             r2_qc = coef_dgam * (rf_qc * rf_qc) * ( Nc(k,ij) * 1.0E+6 );  // r^2 moment  of qc
//             r2_qr = 0.25 * N0r * GAM_3 * RLMDr_3;                          // r^2 moment  of qr
//             r3_qc = coef_xf * dens * qc;                                   // r^3 moment  of qc
//             r3_qr = coef_xf * dens * qr;                                   // r^3 moment  of qr

//             zerosw = 0.5 - std::copysign(0.5, (r2_qc + r2_qr) - r2_min );
//             rceff(k,ij) = ( r3_qc + r3_qr ) / ( r2_qc + r2_qr + zerosw ) * ( 1.0 - zerosw );

//             sw = ( 0.5 + std::copysign(0.5, (qc + qr) - q_min ) ) * zerosw; // if qc+qr > 0.01[g/kg], sw=1

//             rctop(0,ij) = (       sw ) * rceff(k,ij)
//                          + ( 1.0 - sw ) * UNDEF;
//             tctop(0,ij) = (       sw ) * temp
//                          + ( 1.0 - sw ) * UNDEF;

//             zerosw = 0.5 - std::copysign( 0.5, r2_qc - r2_min );
//             rceff_cld(k,ij) = r3_qc / ( r2_qc + zerosw ) * ( 1.0 - zerosw );

//             sw = ( 0.5 + std::copysign( 0.5, qc - q_min ) ) * zerosw; // if qc > 0.01[g/kg], sw=1

//             rctop_cld(0,ij) = (       sw ) * rceff_cld(k,ij)
//                              + ( 1.0 - sw ) * UNDEF;
//             tctop_cld(0,ij) = (       sw ) * temp
//                              + ( 1.0 - sw ) * UNDEF;
//         }
//     }
//     // ----------------------- Big Loop End Here ------------------------

//     // update rhogq
//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), 
//     KOKKOS_LAMBDA(const size_t k, const size_t ij){
//         rhogq(I_QV,k,ij) += drhogqv(k,ij);
//         rhogq(I_QC,k,ij) += drhogqc(k,ij);
//         rhogq(I_QR,k,ij) += drhogqr(k,ij);
//         rhogq(I_QI,k,ij) += drhogqi(k,ij);
//         rhogq(I_QS,k,ij) += drhogqs(k,ij);
//         rhogq(I_QG,k,ij) += drhogqg(k,ij);
//     });

//     // update rhoge
//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), 
//     KOKKOS_LAMBDA(const size_t k, const size_t ij){
//         rhoge(k,ij) = rhoge(k,ij) - LHV * drhogqv(k,ij)
//                                   + LHF * drhogqi(k,ij)
//                                   + LHF * drhogqs(k,ij)
//                                   + LHF * drhogqg(k,ij);
//     });
//     Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
//         rceff    (kmin-1,ij) = 0.0;
//         rceff_cld(kmin-1,ij) = 0.0;
//         rceff    (kmax+1,ij) = 0.0;
//         rceff_cld(kmax+1,ij) = 0.0;

//         double sw = ( 0.5 + std::copysign(0.5, tctop(0,ij) - TEM00 ) ); // if T > 0C, sw=1

//         rwtop(0,ij) =  (       sw ) * rctop(0,ij)
//                      + ( 1.0 - sw ) * UNDEF;

//         sw = ( 0.5 + std::copysign(0.5, rctop_cld(0,ij) - TEM00 ) ); // if T > 0C, sw=1

//         rwtop_cld(0,ij) =  (       sw ) * rctop_cld(0,ij)
//                          + ( 1.0 - sw ) * UNDEF;
//     });

//     //---< preciptation（降水量） >---
//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{nqmax,ijdim}), 
//     KOKKOS_LAMBDA(const size_t nq, const size_t ij){
//         Vt(nq,kmin-1,ij) = 0.0;
//         Vt(nq,kmax+1,ij) = 0.0;
//     });
//     //--- update mass concentration（質量濃度）
//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<3>>({NQW_STR,kmin,IDX_ZERO},{NQW_END+1,kmax+1,ijdim}), 
//     KOKKOS_LAMBDA(const size_t nq, const size_t k, const size_t ij){
//         q(nq,k,ij) = rhogq(nq,k,ij) / rhog(k,ij);
//     });

//     THRMDYN_qd(q, qd);
//     THRMDYN_cv(qd, q, cva);

//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,0},{kmax+1,ijdim}), 
//     KOKKOS_LAMBDA(const size_t k, const size_t ij){
//             tem (k,ij) = rhoge(k,ij) / ( rhog(k,ij) * cva(k,ij) );
//             rgs (k,ij) = gam2 (k,ij) / gsgam2 (k,ij);
//             rgsh(k,ij) = gam2h(k,ij) / gsgam2h(k,ij);
//     });

//     // for(int nq = 0; nq < nqmax; nq++)
//     //     preciptation_flag[nq] = false;
//     preciptation_flag[I_QR] = true;
//     preciptation_flag[I_QI] = true;
//     preciptation_flag[I_QS] = true;
//     preciptation_flag[I_QG] = true;
//     if(precip_transport_type == "3WATER")
//     {
//         preciptation_flag[I_QI] = false;
//     }

//     precip_transport_new( rhog,
//                           rhogvx,
//                           rhogvy,
//                           rhogvz,
//                           rhogw,
//                           rhoge,
//                           rhogq,
//                           rho,
//                           tem,
//                           pre,
//                           vx,
//                           vy,
//                           vz,
//                           w,
//                           q,
//                           qd,
//                           z,
//                           Vt,
//                           preciptation_flag,
//                           precip,precip_rhoe,
//                           precip_lh_heat,
//                           precip_rhophi,
//                           precip_rhokin,
//                           gprec,
//                           gsgam2,
//                           gsgam2h,
//                           rgs,
//                           rgsh,
//                           ix,
//                           iy,
//                           iz,
//                           jx,
//                           jy,
//                           jz,
//                           dt,
//                           nullptr); // precip_trc**

//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
//     KOKKOS_LAMBDA(const size_t k, const size_t ij){
//             ml_Pconv(k,ij) = q(I_QV,k,ij);
//             ml_Pconw(k,ij) = q(I_QC,k,ij);
//             ml_Pconi(k,ij) = q(I_QI,k,ij);
//     });

//     bool ice_adjust;
//     if(OPT_EXPLICIT_ICEGEN)
//     {
//         ice_adjust = false;
//         SATURATION_Setrange(100.0, 90.0);

//         SATURATION_adjustment(rhog, rhoge, rhogq, tem, q, qd, gsgam2, ice_adjust);
//     }
//     else
//     {
//         ice_adjust = true;
//         SATURATION_Setrange(273.16, 233.16);

//         SATURATION_adjustment(rhog, rhoge, rhogq, tem , q, qd, gsgam2, ice_adjust);
//     }

//     Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
//     KOKKOS_LAMBDA(const size_t k, const size_t ij){
//         ml_Pconv(k,ij) = q(I_QV,k,ij) - ml_Pconv(k,ij);
//         ml_Pconw(k,ij) = q(I_QC,k,ij) - ml_Pconw(k,ij);
//         ml_Pconi(k,ij) = q(I_QI,k,ij) - ml_Pconi(k,ij);
//     });

//     negative_filter(rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2);

//     // End of nsw6
// }

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
//             )
// {
// #ifdef ENABLE_DEBUG
//     std::cout << __PRETTY_FUNCTION__ << std::endl;
// #endif

//     // working
//     double drhogqv[kdim][ijdim];
//     double drhogqc[kdim][ijdim];
//     double drhogqi[kdim][ijdim];
//     double drhogqr[kdim][ijdim];
//     double drhogqs[kdim][ijdim];
//     double drhogqg[kdim][ijdim];

//     double psatl[kdim][ijdim];
//     double qsatl[kdim][ijdim];                //< saturated water vapor for liquid water [kg/kg]
//     double psati[kdim][ijdim];
//     double qsati[kdim][ijdim];                //< saturated water vapor for ice water    [kg/kg]
//     double Nc   [kdim][ijdim];                //< Number concentration of cloud water [1/cc]

//     double dens    ;                         //< density
//     double temp    ;                         //< T [K]
//     double qv      ;                         //< mixing ratio of water vapor  [kg/kg]
//     double qc      ;                         //< mixing ratio of liquid water [kg/kg]
//     double qr      ;                         //< mixing ratio of rain         [kg/kg]
//     double qi      ;                         //< mixing ratio of ice water    [kg/kg]
//     double qs      ;                         //< mixing ratio of snow         [kg/kg]
//     double qg      ;                         //< mixing ratio of graupel      [kg/kg]
//     double qv_t    ;                         //< tendency     of water vapor  [kg/kg/s]
//     double qc_t    ;                         //< tendency     of liquid water [kg/kg/s]
//     double qr_t    ;                         //< tendency     of rain         [kg/kg/s]
//     double qi_t    ;                         //< tendency     of ice water    [kg/kg/s]
//     double qs_t    ;                         //< tendency     of snow         [kg/kg/s]
//     double qg_t    ;                         //< tendency     of graupel      [kg/kg/s]
//     double Sliq    ;                         //< saturated ratio S for liquid water [0-1]
//     double Sice    ;                         //< saturated ratio S for ice water    [0-1]
//     double Rdens   ;                         //< 1 / density
//     double rho_fact;                         //< density factor
//     double temc    ;                         //< T - T0 [K]

//     double RLMDr, RLMDr_2, RLMDr_3        ;
//     double RLMDs, RLMDs_2, RLMDs_3        ;
//     double RLMDg, RLMDg_2, RLMDg_3        ;
//     double RLMDr_1br, RLMDr_2br, RLMDr_3br;
//     double RLMDs_1bs, RLMDs_2bs, RLMDs_3bs;
//     double RLMDr_dr, RLMDr_3dr, RLMDr_5dr ;
//     double RLMDs_ds, RLMDs_3ds, RLMDs_5ds ;
//     double RLMDg_dg, RLMDg_3dg, RLMDg_5dg ;
//     double RLMDr_7                        ;
//     double RLMDr_6dr                      ;

//     //---< Roh and Satoh (2014) >---
//     double tems, Xs2                     ;
//     double MOMs_0, MOMs_1, MOMs_2        ;
//     double MOMs_0bs, MOMs_1bs, MOMs_2bs  ;
//     double MOMs_2ds, MOMs_5ds_h, RMOMs_Vt;
//     double coef_at[4]                    ;
//     double coef_bt[4]                    ;
//     double loga_, b_, nm                 ;

//     double Vti, Vtr, Vts, Vtg   ;           //< terminal velocity
//     double Esi_mod, Egs_mod     ;           //< modified accretion efficiency
//     double rhoqc                ;           //< rho * qc
//     double Pracw_orig,  Pracw_kk;           //< accretion       term by orig  & k-k scheme
//     double Praut_berry, Praut_kk;           //< auto-conversion term by berry & k-k scheme
//     double Dc                   ;           //< relative variance
//     double betai, betas         ;           //< sticky parameter for auto-conversion
//     double Ka                   ;           //< thermal diffusion coefficient of air
//     double Kd                   ;           //< diffusion coefficient of water vapor in air
//     double Nu                   ;           //< kinematic viscosity of air
//     double Glv, Giv, Gil        ;           //< thermodynamic function
//     double ventr, vents, ventg  ;           //< ventilation factor
//     double net, fac, fac_sw     ;
//     double zerosw, tmp          ;

//     //---< Bergeron process >---
//     double sw_bergeron     ;                 //< if 0C<T<30C, sw=1
//     double a1 [kdim][ijdim];                 //<
//     double a2 [kdim][ijdim];                 //<
//     double ma2[kdim][ijdim];                 //< 1-a2
//     double dt1             ;                 //< time during which the an ice particle of 40um grows to 50um
//     double Ni50            ;                 //< number concentration of ice particle of 50um

//     //---< Explicit ice generation >---
//     double sw, rhoqi, XNi, XMi, Di, Ni0, Qi0;

//     //---< Effective radius >---
//     constexpr double r2_min = 1.0E-10;
//     constexpr double q_min  = 1.0E-5;
//     double xf_qc;                 //< mean mass   of qc [kg]
//     double rf_qc;                 //< mean radius of qc [m]
//     double r2_qc;                 //< r^2 moment  of qc
//     double r2_qr;                 //< r^2 moment  of qr
//     double r3_qc;                 //< r^3 moment  of qc
//     double r3_qr;                 //< r^3 moment  of qr
//     double dgamma_a, GAM_dgam23, GAM_dgam;
//     double coef_dgam, coef_xf;

//     //---< Precipitation >---
//     bool preciptation_flag[nqmax] = {false};

//     double Vt    [nqmax][kdim][ijdim];
//     double cva   [kdim][ijdim];
//     double rgs   [kdim][ijdim];
//     double rgsh  [kdim][ijdim];

//     double wk    [wk_nmax];
//     //double ml_wk   [wk_nmax][kdim][ijdim]; tentatively remove due to the stack overflow
//     double ml_Pconv[kdim][ijdim];
//     double ml_Pconw[kdim][ijdim];
//     double ml_Pconi[kdim][ijdim];

//     double UNDEF, EPS, PI, Rvap, LHV0, LHS0, LHF0, PRE00;

//     // constexpr int simdlen = 8;
//     // int blk, vec, veclen;

//     //---------------------------------------------------------------------------
//     UNDEF = CONST_UNDEF;
//     EPS   = CONST_EPS  ;
//     PI    = CONST_PI   ;
//     Rvap  = CONST_Rvap ;
//     LHV0  = CONST_LHV0 ;
//     LHS0  = CONST_LHS0 ;
//     LHF0  = CONST_LHF0 ;
//     PRE00 = CONST_Pstd ;

//     negative_filter(rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2);

//     if(OPT_INDIR) // aerosol indirect effect
//     {
//         for(int k = 0; k < kdim; k++)
//         {
//             for(int ij = 0; ij < ijdim; ij++)
//             {
//                 Nc[k][ij] = std::max( UNCCN[k][ij] * 1.0E-6, Nc_def );
//             }
//         }
//     }
//     else
//     {
//         for(int k = 0; k < kdim; k++)
//         {
//             for(int ij = 0; ij < ijdim; ij++)
//             {
//                 Nc[k][ij] = Nc_def;
//             }
//         }
//     }

//     // saturation water contensts
//     SATURATION_psat_liq(tem, psatl);
//     SATURATION_psat_ice(tem, psati);

//     for(int k = 0; k < kdim; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             qsatl[k][ij] = psatl[k][ij] / ( rho[k][ij] * Rvap * tem[k][ij] );
//             qsati[k][ij] = psati[k][ij] / ( rho[k][ij] * Rvap * tem[k][ij] );
//         }
//     }

//     // Bergeron process parameters
//     Bergeron_param(tem, a1, a2, ma2);

//     // work for Effective Radius of Liquid Water
//     // dgamma_a is parameter of Gamma-Dist.
//     // "Berry and Reinhardt(1974-a) eq.(16-a,b)"
//     // 0.38 * {var1/2x} =~ {var1/2 r}
//     //         dgamma_a = 1.0_RP /{var x}
//     Dc = 0.146 - 5.964E-2 * std::log(Nc_def / 2000.0);
//     dgamma_a = 0.1444 / (Dc * Dc); // Dc**2 in fortran
//     GAM_dgam   = MISC_gammafunc(dgamma_a);
//     GAM_dgam23 = MISC_gammafunc(dgamma_a + 2.0/3.0);
//     coef_dgam  = GAM_dgam23 / GAM_dgam * std::pow(dgamma_a, (-2.0/3.0));
//     coef_xf    = 3.0 / 4.0 / PI / rho_w;

//     // !! Big loop start here !!
//     for(int k = kmin; k <= kmax; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             dens = rho[k][ij];
//             temp = tem[k][ij];
//             qv   = std::max( q[I_QV][k][ij], 0.0 );
//             qc   = std::max( q[I_QC][k][ij], 0.0 );
//             qr   = std::max( q[I_QR][k][ij], 0.0 );
//             qi   = std::max( q[I_QI][k][ij], 0.0 );
//             qs   = std::max( q[I_QS][k][ij], 0.0 );
//             qg   = std::max( q[I_QG][k][ij], 0.0 );

//             // saturation ration S
//             Sliq = qv / (std::max(qsatl[k][ij], EPS));
//             Sice = qv / (std::max(qsati[k][ij], EPS));

//             Rdens = 1.0 / dens;
//             rho_fact = std::sqrt(dens00 * Rdens);
//             temc = temp - TEM00;

//             wk[I_delta1] = ( 0.5 + std::copysign(0.5, qr - 1.0E-4) );

//             wk[I_delta2] = ( 0.5 + std::copysign(0.5, 1.0E-4 - qr) )
//                             * ( 0.5 + std::copysign(0.5, 1.0E-4 - qs) );
            
//             wk[I_spsati] = 0.5 + std::copysign(0.5, Sice - 1.0);

//             wk[I_iceflg] = 0.5 - std::copysign(0.5, temc); // 0: warm, 1: ice

//             wk[I_dqv_dt] = qv / dt;
//             wk[I_dqc_dt] = qc / dt;
//             wk[I_dqr_dt] = qr / dt;
//             wk[I_dqi_dt] = qi / dt;
//             wk[I_dqs_dt] = qs / dt;
//             wk[I_dqg_dt] = qg / dt;

//             sw_bergeron = ( 0.5 + std::copysign(0.5, temc + 30.0) )
//                             * ( 0.5 + std::copysign(0.5, 0.0 - temc) )
//                             * ( 1.0 - sw_expice );
            
//             // slope parameter lambda (Rain)
//             zerosw = 0.5 - std::copysign(0.5, qr - 1.0E-12);
//             RLMDr = std::sqrt( std::sqrt( dens * qr / (Ar * N0r * GAM_1br) + zerosw ) ) * (1.0 - zerosw);

//             RLMDr_dr  = std::sqrt(RLMDr);   // **Dr
//             RLMDr_2   = std::pow(RLMDr, 2);
//             RLMDr_3   = std::pow(RLMDr, 3);
//             RLMDr_7   = std::pow(RLMDr, 7);
//             RLMDr_1br = std::pow(RLMDr, 4); // (1+Br)
//             RLMDr_2br = std::pow(RLMDr, 5); // (2+Br)
//             RLMDr_3br = std::pow(RLMDr, 6); // (3+Br)
//             RLMDr_3dr = std::pow(RLMDr, 3) * RLMDr_dr;
//             RLMDr_5dr = std::pow(RLMDr, 5) * RLMDr_dr;
//             RLMDr_6dr = std::pow(RLMDr, 6) * RLMDr_dr;

//             // slope parameter lambda (Snow)
//             zerosw = 0.5 - std::copysign(0.5, qs - 1.0E-12);
//             RLMDs  = std::sqrt( std::sqrt( dens * qs / ( As * N0s * GAM_1bs ) + zerosw ) ) * ( 1.0 - zerosw );

//             RLMDs_ds  = std::sqrt( std::sqrt(RLMDs) ); // **Ds
//             RLMDs_2   = std::pow(RLMDs, 2);
//             RLMDs_3   = std::pow(RLMDs, 3);
//             RLMDs_1bs = std::pow(RLMDs, 4); // (1+Bs)
//             RLMDs_2bs = std::pow(RLMDs, 5); // (2+Bs)
//             RLMDs_3bs = std::pow(RLMDs, 6); // (3+Bs)
//             RLMDs_3ds = std::pow(RLMDs, 3) * RLMDs_ds;
//             RLMDs_5ds = std::pow(RLMDs, 5) * RLMDs_ds;

//             MOMs_0     = N0s * GAM       * RLMDs    ;       // Ns * 0th moment
//             MOMs_1     = N0s * GAM_2     * RLMDs_2  ;       // Ns * 1st moment
//             MOMs_2     = N0s * GAM_3     * RLMDs_3  ;       // Ns * 2nd moment
//             MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs;       // Ns * 0+bs
//             MOMs_1bs   = N0s * GAM_2bs   * RLMDs_2bs;       // Ns * 1+bs
//             MOMs_2bs   = N0s * GAM_3bs   * RLMDs_3bs;       // Ns * 2+bs
//             MOMs_2ds   = N0s * GAM_3ds   * RLMDs_3ds;       // Ns * 2+ds
//             MOMs_5ds_h = N0s * GAM_5ds_h * std::sqrt(RLMDs_5ds); // Ns * (5+ds)/2
//             RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds;

//             //---< modification by Roh and Satoh (2014) >---

//             // bimodal size distribution of snow
//             Xs2 = dens * qs / As;
//             zerosw = 0.5 - std::copysign(0.5, Xs2 - 1.0E-12);

//             tems = std::min(-0.1, temc);
//             coef_at[0] = coef_a[0] + tems * ( coef_a[1] + tems * ( coef_a[4] + tems * coef_a[8] ) );
//             coef_at[1] = coef_a[2] + tems * ( coef_a[3] + tems *   coef_a[6] );
//             coef_at[2] = coef_a[5] + tems *   coef_a[7];
//             coef_at[3] = coef_a[9];
//             coef_bt[0] = coef_b[0] + tems * ( coef_b[1] + tems * ( coef_b[4] + tems * coef_b[8] ) );
//             coef_bt[1] = coef_b[2] + tems * ( coef_b[3] + tems *   coef_b[6] );
//             coef_bt[2] = coef_b[5] + tems *   coef_b[7];
//             coef_bt[3] = coef_b[9];

//             // 0th moment
//             loga_ = coef_at[0];
//             b_    = coef_bt[0];
//             MOMs_0 = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                      + ( 1.0 - sw_roh2014 ) * MOMs_0;
//             // 1st moment
//             nm = 1.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_1 = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw ) * b_) * ( 1.0 - zerosw )
//                      + ( 1.0 - sw_roh2014 ) * MOMs_1;
//             // 2nd moment
//             MOMs_2 = sw_roh2014 * Xs2 + (1.0 - sw_roh2014) * MOMs_2;
//             // 0 + Bs(=2) moment
//             nm = 2.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_0bs = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_0bs;
//             // 1 + Bs(=2) moment
//             nm = 3.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_1bs = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_1bs;
//             // 2 + Bs(=2) moment
//             nm = 4.0;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_2bs = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_2bs;
//             // 2 + Ds(=0.25) moment
//             nm = 2.25;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_2ds = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                        + ( 1.0 - sw_roh2014 ) * MOMs_2ds;
//             // ( 3 + Ds(=0.25) ) / 2  moment
//             nm = 1.625;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             MOMs_5ds_h = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw )
//                          + ( 1.0 - sw_roh2014 ) * MOMs_5ds_h;
//             // Bs(=2) + Ds(=0.25) moment
//             nm = 2.25;
//             loga_ = coef_at[0] + nm * ( coef_at[1] + nm * ( coef_at[2] + nm * coef_at[3] ) );
//             b_    = coef_bt[0] + nm * ( coef_bt[1] + nm * ( coef_bt[2] + nm * coef_bt[3] ) );
//             RMOMs_Vt = sw_roh2014 * std::exp(ln10 * loga_) * std::exp(std::log(Xs2 + zerosw) * b_) * ( 1.0 - zerosw ) / (MOMs_0bs + zerosw)
//                        + ( 1.0 - sw_roh2014 ) * RMOMs_Vt;

//             // slope parameter lambda (Graupel)
//             zerosw = 0.5 - std::copysign(0.5, qg - 1.0E-12);
//             RLMDg  = std::sqrt( std::sqrt( dens * qg / ( Ag * N0g * GAM_1bg ) + zerosw ) ) * ( 1.0 - zerosw );
//             RLMDg_dg  = std::sqrt( RLMDg );       // **Dg
//             RLMDg_2   = std::pow(RLMDg, 2);
//             RLMDg_3   = std::pow(RLMDg, 3);
//             RLMDg_3dg = std::pow(RLMDg, 3) * RLMDg_dg;
//             RLMDg_5dg = std::pow(RLMDg, 5) * RLMDg_dg;

//             wk[I_RLMDr] = RLMDr;
//             wk[I_RLMDs] = RLMDs;
//             wk[I_RLMDg] = RLMDg;

//             //---< terminal velocity >---
//             zerosw = 0.5 - std::copysign(0.5, qi - 1.0E-8 );
//             Vti = ( 1.0 - sw_constVti ) * (-3.29) * std::exp( std::log( dens * qi + zerosw ) * 0.16 ) * ( 1.0 - zerosw ) 
//                   + ( sw_constVti ) * (-CONST_Vti);
//             Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr;
//             Vts = -Cs * rho_fact * RMOMs_Vt;
//             Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg;

//             //---< Nucleation >---
//             // [Pigen] ice nucleation
//             Ni0 = std::max( std::exp(-0.1 * temc), 1.0 ) * 1000.0;
//             Qi0 = 4.92E-11 * std::exp( std::log(Ni0) * 1.33 ) * Rdens;

//             wk[I_Pigen] = std::max( std::min( Qi0 - qi, qv - qsati[k][ij] ), 0.0 ) / dt;

//             //---< Accretion >---
//             Esi_mod = std::min( Esi, Esi * std::exp( gamma_sacr * temc ) );
//             Egs_mod = std::min( Egs, Egs * std::exp( gamma_gacs * temc ) );

//             // [Pracw] accretion rate of cloud water by rain
//             Pracw_orig = qc * 0.25 * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact;

//             zerosw     = 0.5 - std::copysign(0.5, qc * qr - 1.0E-12 );
//             Pracw_kk   = 67.0 * std::exp( std::log( qc * qr + zerosw ) * 1.15 ) * ( 1.0 - zerosw ); // eq.(33) in KK(2000)

//             // switch orig / k-k scheme
//             wk[I_Pracw] = ( 1.0 - sw_kk2000 ) * Pracw_orig
//                         + (       sw_kk2000 ) * Pracw_kk;

//             // [Psacw] accretion rate of cloud water by snow
//             wk[I_Psacw] = qc * 0.25 * PI * Esw * Cs * MOMs_2ds * rho_fact;

//             // [Pgacw] accretion rate of cloud water by graupel
//             wk[I_Pgacw] = qc * 0.25 * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact;

//             // [Praci] accretion rate of cloud ice by rain
//             wk[I_Praci] = qi * 0.25 * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact;

//             // [Psaci] accretion rate of cloud ice by snow
//             wk[I_Psaci] = qi * 0.25 * PI * Esi_mod * Cs * MOMs_2ds * rho_fact;

//             // [Pgaci] accretion rate of cloud ice by grupel
//             wk[I_Pgaci] = qi * 0.25 * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact;

//             // [Piacr] accretion rate of rain by cloud ice
//             wk[I_Piacr] = qi * Ar / mi * 0.25 * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact;

//             // [Psacr] accretion rate of rain by snow
//             wk[I_Psacr] = Ar * 0.25 * PI * Rdens * Esr * N0r * std::abs(Vtr - Vts)
//                              * ( GAM_1br * RLMDr_1br * MOMs_2          
//                                  + 2.0 * GAM_2br * RLMDr_2br * MOMs_1          
//                                  + GAM_3br * RLMDr_3br * MOMs_0 );

//             // [Pgacr] accretion rate of rain by graupel
//             wk[I_Pgacr] = Ar * 0.25 * PI * Rdens * Egr * N0g * N0r * std::abs(Vtg - Vtr)
//                              * ( GAM_1br * RLMDr_1br * GAM_3 * RLMDg_3
//                                  + 2.0 * GAM_2br * RLMDr_2br * GAM_2 * RLMDg_2
//                                  + GAM_3br * RLMDr_3br * GAM * RLMDg );

//             // [Pracs] accretion rate of snow by rain
//             wk[I_Pracs] = As * 0.25 * PI * Rdens * Esr *  N0r * std::abs(Vtr - Vts)
//                              * ( MOMs_0bs * GAM_3 * RLMDr_3
//                                  + 2.0 * MOMs_1bs * GAM_2 * RLMDr_2
//                                  + MOMs_2bs * GAM * RLMDr );

//             // [Pgacs] accretion rate of snow by graupel
//             wk[I_Pgacs] = As * 0.25 * PI * Rdens * Egs_mod * N0g * std::abs(Vtg - Vts)
//                              * ( MOMs_0bs * GAM_3 * RLMDg_3
//                                  + 2.0 * MOMs_1bs * GAM_2 * RLMDg_2
//                                  + MOMs_2bs * GAM * RLMDg );

//             //---< Autoconversion >---            
//             // [Praut] auto-conversion rate from cloud water to rain
//             rhoqc = dens * qc * 1000.0; // [g/m3]
//             Dc    = 0.146 - 5.964E-2 * std::log( Nc[k][ij] / 2000.0 );
//             Praut_berry = Rdens * 1.67E-5 * rhoqc * rhoqc / ( 5.0 + 3.66E-2 * Nc[k][ij] / ( Dc * rhoqc + EPS ) );

//             zerosw   = 0.5 - std::copysign(0.5, qc - 1.0E-12 );
//             Praut_kk = 1350.0                                           
//                        * std::exp( std::log( qc + zerosw ) * 2.47 ) * ( 1.0 - zerosw ) 
//                        * std::exp( std::log( Nc[k][ij] ) * (-1.79) );                     // eq.(29) in KK(2000)

//             Praut_kk = 1350.0 * std::pow(qc, 2.47) * std::pow(Nc[k][ij], -1.79);


//             // switch berry / k-k scheme
//             wk[I_Praut] = ( 1.0 - sw_kk2000 ) * Praut_berry
//                           + sw_kk2000 * Praut_kk;

//             // [Psaut] auto-conversion rate from cloud ice to snow
//             betai = std::min( beta_saut, beta_saut * std::exp( gamma_saut * temc ) );
//             wk[I_Psaut] = std::max( betai * (qi - qicrt_saut), 0.0 );

//             // [Pgaut] auto-conversion rate from snow to graupel
//             betas = std::min( beta_gaut, beta_gaut * std::exp( gamma_gaut * temc ) );
//             wk[I_Pgaut] = std::max( betas * (qs - qscrt_gaut), 0.0 );

//             //---< Evaporation, Sublimation, Melting, and Freezing >---
//             Ka  = ( Ka0 + dKa_dT * temc );
//             Kd  = ( Kd0 + dKd_dT * temc ) * PRE00 / pre[k][ij];
//             Nu  = ( nu0 + dnu_dT * temc ) * Rdens;

//             Glv = 1.0 / ( LHV0/(Ka * temp) * ( LHV0/(Rvap * temp) - 1.0 ) + 1.0 / (Kd * dens * qsatl[k][ij]) );
//             Giv = 1.0 / ( LHS0/(Ka * temp) * ( LHS0/(Rvap * temp) - 1.0 ) + 1.0 / (Kd * dens * qsati[k][ij]) );
//             Gil = 1.0 / ( LHF0/(Ka * temc) );

//             // [Prevp] evaporation rate of rain
//             ventr = f1r * GAM_2 * RLMDr_2 + f2r * std::sqrt( Cr * rho_fact / Nu * RLMDr_5dr ) * GAM_5dr_h;

//             wk[I_Prevp] = 2.0 * PI * Rdens * N0r * ( 1.0 - std::min(Sliq, 1.0) ) * Glv * ventr;

//             // [Pidep,Pisub] deposition/sublimation rate for ice
//             rhoqi = std::max(dens * qi, EPS);
//             XNi   = std::min( std::max( 5.38E+7 * std::exp( std::log(rhoqi) * 0.75 ), 1.0E+3 ), 1.0E+6 );
//             XMi   = rhoqi / XNi;
//             Di    = std::min( Di_a * std::sqrt(XMi), Di_max );

//             tmp = 4.0 * Di * XNi * Rdens * ( Sice - 1.0 ) * Giv;

//             wk[I_Pidep] = (       wk[I_spsati] ) * ( tmp); // Sice > 2
//             wk[I_Pisub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1

//             // [Pihom] homogenious freezing at T < -40C
//             sw = ( 0.5 - std::copysign(0.5, temc + 40.0) ); // if T < -40C, sw=1

//             wk[I_Pihom] = sw * qc / dt;

//             // [Pihtr] heteroginous freezing at -40C < T < 0C
//             sw =   ( 0.5 + std::copysign(0.5, temc + 40.0) )
//                  * ( 0.5 - std::copysign(0.5, temc       ) ); // if -40C < T < 0C, sw=1

//             wk[I_Pihtr] = sw * ( dens / rho_w * (qc*qc) / ( Nc_ihtr * 1.0E+6 ) )
//                              * B_frz * ( std::exp(-A_frz * temc) - 1.0 );

//             // [Pimlt] ice melting at T > 0C
//             sw = ( 0.5 + std::copysign(0.5, temc) ); // if T > 0C, sw=1

//             wk[I_Pimlt] = sw * qi / dt;

//             // [Psdep,Pssub] deposition/sublimation rate for snow
//             vents = f1s * MOMs_1 + f2s * std::sqrt( Cs * rho_fact / Nu ) * MOMs_5ds_h;

//             tmp = 2.0 * PI * Rdens * ( Sice - 1.0 ) * Giv * vents;

//             wk[I_Psdep] = (       wk[I_spsati] ) * ( tmp); // Sice > 1
//             wk[I_Pssub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1            

//             // [Psmlt] melting rate of snow
//             wk[I_Psmlt] = 2.0 * PI * Rdens * Gil * vents
//                           + CL * temc / LHF0 * ( wk[I_Psacw] + wk[I_Psacr] );
//             wk[I_Psmlt] = std::max( wk[I_Psmlt], 0.0 );

//             // [Pgdep/pgsub] deposition/sublimation rate for graupel
//             ventg = f1g * GAM_2 * RLMDg_2 + f2g * std::sqrt( Cg * rho_fact / Nu * RLMDg_5dg ) * GAM_5dg_h;

//             tmp = 2.0 * PI * Rdens * N0g * ( Sice - 1.0 ) * Giv * ventg;

//             wk[I_Pgdep] = (       wk[I_spsati] ) * ( tmp); // Sice > 1
//             wk[I_Pgsub] = ( 1.0 - wk[I_spsati] ) * (-tmp); // Sice < 1

//             // [Pgmlt] melting rate of graupel
//             wk[I_Pgmlt] = 2.0 * PI * Rdens * N0g * Gil * ventg
//                           + CL * temc / LHF0 * ( wk[I_Pgacw] + wk[I_Pgacr] );
//             wk[I_Pgmlt] = std::max( wk[I_Pgmlt], 0.0 );

//             // [Pgfrz] freezing rate of graupel
//             wk[I_Pgfrz] = 2.0 * PI * Rdens * N0r * 60.0 * B_frz * Ar * ( std::exp(-A_frz * temc) - 1.0 ) * RLMDr_7;

//             // [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
//             dt1  = ( std::exp( std::log(mi50) * ma2[k][ij] )
//                    - std::exp( std::log(mi40) * ma2[k][ij] ) ) / ( a1[k][ij] * ma2[k][ij] );
//             Ni50 = qi * dt / ( mi50 * dt1 );

//             wk[I_Psfw] = Ni50 * ( a1[k][ij] * std::pow(mi50, a2[k][ij])
//                                   + PI * Eiw * dens * qc * Ri50 * Ri50 * vti50 );
//             wk[I_Psfi] = qi / dt1;

//             //---< limiter >---
//             wk[I_Pigen] = std::min( wk[I_Pigen], wk[I_dqv_dt] ) * ( wk[I_iceflg] ) * sw_expice;
//             wk[I_Pidep] = std::min( wk[I_Pidep], wk[I_dqv_dt] ) * ( wk[I_iceflg] ) * sw_expice;
//             wk[I_Psdep] = std::min( wk[I_Psdep], wk[I_dqv_dt] ) * ( wk[I_iceflg] )            ;
//             wk[I_Pgdep] = std::min( wk[I_Pgdep], wk[I_dqv_dt] ) * ( wk[I_iceflg] )            ;

//             wk[I_Pracw] = wk[I_Pracw]                           
//                         + wk[I_Psacw] * ( 1.0 - wk[I_iceflg] )   // c->r by s
//                         + wk[I_Pgacw] * ( 1.0 - wk[I_iceflg] );  // c->r by g

//             wk[I_Praut] = std::min( wk[I_Praut], wk[I_dqc_dt] ); 
//             wk[I_Pracw] = std::min( wk[I_Pracw], wk[I_dqc_dt] ); 
//             wk[I_Pihom] = std::min( wk[I_Pihom], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_expice;
//             wk[I_Pihtr] = std::min( wk[I_Pihtr], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_expice;
//             wk[I_Psacw] = std::min( wk[I_Psacw], wk[I_dqc_dt] ) * (        wk[I_iceflg] )            ;
//             wk[I_Psfw ] = std::min( wk[I_Psfw ], wk[I_dqc_dt] ) * (        wk[I_iceflg] ) * sw_bergeron;
//             wk[I_Pgacw] = std::min( wk[I_Pgacw], wk[I_dqc_dt] ) * (        wk[I_iceflg] )            ;

//             wk[I_Prevp] = std::min( wk[I_Prevp], wk[I_dqr_dt] );
//             wk[I_Piacr] = std::min( wk[I_Piacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
//             wk[I_Psacr] = std::min( wk[I_Psacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
//             wk[I_Pgacr] = std::min( wk[I_Pgacr], wk[I_dqr_dt] ) * (        wk[I_iceflg] );
//             wk[I_Pgfrz] = std::min( wk[I_Pgfrz], wk[I_dqr_dt] ) * (        wk[I_iceflg] );

//             wk[I_Pisub] = std::min( wk[I_Pisub], wk[I_dqi_dt] ) * (       wk[I_iceflg] ) * sw_expice;
//             wk[I_Pimlt] = std::min( wk[I_Pimlt], wk[I_dqi_dt] ) * ( 1.0 - wk[I_iceflg] ) * sw_expice;
//             wk[I_Psaut] = std::min( wk[I_Psaut], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
//             wk[I_Praci] = std::min( wk[I_Praci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
//             wk[I_Psaci] = std::min( wk[I_Psaci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;
//             wk[I_Psfi ] = std::min( wk[I_Psfi ], wk[I_dqi_dt] ) * (       wk[I_iceflg] ) * sw_bergeron;
//             wk[I_Pgaci] = std::min( wk[I_Pgaci], wk[I_dqi_dt] ) * (       wk[I_iceflg] )            ;

//             wk[I_Pssub] = std::min( wk[I_Pssub], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
//             wk[I_Psmlt] = std::min( wk[I_Psmlt], wk[I_dqs_dt] ) * ( 1.0 - wk[I_iceflg] );
//             wk[I_Pgaut] = std::min( wk[I_Pgaut], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
//             wk[I_Pracs] = std::min( wk[I_Pracs], wk[I_dqs_dt] ) * (       wk[I_iceflg] );
//             wk[I_Pgacs] = std::min( wk[I_Pgacs], wk[I_dqs_dt] );

//             wk[I_Pgsub] = std::min( wk[I_Pgsub], wk[I_dqg_dt] ) * (       wk[I_iceflg] );
//             wk[I_Pgmlt] = std::min( wk[I_Pgmlt], wk[I_dqg_dt] ) * ( 1.0 - wk[I_iceflg] );

//             wk[I_Piacr_s] = ( 1.0 - wk[I_delta1] ) * wk[I_Piacr];
//             wk[I_Piacr_g] = (       wk[I_delta1] ) * wk[I_Piacr];
//             wk[I_Praci_s] = ( 1.0 - wk[I_delta1] ) * wk[I_Praci];
//             wk[I_Praci_g] = (       wk[I_delta1] ) * wk[I_Praci];
//             wk[I_Psacr_s] = (       wk[I_delta2] ) * wk[I_Psacr];
//             wk[I_Psacr_g] = ( 1.0 - wk[I_delta2] ) * wk[I_Psacr];
//             wk[I_Pracs  ] = ( 1.0 - wk[I_delta2] ) * wk[I_Pracs];

//             // [QC]
//             net = 
//                + wk[I_Pimlt]  // [prod] i->c
//                - wk[I_Praut]  // [loss] c->r
//                - wk[I_Pracw]  // [loss] c->r
//                - wk[I_Pihom]  // [loss] c->i
//                - wk[I_Pihtr]  // [loss] c->i
//                - wk[I_Psacw]  // [loss] c->s
//                - wk[I_Psfw ]  // [loss] c->s
//                - wk[I_Pgacw]; // [loss] c->g

//             fac_sw = 0.5 + std::copysign( 0.5, net + EPS ); // if production > loss , fac_sw=1
//             fac    = fac_sw + ( 1.0 - fac_sw ) 
//                      * std::min( -wk[I_dqc_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Pimlt] = wk[I_Pimlt] * fac;
//             wk[I_Praut] = wk[I_Praut] * fac;
//             wk[I_Pracw] = wk[I_Pracw] * fac;
//             wk[I_Pihom] = wk[I_Pihom] * fac;
//             wk[I_Pihtr] = wk[I_Pihtr] * fac;
//             wk[I_Psacw] = wk[I_Psacw] * fac;
//             wk[I_Psfw ] = wk[I_Psfw ] * fac;
//             wk[I_Pgacw] = wk[I_Pgacw] * fac;

//             // [QI]
//             net = 
//                + wk[I_Pigen  ]  // [prod] v->i
//                + wk[I_Pidep  ]  // [prod] v->i
//                + wk[I_Pihom  ]  // [prod] c->i
//                + wk[I_Pihtr  ]  // [prod] c->i
//                - wk[I_Pisub  ]  // [loss] i->v
//                - wk[I_Pimlt  ]  // [loss] i->c
//                - wk[I_Psaut  ]  // [loss] i->s
//                - wk[I_Praci_s]  // [loss] i->s
//                - wk[I_Psaci  ]  // [loss] i->s
//                - wk[I_Psfi   ]  // [loss] i->s
//                - wk[I_Praci_g]  // [loss] i->g
//                - wk[I_Pgaci  ]; // [loss] i->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqi_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Pigen  ] = wk[I_Pigen  ] * fac;
//             wk[I_Pidep  ] = wk[I_Pidep  ] * fac;
//             wk[I_Pihom  ] = wk[I_Pihom  ] * fac;
//             wk[I_Pihtr  ] = wk[I_Pihtr  ] * fac;
//             wk[I_Pisub  ] = wk[I_Pisub  ] * fac;
//             wk[I_Pimlt  ] = wk[I_Pimlt  ] * fac;
//             wk[I_Psaut  ] = wk[I_Psaut  ] * fac;
//             wk[I_Praci_s] = wk[I_Praci_s] * fac;
//             wk[I_Psaci  ] = wk[I_Psaci  ] * fac;
//             wk[I_Psfi   ] = wk[I_Psfi   ] * fac;
//             wk[I_Praci_g] = wk[I_Praci_g] * fac;
//             wk[I_Pgaci  ] = wk[I_Pgaci  ] * fac;


//             // [QR]
//             net = 
//                + wk[I_Praut  ]  // [prod] c->r
//                + wk[I_Pracw  ]  // [prod] c->r
//                + wk[I_Psmlt  ]  // [prod] s->r
//                + wk[I_Pgmlt  ]  // [prod] g->r
//                - wk[I_Prevp  ]  // [loss] r->v
//                - wk[I_Piacr_s]  // [loss] r->s
//                - wk[I_Psacr_s]  // [loss] r->s
//                - wk[I_Piacr_g]  // [loss] r->g
//                - wk[I_Psacr_g]  // [loss] r->g
//                - wk[I_Pgacr  ]  // [loss] r->g
//                - wk[I_Pgfrz  ]; // [loss] r->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac    = fac_sw + ( 1.0 - fac_sw ) 
//                      * std::min( -wk[I_dqr_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Praut  ] = wk[I_Praut  ] * fac;
//             wk[I_Pracw  ] = wk[I_Pracw  ] * fac;
//             wk[I_Psmlt  ] = wk[I_Psmlt  ] * fac;
//             wk[I_Pgmlt  ] = wk[I_Pgmlt  ] * fac;
//             wk[I_Prevp  ] = wk[I_Prevp  ] * fac;
//             wk[I_Piacr_s] = wk[I_Piacr_s] * fac;
//             wk[I_Psacr_s] = wk[I_Psacr_s] * fac;
//             wk[I_Piacr_g] = wk[I_Piacr_g] * fac;
//             wk[I_Psacr_g] = wk[I_Psacr_g] * fac;
//             wk[I_Pgacr  ] = wk[I_Pgacr  ] * fac;
//             wk[I_Pgfrz  ] = wk[I_Pgfrz  ] * fac;

//             // [QV]
//             net =
//                + wk[I_Prevp]  // [prod] r->v
//                + wk[I_Pisub]  // [prod] i->v
//                + wk[I_Pssub]  // [prod] s->v
//                + wk[I_Pgsub]  // [prod] g->v
//                - wk[I_Pigen]  // [loss] v->i
//                - wk[I_Pidep]  // [loss] v->i
//                - wk[I_Psdep]  // [loss] v->s
//                - wk[I_Pgdep]; // [loss] v->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqv_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Prevp] = wk[I_Prevp] * fac;
//             wk[I_Pisub] = wk[I_Pisub] * fac;
//             wk[I_Pssub] = wk[I_Pssub] * fac;
//             wk[I_Pgsub] = wk[I_Pgsub] * fac;
//             wk[I_Pigen] = wk[I_Pigen] * fac;
//             wk[I_Pidep] = wk[I_Pidep] * fac;
//             wk[I_Psdep] = wk[I_Psdep] * fac;
//             wk[I_Pgdep] = wk[I_Pgdep] * fac;

//             // [QS]
//             net =
//                + wk[I_Psdep  ]  // [prod] v->s
//                + wk[I_Psacw  ]  // [prod] c->s
//                + wk[I_Psfw   ]  // [prod] c->s
//                + wk[I_Piacr_s]  // [prod] r->s
//                + wk[I_Psacr_s]  // [prod] r->s
//                + wk[I_Psaut  ]  // [prod] i->s
//                + wk[I_Praci_s]  // [prod] i->s
//                + wk[I_Psaci  ]  // [prod] i->s
//                + wk[I_Psfi   ]  // [prod] i->s
//                - wk[I_Pssub  ]  // [loss] s->v
//                - wk[I_Psmlt  ]  // [loss] s->r
//                - wk[I_Pgaut  ]  // [loss] s->g
//                - wk[I_Pracs  ]  // [loss] s->g
//                - wk[I_Pgacs  ]; // [loss] s->g

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqs_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Psdep  ] = wk[I_Psdep  ] * fac;
//             wk[I_Psacw  ] = wk[I_Psacw  ] * fac;
//             wk[I_Psfw   ] = wk[I_Psfw   ] * fac;
//             wk[I_Piacr_s] = wk[I_Piacr_s] * fac;
//             wk[I_Psacr_s] = wk[I_Psacr_s] * fac;
//             wk[I_Psaut  ] = wk[I_Psaut  ] * fac;
//             wk[I_Praci_s] = wk[I_Praci_s] * fac;
//             wk[I_Psaci  ] = wk[I_Psaci  ] * fac;
//             wk[I_Psfi   ] = wk[I_Psfi   ] * fac;
//             wk[I_Pssub  ] = wk[I_Pssub  ] * fac;
//             wk[I_Psmlt  ] = wk[I_Psmlt  ] * fac;
//             wk[I_Pgaut  ] = wk[I_Pgaut  ] * fac;
//             wk[I_Pracs  ] = wk[I_Pracs  ] * fac;
//             wk[I_Pgacs  ] = wk[I_Pgacs  ] * fac;

//             // [QG]
//             net =
//                + wk[I_Pgdep  ]  // [prod] v->g
//                + wk[I_Pgacw  ]  // [prod] c->g
//                + wk[I_Piacr_g]  // [prod] r->g
//                + wk[I_Psacr_g]  // [prod] r->g
//                + wk[I_Pgacr  ]  // [prod] r->g
//                + wk[I_Pgfrz  ]  // [prod] r->g
//                + wk[I_Praci_g]  // [prod] i->g
//                + wk[I_Pgaci  ]  // [prod] i->g
//                + wk[I_Pgaut  ]  // [prod] s->g
//                + wk[I_Pracs  ]  // [prod] s->g
//                + wk[I_Pgacs  ]  // [prod] s->g
//                - wk[I_Pgsub  ]  // [loss] g->v
//                - wk[I_Pgmlt  ]; // [loss] g->r

//             fac_sw = 0.5 + std::copysign( 0.5, net+EPS ); // if production > loss , fac_sw=1
//             fac = fac_sw + ( 1.0 - fac_sw ) 
//                   * std::min( -wk[I_dqg_dt]/(net - fac_sw), 1.0 ); // loss limiter

//             wk[I_Pgdep  ] = wk[I_Pgdep  ] * fac;
//             wk[I_Pgacw  ] = wk[I_Pgacw  ] * fac;
//             wk[I_Piacr_g] = wk[I_Piacr_g] * fac;
//             wk[I_Psacr_g] = wk[I_Psacr_g] * fac;
//             wk[I_Pgacr  ] = wk[I_Pgacr  ] * fac;
//             wk[I_Pgfrz  ] = wk[I_Pgfrz  ] * fac;
//             wk[I_Praci_g] = wk[I_Praci_g] * fac;
//             wk[I_Pgaci  ] = wk[I_Pgaci  ] * fac;
//             wk[I_Pgaut  ] = wk[I_Pgaut  ] * fac;
//             wk[I_Pracs  ] = wk[I_Pracs  ] * fac;
//             wk[I_Pgacs  ] = wk[I_Pgacs  ] * fac;
//             wk[I_Pgsub  ] = wk[I_Pgsub  ] * fac;
//             wk[I_Pgmlt  ] = wk[I_Pgmlt  ] * fac;

//             qc_t = + wk[I_Pimlt]  // [prod] i->c
//                    - wk[I_Praut]  // [loss] c->r
//                    - wk[I_Pracw]  // [loss] c->r
//                    - wk[I_Pihom]  // [loss] c->i
//                    - wk[I_Pihtr]  // [loss] c->i
//                    - wk[I_Psacw]  // [loss] c->s
//                    - wk[I_Psfw ]  // [loss] c->s
//                    - wk[I_Pgacw]; // [loss] c->g
 
//             qr_t = + wk[I_Praut  ]  // [prod] c->r
//                    + wk[I_Pracw  ]  // [prod] c->r
//                    + wk[I_Psmlt  ]  // [prod] s->r
//                    + wk[I_Pgmlt  ]  // [prod] g->r
//                    - wk[I_Prevp  ]  // [loss] r->v
//                    - wk[I_Piacr_s]  // [loss] r->s
//                    - wk[I_Psacr_s]  // [loss] r->s
//                    - wk[I_Piacr_g]  // [loss] r->g
//                    - wk[I_Psacr_g]  // [loss] r->g
//                    - wk[I_Pgacr  ]  // [loss] r->g
//                    - wk[I_Pgfrz  ]; // [loss] r->g

//             qi_t = + wk[I_Pigen  ]  // [prod] v->i
//                    + wk[I_Pidep  ]  // [prod] v->i
//                    + wk[I_Pihom  ]  // [prod] c->i
//                    + wk[I_Pihtr  ]  // [prod] c->i
//                    - wk[I_Pisub  ]  // [loss] i->v
//                    - wk[I_Pimlt  ]  // [loss] i->c
//                    - wk[I_Psaut  ]  // [loss] i->s
//                    - wk[I_Praci_s]  // [loss] i->s
//                    - wk[I_Psaci  ]  // [loss] i->s
//                    - wk[I_Psfi   ]  // [loss] i->s
//                    - wk[I_Praci_g]  // [loss] i->g
//                    - wk[I_Pgaci  ]; // [loss] i->g

//             qs_t = + wk[I_Psdep  ]  // [prod] v->s
//                    + wk[I_Psacw  ]  // [prod] c->s
//                    + wk[I_Psfw   ]  // [prod] c->s
//                    + wk[I_Piacr_s]  // [prod] r->s
//                    + wk[I_Psacr_s]  // [prod] r->s
//                    + wk[I_Psaut  ]  // [prod] i->s
//                    + wk[I_Praci_s]  // [prod] i->s
//                    + wk[I_Psaci  ]  // [prod] i->s
//                    + wk[I_Psfi   ]  // [prod] i->s
//                    - wk[I_Pssub  ]  // [loss] s->v
//                    - wk[I_Psmlt  ]  // [loss] s->r
//                    - wk[I_Pgaut  ]  // [loss] s->g
//                    - wk[I_Pracs  ]  // [loss] s->g
//                    - wk[I_Pgacs  ]; // [loss] s->g

//             qg_t = + wk[I_Pgdep  ]  // [prod] v->g
//                    + wk[I_Pgacw  ]  // [prod] c->g
//                    + wk[I_Piacr_g]  // [prod] r->g
//                    + wk[I_Psacr_g]  // [prod] r->g
//                    + wk[I_Pgacr  ]  // [prod] r->g
//                    + wk[I_Pgfrz  ]  // [prod] r->g
//                    + wk[I_Praci_g]  // [prod] i->g
//                    + wk[I_Pgaci  ]  // [prod] i->g
//                    + wk[I_Pgaut  ]  // [prod] s->g
//                    + wk[I_Pracs  ]  // [prod] s->g
//                    + wk[I_Pgacs  ]  // [prod] s->g
//                    - wk[I_Pgsub  ]  // [loss] g->v
//                    - wk[I_Pgmlt  ]; // [loss] g->r

//             qc_t = std::max( qc_t, -wk[I_dqc_dt] );
//             qr_t = std::max( qr_t, -wk[I_dqr_dt] );
//             qi_t = std::max( qi_t, -wk[I_dqi_dt] );
//             qs_t = std::max( qs_t, -wk[I_dqs_dt] );
//             qg_t = std::max( qg_t, -wk[I_dqg_dt] );

//             qv_t = - ( qc_t 
//                      + qr_t 
//                      + qi_t 
//                      + qs_t 
//                      + qg_t );

//             drhogqv[k][ij] = rhog[k][ij] * qv_t * dt;
//             drhogqc[k][ij] = rhog[k][ij] * qc_t * dt;
//             drhogqr[k][ij] = rhog[k][ij] * qr_t * dt;
//             drhogqi[k][ij] = rhog[k][ij] * qi_t * dt;
//             drhogqs[k][ij] = rhog[k][ij] * qs_t * dt;
//             drhogqg[k][ij] = rhog[k][ij] * qg_t * dt;

//             Vt[I_QR][k][ij] = Vtr;
//             Vt[I_QI][k][ij] = Vti;
//             Vt[I_QS][k][ij] = Vts;
//             Vt[I_QG][k][ij] = Vtg;

//             //--- Effective Radius of Liquid Water
//             xf_qc = rhoqc / Nc[k][ij] * 1.0E-9;                            // mean mass   of qc [kg]
//             rf_qc = std::pow(( coef_xf * xf_qc + EPS ), 0.33333333);       // mean radius of qc [m]
//             r2_qc = coef_dgam * (rf_qc * rf_qc) * ( Nc[k][ij] * 1.0E+6 );  // r^2 moment  of qc
//             r2_qr = 0.25 * N0r * GAM_3 * RLMDr_3;                          // r^2 moment  of qr
//             r3_qc = coef_xf * dens * qc;                                   // r^3 moment  of qc
//             r3_qr = coef_xf * dens * qr;                                   // r^3 moment  of qr

//             zerosw = 0.5 - std::copysign(0.5, (r2_qc + r2_qr) - r2_min );
//             rceff[k][ij] = ( r3_qc + r3_qr ) / ( r2_qc + r2_qr + zerosw ) * ( 1.0 - zerosw );

//             sw = ( 0.5 + std::copysign(0.5, (qc + qr) - q_min ) ) * zerosw; // if qc+qr > 0.01[g/kg], sw=1

//             rctop[0][ij] = (       sw ) * rceff[k][ij]
//                          + ( 1.0 - sw ) * UNDEF;
//             tctop[0][ij] = (       sw ) * temp
//                          + ( 1.0 - sw ) * UNDEF;

//             zerosw = 0.5 - std::copysign( 0.5, r2_qc - r2_min );
//             rceff_cld[k][ij] = r3_qc / ( r2_qc + zerosw ) * ( 1.0 - zerosw );

//             sw = ( 0.5 + std::copysign( 0.5, qc - q_min ) ) * zerosw; // if qc > 0.01[g/kg], sw=1

//             rctop_cld[0][ij] = (       sw ) * rceff_cld[k][ij]
//                              + ( 1.0 - sw ) * UNDEF;
//             tctop_cld[0][ij] = (       sw ) * temp
//                              + ( 1.0 - sw ) * UNDEF;
//         }
//     }


//     // update rhogq
//     for(int k = kmin; k <= kmax; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             rhogq[I_QV][k][ij] += drhogqv[k][ij];
//             rhogq[I_QC][k][ij] += drhogqc[k][ij];
//             rhogq[I_QR][k][ij] += drhogqr[k][ij];
//             rhogq[I_QI][k][ij] += drhogqi[k][ij];
//             rhogq[I_QS][k][ij] += drhogqs[k][ij];
//             rhogq[I_QG][k][ij] += drhogqg[k][ij];
//         }
//     }

//     // update rhoge
//     for(int k = kmin; k <= kmax; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             rhoge[k][ij] = rhoge[k][ij] - LHV * drhogqv[k][ij]
//                                         + LHF * drhogqi[k][ij]
//                                         + LHF * drhogqs[k][ij]
//                                         + LHF * drhogqg[k][ij];
//         }
//     }
//     for(int ij = 0; ij < ijdim; ij++)
//     {
//          rceff    [kmin-1][ij] = 0.0;
//          rceff_cld[kmin-1][ij] = 0.0;
//          rceff    [kmax+1][ij] = 0.0;
//          rceff_cld[kmax+1][ij] = 0.0;

//          sw = ( 0.5 + std::copysign(0.5, tctop[0][ij] - TEM00 ) ); // if T > 0C, sw=1

//          rwtop[0][ij] =  (       sw ) * rctop[0][ij]
//                        + ( 1.0 - sw ) * UNDEF;

//          sw = ( 0.5 + std::copysign(0.5, rctop_cld[0][ij] - TEM00 ) ); // if T > 0C, sw=1

//          rwtop_cld[0][ij] =  (       sw ) * rctop_cld[0][ij]
//                            + ( 1.0 - sw ) * UNDEF;
//     }
    

//     //---< preciptation（降水量） >---
//     for(int nq = 0; nq < nqmax; nq++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             Vt[nq][kmin-1][ij] = 0.0;
//             Vt[nq][kmax+1][ij] = 0.0;
//         }
//     }
//     //--- update mass concentration（質量濃度）
//     for(int nq = NQW_STR; nq <= NQW_END; nq++)
//     {
//         for(int k = kmin; k <= kmax; k++)
//         {
//             for(int ij = 0; ij < ijdim; ij++)
//             {
//                 q[nq][k][ij] = rhogq[nq][k][ij] / rhog[k][ij];
//             }
//         }
//     }

//     THRMDYN_qd(q, qd);
//     THRMDYN_cv(qd, q, cva);

//     for(int k = kmin; k <= kmax; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             tem [k][ij] = rhoge[k][ij] / ( rhog[k][ij] * cva[k][ij] );
//             rgs [k][ij] = gam2 [k][ij] / gsgam2 [k][ij];
//             rgsh[k][ij] = gam2h[k][ij] / gsgam2h[k][ij];
//         }
//     }

//     // for(int nq = 0; nq < nqmax; nq++)
//     //     preciptation_flag[nq] = false;
//     preciptation_flag[I_QR] = true;
//     preciptation_flag[I_QI] = true;
//     preciptation_flag[I_QS] = true;
//     preciptation_flag[I_QG] = true;
//     if(precip_transport_type == "3WATER")
//     {
//         preciptation_flag[I_QI] = false;
//     }

//     precip_transport_new( rhog,
//                           rhogvx,
//                           rhogvy,
//                           rhogvz,
//                           rhogw,
//                           rhoge,
//                           rhogq,
//                           rho,
//                           tem,
//                           pre,
//                           vx,
//                           vy,
//                           vz,
//                           w,
//                           q,
//                           qd,
//                           z,
//                           Vt,
//                           preciptation_flag,
//                           precip,precip_rhoe,
//                           precip_lh_heat,
//                           precip_rhophi,
//                           precip_rhokin,
//                           gprec,
//                           gsgam2,
//                           gsgam2h,
//                           rgs,
//                           rgsh,
//                           ix,
//                           iy,
//                           iz,
//                           jx,
//                           jy,
//                           jz,
//                           dt,
//                           nullptr); // precip_trc**
    
//     for(int k = 0; k < kdim; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             ml_Pconv[k][ij] = q[I_QV][k][ij];
//             ml_Pconw[k][ij] = q[I_QC][k][ij];
//             ml_Pconi[k][ij] = q[I_QI][k][ij];
//         }
//     }

//     bool ice_adjust;
//     if(OPT_EXPLICIT_ICEGEN)
//     {
//         ice_adjust = false;
//         SATURATION_Setrange(100.0, 90.0);

//         SATURATION_adjustment(rhog, rhoge, rhogq, tem, q, qd, gsgam2, ice_adjust);
//     }
//     else
//     {
//         ice_adjust = true;
//         SATURATION_Setrange(273.16, 233.16);

//         SATURATION_adjustment(rhog, rhoge, rhogq, tem , q, qd, gsgam2, ice_adjust);
//     }

//     for(int k = 0; k < kdim; k++)
//     {
//         for(int ij = 0; ij < ijdim; ij++)
//         {
//             ml_Pconv[k][ij] = q[I_QV][k][ij] - ml_Pconv[k][ij];
//             ml_Pconw[k][ij] = q[I_QC][k][ij] - ml_Pconw[k][ij];
//             ml_Pconi[k][ij] = q[I_QI][k][ij] - ml_Pconi[k][ij];
//         }
//     }

//     negative_filter(rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2);

//     // End of nsw6
// }