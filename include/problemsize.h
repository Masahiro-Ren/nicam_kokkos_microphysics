#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

template <typename T>
using Vec1d = std::vector<T>;

template <typename T>
using Vec2d = std::vector<std::vector<T>>;

namespace PROBLEM_SIZE 
{

constexpr int H_SHORT = 32;
constexpr int H_LONG = 1024;
constexpr int IO_FID_LOG = 6;
constexpr int ADM_prc_all = 1;

constexpr int IO_FREAD = 0;
constexpr int ADM_KNONE = 1;

constexpr int ADM_gall_in_orig = 16641; // horizontal grid number (for physics)
constexpr int ADM_gall_in      = 16641; // horizontal grid number (for physics)
constexpr int ADM_kall     = 96;   // vertical   grid number
constexpr int ADM_vlayer   = 94;
constexpr int ADM_kmin     = 1 ; // 2 in fortran
constexpr int ADM_kmax     = 94; // 95 in fortran
constexpr int ADM_lall     = 1 ;
constexpr int TRC_VMAX     = 6 ;

constexpr int ijdim        = ADM_gall_in;
constexpr int kdim         = ADM_kall;

constexpr int kmin         = 1 ; // 2 in fortran
constexpr int kmax         = 94; // 95 in fortran
constexpr int knone        = 1 ;
constexpr int nqmax        = 6 ;

constexpr double TIME_DTL     = 6.0 ;// [sec]
constexpr double TIME_CTIME   = 0.0 ;
constexpr int TIME_CSTEP   = 49     ;
constexpr int SET_l        = 1      ;
constexpr int SET_rgnid    = 4      ;

constexpr double CONST_PI     = 3.141592653589793; //< pi
constexpr double CONST_EPS    = 2.220446E-16; //< small number

constexpr double CONST_UNDEF  = -9.9999E+30; 
constexpr double CONST_GRAV   =   9.79764;   
constexpr double CONST_STB    =   5.670373E-8;

constexpr double CONST_Rdry   =  287.04;
constexpr double CONST_CPdry  = 1004.64;
constexpr double CONST_CVdry  =  717.56;
constexpr double CONST_Rvap   =  461.46;
constexpr double CONST_CPvap  = 1846.0; 
constexpr double CONST_CVvap  = 1384.54;
constexpr double CONST_CL     = 4218.0; 
constexpr double CONST_CI     = 2006.0; 

constexpr double CONST_KAPPA  = 0.2857256619550070;
constexpr double CONST_EPSvap = 0.6219718309859156;
constexpr double CONST_PSAT0  =     610.7;
constexpr double CONST_LHV0   = 2501000.0;
constexpr double CONST_LHS0   = 2834000.0;
constexpr double CONST_LHF0   =  333000.0;
constexpr double CONST_LHV00  = 3148911.8;
constexpr double CONST_LHS00  = 2877704.0;
constexpr double CONST_LHF00  = -271207.8;
constexpr double CONST_Pstd   =  101325.0;
constexpr double CONST_PRE00  =  100000.0;
constexpr double CONST_Tstd   =    288.15;
constexpr double CONST_TEM00  =    273.15;
constexpr double CONST_TMELT  =    273.15;
constexpr double CONST_DWATR  =    1000.0;
constexpr double CONST_DICE   =     916.8;
constexpr double CONST_LHV    = CONST_LHV0;
constexpr double CONST_LHS    = CONST_LHS0;
constexpr double CONST_LHF    = CONST_LHF0;

constexpr double PI    = 3.141592653589793; //< pi
constexpr double EPS   = 2.220446E-16;      //< small number
constexpr double UNDEF = CONST_UNDEF ;
constexpr double GRAV  = CONST_GRAV  ;
constexpr double Rdry  = CONST_Rdry  ;
constexpr double Rvap  = CONST_Rvap  ;
constexpr double CL    = CONST_CL    ;
constexpr double LHV   = CONST_LHV0  ;
constexpr double LHS   = CONST_LHS0  ;
constexpr double LHF   = CONST_LHF0  ;
constexpr double rho_w = CONST_DWATR ;
constexpr double rhow  = CONST_DWATR ;
constexpr double DWATR = CONST_DWATR ;
constexpr double DICE  = CONST_DICE  ;
constexpr double Pstd  = CONST_Pstd  ;
constexpr double Tstd  = CONST_Tstd  ;
constexpr double TEM00 = CONST_TEM00 ;

// RP_PREC = precision(0.D0) in fortran
constexpr int    RP_PREC = 15;

// grd, gmtr, vmtr
extern double GRD_gz   [ADM_kall];
extern double GRD_gzh  [ADM_kall];
extern double GRD_dgz  [ADM_kall];
extern double GRD_dgzh [ADM_kall];
extern double GRD_rdgz [ADM_kall];
extern double GRD_rdgzh[ADM_kall];
extern double GRD_afact[ADM_kall];
extern double GRD_bfact[ADM_kall];
extern double GRD_cfact[ADM_kall];
extern double GRD_dfact[ADM_kall];

// run conf
extern std::string EIN_TYPE  ;
extern std::string RAIN_TYPE ;
extern std::string MP_TYPE   ;
extern std::string vgrid_fname;
extern int MP_DIV_NUM;
extern bool opt_2moment_water   ;
extern bool ISOTOPE             ;
extern bool opt_offline_aerosol ;
extern bool opt_aerosol_forcing ;
extern std::string PRCIP_TRN_ECORRECT;

extern std::string RD_TYPE;
extern std::string AE_TYPE;
extern int KAPCL          ;
extern int NCRF           ;
extern int NRBND          ;
extern int NRBND_VIS      ;
extern int NRBND_NIR      ;
extern int NRBND_IR       ;
extern int NRDIR          ;
extern int NRDIR_DIRECT   ;
extern int NRDIR_DIFFUSE  ;
extern int NTAU           ;
extern int NPRES          ;
extern int HYDRO_MAX      ;

extern int NQW_STR; 
extern int NQW_END; 
extern int I_QV ;
extern int I_QC ;
extern int I_QR ;
extern int I_QI ;
extern int I_QS ;
extern int I_QG ;
extern int I_QH ;
extern int I_NC ;
extern int I_NR ;
extern int I_NI ;
extern int I_NS ;
extern int I_NG ;

extern std::vector<std::string> TRC_name;

extern std::vector<double> CVW;
extern std::vector<double> CPW;

extern int SET_iteration;
extern bool SET_check;

}
