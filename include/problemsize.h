#pragma once

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <Kokkos_Core.hpp>
#include <Kokkos_Abort.hpp>
#include <Kokkos_Timer.hpp>
// #include <Kokkos_StdAlgorithms.hpp>

using Kokkos::View;
using Kokkos::subview;
using Kokkos::RangePolicy;
using Kokkos::MDRangePolicy;
using Kokkos::Schedule;

// Rename execution space
using HOST_SPACE = Kokkos::OpenMP;
// Rename memory space
using HOST_MEM = Kokkos::HostSpace;

#if defined(USE_OPENMP)
using DEVICE_SPACE = Kokkos::OpenMP;
using DEVICE_MEM = Kokkos::HostSpace;
using DEFAULT_MEM = HOST_MEM;
#elif defined(USE_CUDA)
using DEVICE_SPACE = Kokkos::Cuda;
using DEVICE_MEM = Kokkos::CudaSpace;
using DEFAULT_MEM = DEVICE_MEM;
#else
using DEVICE_SPACE = Kokkos::OpenMP;
using DEVICE_MEM = Kokkos::HostSpace;
using DEFAULT_MEM = HOST_MEM;
#endif

template<typename T, typename S>
using View1D = Kokkos::View<T*, Kokkos::LayoutRight, S>;
template<typename T, typename S>
using View2D = Kokkos::View<T**, Kokkos::LayoutRight, S>;
template<typename T, typename S>
using View3D = Kokkos::View<T***, Kokkos::LayoutRight, S>;
template<typename T, typename S>
using View4D = Kokkos::View<T****, Kokkos::LayoutRight, S>;

namespace PROBLEM_SIZE 
{

constexpr size_t IDX_ZERO = 0;
constexpr size_t H_SHORT = 32;
constexpr size_t H_LONG = 1024;
constexpr size_t IO_FID_LOG = 6;
constexpr size_t ADM_prc_all = 1;

constexpr size_t IO_FREAD = 0;
constexpr size_t ADM_KNONE = 1;

constexpr size_t ADM_gall_in_orig = 16641; // horizontal grid number (for physics)
constexpr size_t ADM_gall_in      = 16641; // horizontal grid number (for physics)
constexpr size_t ADM_kall     = 96;   // vertical   grid number
constexpr size_t ADM_vlayer   = 94;
constexpr size_t ADM_kmin     = 1 ; // 2 in fortran
constexpr size_t ADM_kmax     = 94; // 95 in fortran
constexpr size_t ADM_lall     = 1 ;
constexpr size_t TRC_VMAX     = 6 ;

constexpr size_t ijdim        = ADM_gall_in;
constexpr size_t kdim         = ADM_kall;

constexpr size_t kmin         = 1 ; // 2 in fortran
constexpr size_t kmax         = 94; // 95 in fortran
constexpr size_t knone        = 1 ;
constexpr size_t nqmax        = 6 ;

constexpr double TIME_DTL     = 6.0 ;// [sec]
constexpr double TIME_CTIME   = 0.0 ;
constexpr int TIME_CSTEP   = 49     ;
constexpr int SET_l        = 1      ;
constexpr int SET_rgnid    = 4      ;

constexpr double CONST_PI     = 3.141592653589793; //< pi
constexpr double CONST_EPS    = 2.220446E-16; //< small number

constexpr double CONST_UNDEF  = -9.9999E+30; 
constexpr double CONST_GRAV   =  9.79764;   
constexpr double CONST_STB    =  5.670373E-8;

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
extern double arrGRD_gz   [ADM_kall];
extern double arrGRD_gzh  [ADM_kall];
extern double arrGRD_dgz  [ADM_kall];
extern double arrGRD_dgzh [ADM_kall];
extern double arrGRD_rdgz [ADM_kall];
extern double arrGRD_rdgzh[ADM_kall];
extern double arrGRD_afact[ADM_kall];
extern double arrGRD_bfact[ADM_kall];
extern double arrGRD_cfact[ADM_kall];
extern double arrGRD_dfact[ADM_kall];
// view version of grd, gmtr, vmtr
extern View1D<double, DEFAULT_MEM> GRD_gz   ;
extern View1D<double, DEFAULT_MEM> GRD_gzh  ;
extern View1D<double, DEFAULT_MEM> GRD_dgz  ;
extern View1D<double, DEFAULT_MEM> GRD_dgzh ;
extern View1D<double, DEFAULT_MEM> GRD_rdgz ;
extern View1D<double, DEFAULT_MEM> GRD_rdgzh;
extern View1D<double, DEFAULT_MEM> GRD_afact;
extern View1D<double, DEFAULT_MEM> GRD_bfact;
extern View1D<double, DEFAULT_MEM> GRD_cfact;
extern View1D<double, DEFAULT_MEM> GRD_dfact;

// run conf
extern std::string EIN_TYPE  ;
extern std::string RAIN_TYPE ;
extern std::string MP_TYPE   ;
extern std::string vgrid_fname;
// extern int MP_DIV_NUM;
// extern bool opt_2moment_water   ;
// extern bool ISOTOPE             ;
// extern bool opt_offline_aerosol ;
// extern bool opt_aerosol_forcing ;
constexpr int MP_DIV_NUM           = 1;
constexpr bool opt_2moment_water   = false;
constexpr bool ISOTOPE             = false;
constexpr bool opt_offline_aerosol = false;
constexpr bool opt_aerosol_forcing = false;
extern std::string PRCIP_TRN_ECORRECT;

extern std::string RD_TYPE;
extern std::string AE_TYPE;
// extern int KAPCL          ;
// extern int NCRF           ;
// extern int NRBND          ;
// extern int NRBND_VIS      ;
// extern int NRBND_NIR      ;
// extern int NRBND_IR       ;
// extern int NRDIR          ;
// extern int NRDIR_DIRECT   ;
// extern int NRDIR_DIFFUSE  ;
// extern int NTAU           ;
// extern int NPRES          ;
// extern int HYDRO_MAX      ;
constexpr int KAPCL               = 7;
constexpr int NCRF                = 2;
constexpr int NRBND               = 3;
constexpr int NRBND_VIS           = 1;
constexpr int NRBND_NIR           = 2;
constexpr int NRBND_IR            = 3;
constexpr int NRDIR               = 2;
constexpr int NRDIR_DIRECT        = 1;
constexpr int NRDIR_DIFFUSE       = 2;
constexpr int NTAU                = 7;
constexpr int NPRES               = 7;
constexpr int HYDRO_MAX           = 7;

// extern size_t NQW_STR; 
// extern size_t NQW_END; 
// extern int I_QV ;
// extern int I_QC ;
// extern int I_QR ;
// extern int I_QI ;
// extern int I_QS ;
// extern int I_QG ;
// extern int I_QH ;
// extern int I_NC ;
// extern int I_NR ;
// extern int I_NI ;
// extern int I_NS ;
// extern int I_NG ;
constexpr size_t NQW_STR = 1; // 2 in fortran
constexpr size_t NQW_END = 5; // 6 in fortran
// values below are all +1 in fortran impl.
constexpr int I_QV =  0;
constexpr int I_QC =  1;
constexpr int I_QR =  2;
constexpr int I_QI =  3;
constexpr int I_QS =  4;
constexpr int I_QG =  5;
constexpr int I_QH =  6;
constexpr int I_NC =  6;
constexpr int I_NR =  7;
constexpr int I_NI =  8;
constexpr int I_NS =  9;
constexpr int I_NG = 10;

extern std::vector<std::string> TRC_name;

// extern double CVW[];
// extern double CPW[];
extern View1D<double, DEFAULT_MEM> CVW;
extern View1D<double, DEFAULT_MEM> CPW;

extern int SET_iteration;
extern bool SET_check;

constexpr int wk_nmax = 49;

}
