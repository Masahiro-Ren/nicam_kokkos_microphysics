#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

#include <string>

const int H_SHORT   = 32;
const int H_LONG    = 1024;
const int ADM_prc_all = 1;
const int IO_FREAD = 0;
const int ADM_KNONE = 1;

int IO_FID_LOG = 6;

int ADM_gall_in_orig = 16641;
int ADM_gall_in = 16641;
int ADM_vlayer = 94;
int ADM_kall = 96;
int ADM_kmin = 2;
int ADM_kmax = 95;
int ADM_lall = 1;
int TRC_VMAX = 6;

int kdim = 96;
int kmin = 2;
int kmax = 95;
int knone = 1;
int nqmax = 6;

double TIME_DTL = 6.0;
double TIME_CTIME = 0.0;
int TIME_CSTEP = 49;
int SET_l = 1;
int SET_rgnid = 4;

const double CONST_PI = 3.141592653589793;
const double CONST_EPS = 2.220446E-16;
const double CONST_UNDEF = -9.9999E+30;
const double CONST_GRAV = 9.79764;
const double CONST_STB = 5.670373E-8;

const double CONST_Rdry = 287.04;
const double CONST_CPdry = 1004.64;
const double CONST_CVdry = 717.56;
const double CONST_Rvap = 461.46;
const double CONST_CPvap = 1846.0;
const double CONST_CVvap = 1384.54;
const double CONST_CL = 4218.0;
const double CONST_CI = 2006.0;

const double CONST_KAPPA  = 0.2857256619550070;
const double CONST_EPSvap = 0.6219718309859156;
const double CONST_PSAT0  =     610.7;
const double CONST_LHV0   = 2501000.0;
const double CONST_LHS0   = 2834000.0;
const double CONST_LHF0   =  333000.0;
const double CONST_LHV00  = 3148911.8;
const double CONST_LHS00  = 2877704.0;
const double CONST_LHF00  = -271207.8;
const double CONST_Pstd   =  101325.0;
const double CONST_PRE00  =  100000.0;
const double CONST_Tstd   =    288.15;
const double CONST_TEM00  =    273.15;
const double CONST_TMELT  =    273.15;
const double CONST_DWATR  =    1000.0;
const double CONST_DICE   =     916.8;
const double CONST_LHV    = CONST_LHV0;
const double CONST_LHS    = CONST_LHS0;
const double CONST_LHF    = CONST_LHF0;

double PI    = 3.141592653589793; //< pi
double EPS   = 2.220446E-16;      //< small number
const double UNDEF = CONST_UNDEF;
const double GRAV  = CONST_GRAV;
const double Rdry  = CONST_Rdry;
const double Rvap  = CONST_Rvap;
const double CL    = CONST_CL;
double LHV   = CONST_LHV0;
double LHS   = CONST_LHS0;
double LHF   = CONST_LHF0;
const double rho_w = CONST_DWATR;
const double rhow  = CONST_DWATR;
const double DWATR = CONST_DWATR;
const double DICE  = CONST_DICE;
const double Pstd  = CONST_Pstd;
const double Tstd  = CONST_Tstd;
const double TEM00 = CONST_TEM00;

// grd, gmtr, vmtr
const char* vgrid_fname = "./vgrid94.dat";

double* GRD_gz   ;
double* GRD_gzh  ;
double* GRD_dgz  ;
double* GRD_dgzh ;
double* GRD_rdgz ;
double* GRD_rdgzh;
double* GRD_afact;
double* GRD_bfact;
double* GRD_cfact;
double* GRD_dfact;

// run conf
const char* EIN_TYPE            = "SIMPLE";
const char* RAIN_TYPE           = "COLD";
const char* MP_TYPE             = "NSW6";
const int MP_DIV_NUM          = 1;
const bool opt_2moment_water   = false;
const bool ISOTOPE             = false;
const bool opt_offline_aerosol = false;
const bool opt_aerosol_forcing = false;
const char* PRCIP_TRN_ECORRECT  = "KIN2EIN";

const char* RD_TYPE             = "MSTRNX";
const char* AE_TYPE             = "NONE";
const int KAPCL               = 7;
const int NCRF                = 2;
const int NRBND               = 3;
const int NRBND_VIS           = 1;
const int NRBND_NIR           = 2;
const int NRBND_IR            = 3;
const int NRDIR               = 2;
const int NRDIR_DIRECT        = 1;
const int NRDIR_DIFFUSE       = 2;
const int NTAU                = 7;
const int NPRES               = 7;
const int HYDRO_MAX           = 7;

int NQW_STR = 2;
int NQW_END = 6;
int I_QV =  1;
int I_QC =  2;
int I_QR =  3;
int I_QI =  4;
int I_QS =  5;
int I_QG =  6;
int I_QH =  7;
int I_NC =  7;
int I_NR =  8;
int I_NI =  9;
int I_NS = 10;
int I_NG = 11;

// not sure ??
const std::string TRC_name[6] = {
    "QV",
    "QC",
    "QR",
    "QI",
    "QS",
    "QG"
};

double CVW[6] = {
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry
};

double CPW[6] = {
    CONST_CPdry,
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry,
    CONST_CVdry,
};

int SET_iteration = 1;
bool SET_check     = true;


#endif