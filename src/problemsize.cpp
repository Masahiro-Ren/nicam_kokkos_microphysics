#include "problemsize.h"

namespace PROBLEM_SIZE
{

double GRD_gz   [ADM_kall];
double GRD_gzh  [ADM_kall];
double GRD_dgz  [ADM_kall];
double GRD_dgzh [ADM_kall];
double GRD_rdgz [ADM_kall];
double GRD_rdgzh[ADM_kall];
double GRD_afact[ADM_kall];
double GRD_bfact[ADM_kall];
double GRD_cfact[ADM_kall];
double GRD_dfact[ADM_kall];

// define run conf
std::string EIN_TYPE     = "SIMPLE";
std::string RAIN_TYPE    = "COLD";
std::string MP_TYPE      = "NSW6";
// std::string vgrid_fname  = "./vgrid94.dat";
int MP_DIV_NUM           = 1;
bool opt_2moment_water   = false;
bool ISOTOPE             = false;
bool opt_offline_aerosol = false;
bool opt_aerosol_forcing = false;
std::string PRCIP_TRN_ECORRECT  = "KIN2EIN";

std::string RD_TYPE     = "MSTRNX";
std::string AE_TYPE     = "NONE";
int KAPCL               = 7;
int NCRF                = 2;
int NRBND               = 3;
int NRBND_VIS           = 1;
int NRBND_NIR           = 2;
int NRBND_IR            = 3;
int NRDIR               = 2;
int NRDIR_DIRECT        = 1;
int NRDIR_DIFFUSE       = 2;
int NTAU                = 7;
int NPRES               = 7;
int HYDRO_MAX           = 7;

int NQW_STR = 1; // 2 in fortran
int NQW_END = 5; // 6 in fortran
// values below are all 1 less than fortran impl.
int I_QV =  0;
int I_QC =  1;
int I_QR =  2;
int I_QI =  3;
int I_QS =  4;
int I_QG =  5;
int I_QH =  6;
int I_NC =  6;
int I_NR =  7;
int I_NI =  8;
int I_NS =  9;
int I_NG = 10;

std::vector<std::string> TRC_name{ "QV",
                                   "QC",
                                   "QR",
                                   "QI",
                                   "QS",
                                   "QG"  };

std::vector<double> CVW{    CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry, };

std::vector<double> CPW{    CONST_CPdry,
                            CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry,
                            CONST_CVdry, };

int SET_iteration = 1;
bool SET_check = true;

}
