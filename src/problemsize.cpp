#include "problemsize.h"

// define run conf
std::string EIN_TYPE     = "SIMPLE";
std::string RAIN_TYPE    = "COLD";
std::string MP_TYPE      = "NSW6";
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

std::vector<std::string> TRC_name{ "QV",
                                   "QC",
                                   "QR",
                                   "QI",
                                   "QS",
                                   "QG"  };

std::vector<double> CVM{    CONST_CVdry,
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