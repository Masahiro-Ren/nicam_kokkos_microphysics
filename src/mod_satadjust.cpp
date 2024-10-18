#include "mod_satadjust.h"

/**
 * Some private parameters
 */
double TEM_MIN = 10.0;
double DTEM_EPS0 = 1.0E-8;

std::string ALPHA_TYPE = "LINEAR";

double SATURATION_ULIMIT_TEMP = 273.15;
double SATURATION_LLIMIT_TEMP = 233.15;

namespace SATADJUST{

double CPovR_liq;
double CPovR_ice;
double CVovR_liq;
double CVovR_ice;
double LovR_liq;
double LovR_ice;


void SATURATION_Setup()
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    if(EIN_TYPE == "EXACT")
    {
        CPovR_liq = (CONST_CPvap - CONST_CL) / CONST_Rvap;
        CPovR_ice = (CONST_CPvap - CONST_CI) / CONST_Rvap;
        CVovR_liq = (CONST_CVvap - CONST_CL) / CONST_Rvap;
        CVovR_ice = (CONST_CVvap - CONST_CI) / CONST_Rvap;
    }
    else if (EIN_TYPE == "SIMPLE" || EIN_TYPE == "SIMPLE2")
    {
        CPovR_liq = 0.0;
        CPovR_ice = 0.0;
        CVovR_liq = 0.0;
        CVovR_ice = 0.0;
    }

    LovR_liq = LHV / CONST_Rvap;
    LovR_ice = LHS / CONST_Rvap;
}

};
