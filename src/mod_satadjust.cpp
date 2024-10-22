#include "mod_satadjust.h"

/**
 * Some private parameters
 */
double TEM_MIN = 10.0;
double DTEM_EPS0 = 1.0E-8;

std::string ALPHA_TYPE = "LINEAR";

double SATURATION_ULIMIT_TEMP = 273.15;
double SATURATION_LLIMIT_TEMP = 233.15;

// private procedures
void satadjust_all( double rho    [kdim][ijdim],
                    double Emoist [kdim][ijdim],
                    double qsum   [kdim][ijdim],
                    double tem    [kdim][ijdim],
                    double q      [nqmax][kdim][ijdim] );

void satadjust_liq( double rho    [kdim][ijdim],
                    double Emoist [kdim][ijdim],
                    double qsum   [kdim][ijdim],
                    double tem    [kdim][ijdim],
                    double q      [nqmax][kdim][ijdim] );

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

void SATURATION_Setrange(double Tw, double Ti)
{
    SATURATION_ULIMIT_TEMP = Tw;
    SATURATION_LLIMIT_TEMP = Ti;
}

void SATURATION_psat_liq(double tem[kdim][ijdim], double psat[kdim][ijdim])
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    double rtem;

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0 = CONST_PSAT0;

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

            psat[k][ij] = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_liq) * 
                                  std::exp(LovR_liq * (RTEM00 - rtem));
        }
    }
}

void SATURATION_psat_ice(double tem[kdim][ijdim], double psat[kdim][ijdim])
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    double rtem;

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0 = CONST_PSAT0;

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

            psat[k][ij] = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_ice) * 
                                  std::exp(LovR_ice * (RTEM00 - rtem));
        }
    }
}

void SATURATION_adjustment( double rhog   [kdim][ijdim],
                            double rhoge  [kdim][ijdim],
                            double rhogq  [nqmax][kdim][ijdim],
                            double tem    [kdim][ijdim],
                            double q      [nqmax][kdim][ijdim],
                            double qd     [kdim][ijdim],
                            double gsgam2 [kdim][ijdim],
                            bool   ice_adjust )
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    double ein_mosit[kdim][ijdim];
    double qsum     [kdim][ijdim];
    double CVtot    [kdim][ijdim];
    double rho      [kdim][ijdim];

    // ein_moist = U1(rho,qsum,T1) : "unsaturated temperature"
    if(I_QI > 0 && ice_adjust)
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                ein_mosit[k][ij] = rhoge[k][ij] / rhog[k][ij] 
                                    + q[I_QV][k][ij] * LHV 
                                    - q[I_QI][k][ij] * LHF;
                qsum[k][ij] = q[I_QV][k][ij] 
                              + q[I_QC][k][ij]
                              + q[I_QI][k][ij];

                q[I_QV][k][ij] = qsum[k][ij];
                q[I_QC][k][ij] = 0.0;
                q[I_QI][k][ij] = 0.0;
            }
        }
    }
    else
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                ein_mosit[k][ij] = rhoge[k][ij] / rhog[k][ij] 
                                    + q[I_QV][k][ij] * LHV;

                qsum[k][ij] = q[I_QV][k][ij] 
                              + q[I_QC][k][ij];

                q[I_QV][k][ij] = qsum[k][ij];
                q[I_QC][k][ij] = 0.0;
            }
        }
    }

    THRMDYN_cv(qd, q, CVtot);

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rho[k][ij] = rhog[k][ij] /  gsgam2[k][ij];
            tem[k][ij] = ( ein_mosit[k][ij] - q[I_QV][k][ij] * LHV ) / CVtot[k][ij];
        }
    }

    if(I_QI > 0 && ice_adjust)
    {
        satadjust_all(rho, ein_mosit, qsum, tem, q);
    }
    else
    {
        satadjust_liq(rho, ein_mosit, qsum, tem, q);
    }


    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rhogq[I_QV][k][ij] = rhog[k][ij] * q[I_QV][k][ij];
            rhogq[I_QC][k][ij] = rhog[k][ij] * q[I_QC][k][ij];
        }
    }

    if(I_QI > 0 && ice_adjust)
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhogq[I_QI][k][ij] = rhog[k][ij] * q[I_QI][k][ij];
            }
        }
    }

    THRMDYN_cv(qd, q, CVtot);

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rhoge[k][ij] = rhog[k][ij] * tem[k][ij] * CVtot[k][ij];
        }
    }

}

}

void satadjust_all( double rho    [kdim][ijdim],
                    double Emoist [kdim][ijdim],
                    double qsum   [kdim][ijdim],
                    double tem    [kdim][ijdim],
                    double q      [nqmax][kdim][ijdim] )
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void satadjust_liq( double rho    [kdim][ijdim],
                    double Emoist [kdim][ijdim],
                    double qsum   [kdim][ijdim],
                    double tem    [kdim][ijdim],
                    double q      [nqmax][kdim][ijdim] )
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}
