#pragma once

#include "problemsize.h"
#include "mod_thrmdyn.h"

using namespace PROBLEM_SIZE;
using namespace THRMDYN;

namespace SATADJUST{

extern double CPovR_liq;
extern double CPovR_ice;
extern double CVovR_liq;
extern double CVovR_ice;
extern double LovR_liq;
extern double LovR_ice;

void SATURATION_Setup();

void SATURATION_Setrange(double Tw, double Ti);

/**
 * > calc saturation vapor pressure from Clausius-Clapeyron equation (2D) 
 * tem  => temperature [K]
 * psat => saturation vapor pressure [Pa]
 */
void SATURATION_psat_liq(double tem[kdim][ijdim], double psat[kdim][ijdim]);

/**
 * > calc saturation vapor pressure from Clausius-Clapeyron equation (2D) 
 * tem  => temperature [K]
 * psat => saturation vapor pressure [Pa]
 */
void SATURATION_psat_ice(double tem[kdim][ijdim], double psat[kdim][ijdim]);

void SATURATION_adjustment( double rhog   [kdim][ijdim],
                            double rhoge  [kdim][ijdim],
                            double rhogq  [nqmax][kdim][ijdim],
                            double tem    [kdim][ijdim],
                            double q      [nqmax][kdim][ijdim],
                            double qd     [kdim][ijdim],
                            double gsgam2 [kdim][ijdim],
                            bool   ice_adjust );

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
}
