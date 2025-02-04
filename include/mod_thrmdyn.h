#pragma once

#include "problemsize.h"
#include "mod_debug.h"

using namespace PROBLEM_SIZE;
using namespace DEBUG;

namespace THRMDYN {

    void THRMDYN_qd(double q[nqmax][kdim][ijdim], double qd[kdim][ijdim]);

    void THRMDYN_cv(double qd[kdim][ijdim], double q[nqmax][kdim][ijdim], double cv[kdim][ijdim]);

    void THRMDYN_tempre( double ein[kdim][ijdim], 
                         double rho[kdim][ijdim],
                         double q[nqmax][kdim][ijdim],
                         double tem[kdim][ijdim],
                         double pre[kdim][ijdim] );
    
    /**
     * Kokkos ver.
     */
    void THRMDYN_qd(View<double***>& q, View<double**>& qd);

    void THRMDYN_cv(View<double**>& qd, View<double***>& q, View<double**>& cv);

    void THRMDYN_tempre( View<double**>&  ein, 
                         View<double**>&  rho,
                         View<double***>& q  ,
                         View<double**>&  tem,
                         View<double**>&  pre );
}
