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
    void THRMDYN_qd(const View<double***>& q, const View<double**>& qd);

    void THRMDYN_cv(const View<double**>& qd, const View<double***>& q, const View<double**>& cv);

    void THRMDYN_tempre( const View<double**>&  ein, 
                         const View<double**>&  rho,
                         const View<double***>& q  ,
                         const View<double**>&  tem,
                         const View<double**>&  pre );
}
