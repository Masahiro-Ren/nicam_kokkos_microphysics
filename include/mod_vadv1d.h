#pragma once

#include "problemsize.h"
#include "mod_debug.h"

using namespace PROBLEM_SIZE;
using namespace DEBUG;

namespace VADV1D{

/**
 * ijdim, kdim, are all in the problemsize
 */
void vadv1d_prep( int    mkmin,
                  int    mkmax,
                  double dz       [kdim],
                  double zh       [kdim],
                  double wp       [kdim][ijdim],
                  double zdis     [kdim][ijdim],
                  int    kcell    [kdim][ijdim],
                  int    kcell_max[kdim],
                  int    kcell_min[kdim],
                  double dt);

void vadv1d_getflux_new( int    mkmin,
                         int    mkmax,
                         double dz       [kdim],
                         double rhof     [kdim][ijdim],
                         double zdis0    [kdim][ijdim],
                         int    kcell    [kdim][ijdim],
                         int    kcell_max[kdim],
                         int    kcell_min[kdim],
                         double frhof    [kdim][ijdim] );

/**
 * Kokkos ver.
 */
void vadv1d_prep( int    mkmin,
                  int    mkmax,
                  View<double*>&  dz       ,
                  View<double*>&  zh       ,
                  const View<double**>& wp       ,
                  View<double**>& zdis     ,
                  View<int**>&    kcell    ,
                  View<int*>&     kcell_max,
                  View<int*>&     kcell_min,
                  double dt);

void vadv1d_getflux_new( int    mkmin,
                         int    mkmax,
                         View<double*>&  dz       ,
                         View<double**>& rhof     ,
                         View<double**>& zdis0    ,
                         View<int**>&    kcell    ,
                         View<int*>&     kcell_max,
                         View<int*>&     kcell_min,
                         View<double**>&  frhof     );
}
