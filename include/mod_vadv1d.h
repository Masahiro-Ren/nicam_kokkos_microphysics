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
                  const View1D<double, DEFAULT_MEM>&  dz       ,
                  const View1D<double, DEFAULT_MEM>&  zh       ,
                  const View2D<double, DEFAULT_MEM>&  wp       ,
                  const View2D<double, DEFAULT_MEM>&  zdis     ,
                  const View2D<int,    DEFAULT_MEM>&  kcell    ,
                  const View1D<int,    DEFAULT_MEM>&  kcell_max,
                  const View1D<int,    DEFAULT_MEM>&  kcell_min,
                  double dt);

void vadv1d_getflux_new( int    mkmin,
                         int    mkmax,
                         View1D<double, DEFAULT_MEM>& dz       ,
                         View2D<double, DEFAULT_MEM>& rhof     ,
                         View2D<double, DEFAULT_MEM>& zdis0    ,
                         View2D<int,    DEFAULT_MEM>& kcell    ,
                         View1D<int,    DEFAULT_MEM>& kcell_max,
                         View1D<int,    DEFAULT_MEM>& kcell_min,
                         View2D<double, DEFAULT_MEM>& frhof     );
}
