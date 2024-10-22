#pragma once

#include "problemsize.h"
#include "mod_debug.h"

using namespace PROBLEM_SIZE;
using namespace DEBUG;

namespace VADV1D{

/**
 * ijdim, kdim, kmin, kmax are all in the problemsize
 */
void vadv1d_prep( double dz       [kdim],
                  double zh       [kdim],
                  double wp       [kdim][ijdim],
                  double zdis     [kdim][ijdim],
                  double kcell    [kdim][ijdim],
                  double kcell_max[kdim],
                  double kcell_min[kdim],
                  double dt);

void vadv1d_getflux_new( double dz       [kdim],
                         double rhof     [kdim][ijdim],
                         double zdis0    [kdim][ijdim],
                         double kcell    [kdim][ijdim],
                         double kcell_max[kdim],
                         double kcell_min[kdim],
                         double frhof    [kdim][ijdim] );

}
