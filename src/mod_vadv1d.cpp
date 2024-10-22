#include "mod_vadv1d.h"

void find_max_min(double* arr, size_t len, double& max_val, double& min_val);

namespace VADV1D {

void vadv1d_prep( double dz       [kdim],
                  double zh       [kdim],
                  double wp       [kdim][ijdim],
                  double zdis     [kdim][ijdim],
                  int kcell       [kdim][ijdim],
                  int kcell_max   [kdim],
                  int kcell_min   [kdim],
                  double dt)
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    double wh[kdim][ijdim];
    double dt2 = std::pow(dt, 2);
    double dt3 = std::pow(dt, 3);
    // double zzmax, zzmin;

    // vetical velocity at the half level
    for(int k = kmin + 1; k <= kmax; k++)
        for(int ij = 0; ij < ijdim; ij++)
            wh[k][ij] = 0.5 * (wp[k-1][ij] + wp[k][ij]);
    // bottom boundary for wh
    // top    boundary for wh : same as inner region
    for(int ij = 0; ij < ijdim; ij++)
    {
        wh[kmin - 1][ij] = wp[kmin - 1][ij];
        wh[kmin][ij]     = wp[kmin][ij];
        wh[kmax + 1][ij] = wp[kmax + 1][ij];
    }

    // calculation of distance of cell wall during dt
    for(int k = kmin + 1; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            zdis[k][ij] = dt  * wh[k][ij] -
                          dt2 * wh[k][ij] * ( wh[k+1][ij] - wh[k-1][ij] ) / ( dz[k-1] + dz[k] ) / 2.0 +
                          dt3 * wh[k][ij] * ( std::pow( ( wh[k+1][ij] - wh[k-1][ij] ) / ( dz[k-1] + dz[k] ), 2) +
                                              wh[k][ij] * ( ( ( wh[k+1][ij] - wh[k][ij] ) / dz[k] - 
                                                              ( wh[k][ij] - wh[k-1][ij] ) / dz[k-1] ) / (dz[k-1] + dz[k]) * 2.0 ) ) / 6.0; 
        }
    }

    // bottom and top boundary for zdis
    for(int ij = 0; ij < ijdim; ij++)
    {
        zdis[kmin - 1][ij] = 0.0;
        zdis[kmin][ij] = dt  * wh[kmin][ij] -
                         dt2 * wh[kmin][ij] * ( wh[kmin+1][ij] - wh[kmin][ij] ) / dz[kmin] / 2.0;
        zdis[kmax][ij] = dt  * wh[kmax][ij] -
                         dt2 * wh[kmax][ij] * ( wh[kmax+1][ij] - wh[kmax][ij] ) / dz[kmax] / 2.0;
    }

    // calculation of kcell
    // top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            kcell[k][ij] = k;
            kcell_min[k] = k;
            kcell_max[k] = k;
        }
    }

    // setup limiter of max and min of kcell
    for(int k = kmin; k <= kmax; k++)
    {
        // find_max_min(zdis[k], ijdim, zzmax, zzmin);
        double zzmax = *std::max_element(zdis[k], zdis[k] + ijdim);
        double zzmin = *std::min_element(zdis[k], zdis[k] + ijdim);

        if(zzmax > 0.0)
        {
            for(int k2 = k; k2 >= kmin; k2--)
            {
                if( (zh[k2] <= zh[k] - zzmax) && (zh[k2+1] > zh[k] - zzmax) )
                {
                    kcell_min[k] = k2;
                    break;
                } 
            }
        }

        if(zzmin < 0.0)
        {
            for(int k2 = k; k <= kmax; k2++)
            {
                if( (zh[k2] <= zh[k] - zzmin) && (zh[k2+1] > zh[k] - zzmin) )
                {
                    kcell_max[k] = k2;
                }
            } 
        }
    }

    // determine the kcell at each point.
    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            if(kcell_min[k] == k && kcell_max[k] == k)
            {
                kcell[k][ij] = k;
            }
            else
            {
                kcell[k][ij] = 0;
                for(int k2 = kcell_min[k]; k2 <= kcell_max[k]; k2++)
                {
                    int tmp = int(k2 * std::copysign(1.0, (zh[k]-zdis[k][ij])-zh[k2]) * std::copysign(1.0, zh[k2+1] - (zh[k] - zdis[k][ij])) );
                    kcell[k][ij] = std::max(kcell[k][ij], tmp);
                }
            }
        }
    }

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            if(kcell[k][ij] == 0)
                kcell[k][ij] = kmin;
            if(kcell[k][ij] == kmax + 1)
                kcell[k][ij] = kmax;
        }
    }
}

void vadv1d_getflux_new( double dz       [kdim],
                         double rhof     [kdim][ijdim],
                         double zdis0    [kdim][ijdim],
                         double kcell    [kdim][ijdim],
                         double kcell_max[kdim],
                         double kcell_min[kdim],
                         double frhof    [kdim][ijdim] )
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    /* TO DO */
}

}

void find_max_min(double* arr, size_t len, double& max_val, double& min_val)
{
    max_val = arr[0];
    min_val = arr[0];

    if(len == 1)
    {
        return;
    }

    for(int i = 0; i < len; i++)
    {
        if(arr[i] > max_val) max_val = arr[i];
        if(arr[i] < min_val) min_val = arr[i];
    }
}