#include "mod_vadv1d.h"

// void find_max_min(double* arr, size_t len, double& max_val, double& min_val);

namespace VADV1D {

void vadv1d_prep( int    mkmin,
                  int    mkmax,
                  double dz       [kdim],
                  double zh       [kdim],
                  double wp       [kdim][ijdim],
                  double zdis     [kdim][ijdim],
                  int kcell       [kdim][ijdim],
                  int kcell_max   [kdim],
                  int kcell_min   [kdim],
                  double dt)
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double wh[kdim][ijdim];
    
    double dt2 = std::pow(dt, 2);
    double dt3 = std::pow(dt, 3);
    // double zzmax, zzmin;

    // vetical velocity at the half level
    for(int k = mkmin + 1; k <= mkmax; k++)
        for(int ij = 0; ij < ijdim; ij++)
            wh[k][ij] = 0.5 * (wp[k-1][ij] + wp[k][ij]);
    

    // bottom boundary for wh
    // top    boundary for wh : same as inner region
    for(int ij = 0; ij < ijdim; ij++)
    {
        wh[mkmin  ][ij] = wp[mkmin  ][ij];
        wh[mkmin-1][ij] = wp[mkmin-1][ij];
        wh[mkmax+1][ij] = wp[mkmax  ][ij]; 
    }

    // calculation of distance of cell wall during dt
    for(int k = mkmin + 1; k <= mkmax; k++)
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
       zdis[mkmin-1][ij] = 0.0;
       zdis[mkmin  ][ij] = dt * wh[mkmin][ij] - dt2 * wh[mkmin][ij] * ( wh[mkmin+1][ij] - wh[mkmin][ij] ) / dz[mkmin] / 2.0;
       zdis[mkmax+1][ij] = dt * wh[mkmax+1][ij] - dt2 * wh[mkmax+1][ij] * ( wh[mkmax+1][ij] - wh[mkmax][ij] ) / dz[mkmax] / 2.0;
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
    for(int k = mkmin; k <= mkmax; k++)
    {
        // find_max_min(zdis[k], ijdim, zzmax, zzmin);
        double zzmax = *std::max_element(zdis[k], zdis[k] + ijdim);
        double zzmin = *std::min_element(zdis[k], zdis[k] + ijdim);

        if(zzmax > 0.0)
        {
            for(int k2 = k; k2 >= mkmin; k2--)
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
            for(int k2 = k; k2 <= mkmax; k2++)
            {
                if( (zh[k2] <= zh[k] - zzmin) && (zh[k2+1] > zh[k] - zzmin) )
                {
                    kcell_max[k] = k2;
                    break;
                }
            } 
        }
    }

    // determine the kcell at each point.
    for(int k = mkmin; k <= mkmax; k++)
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
                kcell[k][ij] = mkmin;
            if(kcell[k][ij] == mkmax + 1)
                kcell[k][ij] = mkmax;
        }
    }
}

void vadv1d_getflux_new( int    mkmin,
                         int    mkmax,
                         double dz       [kdim],
                         double rhof     [kdim][ijdim],
                         double zdis0    [kdim][ijdim],
                         int    kcell    [kdim][ijdim],
                         int    kcell_max[kdim],
                         int    kcell_min[kdim],
                         double frhof    [kdim][ijdim] )
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double zdis[kdim][ijdim];
    double fact;

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            frhof[k][ij] = 0.0;
        }
    }

    for(int k = mkmin; k <= mkmax; k++)
    {
        if( kcell_min[k] == k && kcell_max[k] == k )
        {
            for(int ij = 0; ij < ijdim; ij++)
                zdis[k][ij] = zdis0[k][ij];
        }
        else
        {
            for(int k2 = kcell_min[k]; k2 <= kcell_max[k]; k2++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    // int kc1 = kcell[k][ij] + 1;
                    // int kc2 = kcell[k][ij] - 1;
                    // fact = dz[k2] * 0.25 
                    //         * ( (std::copysign(1, k2 - kc1 + 1) + 1.0) * ( std::copysign(1, (k - 1) - k2) + 1.0) 
                    //             - (std::copysign(1, k2 - k) + 1.0) * (std::copysign(1, (kc2 - 1) - k2) + 1.0) );
                    fact = dz[k2] * 0.25
                                  * (  ( std::copysign(1, k2 - (kcell[k][ij]+1)) + 1.0 ) * ( std::copysign(1, (k-1) - k2) + 1.0 )
                                     - ( std::copysign(1, k2 - k) + 1.0 ) * ( std::copysign(1, (kcell[k][ij] - 1) - k2) + 1.0 ) );
                    
                    frhof[k][ij] = frhof[k][ij] + rhof[k2][ij] * fact;
                    zdis[k][ij] = zdis0[k][ij] - fact;
                }

            }
        }

        for(int ij = 0; ij < ijdim; ij++)
        {
            int kc = kcell[k][ij];
            frhof[k][ij] = frhof[k][ij] + rhof[kc][ij] * zdis[k][ij];
        }
    }

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            // double val_abs = std::abs(frhof[k][ij]) - CONST_EPS;
            // frhof[k][ij] = frhof[k][ij] * ( 0.5 + std::copysign(0.5, val_abs));
            frhof[k][ij] = frhof[k][ij] * ( 0.5 + std::copysign(0.5, std::abs(frhof[k][ij]) - CONST_EPS ) ); // small negative filter
        }
    }
}


// void find_max_min(double* arr, size_t len, double& max_val, double& min_val)
// {
//     max_val = arr[0];
//     min_val = arr[0];

//     if(len == 1)
//     {
//         return;
//     }

//     for(int i = 0; i < len; i++)
//     {
//         if(arr[i] > max_val) max_val = arr[i];
//         if(arr[i] < min_val) min_val = arr[i];
//     }
// }

void vadv1d_prep( int    mkmin,
                  int    mkmax,
                  const View1D<double, DEFAULT_MEM>&  dz       ,
                  const View1D<double, DEFAULT_MEM>&  zh       ,
                  const View2D<double, DEFAULT_MEM>&  wp       ,
                  const View2D<double, DEFAULT_MEM>&  zdis     ,
                  const View2D<int,    DEFAULT_MEM>&  kcell    ,
                  const View1D<int,    DEFAULT_MEM>&  kcell_max,
                  const View1D<int,    DEFAULT_MEM>&  kcell_min,
                  double dt)
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    // Copy data to host
    auto h_dz = Kokkos::create_mirror_view(dz);
    auto h_zh = Kokkos::create_mirror_view(zh);
    auto h_wp = Kokkos::create_mirror_view(wp);
    auto h_zdis = Kokkos::create_mirror_view(zdis);
    auto h_kcell = Kokkos::create_mirror_view(kcell);
    auto h_kcell_max = Kokkos::create_mirror_view(kcell_max);
    auto h_kcell_min = Kokkos::create_mirror_view(kcell_min);

    Kokkos::deep_copy(h_dz, dz);
    Kokkos::deep_copy(h_zh, zh);
    Kokkos::deep_copy(h_wp, wp);
    Kokkos::deep_copy(h_zdis, zdis);
    Kokkos::deep_copy(h_kcell, kcell);
    Kokkos::deep_copy(h_kcell_max, kcell_max);
    Kokkos::deep_copy(h_kcell_min, kcell_min);

    View2D<double, HOST_MEM> wh("wh", kdim, ijdim);
    
    double dt2 = std::pow(dt, 2);
    double dt3 = std::pow(dt, 3);
    // double zzmax, zzmin;

    // vetical velocity at the half level
    for(size_t k = mkmin + 1; k <= mkmax; k++)
        for(size_t ij = 0; ij < ijdim; ij++)
            wh(k,ij) = 0.5 * (h_wp(k-1,ij) + h_wp(k,ij));

    // bottom boundary for wh
    // top    boundary for wh : same as inner region
    for(size_t ij = 0; ij < ijdim; ij++)
    {
        wh(mkmin  ,ij) = h_wp(mkmin  ,ij);
        wh(mkmin-1,ij) = h_wp(mkmin-1,ij);
        wh(mkmax+1,ij) = h_wp(mkmax  ,ij); 
    }

    // calculation of distance of cell wall during dt
    for(size_t k = mkmin + 1; k <= mkmax; k++)
    {
        for(size_t ij = 0; ij < ijdim; ij++)
        {
            h_zdis(k,ij) = dt  * wh(k,ij) -
                           dt2 * wh(k,ij) * ( wh(k+1,ij) - wh(k-1,ij) ) / ( h_dz(k-1) + h_dz(k) ) / 2.0 +
                           dt3 * wh(k,ij) * ( std::pow( ( wh(k+1,ij) - wh(k-1,ij) ) / ( h_dz(k-1) + h_dz(k) ), 2) +
                                              wh(k,ij) * ( ( ( wh(k+1,ij) - wh(k,ij) ) / h_dz(k) - 
                                                             ( wh(k,ij) - wh(k-1,ij) ) / h_dz(k-1) ) / (h_dz(k-1) + h_dz(k)) * 2.0 ) ) / 6.0; 
        }
    }

    // bottom and top boundary for zdis
    for(size_t ij = 0; ij < ijdim; ij++)
    {
       h_zdis(mkmin-1,ij) = 0.0;
       h_zdis(mkmin  ,ij) = dt * wh(mkmin,ij) - dt2 * wh(mkmin,ij) * ( wh(mkmin+1,ij) - wh(mkmin,ij) ) / h_dz(mkmin) / 2.0;
       h_zdis(mkmax+1,ij) = dt * wh(mkmax+1,ij) - dt2 * wh(mkmax+1,ij) * ( wh(mkmax+1,ij) - wh(mkmax,ij) ) / h_dz(mkmax) / 2.0;
    }

    // calculation of kcell
    // top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    for(size_t k = 0; k < kdim; k++)
    {
        for(size_t ij = 0; ij < ijdim; ij++)
        {
            h_kcell(k,ij) = k;
            h_kcell_min(k) = k;
            h_kcell_max(k) = k;
        }
    }

    // setup limiter of max and min of kcell
    for(size_t k = mkmin; k <= mkmax; k++)
    {
        double zzmax = std::numeric_limits<double>::min();
        double zzmin = std::numeric_limits<double>::min();

        for(size_t ij = 0; ij < ijdim; ij++)
        {
            double val = h_zdis(k,ij);
            zzmax = val > zzmax ? val : zzmax;
            zzmin = val < zzmin ? val : zzmin;
        }

        if(zzmax > 0.0)
        {
            for(size_t k2 = k; k2 >= mkmin; k2--)
            {
                if( (h_zh(k2) <= h_zh(k) - zzmax) && (h_zh(k2+1) > h_zh(k) - zzmax) )
                {
                    h_kcell_min(k) = k2;
                    break;
                } 
            }
        }

        if(zzmin < 0.0)
        {
            for(size_t k2 = k; k2 <= mkmax; k2++)
            {
                if( (h_zh(k2) <= h_zh(k) - zzmin) && (h_zh(k2+1) > h_zh(k) - zzmin) )
                {
                    h_kcell_max(k) = k2;
                    break;
                }
            } 
        }
    }

    // determine the kcell at each point.
    for(size_t k = mkmin; k <= mkmax; k++)
    {
        for(size_t ij = 0; ij < ijdim; ij++)
        {
            if(h_kcell_min(k) == k && h_kcell_max(k) == k)
            {
                h_kcell(k,ij) = k;
            }
            else
            {
                h_kcell(k,ij) = 0;
                for(size_t k2 = h_kcell_min(k); k2 <= h_kcell_max(k); k2++)
                {
                    int tmp = int(k2 * std::copysign(1.0, (h_zh(k)-h_zdis(k,ij))-h_zh(k2)) * std::copysign(1.0, h_zh(k2+1) - (h_zh(k) - h_zdis(k,ij))) );
                    h_kcell(k,ij) = std::max(h_kcell(k,ij), tmp);
                }
            }
        }
    }

    for(size_t k = 0; k < kdim; k++)
    {
        for(size_t ij = 0; ij < ijdim; ij++)
        {
            if(h_kcell(k,ij) == 0)
                h_kcell(k,ij) = mkmin;
            if(h_kcell(k,ij) == mkmax + 1)
                h_kcell(k,ij) = mkmax;
        }
    }

    // Copy editted data to device
    Kokkos::deep_copy(zdis, h_zdis);
    Kokkos::deep_copy(kcell, h_kcell);
    Kokkos::deep_copy(kcell_max, h_kcell_max);
    Kokkos::deep_copy(kcell_min, h_kcell_min);
}

void vadv1d_getflux_new( int    mkmin,
                         int    mkmax,
                         View1D<double, DEFAULT_MEM>& dz       ,
                         View2D<double, DEFAULT_MEM>& rhof     ,
                         View2D<double, DEFAULT_MEM>& zdis0    ,
                         View2D<int,    DEFAULT_MEM>& kcell    ,
                         View1D<int,    DEFAULT_MEM>& kcell_max,
                         View1D<int,    DEFAULT_MEM>& kcell_min,
                         View2D<double, DEFAULT_MEM>& frhof     )
{
#ifdef ENABLE_DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    View2D<double, DEFAULT_MEM> zdis("zdis",kdim,ijdim);
    // double fact;

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}),
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        frhof(k,ij) = 0.0;
    });

    // for(int k = mkmin; k <= mkmax; k++)
    // {
    //     if( kcell_min(k) == k && kcell_max(k) == k )
    //     {
    //         Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
    //             zdis(k,ij) = zdis0(k,ij);
    //         });
    //     }
    //     else
    //     {
    //         for(int k2 = kcell_min(k); k2 <= kcell_max(k); k2++)
    //         {
    //             Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
    //                 double fact = dz(k2) * 0.25
    //                                      * (  ( copysign(1, k2 - (kcell(k,ij)+1)) + 1.0 ) * ( copysign(1, (k-1) - k2) + 1.0 )
    //                                         - ( copysign(1, k2 - k) + 1.0 ) * ( copysign(1, (kcell(k,ij) - 1) - k2) + 1.0 ) );
                    
    //                 frhof(k,ij) = frhof(k,ij) + rhof(k2,ij) * fact;
    //                 zdis(k,ij) = zdis0(k,ij) - fact;
    //             });

    //         }
    //     }

    //     Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
    //         int kc = kcell(k,ij);
    //         frhof(k,ij) = frhof(k,ij) + rhof(kc,ij) * zdis(k,ij);
    //     });
    // }

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({mkmin,0},{mkmax+1, ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        if( kcell_min(k) == k && kcell_max(k) == k )
        {
            zdis(k,ij) = zdis0(k,ij);
        }
        else
        {
            for(int k2 = kcell_min(k); k2 <= kcell_max(k); k2++)
            {
                double f1 = static_cast<double>(k2 - (kcell(k,ij)+1));
                double f2 = static_cast<double>((k-1) - k2);
                double f3 = static_cast<double>(k2 - k);
                double f4 = static_cast<double>((kcell(k,ij) - 1) - k2);
                // double fact = dz(k2) * 0.25
                //                      * (  ( copysign(1, k2 - (kcell(k,ij)+1)) + 1.0 ) * ( copysign(1, (k-1) - k2) + 1.0 )
                //                      - ( copysign(1, k2 - k) + 1.0 ) * ( copysign(1, (kcell(k,ij) - 1) - k2) + 1.0 ) );
                double fact = dz(k2) * 0.25
                                     * (  ( Kokkos::copysign(1.0, f1) + 1.0 ) * ( Kokkos::copysign(1.0, f2) + 1.0 )
                                     - ( Kokkos::copysign(1.0, f3) + 1.0 ) * ( Kokkos::copysign(1.0, f4) + 1.0 ) );
                
                frhof(k,ij) = frhof(k,ij) + rhof(k2,ij) * fact;
                zdis(k,ij) = zdis0(k,ij) - fact;
            }
        }

        int kc = kcell(k,ij);
        frhof(k,ij) = frhof(k,ij) + rhof(kc,ij) * zdis(k,ij);
    });

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        frhof(k,ij) = frhof(k,ij) * ( 0.5 + Kokkos::copysign(0.5, abs(frhof(k,ij)) - CONST_EPS ) ); // small negative filter
    });
}

}