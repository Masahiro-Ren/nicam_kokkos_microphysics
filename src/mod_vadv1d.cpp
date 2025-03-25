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
    std::cout << __PRETTY_FUNCTION__ << std::endl;

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
    std::cout << __PRETTY_FUNCTION__ << std::endl;

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
                  View<double*>&  dz       ,
                  View<double*>&  zh       ,
                  const View<double**>& wp       ,
                  View<double**>& zdis     ,
                  View<int**>& kcell       ,
                  View<int*>&  kcell_max   ,
                  View<int*>&  kcell_min   ,
                  double dt)
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    View<double**> wh("wh", kdim, ijdim);
    
    double dt2 = std::pow(dt, 2);
    double dt3 = std::pow(dt, 3);
    // double zzmax, zzmin;

    // vetical velocity at the half level
    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({mkmin+1,IDX_ZERO},{mkmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        wh(k,ij) = 0.5 * (wp(k-1,ij) + wp(k,ij));
    });

    // bottom boundary for wh
    // top    boundary for wh : same as inner region
    Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
        wh(mkmin  ,ij) = wp(mkmin  ,ij);
        wh(mkmin-1,ij) = wp(mkmin-1,ij);
        wh(mkmax+1,ij) = wp(mkmax  ,ij); 
    });

    // calculation of distance of cell wall during dt
    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({mkmin+1,IDX_ZERO},{mkmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
            zdis(k,ij) = dt  * wh(k,ij) -
                         dt2 * wh(k,ij) * ( wh(k+1,ij) - wh(k-1,ij) ) / ( dz(k-1) + dz(k) ) / 2.0 +
                         dt3 * wh(k,ij) * ( std::pow( ( wh(k+1,ij) - wh(k-1,ij) ) / ( dz(k-1) + dz(k) ), 2) +
                                            wh(k,ij) * ( ( ( wh(k+1,ij) - wh(k,ij) ) / dz(k) - 
                                                           ( wh(k,ij) - wh(k-1,ij) ) / dz(k-1) ) / (dz(k-1) + dz(k)) * 2.0 ) ) / 6.0; 
    });

    // bottom and top boundary for zdis
    Kokkos::parallel_for(RangePolicy<>(0,ijdim), KOKKOS_LAMBDA(const size_t ij){
       zdis(mkmin-1,ij) = 0.0;
       zdis(mkmin  ,ij) = dt * wh(mkmin,ij) - dt2 * wh(mkmin,ij) * ( wh(mkmin+1,ij) - wh(mkmin,ij) ) / dz(mkmin) / 2.0;
       zdis(mkmax+1,ij) = dt * wh(mkmax+1,ij) - dt2 * wh(mkmax+1,ij) * ( wh(mkmax+1,ij) - wh(mkmax,ij) ) / dz(mkmax) / 2.0;
    });

    // calculation of kcell
    // top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}),
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        kcell(k,ij) = k;
    });

    Kokkos::parallel_for(RangePolicy<>(0,kdim), KOKKOS_LAMBDA(const size_t k){
        kcell_min(k) = k;
        kcell_max(k) = k;
    });

    // setup limiter of max and min of kcell
    Kokkos::parallel_for(RangePolicy<Kokkos::Serial>(mkmin,mkmax+1), KOKKOS_LAMBDA(const size_t k){
        double zzmax = std::numeric_limits<double>::min();
        double zzmin = std::numeric_limits<double>::min();

        for(size_t ij = 0; ij < ijdim; ij++)
        {
            double val = zdis(k,ij);
            zzmax = val > zzmax ? val : zzmax;
            zzmin = val < zzmin ? val : zzmin;
        }

        if(zzmax > 0.0)
        {
            for(int k2 = k; k2 >= mkmin; k2--)
            {
                if( (zh(k2) <= zh(k) - zzmax) && (zh(k2+1) > zh(k) - zzmax) )
                {
                    kcell_min(k) = k2;
                    break;
                } 
            }
        }

        if(zzmin < 0.0)
        {
            for(int k2 = k; k2 <= mkmax; k2++)
            {
                if( (zh(k2) <= zh(k) - zzmin) && (zh(k2+1) > zh(k) - zzmin) )
                {
                    kcell_max(k) = k2;
                    break;
                }
            } 
        }
    });

    // determine the kcell at each point.
    Kokkos::parallel_for(MDRangePolicy<Kokkos::Serial, Kokkos::Rank<2>>({mkmin,IDX_ZERO},{mkmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        if(kcell_min(k) == k && kcell_max(k) == k)
        {
            kcell(k,ij) = k;
        }
        else
        {
            kcell(k,ij) = 0;
            for(int k2 = kcell_min(k); k2 <= kcell_max(k); k2++)
            {
                int tmp = int(k2 * std::copysign(1.0, (zh(k)-zdis(k,ij))-zh(k2)) * std::copysign(1.0, zh(k2+1) - (zh(k) - zdis(k,ij))) );
                kcell(k,ij) = std::max(kcell(k,ij), tmp);
            }
        }
    });

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}),
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        if(kcell(k,ij) == 0)
            kcell(k,ij) = mkmin;
        if(kcell(k,ij) == mkmax + 1)
            kcell(k,ij) = mkmax;
    });
}

void vadv1d_getflux_new( int    mkmin,
                         int    mkmax,
                         View<double*>&  dz       ,
                         View<double**>& rhof     ,
                         View<double**>& zdis0    ,
                         View<int**>&    kcell    ,
                         View<int*>&     kcell_max,
                         View<int*>&     kcell_min,
                         View<double**>&  frhof     )
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    View<double**> zdis("zdis",kdim,ijdim);
    // double fact;

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}),
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        frhof(k,ij) = 0.0;
    });

    Kokkos::parallel_for(RangePolicy<Kokkos::Serial>(mkmin,mkmax+1), KOKKOS_LAMBDA(const size_t k){
        if( kcell_min(k) == k && kcell_max(k) == k )
        {
            for(int ij = 0; ij < ijdim; ij++)
                zdis(k,ij) = zdis0(k,ij);
        }
        else
        {
            for(int k2 = kcell_min(k); k2 <= kcell_max(k); k2++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    // int kc1 = kcell[k][ij] + 1;
                    // int kc2 = kcell[k][ij] - 1;
                    // fact = dz[k2] * 0.25 
                    //         * ( (std::copysign(1, k2 - kc1 + 1) + 1.0) * ( std::copysign(1, (k - 1) - k2) + 1.0) 
                    //             - (std::copysign(1, k2 - k) + 1.0) * (std::copysign(1, (kc2 - 1) - k2) + 1.0) );
                    double fact = dz(k2) * 0.25
                                         * (  ( std::copysign(1, k2 - (kcell(k,ij)+1)) + 1.0 ) * ( std::copysign(1, (k-1) - k2) + 1.0 )
                                            - ( std::copysign(1, k2 - k) + 1.0 ) * ( std::copysign(1, (kcell(k,ij) - 1) - k2) + 1.0 ) );
                    
                    frhof(k,ij) = frhof(k,ij) + rhof(k2,ij) * fact;
                    zdis(k,ij) = zdis0(k,ij) - fact;
                }

            }
        }

        for(int ij = 0; ij < ijdim; ij++)
        {
            int kc = kcell(k,ij);
            frhof(k,ij) = frhof(k,ij) + rhof(kc,ij) * zdis(k,ij);
        }
    });

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}),
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        frhof(k,ij) = frhof(k,ij) * ( 0.5 + std::copysign(0.5, std::abs(frhof(k,ij)) - CONST_EPS ) ); // small negative filter
    });
}
