#include "mod_debug.h"
#include "limits"

using namespace PROBLEM_SIZE;
using namespace DATA_IO;

namespace DEBUG {

    int EX_STEP = 49;
    int EX_rgnid;

    int sat_ite_sum = 0;
    int sat_ite_count = 0;
    int sat_ite_max = 0;
    int sat_ite_min = -1;

    void ADM_Proc_stop()
    {
        std::cerr << "Process is stopped unpeacefully." << std::endl;
        exit(1);
    }

    void PROF_val_check(const std::string& val_name, double arr2d[ADM_kall][ADM_gall_in], double CHECK_arr2d[ADM_kall][ADM_gall_in])
    {
        // std::cout << __PRETTY_FUNCTION__ << std::endl;

        double err_sum = 0.0;
        double err_max = std::numeric_limits<double>::min();
        double err_min = std::numeric_limits<double>::max();

        for(int k = 0; k < ADM_kall; k++)
        {
            for(int ij = 0; ij < ADM_gall_in; ij++)
            {
                double err = std::abs( ( arr2d[k][ij] - CHECK_arr2d[k][ij] ) );
                err_sum += err;
                err_max = std::max(err_max, err);
                err_min = std::min(err_min, err);
            }
        }

        std::cout << "Checking [" << val_name << "] ";
        std::cout << "Max = " << std::setprecision(16) << std::scientific << err_max << "; ";
        std::cout << "Min = " << std::setprecision(16) << std::scientific << err_min << "; ";
        std::cout << "Sum = " << std::setprecision(16) << std::scientific << err_sum << "; " << std::endl;
    }

    void PROF_val_check(const std::string& val_name, double arr3d[ADM_lall][ADM_kall][ADM_gall_in], double CHECK_arr3d[ADM_lall][ADM_kall][ADM_gall_in])
    {
        // std::cout << __PRETTY_FUNCTION__ << std::endl;

        double err_sum = 0.0;
        double err_max = std::numeric_limits<double>::min();
        double err_min = std::numeric_limits<double>::max();

        for(int l = 0; l < ADM_lall; l++)
        {
            for(int k = 0; k < ADM_kall; k++)
            {
                for(int ij = 0; ij < ADM_gall_in; ij++)
                {
                    double err = std::abs( ( arr3d[l][k][ij] - CHECK_arr3d[l][k][ij] ) );
                    err_sum += err;
                    err_max = std::max(err_max, err);
                    err_min = std::min(err_min, err);
                }
            }
        }

        std::cout << "Checking [" << val_name << "] ";
        std::cout << "Max = " << std::setprecision(16) << std::scientific << err_max << "; ";
        std::cout << "Min = " << std::setprecision(16) << std::scientific << err_min << "; ";
        std::cout << "Sum = " << std::setprecision(16) << std::scientific << err_sum << "; " << std::endl;
    }

    void CVW_CPW_Setup()
    {
        auto h_CVW = Kokkos::create_mirror_view(CVW);
        auto h_CPW = Kokkos::create_mirror_view(CPW);

        for(int i = 0; i < 6; i++)
        {
            h_CVW(i) = CONST_CVdry;
            h_CPW(i) = CONST_CVdry;
        }
        h_CPW(0) = CONST_CPdry;

        Kokkos::deep_copy(CVW, h_CVW);
        Kokkos::deep_copy(CPW, h_CPW);
    }

    void GRD_Setup()
    {
        std::cout << __PRETTY_FUNCTION__ << std::endl;

        // host mirror of GRD
        auto h_GRD_gz    = Kokkos::create_mirror_view(GRD_gz   );
        auto h_GRD_gzh   = Kokkos::create_mirror_view(GRD_gzh  );
        auto h_GRD_dgz   = Kokkos::create_mirror_view(GRD_dgz  );
        auto h_GRD_dgzh  = Kokkos::create_mirror_view(GRD_dgzh );
        auto h_GRD_rdgz  = Kokkos::create_mirror_view(GRD_rdgz );
        auto h_GRD_rdgzh = Kokkos::create_mirror_view(GRD_rdgzh);
        auto h_GRD_afact = Kokkos::create_mirror_view(GRD_afact);
        auto h_GRD_bfact = Kokkos::create_mirror_view(GRD_bfact);
        auto h_GRD_cfact = Kokkos::create_mirror_view(GRD_cfact);
        auto h_GRD_dfact = Kokkos::create_mirror_view(GRD_dfact);

        // Reading data from vgrid94.dat
        // GRD_Input_vgrid(vgrid_fname);
        GRD_Input_vgrid();

        // calculation of grid intervals ( cell center )
        for(int k = ADM_kmin - 1; k <= ADM_kmax; k++)
        {
            arrGRD_dgz[k] = arrGRD_gzh[k+1] - arrGRD_gzh[k];
        }
        arrGRD_dgz[ADM_kmax + 1] = arrGRD_dgz[ADM_kmax];

        // calculation of grid intervals ( cell wall )
        for(int k = ADM_kmin; k <= ADM_kmax + 1; k++)
        {
            arrGRD_dgzh[k] = arrGRD_gz[k] - arrGRD_gz[k-1];
        }
        arrGRD_dgzh[ADM_kmin - 1] = arrGRD_dgzh[ADM_kmin];

        // calculation of 1/dgz and 1/dgzh
        for(int k = 0; k < ADM_kall; k++)
        {
            arrGRD_rdgz[k]  = 1.0 / arrGRD_dgz [k];
            arrGRD_rdgzh[k] = 1.0 / arrGRD_dgzh[k];
        }

        //---< vertical interpolation factor >---

        // vertical interpolation factor
        for(int k = ADM_kmin; k <= ADM_kmax + 1; k++)
        {
            arrGRD_afact[k] = ( arrGRD_gzh[k] - arrGRD_gz[k-1] ) /
                              ( arrGRD_gz[k] -  arrGRD_gz[k-1] );
        }
        arrGRD_afact[ADM_kmin-1] = 1.0;

        for(int k = 0; k < ADM_kall; k++)
        {
            arrGRD_bfact[k] = 1.0 - arrGRD_afact[k];
        }

        for(int k = ADM_kmin; k <= ADM_kmax; k++)
        {
            arrGRD_cfact[k] =   ( arrGRD_gz[k] - arrGRD_gzh[k] )
                              / ( arrGRD_gzh[k+1] - arrGRD_gzh[k] );
        }
        arrGRD_cfact[ADM_kmin - 1] = 1.0;
        arrGRD_cfact[ADM_kmax + 1] = 0.0;

        for(int k = 0; k < ADM_kall; k++)
        {
            arrGRD_dfact[k] = 1.0 - arrGRD_cfact[k];
        }

        // copy data to host views
        for(size_t k = 0; k < ADM_kall; k++)
        {
            h_GRD_gz   (k) = arrGRD_gz   [k];
            h_GRD_gzh  (k) = arrGRD_gzh  [k];
            h_GRD_dgz  (k) = arrGRD_dgz  [k];
            h_GRD_dgzh (k) = arrGRD_dgzh [k];
            h_GRD_rdgz (k) = arrGRD_rdgz [k];
            h_GRD_rdgzh(k) = arrGRD_rdgzh[k];
            h_GRD_afact(k) = arrGRD_afact[k];
            h_GRD_bfact(k) = arrGRD_bfact[k];
            h_GRD_cfact(k) = arrGRD_cfact[k];
            h_GRD_dfact(k) = arrGRD_dfact[k];
        }

        Kokkos::deep_copy(GRD_gz   , h_GRD_gz   );
        Kokkos::deep_copy(GRD_gzh  , h_GRD_gzh  );
        Kokkos::deep_copy(GRD_dgz  , h_GRD_dgz  );
        Kokkos::deep_copy(GRD_dgzh , h_GRD_dgzh );
        Kokkos::deep_copy(GRD_rdgz , h_GRD_rdgz );
        Kokkos::deep_copy(GRD_rdgzh, h_GRD_rdgzh);
        Kokkos::deep_copy(GRD_afact, h_GRD_afact);
        Kokkos::deep_copy(GRD_bfact, h_GRD_bfact);
        Kokkos::deep_copy(GRD_cfact, h_GRD_cfact);
        Kokkos::deep_copy(GRD_dfact, h_GRD_dfact);
        Kokkos::fence();
    }

    void GRD_Input_vgrid()
    {
        /** Waiting for implementing */
        std::cout << __PRETTY_FUNCTION__ << " Reading vgrid data " << std::endl;
        read_data_1d("data/vgrid/GRD_gz.dat", arrGRD_gz);
        read_data_1d("data/vgrid/GRD_gzh.dat", arrGRD_gzh);
    }

    void cnvvar_rhogkin_in(
        double rhog     [kdim][ijdim],
        double rhogvx   [kdim][ijdim],
        double rhogvy   [kdim][ijdim],
        double rhogvz   [kdim][ijdim],
        double rhogw    [kdim][ijdim],
        double C2Wfact  [2][kdim][ijdim],
        double W2Cfact  [2][kdim][ijdim],
        double rhogkin  [kdim][ijdim],
        double rhogkin_h[kdim][ijdim],
        double rhogkin_v[kdim][ijdim] )
    {
#ifdef DEBUG
        std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

        size_t gall = ijdim;

        // --- horizontal kinetic energy
        for(int k = kmin; k <= kmax; k++)
        {
            for(int g = 0; g < gall; g++)
            {
                rhogkin_h[k][g] = 0.5 * ( rhogvx[k][g] * rhogvx[k][g] +
                                          rhogvy[k][g] * rhogvy[k][g] +
                                          rhogvz[k][g] * rhogvz[k][g] ) / rhog[k][g];
            }
        }

        for(int g = 0; g < gall; g++)
        {
            rhogkin_h[kmin - 1][g] = 0.0;
            rhogkin_h[kmax + 1][g] = 0.0;
        }

        // --- vertical kinetic energy
        for(int k = kmin + 1; k <= kmax; k++)
        {
            for(int g = 0; g < gall; g++)
            {
                rhogkin_v[k][g] = 0.5 * ( rhogw[k][g] * rhogw[k][g] ) /
                                        ( C2Wfact[0][k][g] * rhog[k][g] +
                                          C2Wfact[1][k][g] * rhog[k-1][g] );
            }
        }

        for(int g = 0; g < gall; g++)
        {
            rhogkin_v[kmin - 1][g] = 0.0;
            rhogkin_v[kmin][g]     = 0.0;
            rhogkin_v[kmax + 1][g] = 0.0;
        }

        // -- total kinetic energy
        for(int k = kmin; k <= kmax; k++)
        {
            for(int g = 0; g < gall; g++)
            {
                rhogkin[k][g] = rhogkin_h[k][g] +                            // horizontal
                                ( W2Cfact[0][k][g] * rhogkin_v[k + 1][g] +   // vertical
                                  W2Cfact[1][k][g] * rhogkin_v[k][g] );
            }
        }

        for(int g = 0; g < gall; g++)
        {
            rhogkin[kmin - 1][g] = 0.0;
            rhogkin[kmax + 1][g] = 0.0;
        }

    }

    void cnvvar_rhogkin_in(
        const View<double**>&  rhog     ,
        const View<double**>&  rhogvx   ,
        const View<double**>&  rhogvy   ,
        const View<double**>&  rhogvz   ,
        const View<double**>&  rhogw    ,
        const View<double***>& C2Wfact  ,
        const View<double***>& W2Cfact  ,
        const View<double**>&  rhogkin  ,
        const View<double**>&  rhogkin_h,
        const View<double**>&  rhogkin_v )
    {
#ifdef DEBUG
        std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

        size_t gmin = 0;
        size_t gall = ijdim;

        // --- horizontal kinetic energy
        Kokkos::parallel_for("horizontal kinetic energy", 
                MDRangePolicy<Kokkos::Rank<2>>({kmin, gmin},{kmax+1, gall}), 
                KOKKOS_LAMBDA(const size_t k, const size_t g){
                    rhogkin_h(k,g) = 0.5 * ( rhogvx(k,g) * rhogvx(k,g) +
                                             rhogvy(k,g) * rhogvy(k,g) +
                                             rhogvz(k,g) * rhogvz(k,g) ) / rhog(k,g);
                });

        Kokkos::parallel_for(RangePolicy<>(0,gall), KOKKOS_LAMBDA(const size_t g){
            rhogkin_h(kmin - 1,g) = 0.0;
            rhogkin_h(kmax + 1,g) = 0.0;
        });

        // // --- vertical kinetic energy
        Kokkos::parallel_for("horizontal kinetic energy", 
                MDRangePolicy<Kokkos::Rank<2>>({kmin+1, gmin},{kmax+1, gall}), 
                KOKKOS_LAMBDA(const size_t k, const size_t g){
                    rhogkin_v(k,g) = 0.5 * ( rhogw(k,g) * rhogw(k,g) ) /
                                            ( C2Wfact(0,k,g) * rhog(k,g) +
                                              C2Wfact(1,k,g) * rhog(k-1,g) );
                });

        Kokkos::parallel_for(RangePolicy<>(0,gall), KOKKOS_LAMBDA(const size_t g){
            rhogkin_v(kmin - 1,g) = 0.0;
            rhogkin_v(kmin,g)     = 0.0;
            rhogkin_v(kmax + 1,g) = 0.0;
        });

        // // -- total kinetic energy
        Kokkos::parallel_for("horizontal kinetic energy", 
                MDRangePolicy<Kokkos::Rank<2>>({kmin, gmin},{kmax+1, gall}), 
                KOKKOS_LAMBDA(const size_t k, const size_t g){
                    rhogkin(k,g) = rhogkin_h(k,g) +                            // horizontal
                                   ( W2Cfact(0,k,g) * rhogkin_v(k + 1,g) +   // vertical
                                     W2Cfact(1,k,g) * rhogkin_v(k,g) );
                });

        Kokkos::parallel_for(RangePolicy<>(0,gall), KOKKOS_LAMBDA(const size_t g){
            rhogkin(kmin - 1,g) = 0.0;
            rhogkin(kmax + 1,g) = 0.0;
        });
    }

    double MISC_gammafunc(double xx)
    {
        double f;
        double coef[6] = {
            +76.18009172947146,
            -86.50532032941677,
            +24.01409824083091,
            -1.231739572450155,
            +0.1208650973866179E-2,
            -0.5395239384953E-5
        };

        double x, y, tmp, ser;

        x = xx;
        y = x;
        tmp = x + 5.5;
        tmp = tmp - ( x + 0.5 ) * std::log(tmp);

        ser = 1.000000000190015;
        for(int j = 0; j < 6; j++)
        {
            y += 1;
            ser += coef[j] / y;
        }

        f = std::exp(-tmp + std::log(2.5066282746310005 * ser / x));

        return f;
    }

    void PROF_val_check(const std::string& val_name, const View<double****>& arr4d, const size_t idx_arr2d, const View<double***>& CHECK_arr3d)
    {
        auto sub_arr2d = subview(arr4d, 0, idx_arr2d, Kokkos::ALL(), Kokkos::ALL());
        auto sub_check_arr2d = subview(CHECK_arr3d, 0, Kokkos::ALL(), Kokkos::ALL());

        double err_sum;
        double err_max;
        double err_min;

        Kokkos::parallel_reduce(MDRangePolicy<Kokkos::Rank<2>>({0,0},{ADM_kall,ADM_gall_in}),
        KOKKOS_LAMBDA(const size_t k, const size_t ij, double& local_sum, double& local_min, double& local_max){
            double err = std::abs(sub_arr2d(k,ij) - sub_check_arr2d(k,ij));
            local_sum += err;
            local_min = std::min(local_min, err);
            local_max = std::max(local_max, err);
        }, err_sum, Kokkos::Min<double>(err_min), Kokkos::Max<double>(err_max));

        std::cout << "Checking [" << val_name << "] ";
        std::cout << "Max = " << std::setprecision(16) << std::scientific << err_max << "; ";
        std::cout << "Min = " << std::setprecision(16) << std::scientific << err_min << "; ";
        std::cout << "Sum = " << std::setprecision(16) << std::scientific << err_sum << "; " << std::endl;
    }

    void PROF_val_check(const std::string& val_name, const View<double***>& arr3d, const View<double***>& CHECK_arr3d)
    {
        double err_sum;
        double err_max;
        double err_min;

        Kokkos::parallel_reduce(MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{ADM_lall,ADM_kall,ADM_gall_in}),
        KOKKOS_LAMBDA(const size_t l, const size_t k, const size_t ij, double& local_sum, double& local_min, double& local_max){
            double err = std::abs(arr3d(l,k,ij) - CHECK_arr3d(l,k,ij));
            local_sum += err;
            local_min = std::min(local_min, err);
            local_max = std::max(local_max, err);
        }, err_sum, Kokkos::Min<double>(err_min), Kokkos::Max<double>(err_max));

        std::cout << "Checking [" << val_name << "] ";
        std::cout << "Max = " << std::setprecision(16) << std::scientific << err_max << "; ";
        std::cout << "Min = " << std::setprecision(16) << std::scientific << err_min << "; ";
        std::cout << "Sum = " << std::setprecision(16) << std::scientific << err_sum << "; " << std::endl;
    }
}
