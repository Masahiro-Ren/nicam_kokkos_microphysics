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
        std::cout << "Checking " << val_name << " RAE " << std::endl;

        double err_sum = 0.0;
        double err_max = -1.0;
        double err_min = std::numeric_limits<double>::infinity();

        for(int k = 0; k < ADM_kall; k++)
        {
            for(int ij = 0; ij < ADM_gall_in; ij++)
            {
                double err;
                if(std::abs(CHECK_arr2d[k][ij]) > CONST_EPS)
                {
                    err = std::abs( ( arr2d[k][ij] - CHECK_arr2d[k][ij] ) / CHECK_arr2d[k][ij] );
                }
                else
                {
                    err = std::abs(arr2d[k][ij]);
                }

                err_sum += err;
                err_max = std::max(err_max, err);
                err_min = std::min(err_min, err);
            }
        }

        std::cout << "Max = " << std::setprecision(16) << std::scientific << err_max << "; ";
        std::cout << "Min = " << std::setprecision(16) << std::scientific << err_min << "; ";
        std::cout << "Sum = " << std::setprecision(16) << std::scientific << err_sum << "; " << std::endl;
    }

    void PROF_val_check(const std::string& val_name, double arr3d[ADM_lall][ADM_kall][ADM_gall_in], double CHECK_arr3d[ADM_lall][ADM_kall][ADM_gall_in])
    {
        // std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "Checking " << val_name << " RAE " << std::endl;

        double err_sum = 0.0;
        double err_max = -1.0;
        double err_min = std::numeric_limits<double>::infinity();

        for(int l = 0; l < ADM_lall; l++)
        {
            for(int k = 0; k < ADM_kall; k++)
            {
                for(int ij = 0; ij < ADM_gall_in; ij++)
                {
                    double err;
                    if(std::abs(CHECK_arr3d[l][k][ij]) > CONST_EPS)
                    {
                        err = std::abs( ( arr3d[l][k][ij] - CHECK_arr3d[l][k][ij] ) / CHECK_arr3d[l][k][ij] );
                    }
                    else
                    {
                        err = std::abs(arr3d[l][k][ij]);
                    }

                    err_sum += err;
                    err_max = std::max(err_max, err);
                    err_min = std::min(err_min, err);
                }
            }
        }

        std::cout << "Max = " << std::setprecision(16) << std::scientific << err_max << "; ";
        std::cout << "Min = " << std::setprecision(16) << std::scientific << err_min << "; ";
        std::cout << "Sum = " << std::setprecision(16) << std::scientific << err_sum << "; " << std::endl;
    }

    void GRD_Setup()
    {
        std::cout << __PRETTY_FUNCTION__ << std::endl;

        // Setting the vertical coordinate
        // GRD_gz    = new double[ADM_kall];
        // GRD_gzh   = new double[ADM_kall];
        // GRD_dgz   = new double[ADM_kall];
        // GRD_dgzh  = new double[ADM_kall];
        // GRD_rdgz  = new double[ADM_kall];
        // GRD_rdgzh = new double[ADM_kall];
        // Vertical interpolation factor
        // GRD_afact = new double[ADM_kall];
        // GRD_bfact = new double[ADM_kall];
        // GRD_cfact = new double[ADM_kall];
        // GRD_dfact = new double[ADM_kall];

        // Reading data from vgrid94.dat
        // GRD_Input_vgrid(vgrid_fname);
        GRD_Input_vgrid();

        for(int k = 0; k <= ADM_kmax; k++)
        {
            GRD_dgz[k] = GRD_gzh[k+1] - GRD_gzh[k];
        }
        GRD_dgz[ADM_kmax + 1] = GRD_dgz[ADM_kmax];

        for(int k = 0; k <= ADM_kmax; k++)
        {
            GRD_dgzh[k] = GRD_gz[k+1] - GRD_gz[k];
        }
        GRD_dgz[ADM_kmin - 1] = GRD_dgz[ADM_kmin];

        for(int k = 0; k < ADM_kall; k++)
        {
            GRD_rdgz[k] = 1.0 / GRD_dgz[k];
            GRD_rdgzh[k] = 1.0 / GRD_dgzh[k];
        }

        for(int k = ADM_kmin; k <= ADM_kmax + 1; k++)
        {
            GRD_afact[k] = ( GRD_gzh[k] - GRD_gz[k-1] ) /
                            ( GRD_gz[k] - GRD_gz[k-1] );
        }
        GRD_cfact[ADM_kmin - 1] = 1.0;
        GRD_cfact[ADM_kmax + 1] = 0.0;

        for(int k = 0; k < ADM_kall; k++)
        {
            GRD_dfact[k] = 1.0 - GRD_cfact[k];
        }
    }

    void GRD_Input_vgrid()
    {
        /** Waiting for implementing */
        std::cout << __PRETTY_FUNCTION__ << " Reading vgrid data " << std::endl;
        read_data_1d("data/vgrid/GRD_gz.dat", GRD_dgz);
        read_data_1d("data/vgrid/GRD_gzh.dat", GRD_dgzh);
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
        /** TO DO */
        std::cout << __PRETTY_FUNCTION__ << std::endl;

        size_t gall = ijdim;

        // --- horizontal kinetic energy
        for(int k = kmin; k <= kmax; k++)
        {
            for(int g = 0; g < gall; g++)
            {
                rhogkin[k][g] = 0.5 * ( rhogvx[k][g] * rhogvx[k][g] +
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
        for(int k = kmin; k <= kmax; k++)
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
}
