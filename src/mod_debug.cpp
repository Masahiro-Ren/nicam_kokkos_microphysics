#include "mod_debug.h"

using namespace PROBLEM_SIZE;
using namespace DATA_IO;

namespace DEBUG {
    void ADM_Proc_stop()
    {
        exit(1);
    }

    void PROF_val_check()
    {
        /** Waiting for implementing */
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    void GRD_Setup()
    {
        /** Waiting for implementing */
        // std::cout << __PRETTY_FUNCTION__ << std::endl;

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
        // std::cout << __PRETTY_FUNCTION__ << std::endl;
        read_data_1d("data/vgrid/GRD_gz.dat", GRD_dgz);
        read_data_1d("data/vgrid/GRD_gzh.dat", GRD_dgzh);
    }

    template<size_t ijdim, size_t kdim>
    void cnvvar_rhogkin_in(
        double rhog     [kdim][ijdim],
        double rhogvx   [kdim][ijdim],
        double rhogvy   [kdim][ijdim],
        double rhogvz   [kdim][ijdim],
        double rhogw    [kdim][ijdim],
        double C2Wfact  [2][kdim][ijdim],
        double W2Cfaxt  [2][kdim][ijdim],
        double rhogkin  [kdim][ijdim],
        double rhogkin_h[kdim][ijdim],
        double rhogkin_v[kdim][ijdim] )
    {
        /** Waiting for implementing */
        std::cout << __PRETTY_FUNCTION__ << std::endl;
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
};