#include "mod_debug.h"

using namespace PROBLEM_SIZE;

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
        GRD_Input_vgrid(vgrid_fname);

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

    void GRD_Input_vgrid(const std::string& fname)
    {
        /** Waiting for implementing */
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    void cnvvar_rhogkin_in()
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