#include "mod_debug.h"

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
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    void GRD_Input_vgrid()
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