#include "mod_thrmdyn.h"

namespace THRMDYN{
    void THRMDYN_qd(double q[nqmax][kdim][ijdim], double qd[kdim][ijdim])
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                qd[k][ij] = 1.0;
            }
        }

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            for(int k = 0; k < kdim; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    qd[k][ij] = qd[k][ij] - q[nq][k][ij];
                }
            }
        }
    }

    void THRMDYN_cv(double qd[kdim][ijdim], double q[nqmax][kdim][ijdim], double cv[kdim][ijdim])
    {
        double CVdry = CONST_CVdry;

        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                cv[k][ij] = qd[k][ij] * CVdry;
            }
        }

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            for(int k = 0; k < kdim; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    cv[k][ij] = cv[k][ij] + q[nq][k][ij] * CVW[nq];
                }
            }
        }

    }

    void THRMDYN_tempre( double ein[kdim][ijdim], 
                         double rho[kdim][ijdim],
                         double q[nqmax][kdim][ijdim],
                         double tem[kdim][ijdim],
                         double pre[kdim][ijdim] )
    {
        double cv[kdim][ijdim];
        double qd[kdim][ijdim];
        double CVdry = CONST_CVdry;
        double Rdry = CONST_Rdry;
        double Rvap = CONST_Rvap;

        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                cv[k][ij] = 0.0;
                qd[k][ij] = 1.0;
            }
        }

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            for(int k = 0; k < kdim; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    cv[k][ij] = cv[k][ij] + q[nq][k][ij] * CVW[nq];
                    qd[k][ij] = qd[k][ij] - q[nq][k][ij];
                }
            }
        }

        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                cv[k][ij] = cv[k][ij] + qd[k][ij] * CVdry;
                tem[k][ij] = ein[k][ij] / cv[k][ij];
                pre[k][ij] = rho[k][ij] * tem[k][ij] * ( qd[k][ij] * Rdry + q[I_QV][k][ij] * Rvap );
            }
        }
    }
};
