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

    /**
     * Kokkos ver.
     */
    void THRMDYN_qd(const View<double***>& q, const View<double**>& qd)
    {
        // for(int k = 0; k < kdim; k++)
        // {
        //     for(int ij = 0; ij < ijdim; ij++)
        //     {
        //         qd(k,ij) = 1.0;
        //     }
        // }
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            qd(k,ij) = 1.0;
        });

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            // for(int k = 0; k < kdim; k++)
            // {
            //     for(int ij = 0; ij < ijdim; ij++)
            //     {
            //         qd(k,ij) = qd(k,ij) - q(nq,k,ij);
            //     }
            // }
            Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
            KOKKOS_LAMBDA(const size_t k, const size_t ij){
                qd(k,ij) = qd(k,ij) - q(nq,k,ij);
            });
        }

    }

    void THRMDYN_cv(const View<double**>& qd, const View<double***>& q, const View<double**>& cv)
    {
        double CVdry = CONST_CVdry;

        // for(int k = 0; k < kdim; k++)
        // {
        //     for(int ij = 0; ij < ijdim; ij++)
        //     {
        //         cv(k,ij) = qd(k,ij) * CVdry;
        //     }
        // }
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            cv(k,ij) = qd(k,ij) * CVdry;
        });

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            // for(int k = 0; k < kdim; k++)
            // {
            //     for(int ij = 0; ij < ijdim; ij++)
            //     {
            //         cv(k,ij) = cv(k,ij) + q(nq,k,ij) * CVW[nq];
            //     }
            // }
            Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
            KOKKOS_LAMBDA(const size_t k, const size_t ij){
                cv(k,ij) = cv(k,ij) + q(nq,k,ij) * CVW[nq];
            });
        }
    }

    void THRMDYN_tempre( const View<double**>&  ein, 
                         const View<double**>&  rho,
                         const View<double***>& q  ,
                         const View<double**>&  tem,
                         const View<double**>&  pre )
    {
        View<double**> cv("cv", kdim, ijdim);
        View<double**> qd("qd", kdim, ijdim);
        double CVdry = CONST_CVdry;
        double Rdry = CONST_Rdry;
        double Rvap = CONST_Rvap;

        // for(int k = 0; k < kdim; k++)
        // {
        //     for(int ij = 0; ij < ijdim; ij++)
        //     {
        //         cv(k,ij) = 0.0;
        //         qd(k,ij) = 1.0;
        //     }
        // }
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            cv(k,ij) = 0.0;
            qd(k,ij) = 1.0;
        });

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            // for(int k = 0; k < kdim; k++)
            // {
            //     for(int ij = 0; ij < ijdim; ij++)
            //     {
            //         cv(k,ij) = cv(k,ij) + q(nq,k,ij) * CVW[nq];
            //         qd(k,ij) = qd(k,ij) - q(nq,k,ij);
            //     }
            // }
            Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
            KOKKOS_LAMBDA(const size_t k, const size_t ij){
                cv(k,ij) = cv(k,ij) + q(nq,k,ij) * CVW[nq];
                qd(k,ij) = qd(k,ij) - q(nq,k,ij);
            });
        }

        // for(int k = 0; k < kdim; k++)
        // {
        //     for(int ij = 0; ij < ijdim; ij++)
        //     {
        //         cv (k,ij) = cv (k,ij) + qd (k,ij) * CVdry;
        //         tem(k,ij) = ein(k,ij) / cv (k,ij);
        //         pre(k,ij) = rho(k,ij) * tem(k,ij) * ( qd(k,ij) * Rdry + q(I_QV,k,ij) * Rvap );

        //     }
        // }
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
                cv (k,ij) = cv (k,ij) + qd (k,ij) * CVdry;
                tem(k,ij) = ein(k,ij) / cv (k,ij);
                pre(k,ij) = rho(k,ij) * tem(k,ij) * ( qd(k,ij) * Rdry + q(I_QV,k,ij) * Rvap );
        });

    }
};
