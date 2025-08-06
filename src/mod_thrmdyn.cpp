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
    void THRMDYN_qd(const View3D<double, DEFAULT_MEM>& q, const View2D<double, DEFAULT_MEM>& qd)
    {
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            qd(k,ij) = 1.0;
        });

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
            KOKKOS_LAMBDA(const size_t k, const size_t ij){
                qd(k,ij) = qd(k,ij) - q(nq,k,ij);
            });
        }

    }

    void THRMDYN_cv(const View2D<double, DEFAULT_MEM>& qd, const View3D<double, DEFAULT_MEM>& q, const View2D<double, DEFAULT_MEM>& cv)
    {
        // double CVdry = CONST_CVdry;
        auto CVW = PROBLEM_SIZE::CVW;
        // View<double> CVdry("CVdry");
        // Kokkos::deep_copy(CVdry, CONST_CVdry);

        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            const double CVdry = CONST_CVdry;
            cv(k,ij) = qd(k,ij) * CVdry;
        });

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
            KOKKOS_LAMBDA(const size_t k, const size_t ij){
                cv(k,ij) = cv(k,ij) + q(nq,k,ij) * CVW(nq);
            });
        }
    }

    void THRMDYN_tempre( const View2D<double, DEFAULT_MEM>&  ein, 
                         const View2D<double, DEFAULT_MEM>&  rho,
                         const View3D<double, DEFAULT_MEM>&  q  ,
                         const View2D<double, DEFAULT_MEM>&  tem,
                         const View2D<double, DEFAULT_MEM>&  pre )
    {
        auto CVW = PROBLEM_SIZE::CVW;
        View2D<double, DEFAULT_MEM> cv("cv", kdim, ijdim);
        View2D<double, DEFAULT_MEM> qd("qd", kdim, ijdim);
        // View<double> CVdry("CVdry"); 
        // View<double> Rdry ("Rdry "); 
        // View<double> Rvap ("Rvap "); 

        // Kokkos::deep_copy(CVdry, CONST_CVdry);
        // Kokkos::deep_copy(Rdry, CONST_Rdry);
        // Kokkos::deep_copy(Rvap, CONST_Rvap);
        // Kokkos::deep_copy(d_CVW, CVW);

        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            cv(k,ij) = 0.0;
            qd(k,ij) = 1.0;
        });

        for(int nq = NQW_STR; nq <= NQW_END; nq++)
        {
            Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
            KOKKOS_LAMBDA(const size_t k, const size_t ij){
                cv(k,ij) = cv(k,ij) + q(nq,k,ij) * CVW(nq);
                qd(k,ij) = qd(k,ij) - q(nq,k,ij);
            });
        }

        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
                const double CVdry = CONST_CVdry;
                const double Rdry = CONST_Rdry;
                const double Rvap = CONST_Rvap;
                cv (k,ij) = cv (k,ij) + qd (k,ij) * CVdry;
                tem(k,ij) = ein(k,ij) / cv (k,ij);
                pre(k,ij) = rho(k,ij) * tem(k,ij) * ( qd(k,ij) * Rdry + q(I_QV,k,ij) * Rvap );
        });

    }
};
