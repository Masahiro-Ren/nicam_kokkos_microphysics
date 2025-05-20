#include "mod_satadjust.h"

/**
 * Some private parameters
 */
double TEM_MIN = 10.0;
double DTEM_EPS0 = 1.0E-8;

std::string ALPHA_TYPE = "LINEAR";

double SATURATION_ULIMIT_TEMP = 273.15;
double SATURATION_LLIMIT_TEMP = 233.15;



namespace SATADJUST{

double CPovR_liq;
double CPovR_ice;
double CVovR_liq;
double CVovR_ice;
double LovR_liq;
double LovR_ice;


void SATURATION_Setup()
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    if(EIN_TYPE == "EXACT")
    {
        CPovR_liq = (CONST_CPvap - CONST_CL) / CONST_Rvap;
        CPovR_ice = (CONST_CPvap - CONST_CI) / CONST_Rvap;
        CVovR_liq = (CONST_CVvap - CONST_CL) / CONST_Rvap;
        CVovR_ice = (CONST_CVvap - CONST_CI) / CONST_Rvap;
    }
    else if (EIN_TYPE == "SIMPLE" || EIN_TYPE == "SIMPLE2")
    {
        CPovR_liq = 0.0;
        CPovR_ice = 0.0;
        CVovR_liq = 0.0;
        CVovR_ice = 0.0;
    }

    LovR_liq = LHV / CONST_Rvap;
    LovR_ice = LHS / CONST_Rvap;
}

void SATURATION_Setrange(double Tw, double Ti)
{
    SATURATION_ULIMIT_TEMP = Tw;
    SATURATION_LLIMIT_TEMP = Ti;
}

void SATURATION_psat_liq(double tem[kdim][ijdim], double psat[kdim][ijdim])
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double rtem;

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0 = CONST_PSAT0;

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

            psat[k][ij] = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_liq) * 
                                  std::exp(LovR_liq * (RTEM00 - rtem));
        }
    }
}

// Kokkos ver.
void SATURATION_psat_liq(const View<double**>& tem, const View<double**>& psat)
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0 = CONST_PSAT0;

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
            double rtem = 1.0 / ( std::max(tem(k,ij), TEM_MIN) );

            psat(k,ij) = PSAT0 * std::pow((tem(k,ij) * RTEM00), CPovR_liq) * 
                                 std::exp(LovR_liq * (RTEM00 - rtem));
    });
}

void SATURATION_psat_ice(double tem[kdim][ijdim], double psat[kdim][ijdim])
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double rtem;

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0 = CONST_PSAT0;

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

            psat[k][ij] = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_ice) * 
                                  std::exp(LovR_ice * (RTEM00 - rtem));
        }
    }
}

//Kokkos ver.
void SATURATION_psat_ice(const View<double**>& tem, const View<double**>& psat)
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0 = CONST_PSAT0;

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
            double rtem = 1.0 / ( std::max(tem(k,ij), TEM_MIN) );

            psat(k,ij) = PSAT0 * std::pow((tem(k,ij) * RTEM00), CPovR_ice) * 
                                 std::exp(LovR_ice * (RTEM00 - rtem));
    });
}

void SATURATION_adjustment( double rhog   [kdim][ijdim],
                            double rhoge  [kdim][ijdim],
                            double rhogq  [nqmax][kdim][ijdim],
                            double tem    [kdim][ijdim],
                            double q      [nqmax][kdim][ijdim],
                            double qd     [kdim][ijdim],
                            double gsgam2 [kdim][ijdim],
                            bool   ice_adjust )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double ein_mosit[kdim][ijdim];
    double qsum     [kdim][ijdim];
    double CVtot    [kdim][ijdim];
    double rho      [kdim][ijdim];

    // ein_moist = U1(rho,qsum,T1) : "unsaturated temperature"
    if(I_QI > 0 && ice_adjust)
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                ein_mosit[k][ij] = rhoge[k][ij] / rhog[k][ij] 
                                    + q[I_QV][k][ij] * LHV 
                                    - q[I_QI][k][ij] * LHF;
                qsum[k][ij] = q[I_QV][k][ij] 
                              + q[I_QC][k][ij]
                              + q[I_QI][k][ij];

                q[I_QV][k][ij] = qsum[k][ij];
                q[I_QC][k][ij] = 0.0;
                q[I_QI][k][ij] = 0.0;
            }
        }
    }
    else
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                ein_mosit[k][ij] = rhoge[k][ij] / rhog[k][ij] 
                                    + q[I_QV][k][ij] * LHV;

                qsum[k][ij] = q[I_QV][k][ij] 
                              + q[I_QC][k][ij];

                q[I_QV][k][ij] = qsum[k][ij];
                q[I_QC][k][ij] = 0.0;
            }
        }
    }

    THRMDYN_cv(qd, q, CVtot);

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rho[k][ij] = rhog[k][ij] /  gsgam2[k][ij];
            tem[k][ij] = ( ein_mosit[k][ij] - q[I_QV][k][ij] * LHV ) / CVtot[k][ij];
        }
    }

    if(I_QI > 0 && ice_adjust)
    {
        satadjust_all(rho, ein_mosit, qsum, tem, q);
    }
    else
    {
        satadjust_liq(rho, ein_mosit, qsum, tem, q);
    }


    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rhogq[I_QV][k][ij] = rhog[k][ij] * q[I_QV][k][ij];
            rhogq[I_QC][k][ij] = rhog[k][ij] * q[I_QC][k][ij];
        }
    }

    if(I_QI > 0 && ice_adjust)
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhogq[I_QI][k][ij] = rhog[k][ij] * q[I_QI][k][ij];
            }
        }
    }

    THRMDYN_cv(qd, q, CVtot);

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rhoge[k][ij] = rhog[k][ij] * tem[k][ij] * CVtot[k][ij];
        }
    }

}

/**
 * Kokkos ver.
 */
void SATURATION_adjustment( const View<double**>&  rhog   ,
                            const View<double**>&  rhoge  ,
                            const View<double***>& rhogq  ,
                            const View<double**>&  tem    ,
                            const View<double***>& q      ,
                            const View<double**>&  qd     ,
                            const View<double**>&  gsgam2 ,
                            bool   ice_adjust )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    View<double**> ein_mosit("ein_mosit", kdim, ijdim);
    View<double**> qsum     ("qsum     ", kdim, ijdim);
    View<double**> CVtot    ("CVtot    ", kdim, ijdim);
    View<double**> rho      ("rho      ", kdim, ijdim);

    // ein_moist = U1(rho,qsum,T1) : "unsaturated temperature"
    if(I_QI > 0 && ice_adjust)
    {
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            ein_mosit(k,ij) = rhoge(k,ij) / rhog(k,ij) 
                              + q(I_QV,k,ij) * LHV 
                              - q(I_QI,k,ij) * LHF;
            qsum(k,ij) =  q(I_QV,k,ij) 
                        + q(I_QC,k,ij)
                        + q(I_QI,k,ij);

            q(I_QV,k,ij) = qsum(k,ij);
            q(I_QC,k,ij) = 0.0;
            q(I_QI,k,ij) = 0.0;  
        });
    }
    else
    {

        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            ein_mosit(k,ij) = rhoge(k,ij) / rhog(k,ij) 
                                + q(I_QV,k,ij) * LHV;

            qsum(k,ij) =  q(I_QV,k,ij) 
                        + q(I_QC,k,ij);

            q(I_QV,k,ij) = qsum(k,ij);
            q(I_QC,k,ij) = 0.0;
        });
    }

    THRMDYN_cv(qd, q, CVtot);

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        rho(k,ij) = rhog(k,ij) /  gsgam2(k,ij);
        tem(k,ij) = ( ein_mosit(k,ij) - q(I_QV,k,ij) * LHV ) / CVtot(k,ij);
    });

    if(I_QI > 0 && ice_adjust)
    {
        satadjust_all(rho, ein_mosit, qsum, tem, q);
    }
    else
    {
        satadjust_liq(rho, ein_mosit, qsum, tem, q);
    }

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        rhogq(I_QV,k,ij) = rhog(k,ij) * q(I_QV,k,ij);
        rhogq(I_QC,k,ij) = rhog(k,ij) * q(I_QC,k,ij); 
    });

    if(I_QI > 0 && ice_adjust)
    {
        Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
        KOKKOS_LAMBDA(const size_t k, const size_t ij){
            rhogq(I_QI,k,ij) = rhog(k,ij) * q(I_QI,k,ij);
        });
    }

    THRMDYN_cv(qd, q, CVtot);

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({0,0},{kdim,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        rhoge(k,ij) = rhog(k,ij) * tem(k,ij) * CVtot(k,ij);
    });
}

void satadjust_all( double rho    [kdim][ijdim],
                    double Emoist [kdim][ijdim],
                    double qsum   [kdim][ijdim],
                    double tem    [kdim][ijdim],
                    double q      [nqmax][kdim][ijdim] )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double qd[kdim][ijdim];

    double rtem, alpha, psatl, psati, psat, qsatl, qsati,qsat;
    double CVtot, Emoist_new, dtemp, lim1, lim2;
    double dalpha_dT, dqsatl_dT, dqsati_dT, dqsat_dT, dqc_dT, dqi_dT, dCVtot_dT, dEmoist_dT;

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0  = CONST_PSAT0;
    double Rvap   = CONST_Rvap;
    double CVdry  = CONST_CVdry;
    double EPS    = CONST_EPS;

    double dtemp_criteria = std::pow(10.0, (-(RP_PREC + 1)/2));

    const int itelim = 100;
    bool converged;

    THRMDYN_qd(q, qd);

    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            alpha = ( tem[k][ij] - SATURATION_LLIMIT_TEMP ) / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP );
            alpha = std::min( std::max(alpha, 0.0), 1.0 );

            rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

            psatl = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_liq) * std::exp(LovR_liq * (RTEM00 - rtem));
            psati = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_ice) * std::exp(LovR_ice * (RTEM00 - rtem));
            psat  = psatl * alpha + psati * (1.0 - alpha);

            qsat = psat / ( rho[k][ij] * Rvap * tem[k][ij] );

            if(qsum[k][ij] - qsat > EPS)
            {
                converged = false;
                int ite;
                for(ite = 0; ite < itelim; ite++)
                {
                    alpha = ( tem[k][ij] - SATURATION_LLIMIT_TEMP ) / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP );
                    alpha = std::min( std::max(alpha, 0.0), 1.0 );

                    rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

                    psatl = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_liq) * std::exp(LovR_liq * (RTEM00 - rtem));
                    psati = PSAT0 * std::pow((tem[k][ij] * RTEM00), CPovR_ice) * std::exp(LovR_ice * (RTEM00 - rtem));
                    psat  = psatl * alpha + psati * (1.0 - alpha);

                    qsatl = psatl / ( rho[k][ij] * Rvap * tem[k][ij] );
                    qsati = psati / ( rho[k][ij] * Rvap * tem[k][ij] );
                    qsat  = psat  / ( rho[k][ij] * Rvap * tem[k][ij] );

                    // Sepration
                    q[I_QV][k][ij] = qsat;
                    q[I_QC][k][ij] = ( qsum[k][ij] - qsat ) * alpha;
                    q[I_QI][k][ij] = ( qsum[k][ij] - qsat ) * (1.0 - alpha);

                    CVtot = qd[k][ij] * CVdry;
                    for(int nq = NQW_STR; nq <= NQW_END; nq++)
                        CVtot = CVtot + q[nq][k][ij] * CVW[nq];
                    
                    Emoist_new = tem[k][ij] * CVtot + q[I_QV][k][ij] * LHV - q[I_QI][k][ij] * LHF;

                    // dx/dT
                    lim1 = 0.5 + std::copysign(0.5, SATURATION_ULIMIT_TEMP - tem[k][ij]);
                    lim2 = 0.5 + std::copysign(0.5, tem[k][ij] - SATURATION_LLIMIT_TEMP);
                    dalpha_dT = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP );

                    dqsatl_dT = ( LovR_liq / (std::pow(tem[k][ij], 2)) + CVovR_liq / tem[k][ij] ) * qsatl;
                    dqsati_dT = ( LovR_ice / (std::pow(tem[k][ij], 2)) + CVovR_ice / tem[k][ij] ) * qsati;
                    dqsat_dT  = qsatl * dalpha_dT + dqsatl_dT * alpha -
                                qsati * dalpha_dT + dqsati_dT * (1.0 - alpha);
                    
                    dqc_dT =  ( qsum[k][ij] - qsat ) * dalpha_dT - dqsat_dT * alpha;
                    dqi_dT = -( qsum[k][ij] - qsat ) * dalpha_dT - dqsat_dT * (1.0 - alpha);

                    dCVtot_dT = dqsat_dT * CVW[I_QV]
                                + dqc_dT * CVW[I_QC]
                                + dqi_dT * CVW[I_QI];

                    dEmoist_dT = tem[k][ij] * dCVtot_dT
                                 + CVtot
                                 + dqsat_dT * LHV
                                 - dqi_dT   * LHF;
                    
                    dtemp = ( Emoist_new - Emoist[k][ij] ) / dEmoist_dT;

                    tem[k][ij] = tem[k][ij] - dtemp;

                    if(std::abs(dtemp) < dtemp_criteria)
                    {
                        converged = true;
                        break;
                    }

                    if( tem[k][ij] * 0.0 != 0.0 ) break;

                }

                if(!converged)
                {
                    std::cerr << rho[k][ij] << "\t" << tem[k][ij] << "\t" << q[I_QV][k][ij] << "\t" << q[I_QC][k][ij] << "\t" << q[I_QI][k][ij] << std::endl;
                    std::cerr << "xxx [satadjust_all] not converged! dtemp = " << dtemp << " ij= " << ij << " k= " << k << " ite= " << ite << std::endl;
                    ADM_Proc_stop();
                }
            }
        }
    }

}

void satadjust_liq( double rho    [kdim][ijdim],
                    double Emoist [kdim][ijdim],
                    double qsum   [kdim][ijdim],
                    double tem    [kdim][ijdim],
                    double q      [nqmax][kdim][ijdim] )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double qd[kdim][ijdim];

    double rtem, psat, qsat;
    double CVtot, Emosit_new, dtemp;
    double dqsat_dT, dCVtot_dT, dEmosit_dT;
    
    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0  = CONST_PSAT0;
    double Rvap   = CONST_Rvap;
    double CVdry  = CONST_CVdry;
    double EPS    = CONST_EPS;

    double dtemp_criteria = std::pow(10.0, (-(RP_PREC + 1)/2));

    const int itelim = 100;
    bool converged;

    // In Fortran: ite_temp(2:95, 1:16641)
    int ite_temp[kdim][ijdim] = {{0}};

    THRMDYN_qd(q, qd);

    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

            psat = PSAT0 * ( std::pow(tem[k][ij] * RTEM00, CPovR_liq) ) * std::exp(LovR_liq * (RTEM00 - rtem));

            qsat = psat / ( rho[k][ij] * Rvap * tem[k][ij] );

            if( qsum[k][ij] - qsat > EPS )
            {
                converged = false;

                int ite;

                for(ite = 0; ite < itelim; ite++)
                {
                    rtem = 1.0 / ( std::max(tem[k][ij], TEM_MIN) );

                    psat = PSAT0 * ( std::pow(tem[k][ij] * RTEM00, CPovR_liq) ) * std::exp(LovR_liq * (RTEM00 - rtem));

                    qsat = psat / ( rho[k][ij] * Rvap * tem[k][ij] );

                    // Sepration
                    q[I_QV][k][ij] = qsat;
                    q[I_QC][k][ij] = qsum[k][ij] - qsat;

                    CVtot = qd[k][ij] * CVdry;
                    for(int nq = NQW_STR; nq <= NQW_END; nq++)
                        CVtot = CVtot + q[nq][k][ij] * CVW[nq];
                    
                    Emosit_new = tem[k][ij] * CVtot + q[I_QV][k][ij] * LHV;

                    // dx/dT
                    dqsat_dT = ( LovR_liq / (std::pow(tem[k][ij], 2)) + CVovR_liq / tem[k][ij] ) * qsat;

                    dCVtot_dT = dqsat_dT * ( CVW[I_QV] - CVW[I_QC] );

                    dEmosit_dT = tem[k][ij] * dCVtot_dT
                                 + CVtot
                                 + dqsat_dT * LHV;
                    
                    dtemp = ( Emosit_new - Emoist[k][ij] ) / dEmosit_dT;

                    tem[k][ij] = tem[k][ij] - dtemp;

                    if(std::abs(dtemp) < dtemp_criteria)
                    {
                        converged = true;
                        break;
                    }

                    if(tem[k][ij] * 0.0 != 0.0) break;
                }

                ite_temp[k][ij] = ite;

                if(!converged)
                {
                    std::cerr << rho[k][ij] << "\t" << tem[k][ij] << "\t" << q[I_QV][k][ij] << "\t" << q[I_QC][k][ij] << "\t" << q[I_QI][k][ij] << std::endl;
                    std::cerr << "xxx [satadjust_all] not converged! dtemp = " << dtemp << " ij= " << ij << " k= " << k << " ite= " << ite << std::endl;
                    ADM_Proc_stop();
                }
            }

        }
    }

    // find sum, max in one time
    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            sat_ite_sum += ite_temp[k][ij];
            if(sat_ite_max < ite_temp[k][ij]) sat_ite_max = ite_temp[k][ij];
        }
    }

    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            if(ite_temp[k][ij] != 0)
                sat_ite_count += 1;
            else
                ite_temp[k][ij] = itelim;
        }
    }

    // find min value
    for(int k = kmin; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            if(sat_ite_min > ite_temp[k][ij] || sat_ite_min < 0) sat_ite_min = ite_temp[k][ij];
        }
    }

}

// Kokkos private procedures
void satadjust_all( const View<double**>&  rho    ,
                    const View<double**>&  Emoist ,
                    const View<double**>&  qsum   ,
                    const View<double**>&  tem    ,
                    const View<double***>& q       )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    View<double**> qd("qd", kdim, ijdim);

    // double rtem, alpha, psatl, psati, psat, qsatl, qsati, qsat;
    // double CVtot, Emoist_new, dtemp, lim1, lim2;
    // double dalpha_dT, dqsatl_dT, dqsati_dT, dqsat_dT, dqc_dT, dqi_dT, dCVtot_dT, dEmoist_dT;

    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0  = CONST_PSAT0;
    double Rvap   = CONST_Rvap;
    double CVdry  = CONST_CVdry;
    double EPS    = CONST_EPS;

    double dtemp_criteria = std::pow(10.0, (-(RP_PREC + 1)/2));

    const int itelim = 100;
    // bool converged;

    THRMDYN_qd(q, qd);

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k,  const size_t ij){
        double alpha = ( tem(k,ij) - SATURATION_LLIMIT_TEMP ) / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP );
        alpha = std::min( std::max(alpha, 0.0), 1.0 );

        double rtem = 1.0 / ( std::max(tem(k,ij), TEM_MIN) );

        double psatl = PSAT0 * std::pow((tem(k,ij) * RTEM00), CPovR_liq) * std::exp(LovR_liq * (RTEM00 - rtem));
        double psati = PSAT0 * std::pow((tem(k,ij) * RTEM00), CPovR_ice) * std::exp(LovR_ice * (RTEM00 - rtem));
        double psat  = psatl * alpha + psati * (1.0 - alpha);

        double qsat = psat / ( rho(k,ij) * Rvap * tem(k,ij) );

        if(qsum(k,ij) - qsat > EPS)
        {
            bool converged = false;
            int ite;
            double dtemp;
            for(ite = 0; ite < itelim; ite++)
            {
                alpha = ( tem(k,ij) - SATURATION_LLIMIT_TEMP ) / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP );
                alpha = std::min( std::max(alpha, 0.0), 1.0 );

                rtem = 1.0 / ( std::max(tem(k,ij), TEM_MIN) );

                psatl = PSAT0 * std::pow((tem(k,ij) * RTEM00), CPovR_liq) * std::exp(LovR_liq * (RTEM00 - rtem));
                psati = PSAT0 * std::pow((tem(k,ij) * RTEM00), CPovR_ice) * std::exp(LovR_ice * (RTEM00 - rtem));
                psat  = psatl * alpha + psati * (1.0 - alpha);

                double qsatl = psatl / ( rho(k,ij) * Rvap * tem(k,ij) );
                double qsati = psati / ( rho(k,ij) * Rvap * tem(k,ij) );
                qsat  = psat  / ( rho(k,ij) * Rvap * tem(k,ij) );

                // Sepration
                q(I_QV,k,ij) = qsat;
                q(I_QC,k,ij) = ( qsum(k,ij) - qsat ) * alpha;
                q(I_QI,k,ij) = ( qsum(k,ij) - qsat ) * (1.0 - alpha);

                double CVtot = qd(k,ij) * CVdry;
                for(int nq = NQW_STR; nq <= NQW_END; nq++)
                    CVtot = CVtot + q(nq,k,ij) * CVW[nq];
                
                double Emoist_new = tem(k,ij) * CVtot + q(I_QV,k,ij) * LHV - q(I_QI,k,ij) * LHF;

                // dx/dT
                double lim1 = 0.5 + std::copysign(0.5, SATURATION_ULIMIT_TEMP - tem(k,ij));
                double lim2 = 0.5 + std::copysign(0.5, tem(k,ij) - SATURATION_LLIMIT_TEMP);
                double dalpha_dT = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP );

                double dqsatl_dT = ( LovR_liq / (std::pow(tem(k,ij), 2)) + CVovR_liq / tem(k,ij) ) * qsatl;
                double dqsati_dT = ( LovR_ice / (std::pow(tem(k,ij), 2)) + CVovR_ice / tem(k,ij) ) * qsati;
                double dqsat_dT  = qsatl * dalpha_dT + dqsatl_dT * alpha -
                        qsati * dalpha_dT + dqsati_dT * (1.0 - alpha);
                
                double dqc_dT =  ( qsum(k,ij) - qsat ) * dalpha_dT - dqsat_dT * alpha;
                double dqi_dT = -( qsum(k,ij) - qsat ) * dalpha_dT - dqsat_dT * (1.0 - alpha);

                double dCVtot_dT = dqsat_dT * CVW[I_QV]
                                   + dqc_dT * CVW[I_QC]
                                   + dqi_dT * CVW[I_QI];

                double dEmoist_dT = tem(k,ij) * dCVtot_dT
                                    + CVtot
                                    + dqsat_dT * LHV
                                    - dqi_dT   * LHF;
                
                dtemp = ( Emoist_new - Emoist(k,ij) ) / dEmoist_dT;

                tem(k,ij) = tem(k,ij) - dtemp;

                if(std::abs(dtemp) < dtemp_criteria)
                {
                    converged = true;
                    break;
                }

                if( tem(k,ij) * 0.0 != 0.0 ) break;
            }

            if(!converged)
            {
                std::cerr << rho(k,ij) << "\t" << tem(k,ij) << "\t" << q(I_QV,k,ij) << "\t" << q(I_QC,k,ij) << "\t" << q(I_QI,k,ij) << std::endl;
                std::cerr << "xxx [satadjust_all] not converged! dtemp = " << dtemp << " ij= " << ij << " k= " << k << " ite= " << ite << std::endl;
                ADM_Proc_stop();
            }
        }

    });

}

void satadjust_liq( const View<double**>&  rho    ,
                    const View<double**>&  Emoist ,
                    const View<double**>&  qsum   ,
                    const View<double**>&  tem    ,
                    const View<double***>& q       )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    View<double**> qd("qd", kdim, ijdim);

    // double rtem, psat, qsat;
    // double CVtot, Emosit_new, dtemp;
    // double dqsat_dT, dCVtot_dT, dEmosit_dT;
    
    double RTEM00 = 1.0 / CONST_TEM00;
    double PSAT0  = CONST_PSAT0;
    double Rvap   = CONST_Rvap;
    double CVdry  = CONST_CVdry;
    double EPS    = CONST_EPS;

    double dtemp_criteria = std::pow(10.0, (-(RP_PREC + 1)/2));

    const int itelim = 100;
    // bool converged;

    // In Fortran: ite_temp(2:95, 1:16641)
    View<int**> ite_temp("ite_temp", kdim, ijdim);

    THRMDYN_qd(q, qd);

    Kokkos::parallel_for(MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO},{kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij){
        double rtem = 1.0 / ( std::max(tem(k,ij), TEM_MIN) );

        double psat = PSAT0 * ( std::pow(tem(k,ij) * RTEM00, CPovR_liq) ) * std::exp(LovR_liq * (RTEM00 - rtem));

        double qsat = psat / ( rho(k,ij) * Rvap * tem(k,ij) ); 

        if( qsum(k,ij) - qsat > EPS )
        {
            bool converged = false;
            int ite;
            double dtemp;

            for(ite = 0; ite < itelim; ite++)
            {
                rtem = 1.0 / ( std::max(tem(k,ij), TEM_MIN) );

                psat = PSAT0 * ( std::pow(tem(k,ij) * RTEM00, CPovR_liq) ) * std::exp(LovR_liq * (RTEM00 - rtem));

                qsat = psat / ( rho(k,ij) * Rvap * tem(k,ij) );

                // Sepration
                q(I_QV,k,ij) = qsat;
                q(I_QC,k,ij) = qsum(k,ij) - qsat;

                double CVtot = qd(k,ij) * CVdry;
                for(int nq = NQW_STR; nq <= NQW_END; nq++)
                    CVtot = CVtot + q(nq,k,ij) * CVW[nq];
                
                double Emosit_new = tem(k,ij) * CVtot + q(I_QV,k,ij) * LHV;

                // dx/dT
                double dqsat_dT = ( LovR_liq / (std::pow(tem(k,ij), 2)) + CVovR_liq / tem(k,ij) ) * qsat;

                double dCVtot_dT = dqsat_dT * ( CVW[I_QV] - CVW[I_QC] );

                double dEmosit_dT = tem(k,ij) * dCVtot_dT
                                + CVtot
                                + dqsat_dT * LHV;
                
                dtemp = ( Emosit_new - Emoist(k,ij) ) / dEmosit_dT;

                tem(k,ij) = tem(k,ij) - dtemp;

                if(std::abs(dtemp) < dtemp_criteria)
                {
                    converged = true;
                    break;
                }

                if(tem(k,ij) * 0.0 != 0.0) break;
            }

            ite_temp(k,ij) = ite;

            if(!converged)
            {
                std::cerr << rho(k,ij) << "\t" << tem(k,ij) << "\t" << q(I_QV,k,ij) << "\t" << q(I_QC,k,ij) << "\t" << q(I_QI,k,ij) << std::endl;
                std::cerr << "xxx [satadjust_all] not converged! dtemp = " << dtemp << " ij= " << ij << " k= " << k << " ite= " << ite << std::endl;
                ADM_Proc_stop();
            }
        }

    });

    Kokkos::parallel_reduce("sat_ite_summax", MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO}, {kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij, int& local_sum, int& local_max){
        local_sum += ite_temp(k,ij);
        if(local_max < ite_temp(k,ij)) local_max = ite_temp(k,ij);
    }, sat_ite_sum, Kokkos::Max<int>(sat_ite_max));

    Kokkos::parallel_reduce("sat_ite_count", MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO}, {kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij, int& update){
        if(ite_temp(k,ij) != 0)
            update += 1;
        else
            ite_temp(k,ij) = itelim;
    }, sat_ite_count);

    Kokkos::parallel_reduce("sat_ite_min", MDRangePolicy<Kokkos::Rank<2>>({kmin,IDX_ZERO}, {kmax+1,ijdim}), 
    KOKKOS_LAMBDA(const size_t k, const size_t ij, int& local_min){
        if(local_min > ite_temp(k,ij)) local_min = ite_temp(k,ij);
    }, Kokkos::Min<int>(sat_ite_min));

}

}

