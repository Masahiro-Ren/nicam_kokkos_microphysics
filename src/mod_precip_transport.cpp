#include "mod_precip_transport.h"

namespace PRECIP_TRANSPORT {

void precip_transport_new(  double rhog               [kdim][ijdim],
                            double rhogvx             [kdim][ijdim],
                            double rhogvy             [kdim][ijdim],
                            double rhogvz             [kdim][ijdim],
                            double rhogw              [kdim][ijdim],
                            double rhoge              [kdim][ijdim],
                            double rhogq              [nqmax][kdim][ijdim],
                            double rho                [kdim][ijdim],
                            double tem                [kdim][ijdim],
                            double pre                [kdim][ijdim],
                            double vx                 [kdim][ijdim],
                            double vy                 [kdim][ijdim],
                            double vz                 [kdim][ijdim],
                            double w                  [kdim][ijdim],
                            double q                  [nqmax][kdim][ijdim],
                            double qd                 [kdim][ijdim],
                            double z                  [kdim][ijdim],
                            double Vterm              [nqmax][kdim][ijdim],
                            bool   precipitating_flag [nqmax],
                            double precip             [2][ijdim],
                            double precip_rhoe        [ijdim],
                            double precip_lh_heat     [ijdim],
                            double precip_rhophi      [ijdim],
                            double precip_rhokin      [ijdim],
                            double frain              [kdim][ijdim],
                            double gsgam2             [kdim][ijdim],
                            double gsgam2h            [kdim][ijdim],
                            double rgs                [kdim][ijdim],
                            double rgsh               [kdim][ijdim],
                            double ix                 [ijdim],
                            double iy                 [ijdim],
                            double iz                 [ijdim],
                            double jx                 [ijdim],
                            double jy                 [ijdim],
                            double jz                 [ijdim],
                            double dt,
                            double **precip_trc = nullptr      // precip[nqmax][ijdim]
                            )
{
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

    double rhogkin   [kdim][ijdim];
    double rhogkin_h [kdim][ijdim];
    double rhogkin_v [kdim][ijdim];

    double rhoq     [kdim][ijdim];
    double rhoeq    [kdim][ijdim];
    double rhophiq  [kdim][ijdim];
    double rhokin_h [kdim][ijdim];
    double rhokin_v [kdim][ijdim];
    double rhouq    [kdim][ijdim];
    double rhovq    [kdim][ijdim];
    double rhowq    [kdim][ijdim];

    double fprec_q        [kdim][ijdim];
    double fprec_rhoe     [kdim][ijdim];
    double fprec_rhophi   [kdim][ijdim];
    double fprec_rhokin_h [kdim][ijdim];
    double fprec_rhokin_v [kdim][ijdim];
    double fprec_rhou     [kdim][ijdim];
    double fprec_rhov     [kdim][ijdim];
    double fprec_rhow     [kdim][ijdim];

    double drhoq     [nqmax][kdim][ijdim];
    double drhoe     [kdim][ijdim];
    double drhophi   [kdim][ijdim];
    double drhokin_h [kdim][ijdim];
    double drhokin_v [kdim][ijdim];
    double drhogu    [kdim][ijdim];
    double drhogv    [kdim][ijdim];
    double drhogw    [kdim][ijdim];

    double kin_h0 [kdim][ijdim];
    double kin_h  [kdim][ijdim];
    double vx_t   [kdim][ijdim];
    double vy_t   [kdim][ijdim];
    double vz_t   [kdim][ijdim];
    double kin_v0 [kdim][ijdim];
    double kin_v  [kdim][ijdim];
    double w_t    [kdim][ijdim];

    double zdis0       [kdim][ijdim];
    int    kcell       [kdim][ijdim];
    int    kcell_max   [kdim];
    int    kcell_min   [kdim];
    double zdis0h      [kdim][ijdim];
    int    kcellh      [kdim][ijdim];
    int    kcellh_max  [kdim];
    int    kcellh_min  [kdim];

    double Vtermh [kdim][ijdim];
    double qh     [kdim][ijdim];
    double rhogh  [kdim][ijdim];
    double ein    [kdim][ijdim];
    double tmp    [kdim][ijdim];
    double tmp_h  [kdim][ijdim];
    double tmp_v  [kdim][ijdim];
    double tmp2   [kdim][ijdim];

    double C2Wfact [2][kdim][ijdim];
    double W2Cfact [2][kdim][ijdim];

    double GRD_gz_shift [kdim];
    double GRAV = CONST_GRAV;


    for(int k = kmin; k <= kmax + 1; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            C2Wfact[0][k][ij] = GRD_afact[k] / gsgam2[k][ij] * gsgam2h[k][ij];
            C2Wfact[1][k][ij] = GRD_bfact[k] / gsgam2[k-1][ij] * gsgam2h[k][ij];
        }
    }

    for(int ij = 0; ij < ijdim; ij++)
    {
        C2Wfact[0][kmin-1][ij] = 0.0;
        C2Wfact[1][kmin-1][ij] = 0.0;
    }

    for(int k = kmin - 1; k <= kmax; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            W2Cfact[0][k][ij] = GRD_cfact[k] * gsgam2[k][ij] / gsgam2h[k+1][ij];
            W2Cfact[1][k][ij] = GRD_dfact[k] * gsgam2[k][ij] / gsgam2h[k+1][ij];
        }
    }

    for(int ij = 0; ij < ijdim; ij++)
    {
        W2Cfact[0][kmax+1][ij] = 0.0; 
        W2Cfact[1][kmax+1][ij] = 0.0; 
    }


    cnvvar_rhogkin_in(rhog, rhogvx, rhogvy, rhogvz, rhogw, C2Wfact, W2Cfact, rhogkin, rhogkin_h, rhogkin_v);

    for(int k = kmin; k <= kmax + 1; k++)
    {
        GRD_gz_shift[k] = GRD_gz[k-1];
    }

    for(int nq = 0; nq < nqmax; nq++)
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                drhoq[nq][k][ij] = 0.0;
            }
        }
    }

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            drhoe     [k][ij] = 0.0;
            drhophi   [k][ij] = 0.0;
            drhokin_h [k][ij] = 0.0;
            drhokin_v [k][ij] = 0.0;
            drhogu    [k][ij] = 0.0;
            drhogv    [k][ij] = 0.0;
            drhogw    [k][ij] = 0.0;
        }
    }

    for(int ij = 0; ij < ijdim; ij++)
    {
        precip         [0][ij] = 0.0;
        precip         [1][ij] = 0.0;
        precip_rhoe    [ij] = 0.0;
        precip_lh_heat [ij] = 0.0;
        precip_rhophi  [ij] = 0.0;
        precip_rhokin  [ij] = 0.0;
    }

    for(int nq = 0; nq < nqmax; nq++)
    {
        if(!precipitating_flag[nq]) continue;

        vadv1d_prep(kmin, kmax, GRD_dgz, GRD_gzh, Vterm[nq], zdis0, kcell, kcell_max, kcell_min, dt);

        // mass
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhoq[k][ij] = rhogq[nq][k][ij] * rgs[k][ij];
            }
        }

        vadv1d_getflux_new(kmin, kmax, GRD_dgz, rhoq, zdis0, kcell, kcell_max, kcell_min, fprec_q);

        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                drhoq[nq][k][ij] = -( fprec_q[k+1][ij] - fprec_q[k][ij] ) * GRD_rdgz[k];
            }
        }

        // -- Only for mass tracer
        if( nq >= NQW_STR && nq <= NQW_END )
        {
            //--- internal energy
            //--- potential energy
            //--- horizontal kinetic energy
            //--- momentum u
            //--- momentum v
            for(int k = 0; k < kdim; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    rhoeq   [k][ij] = rhogq[nq][k][ij] * rgs[k][ij] * CVW[nq] * tem[k][ij];
                    rhophiq [k][ij] = rhogq[nq][k][ij] * rgs[k][ij] * GRAV * z[k][ij];
                    rhokin_h[k][ij] =     q[nq][k][ij] * rgs[k][ij] * rhogkin_h[k][ij];
                    rhouq   [k][ij] = rhogq[nq][k][ij] * rgs[k][ij] * ( vx[k][ij] * ix[ij]
                                                                      + vy[k][ij] * iy[ij]
                                                                      + vz[k][ij] * iz[ij]);
                    rhovq   [k][ij] = rhogq[nq][k][ij] * rgs[k][ij] * ( vx[k][ij] * jx[ij]
                                                                      + vy[k][ij] * jy[ij]
                                                                      + vz[k][ij] * jz[ij]);
                }

            }
            vadv1d_getflux_new(kmin, kmax, GRD_dgz, rhoeq, zdis0, kcell, kcell_max, kcell_min, fprec_rhoe);

            vadv1d_getflux_new(kmin, kmax, GRD_dgz, rhophiq, zdis0, kcell, kcell_max, kcell_min, fprec_rhophi);

            vadv1d_getflux_new(kmin, kmax, GRD_dgz, rhokin_h, zdis0, kcell, kcell_max, kcell_min, fprec_rhokin_h);

            vadv1d_getflux_new(kmin, kmax, GRD_dgz, rhouq, zdis0, kcell, kcell_max, kcell_min, fprec_rhou);

            vadv1d_getflux_new(kmin, kmax, GRD_dgz, rhovq, zdis0, kcell, kcell_max, kcell_min, fprec_rhov);

            for(int k = kmin; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    drhoe[k][ij] = drhoe[k][ij] - ( fprec_rhoe[k+1][ij] - fprec_rhoe[k][ij] ) * GRD_rdgz[k];
                }
            }

            for(int k = kmin; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    drhophi[k][ij] = drhophi[k][ij] - ( fprec_rhophi[k+1][ij] - fprec_rhophi[k][ij] ) * GRD_rdgz[k];
                }
            }

            for(int k = kmin; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    drhokin_h[k][ij] = drhokin_h[k][ij] - ( fprec_rhokin_h[k+1][ij] - fprec_rhokin_h[k][ij] ) * GRD_rdgz[k];
                }
            }

            for(int k = kmin; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    drhogu[k][ij] = drhogu[k][ij] - ( fprec_rhou[k+1][ij] - fprec_rhou[k][ij] ) * GRD_rdgz[k];
                }
            }

            for(int k = kmin; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    drhogv[k][ij] = drhogv[k][ij] - ( fprec_rhov[k+1][ij] - fprec_rhov[k][ij] ) * GRD_rdgz[k];
                }
            }

            // half level

            for(int k = kmin + 1; k <= kmax - 1; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    Vtermh[k][ij] = 0.5 * ( Vterm[nq][k][ij] + Vterm[nq][k-1][ij]);
                }
            }

            for(int ij = 0; ij < ijdim; ij++)
            {
                Vtermh[kmin-1][ij] = 0.0;
                Vtermh[kmin  ][ij] = 0.0;
                Vtermh[kmax  ][ij] = 0.0;
                Vtermh[kmax+1][ij] = 0.0;
            }

            vadv1d_prep(kmin + 1, kmax, GRD_dgzh, GRD_gz_shift, Vtermh, zdis0h, kcellh, kcellh_max, kcellh_min, dt);

            for(int k = kmin + 1; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    qh[k][ij] = 0.5 * ( q[nq][k][ij] + q[nq][k-1][ij] ) * rgsh[k][ij];
                }
            }

            for(int ij = 0; ij < ijdim; ij++)
            {
                qh[kmin-1][ij] = 0.0;
                qh[kmin  ][ij] = 0.0;
                qh[kmax+1][ij] = 0.0;
            }
            //--- vertical kinetic energy
            //--- moment w
            //--- half level

            for(int k = kmin + 1; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    rhokin_v[k][ij] = qh[k][ij] * rhogkin_v[k][ij];
                    rhowq   [k][ij] = qh[k][ij] * rhogw    [k][ij];
                }
            }

            vadv1d_getflux_new(kmin + 1, kmax, GRD_dgzh, rhokin_v, zdis0h, kcellh, kcell_max, kcellh_min, fprec_rhokin_v);

            vadv1d_getflux_new(kmin + 1, kmax, GRD_dgzh, rhowq, zdis0h, kcellh, kcell_max, kcellh_min, fprec_rhow);

            for(int k = kmin + 1; k <= kmax; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    drhokin_v[k][ij] = drhokin_v[k][ij] - ( fprec_rhokin_v[k][ij] - fprec_rhokin_v[k-1][ij] ) * GRD_rdgzh[k];
                    drhogw   [k][ij] = drhogw   [k][ij] - ( fprec_rhow    [k][ij] - fprec_rhow    [k-1][ij] ) * GRD_rdgzh[k];
                }
            }

            // precipitation on the ground
            if( nq == I_QC )
            {
                for(int ij = 0; ij < ijdim; ij++)
                    precip[0][ij] = precip[0][ij] - fprec_q[kmin][ij] / dt;
            }
            else if(nq == I_QR)
            {
                for(int ij = 0; ij < ijdim; ij++)
                    precip[0][ij] = precip[0][ij] - fprec_q[kmin][ij] / dt;

                for(int k = kmin; k <= kmax; k++)
                {
                    for(int ij = 0; ij < ijdim; ij++)
                    {
                        frain[k][ij] = -fprec_q[k][ij] / dt;
                    }
                }

                for(int ij = 0; ij < ijdim; ij++)
                {
                    frain[kmin-1][ij] = 0.0;
                    frain[kmax+1][ij] = 0.0;
                }
            }
            else if( nq == I_QI || nq == I_QS || nq == I_QG )
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    precip[1][ij] = precip[1][ij] - fprec_q[kmin][ij] / dt;
                    precip_lh_heat[ij] = precip_lh_heat[ij] + fprec_q[kmin][ij] * LHF / dt;
                }
            }

            if(precip_trc != nullptr)
            {
                for(int ij = 0; ij < ijdim; ij++)
                    precip_trc[nq][ij] = precip_trc[nq][ij] - fprec_q[kmin][ij] / dt;
            }

            for(int ij = 0; ij < ijdim; ij++)
            {
                precip_rhoe  [ij] = precip_rhoe[ij] - fprec_rhoe[kmin][ij] / dt;
                precip_rhophi[ij] = precip_rhophi[ij] - fprec_rhophi[kmin][ij] / dt;

                precip_rhokin[ij] = precip_rhokin[ij] - fprec_rhokin_h[kmin][ij] / dt
                                                      - fprec_rhokin_v[kmin][ij] / dt;
            }
        }
    } // tracer loop

    // Change in internal energy comes from precipitation and dissipation of kinetic energy due to drag force.
    // See Ooyama(2001) (3.13)
    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rhoge    [k][ij] = rhoge    [k][ij] + drhoe[k][ij] + drhophi[k][ij];
            rhogkin_h[k][ij] = rhogkin_h[k][ij] + drhokin_h[k][ij];
            rhogkin_v[k][ij] = rhogkin_v[k][ij] + drhokin_v[k][ij];
            rhogvx   [k][ij] = rhogvx   [k][ij] + drhogu[k][ij] * ix[ij] + drhogv[k][ij] * jx[ij];
            rhogvy   [k][ij] = rhogvy   [k][ij] + drhogu[k][ij] * iy[ij] + drhogv[k][ij] * jy[ij];
            rhogvz   [k][ij] = rhogvz   [k][ij] + drhogu[k][ij] * iz[ij] + drhogv[k][ij] * jz[ij];
            rhogw    [k][ij] = rhogw    [k][ij] + drhogw[k][ij];
        }
    }

    for(int nq = 0; nq < nqmax; nq++)
    {
        if( nq >= NQW_STR && nq <= NQW_END )
        {
            for(int k = 0; k < kdim; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    rhogq[nq][k][ij] = rhogq[nq][k][ij] + drhoq[nq][k][ij];
                    rhog [k][ij]     = rhog[k][ij] + drhoq[nq][k][ij];
                    rhoge[k][ij]     = rhoge[k][ij] - drhoq[nq][k][ij] * GRAV * z[k][ij];
                }
            }
        }
        else
        {
            for(int k = 0; k < kdim; k++)
            {
                for(int ij = 0; ij < ijdim; ij++)
                {
                    rhogq[nq][k][ij] = rhogq[nq][k][ij] + drhoq[nq][k][ij];
                }
            }
        }
    }

    if( PRCIP_TRN_ECORRECT == "KIN2EIN" )
    {
        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                tmp2[k][ij] = rhogkin_h[k][ij] + ( W2Cfact[0][k][ij] * rhogkin_v[k+1][ij] 
                                                 + W2Cfact[1][k][ij] * rhogkin_v[k][ij] );
            }
        }

        for(int ij = 0; ij < ijdim; ij++)
        {
            tmp2[kmin-1][ij] = 0.0;
            tmp2[kmax+1][ij] = 0.0;
        }

        cnvvar_rhogkin_in(rhog, rhogvx, rhogvy, rhogvz, rhogw, C2Wfact, W2Cfact, tmp, tmp_h, tmp_v);

        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhoge[k][ij] = rhoge[k][ij] + ( tmp2[k][ij] - tmp[k][ij] );
            }
        }
    }
    else if(PRCIP_TRN_ECORRECT == "KIN2KIN")
    {
        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                kin_h0[k][ij] = rhogkin_h[k][ij] / rhog[k][ij];
                vx_t  [k][ij] = rhogvx   [k][ij] / rhog[k][ij];
                vy_t  [k][ij] = rhogvy   [k][ij] / rhog[k][ij];
                vz_t  [k][ij] = rhogvz   [k][ij] / rhog[k][ij];

                kin_h[k][ij] = 0.5 * ( vx_t[k][ij] * vx_t[k][ij] 
                                     + vy_t[k][ij] * vy_t[k][ij]
                                     + vz_t[k][ij] * vz_t[k][ij] );
            }
        }

        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                if(kin_h[k][ij] > 1.0E-20)
                {
                    vx_t[k][ij] = vx_t[k][ij] * std::sqrt( std::abs(kin_h0[k][ij] / kin_h[k][ij]) );
                    vy_t[k][ij] = vy_t[k][ij] * std::sqrt( std::abs(kin_h0[k][ij] / kin_h[k][ij]) );
                    vz_t[k][ij] = vz_t[k][ij] * std::sqrt( std::abs(kin_h0[k][ij] / kin_h[k][ij]) );
                }
            }
        }

        for(int k = kmin; k <= kmax; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhogvx[k][ij] = vx_t[k][ij] * rhog[k][ij];
                rhogvy[k][ij] = vy_t[k][ij] * rhog[k][ij];
                rhogvz[k][ij] = vz_t[k][ij] * rhog[k][ij];
            }
        }

        for(int k = kmin; k <= kmax + 1; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhogh[k][ij] = ( C2Wfact[0][k][ij] * rhog[k][ij] 
                                + C2Wfact[1][k][ij] * rhog[k-1][ij]);
            }
        }

        for(int k = kmin; k <= kmax + 1; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                kin_h0[k][ij] = rhogkin_v[k][ij] / rhogh[k][ij];
                w_t   [k][ij] = rhogw[k][ij] / rhogh[k][ij];
                kin_v [k][ij] = 0.5 * std::pow(w_t[k][ij], 2);
            }
        }

        for(int k = kmin; k <= kmax + 1; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                if( kin_v[k][ij] > 1.0E-20 )
                {
                    w_t[k][ij] = w_t[k][ij] * std::sqrt( std::abs(kin_v0[k][ij] / kin_v[k][ij]) );
                }
            }
        }

        for(int k = kmin; k <= kmax + 1; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                rhogw[k][ij] = w_t[k][ij] * rhogh[k][ij];
            }
        }
    }
    else
    {
        std::cerr << "Error in PRCIP_TRN_ECORRECT: " << PRCIP_TRN_ECORRECT << std::endl;
    }

    for(int nq = 0; nq < nqmax; nq++)
    {
        for(int k = 0; k < kdim; k++)
        {
            for(int ij = 0; ij < ijdim; ij++)
            {
                q[nq][k][ij] = rhogq[nq][k][ij] / rhog[k][ij];
            }
        }
    }

    for(int k = 0; k < kdim; k++)
    {
        for(int ij = 0; ij < ijdim; ij++)
        {
            rho[k][ij] = rhog[k][ij] / gsgam2[k][ij];
            ein[k][ij] = rhoge[k][ij] / rhog[k][ij];
        }
    }

    THRMDYN_qd(q, qd);
    THRMDYN_tempre(ein, rho, q, tem, pre);
}

}
