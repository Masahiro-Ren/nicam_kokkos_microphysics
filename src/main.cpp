#include "problemsize.h"
#include "data_io.h"
#include "mod_debug.h"
#include "mod_satadjust.h"
#include "mod_mp_driver.h"

using namespace PROBLEM_SIZE;
using namespace DATA_IO;
using namespace DEBUG;
using namespace SATADJUST;
using namespace MP_DRIVER;

double GDCLW  [ADM_lall][ADM_kall][ADM_gall_in];
double GDCFRC [ADM_lall][ADM_kall][ADM_gall_in];

/**
 * For result checking
 */
double CHECK_rhog  [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_rhogvx[ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_rhogvy[ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_rhogvz[ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_rhogw [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_rhoge [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_QV1   [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_QC2   [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_QR3   [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_QI4   [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_QS5   [ADM_lall][ADM_kall][ADM_gall_in];
double CHECK_QG6   [ADM_lall][ADM_kall][ADM_gall_in];

int main(int argc, char* argv[])
{
Kokkos::initialize(argc, argv);
{
    // declare all variables
    View<double***> rhog  ("rhog  ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rhogvx("rhogvx", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rhogvy("rhogvy", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rhogvz("rhogvz", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rhogw ("rhogw ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rhoge ("rhoge ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> vx    ("vx    ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> vy    ("vy    ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> vz    ("vz    ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> w     ("w     ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> unccn ("unccn ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rho   ("rho   ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> pre   ("pre   ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> tem   ("tem   ", ADM_lall, ADM_kall, ADM_gall_in);

    View<double****> rhogq_Lswp ("rhogq_Lswp", ADM_lall, TRC_VMAX, ADM_kall, ADM_gall_in);
    View<double****> q_Lswp     ("q_Lswp    ", ADM_lall, TRC_VMAX, ADM_kall, ADM_gall_in);

    View<double****> precip_mp  ("precip_mp ", 2, ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double****> precip1_mp ("precip1_mp", 2, ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double****> precip2_mp ("precip2_mp", 2, ADM_lall, ADM_KNONE, ADM_gall_in);

    View<double***> rhoein_precip_mp ("rhoein_precip_mp", ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double***> lh_precip_mp     ("lh_precip_mp    ", ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double***> rhophi_precip_mp ("rhophi_precip_mp", ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double***> rhokin_precip_mp ("rhokin_precip_mp", ADM_lall, ADM_KNONE, ADM_gall_in);

    View<double***> frhoge_af ("frhoge_af ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> frhogqv_af("frhogqv_af", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> frhoge_rad("frhoge_rad", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> qke       ("qke       ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> gsgam2    ("gsgam2    ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> gsgam2h   ("gsgam2h   ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> gam2      ("gam2      ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> gam2h     ("gam2h     ", ADM_lall, ADM_kall, ADM_gall_in);

    View<double**> ix ("ix", ADM_lall, ADM_gall_in);
    View<double**> iy ("iy", ADM_lall, ADM_gall_in);
    View<double**> iz ("iz", ADM_lall, ADM_gall_in);
    View<double**> jx ("jx", ADM_lall, ADM_gall_in);
    View<double**> jy ("jy", ADM_lall, ADM_gall_in);
    View<double**> jz ("jz", ADM_lall, ADM_gall_in);

    View<double***> z           ("z          ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> zh          ("zh         ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> GPREC       ("GPREC      ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> CBMFX       ("CBMFX      ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> qd          ("qd         ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rceff       ("rceff      ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rceff_solid ("rceff_solid", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> rceff_cld   ("rceff_cld  ", ADM_lall, ADM_kall, ADM_gall_in);

    View<double***> rctop ("rctop", ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double***> rwtop ("rwtop", ADM_lall, ADM_KNONE, ADM_gall_in);
    View<double***> tctop ("tctop", ADM_lall, ADM_KNONE, ADM_gall_in);

    /**
     * Display Simulation Configuerations
    */
    std::cout << "[KERNEL] physicskernel_microphysics \n";
    std::cout << "ADM_gall_in_orig = " <<  PROBLEM_SIZE::ADM_gall_in_orig << std::endl;
    std::cout << "ADM_kall = " << PROBLEM_SIZE::ADM_kall << std::endl;
    std::cout << "SET_l = " << PROBLEM_SIZE::SET_l << std::endl;
    std::cout << "MP_TYPE = " << PROBLEM_SIZE::MP_TYPE << std::endl;
    std::cout << "PI = " << PROBLEM_SIZE::PI << std::endl;
    std::cout << "EPS = " << PROBLEM_SIZE::EPS << std::endl;

    std::cout << "============= Start Initialize =============== \n";

    read_data_3d("data/rhog.dat", rhog);
    read_data_3d("data/rhogvx.dat", rhogvx);
    read_data_3d("data/rhogvy.dat", rhogvy);
    read_data_3d("data/rhogvz.dat", rhogvz);
    read_data_3d("data/rhogw.dat", rhogw);
    read_data_3d("data/rhoge.dat", rhoge);
    read_data_3d("data/vx.dat", vx);
    read_data_3d("data/vy.dat", vy);
    read_data_3d("data/vz.dat", vz);
    read_data_3d("data/w.dat", w);
    read_data_3d("data/unccn.dat", unccn);
    read_data_3d("data/rho.dat", rho);
    read_data_3d("data/pre.dat", pre);
    read_data_3d("data/tem.dat", tem);

    read_data_4d("data/rhogq_Lswp.dat", rhogq_Lswp);
    read_data_4d("data/q_Lswp.dat", q_Lswp);

    read_data_4d("data/precip_mp.dat", precip_mp);
    read_data_4d("data/precip1_mp.dat", precip1_mp);
    read_data_4d("data/precip2_mp.dat", precip2_mp);

    read_data_3d("data/rhoein_precip_mp.dat", rhoein_precip_mp);
    read_data_3d("data/lh_precip_mp.dat", lh_precip_mp);
    read_data_3d("data/rhophi_precip_mp.dat", rhophi_precip_mp);
    read_data_3d("data/rhokin_precip_mp.dat", rhokin_precip_mp);

    read_data_3d("data/frhoge_af.dat", frhoge_af);
    read_data_3d("data/frhogqv_af.dat", frhogqv_af);
    read_data_3d("data/frhoge_rad.dat", frhoge_rad);
    read_data_3d("data/qke.dat", qke);
    read_data_3d("data/gsgam2.dat", gsgam2);
    read_data_3d("data/gsgam2h.dat", gsgam2h);
    read_data_3d("data/gam2.dat", gam2);
    read_data_3d("data/gam2h.dat", gam2h);

    read_data_2d("data/ix.dat", ix);
    read_data_2d("data/iy.dat", iy);
    read_data_2d("data/iz.dat", iz);
    read_data_2d("data/jx.dat", jx);
    read_data_2d("data/jy.dat", jy);
    read_data_2d("data/jz.dat", jz);

    read_data_3d("data/z.dat", z);
    read_data_3d("data/zh.dat", zh);
    read_data_3d("data/GPREC.dat", GPREC);
    read_data_3d("data/CBMFX.dat", CBMFX);

    // /**
    //  * Vertical grid setup
    //  */
    // GRD_Setup();
    // /**
    //  * Saturation set_up
    //  */
    // SATURATION_Setup();
    // /**
    //  * microphysics initialization
    //  */
    // mp_init(MP_TYPE);

    // int l = SET_l;

    // std::cout << "============= Finish Initialize =============== \n";

    // /**
    //  * Start Simulation
    // */
    // std::cout << "============= Start Kernel =============== \n";

    // /**
    //  * create temp arrays to stroe :
    //  *      precip_mp(:,ADM_KNONE,l,:)
    //  *      precip1_mp(:,ADM_KNONE,l,:)
    //  *      precip2_mp(:,ADM_KNONE,l,:)
    //  */
    // double precip_mp_tmp [2][ADM_gall_in];
    // double precip1_mp_tmp[2][ADM_gall_in];
    // double precip2_mp_tmp[2][ADM_gall_in];

    // for(int k = 0; k < 2; k++)
    // {
    //     for(int ij = 0; ij < ADM_gall_in; ij++)
    //     {
    //         precip_mp_tmp [k][ij] = precip_mp [k][0][0][ij];
    //         precip1_mp_tmp[k][ij] = precip1_mp[k][0][0][ij];
    //         precip2_mp_tmp[k][ij] = precip2_mp[k][0][0][ij];
    //     }
    // }
    // // Copy end

    // for(int i = 0; i < SET_iteration; i++)
    // {
    //     // Call mp_driver
    //     mp_driver(
    //               l,
    //               rhog            [0],
    //               rhogvx          [0],
    //               rhogvy          [0],
    //               rhogvz          [0],
    //               rhogw           [0],
    //               rhoge           [0],
    //               rhogq_Lswp      [0],
    //               vx              [0],
    //               vy              [0],
    //               vz              [0],
    //               w               [0],
    //               unccn           [0],
    //               rho             [0],
    //               tem             [0],
    //               pre             [0],
    //               q_Lswp          [0],
    //               qd              [0],
    //               precip_mp_tmp,
    //               precip1_mp_tmp,
    //               precip2_mp_tmp,
    //               rhoein_precip_mp[0][0],
    //               lh_precip_mp    [0][0],
    //               rhophi_precip_mp[0][0],
    //               rhokin_precip_mp[0][0],
    //               rceff           [0],
    //               rceff_solid     [0],
    //               rceff_cld       [0],
    //               rctop           [0],
    //               rwtop           [0],
    //               tctop           [0],
    //               frhoge_af       [0],
    //               frhogqv_af      [0],
    //               frhoge_rad      [0],
    //               qke             [0],
    //               gsgam2          [0],
    //               gsgam2h         [0],
    //               gam2            [0],
    //               gam2h           [0],
    //               ix              [0],
    //               iy              [0],
    //               iz              [0],
    //               jx              [0],
    //               jy              [0],
    //               jz              [0],
    //               z               [0],
    //               zh              [0],
    //               TIME_DTL,
    //               TIME_CTIME,
    //               GDCLW           [0],
    //               GDCFRC          [0],
    //               GPREC           [0],
    //               CBMFX           [0]
    //               );
    // }
    // std::cout << "============= Finish Kernel =============== \n";

    // if (SET_check)
    // {
    //     std::cout << "Checking Reuslts \n";

    //     read_data_3d("ref_verify/calculated_rhog_DP.dat", CHECK_rhog);
    //     read_data_3d("ref_verify/calculated_rhogvx_DP.dat", CHECK_rhogvx);
    //     read_data_3d("ref_verify/calculated_rhogvy_DP.dat", CHECK_rhogvy);
    //     read_data_3d("ref_verify/calculated_rhogvz_DP.dat", CHECK_rhogvz);
    //     read_data_3d("ref_verify/calculated_rhoge_DP.dat", CHECK_rhoge);
    //     read_data_3d("ref_verify/calculated_rhogw_DP.dat", CHECK_rhogw);
    //     read_data_3d("ref_verify/calculated_QV1_DP.dat", CHECK_QV1);
    //     read_data_3d("ref_verify/calculated_QC2_DP.dat", CHECK_QC2);
    //     read_data_3d("ref_verify/calculated_QR3_DP.dat", CHECK_QR3);
    //     read_data_3d("ref_verify/calculated_QI4_DP.dat", CHECK_QI4);
    //     read_data_3d("ref_verify/calculated_QS5_DP.dat", CHECK_QS5);
    //     read_data_3d("ref_verify/calculated_QG6_DP.dat", CHECK_QG6);

    //     PROF_val_check("rhog",   rhog,   CHECK_rhog);
    //     PROF_val_check("rhogvx", rhogvx, CHECK_rhogvx);
    //     PROF_val_check("rhogvy", rhogvy, CHECK_rhogvy);
    //     PROF_val_check("rhogvz", rhogvz, CHECK_rhogvz);
    //     PROF_val_check("rhoge",  rhoge,  CHECK_rhoge);
    //     PROF_val_check("rhogw",  rhogw,  CHECK_rhogw);
    //     PROF_val_check("QV", rhogq_Lswp[0][I_QV], CHECK_QV1[0]);
    //     PROF_val_check("QC", rhogq_Lswp[0][I_QC], CHECK_QC2[0]);
    //     PROF_val_check("QR", rhogq_Lswp[0][I_QR], CHECK_QR3[0]);
    //     PROF_val_check("QI", rhogq_Lswp[0][I_QI], CHECK_QI4[0]);
    //     PROF_val_check("QS", rhogq_Lswp[0][I_QS], CHECK_QS5[0]);
    //     PROF_val_check("QG", rhogq_Lswp[0][I_QG], CHECK_QG6[0]);
    // }

    std::cout << "============= All process finished =============== \n";
}
Kokkos::finalize();
}
