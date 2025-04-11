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


/**
 * For result checking
 */
// double CHECK_rhog  [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_rhogvx[ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_rhogvy[ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_rhogvz[ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_rhogw [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_rhoge [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_QV1   [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_QC2   [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_QR3   [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_QI4   [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_QS5   [ADM_lall][ADM_kall][ADM_gall_in];
// double CHECK_QG6   [ADM_lall][ADM_kall][ADM_gall_in];

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

    View<double***> GDCLW  ("GDCLW ", ADM_lall, ADM_kall, ADM_gall_in);
    View<double***> GDCFRC ("GDCFRC", ADM_lall, ADM_kall, ADM_gall_in);

    /**
     * For result checking
     */
    View<double***> CHECK_rhog  ("CHECK_rhog  ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_rhogvx("CHECK_rhogvx",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_rhogvy("CHECK_rhogvy",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_rhogvz("CHECK_rhogvz",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_rhogw ("CHECK_rhogw ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_rhoge ("CHECK_rhoge ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_QV1   ("CHECK_QV1   ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_QC2   ("CHECK_QC2   ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_QR3   ("CHECK_QR3   ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_QI4   ("CHECK_QI4   ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_QS5   ("CHECK_QS5   ",ADM_lall,ADM_kall,ADM_gall_in);
    View<double***> CHECK_QG6   ("CHECK_QG6   ",ADM_lall,ADM_kall,ADM_gall_in);
    
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

    /**
     * Vertical grid setup
     */
    GRD_Setup();
    // /**
    //  * Saturation set_up
    //  */
    SATURATION_Setup();
    // /**
    //  * microphysics initialization
    //  */
    mp_init(MP_TYPE);

    int l = SET_l;

    std::cout << "============= Finish Initialize =============== \n";

    // /**
    //  * Start Simulation
    // */
    std::cout << "============= Start Kernel =============== \n";

    // /**
    //  * create sub views:
    //  */
    auto sub_rhog   = subview(rhog  , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rhogvx = subview(rhogvx, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rhogvy = subview(rhogvy, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rhogvz = subview(rhogvz, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rhogw  = subview(rhogw , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rhoge  = subview(rhoge , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_vx     = subview(vx    , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_vy     = subview(vy    , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_vz     = subview(vz    , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_w      = subview(w     , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_unccn  = subview(unccn , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rho    = subview(rho   , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_pre    = subview(pre   , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_tem    = subview(tem   , 0, Kokkos::ALL(), Kokkos::ALL());

    auto sub_rhogq_Lswp = subview(rhogq_Lswp, 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto sub_q_Lswp     = subview(q_Lswp    , 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    auto sub_precip_mp  = subview(precip_mp , Kokkos::ALL(), 0, 0, Kokkos::ALL());
    auto sub_precip1_mp = subview(precip1_mp, Kokkos::ALL(), 0, 0, Kokkos::ALL());
    auto sub_precip2_mp = subview(precip2_mp, Kokkos::ALL(), 0, 0, Kokkos::ALL());

    auto sub_rhoein_precip_mp = subview(rhoein_precip_mp, 0, 0, Kokkos::ALL());
    auto sub_lh_precip_mp     = subview(lh_precip_mp    , 0, 0, Kokkos::ALL());
    auto sub_rhophi_precip_mp = subview(rhophi_precip_mp, 0, 0, Kokkos::ALL());
    auto sub_rhokin_precip_mp = subview(rhokin_precip_mp, 0, 0, Kokkos::ALL());

    auto sub_frhoge_af  = subview(frhoge_af , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_frhogqv_af = subview(frhogqv_af, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_frhoge_rad = subview(frhoge_rad, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_qke        = subview(qke       , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_gsgam2     = subview(gsgam2    , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_gsgam2h    = subview(gsgam2h   , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_gam2       = subview(gam2      , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_gam2h      = subview(gam2h     , 0, Kokkos::ALL(), Kokkos::ALL());

    auto sub_ix = subview(ix, 0, Kokkos::ALL());
    auto sub_iy = subview(iy, 0, Kokkos::ALL());
    auto sub_iz = subview(iz, 0, Kokkos::ALL());
    auto sub_jx = subview(jx, 0, Kokkos::ALL());
    auto sub_jy = subview(jy, 0, Kokkos::ALL());
    auto sub_jz = subview(jz, 0, Kokkos::ALL());

    auto sub_z           = subview(z           , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_zh          = subview(zh          , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_GPREC       = subview(GPREC       , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_CBMFX       = subview(CBMFX       , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_qd          = subview(qd          , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rceff       = subview(rceff       , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rceff_solid = subview(rceff_solid , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rceff_cld   = subview(rceff_cld   , 0, Kokkos::ALL(), Kokkos::ALL());

    auto sub_rctop = subview(rctop, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_rwtop = subview(rwtop, 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_tctop = subview(tctop, 0, Kokkos::ALL(), Kokkos::ALL());

    auto sub_GDCLW  = subview(GDCLW , 0, Kokkos::ALL(), Kokkos::ALL());
    auto sub_GDCFRC = subview(GDCFRC, 0, Kokkos::ALL(), Kokkos::ALL());

    for(int i = 0; i < SET_iteration; i++)
    {
        // Call mp_driver
        mp_driver(
                  l,
                  sub_rhog             ,
                  sub_rhogvx           ,
                  sub_rhogvy           ,
                  sub_rhogvz           ,
                  sub_rhogw            ,
                  sub_rhoge            ,
                  sub_rhogq_Lswp       ,
                  sub_vx               ,
                  sub_vy               ,
                  sub_vz               ,
                  sub_w                ,
                  sub_unccn            ,
                  sub_rho              ,
                  sub_tem              ,
                  sub_pre              ,
                  sub_q_Lswp           ,
                  sub_qd               ,
                  sub_precip_mp        ,  
                  sub_precip1_mp       ,
                  sub_precip2_mp       ,
                  sub_rhoein_precip_mp ,
                  sub_lh_precip_mp     ,
                  sub_rhophi_precip_mp ,
                  sub_rhokin_precip_mp ,
                  sub_rceff            ,
                  sub_rceff_solid      ,
                  sub_rceff_cld        ,
                  sub_rctop            ,
                  sub_rwtop            ,
                  sub_tctop            ,
                  sub_frhoge_af        ,
                  sub_frhogqv_af       ,
                  sub_frhoge_rad       ,
                  sub_qke              ,
                  sub_gsgam2           ,
                  sub_gsgam2h          ,
                  sub_gam2             ,
                  sub_gam2h            ,
                  sub_ix               ,
                  sub_iy               ,
                  sub_iz               ,
                  sub_jx               ,
                  sub_jy               ,
                  sub_jz               ,
                  sub_z                ,
                  sub_zh               ,
                  TIME_DTL,
                  TIME_CTIME,
                  sub_GDCLW            ,
                  sub_GDCFRC           ,
                  sub_GPREC            ,
                  sub_CBMFX           
                  );
    }
    std::cout << "============= Finish Kernel =============== \n";

    if (SET_check)
    {
        std::cout << "Checking Reuslts \n";

        read_data_3d("ref_verify/calculated_rhog_DP.dat", CHECK_rhog);
        read_data_3d("ref_verify/calculated_rhogvx_DP.dat", CHECK_rhogvx);
        read_data_3d("ref_verify/calculated_rhogvy_DP.dat", CHECK_rhogvy);
        read_data_3d("ref_verify/calculated_rhogvz_DP.dat", CHECK_rhogvz);
        read_data_3d("ref_verify/calculated_rhoge_DP.dat", CHECK_rhoge);
        read_data_3d("ref_verify/calculated_rhogw_DP.dat", CHECK_rhogw);
        read_data_3d("ref_verify/calculated_QV1_DP.dat", CHECK_QV1);
        read_data_3d("ref_verify/calculated_QC2_DP.dat", CHECK_QC2);
        read_data_3d("ref_verify/calculated_QR3_DP.dat", CHECK_QR3);
        read_data_3d("ref_verify/calculated_QI4_DP.dat", CHECK_QI4);
        read_data_3d("ref_verify/calculated_QS5_DP.dat", CHECK_QS5);
        read_data_3d("ref_verify/calculated_QG6_DP.dat", CHECK_QG6);

        PROF_val_check("rhog",   rhog,   CHECK_rhog);
        PROF_val_check("rhogvx", rhogvx, CHECK_rhogvx);
        PROF_val_check("rhogvy", rhogvy, CHECK_rhogvy);
        PROF_val_check("rhogvz", rhogvz, CHECK_rhogvz);
        PROF_val_check("rhoge",  rhoge,  CHECK_rhoge);
        PROF_val_check("rhogw",  rhogw,  CHECK_rhogw);
        PROF_val_check("QV", rhogq_Lswp, I_QV, CHECK_QV1);
        PROF_val_check("QC", rhogq_Lswp, I_QC, CHECK_QC2);
        PROF_val_check("QR", rhogq_Lswp, I_QR, CHECK_QR3);
        PROF_val_check("QI", rhogq_Lswp, I_QI, CHECK_QI4);
        PROF_val_check("QS", rhogq_Lswp, I_QS, CHECK_QS5);
        PROF_val_check("QG", rhogq_Lswp, I_QG, CHECK_QG6);
    }

    std::cout << "============= All process finished =============== \n";
}
Kokkos::finalize();
}
