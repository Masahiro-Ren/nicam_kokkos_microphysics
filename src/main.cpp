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


int main(int argc, char* argv[])
{
Kokkos::initialize(argc, argv);
{
    // declare all variables
    View3D<double, Kokkos::CudaSpace> rhog  ("rhog  ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhogvx("rhogvx", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhogvy("rhogvy", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhogvz("rhogvz", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhogw ("rhogw ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhoge ("rhoge ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> vx    ("vx    ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> vy    ("vy    ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> vz    ("vz    ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> w     ("w     ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> unccn ("unccn ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rho   ("rho   ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> pre   ("pre   ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> tem   ("tem   ", ADM_lall, ADM_kall, ADM_gall_in);

    View4D<double, Kokkos::CudaSpace> rhogq_Lswp ("rhogq_Lswp", ADM_lall, TRC_VMAX, ADM_kall, ADM_gall_in);
    View4D<double, Kokkos::CudaSpace> q_Lswp     ("q_Lswp    ", ADM_lall, TRC_VMAX, ADM_kall, ADM_gall_in);

    View4D<double, Kokkos::CudaSpace> precip_mp  ("precip_mp ", 2, ADM_lall, ADM_KNONE, ADM_gall_in);
    View4D<double, Kokkos::CudaSpace> precip1_mp ("precip1_mp", 2, ADM_lall, ADM_KNONE, ADM_gall_in);
    View4D<double, Kokkos::CudaSpace> precip2_mp ("precip2_mp", 2, ADM_lall, ADM_KNONE, ADM_gall_in);

    View3D<double, Kokkos::CudaSpace> rhoein_precip_mp ("rhoein_precip_mp", ADM_lall, ADM_KNONE, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> lh_precip_mp     ("lh_precip_mp    ", ADM_lall, ADM_KNONE, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhophi_precip_mp ("rhophi_precip_mp", ADM_lall, ADM_KNONE, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rhokin_precip_mp ("rhokin_precip_mp", ADM_lall, ADM_KNONE, ADM_gall_in);

    View3D<double, Kokkos::CudaSpace> frhoge_af ("frhoge_af ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> frhogqv_af("frhogqv_af", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> frhoge_rad("frhoge_rad", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> qke       ("qke       ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> gsgam2    ("gsgam2    ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> gsgam2h   ("gsgam2h   ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> gam2      ("gam2      ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> gam2h     ("gam2h     ", ADM_lall, ADM_kall, ADM_gall_in);

    View2D<double, Kokkos::CudaSpace> ix ("ix", ADM_lall, ADM_gall_in);
    View2D<double, Kokkos::CudaSpace> iy ("iy", ADM_lall, ADM_gall_in);
    View2D<double, Kokkos::CudaSpace> iz ("iz", ADM_lall, ADM_gall_in);
    View2D<double, Kokkos::CudaSpace> jx ("jx", ADM_lall, ADM_gall_in);
    View2D<double, Kokkos::CudaSpace> jy ("jy", ADM_lall, ADM_gall_in);
    View2D<double, Kokkos::CudaSpace> jz ("jz", ADM_lall, ADM_gall_in);

    View3D<double, Kokkos::CudaSpace> z           ("z          ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> zh          ("zh         ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> GPREC       ("GPREC      ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> CBMFX       ("CBMFX      ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> qd          ("qd         ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rceff       ("rceff      ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rceff_solid ("rceff_solid", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rceff_cld   ("rceff_cld  ", ADM_lall, ADM_kall, ADM_gall_in);

    View3D<double, Kokkos::CudaSpace> rctop ("rctop", ADM_lall, ADM_KNONE, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> rwtop ("rwtop", ADM_lall, ADM_KNONE, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> tctop ("tctop", ADM_lall, ADM_KNONE, ADM_gall_in);

    View3D<double, Kokkos::CudaSpace> GDCLW  ("GDCLW ", ADM_lall, ADM_kall, ADM_gall_in);
    View3D<double, Kokkos::CudaSpace> GDCFRC ("GDCFRC", ADM_lall, ADM_kall, ADM_gall_in);

    GRD_gz   = View1D<double, Kokkos::SharedSpace> ("GRD_gz   ", ADM_kall);
    GRD_gzh  = View1D<double, Kokkos::SharedSpace> ("GRD_gzh  ", ADM_kall);
    GRD_dgz  = View1D<double, Kokkos::SharedSpace> ("GRD_dgz  ", ADM_kall);
    GRD_dgzh = View1D<double, Kokkos::SharedSpace> ("GRD_dgzh ", ADM_kall);
    GRD_rdgz = View1D<double, Kokkos::SharedSpace> ("GRD_rdgz ", ADM_kall);
    GRD_rdgzh= View1D<double, Kokkos::SharedSpace> ("GRD_rdgzh", ADM_kall);
    GRD_afact= View1D<double, Kokkos::SharedSpace> ("GRD_afact", ADM_kall);
    GRD_bfact= View1D<double, Kokkos::SharedSpace> ("GRD_bfact", ADM_kall);
    GRD_cfact= View1D<double, Kokkos::SharedSpace> ("GRD_cfact", ADM_kall);
    GRD_dfact= View1D<double, Kokkos::SharedSpace> ("GRD_dfact", ADM_kall);

    CVW = View1D<double, Kokkos::SharedSpace>("CVW", 6);
    CPW = View1D<double, Kokkos::SharedSpace>("CPW", 6);

    /**
     * For result checking
     */
    View3D<double, Kokkos::HostSpace> CHECK_rhog  ("CHECK_rhog  ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_rhogvx("CHECK_rhogvx",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_rhogvy("CHECK_rhogvy",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_rhogvz("CHECK_rhogvz",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_rhogw ("CHECK_rhogw ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_rhoge ("CHECK_rhoge ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_QV1   ("CHECK_QV1   ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_QC2   ("CHECK_QC2   ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_QR3   ("CHECK_QR3   ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_QI4   ("CHECK_QI4   ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_QS5   ("CHECK_QS5   ",ADM_lall,ADM_kall,ADM_gall_in);
    View3D<double, Kokkos::HostSpace> CHECK_QG6   ("CHECK_QG6   ",ADM_lall,ADM_kall,ADM_gall_in);
    
    /**
     * Display Execution Configurations
     */
    Kokkos::print_configuration(std::cout);
    std::cout << "Number of threads: "
                  << Kokkos::DefaultExecutionSpace().concurrency()
                  << std::endl;
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

    /**
     * Create Host Mirror
     */
    // Example: View3D<double, Kokkos::CudaSpace>::HostMirror h_rhog   = Kokkos::create_mirror_view(rhog  );
    auto h_rhog   = Kokkos::create_mirror_view(rhog  );
    auto h_rhogvx = Kokkos::create_mirror_view(rhogvx);
    auto h_rhogvy = Kokkos::create_mirror_view(rhogvy);
    auto h_rhogvz = Kokkos::create_mirror_view(rhogvz);
    auto h_rhogw  = Kokkos::create_mirror_view(rhogw );
    auto h_rhoge  = Kokkos::create_mirror_view(rhoge );
    auto h_vx     = Kokkos::create_mirror_view(vx    );
    auto h_vy     = Kokkos::create_mirror_view(vy    );
    auto h_vz     = Kokkos::create_mirror_view(vz    );
    auto h_w      = Kokkos::create_mirror_view(w     );
    auto h_unccn  = Kokkos::create_mirror_view(unccn );
    auto h_rho    = Kokkos::create_mirror_view(rho   );
    auto h_pre    = Kokkos::create_mirror_view(pre   );
    auto h_tem    = Kokkos::create_mirror_view(tem   );

    auto h_rhogq_Lswp = Kokkos::create_mirror_view(rhogq_Lswp);
    auto h_q_Lswp     = Kokkos::create_mirror_view(q_Lswp    );

    auto h_precip_mp  = Kokkos::create_mirror_view(precip_mp );
    auto h_precip1_mp = Kokkos::create_mirror_view(precip1_mp);
    auto h_precip2_mp = Kokkos::create_mirror_view(precip2_mp);

    auto h_rhoein_precip_mp = Kokkos::create_mirror_view(rhoein_precip_mp);
    auto h_lh_precip_mp     = Kokkos::create_mirror_view(lh_precip_mp    );
    auto h_rhophi_precip_mp = Kokkos::create_mirror_view(rhophi_precip_mp);
    auto h_rhokin_precip_mp = Kokkos::create_mirror_view(rhokin_precip_mp);

    auto h_frhoge_af  = Kokkos::create_mirror_view(frhoge_af );
    auto h_frhogqv_af = Kokkos::create_mirror_view(frhogqv_af);
    auto h_frhoge_rad = Kokkos::create_mirror_view(frhoge_rad);
    auto h_qke        = Kokkos::create_mirror_view(qke       );
    auto h_gsgam2     = Kokkos::create_mirror_view(gsgam2    );
    auto h_gsgam2h    = Kokkos::create_mirror_view(gsgam2h   );
    auto h_gam2       = Kokkos::create_mirror_view(gam2      );
    auto h_gam2h      = Kokkos::create_mirror_view(gam2h     );

    auto h_ix = Kokkos::create_mirror_view(ix);
    auto h_iy = Kokkos::create_mirror_view(iy);
    auto h_iz = Kokkos::create_mirror_view(iz);
    auto h_jx = Kokkos::create_mirror_view(jx);
    auto h_jy = Kokkos::create_mirror_view(jy);
    auto h_jz = Kokkos::create_mirror_view(jz);

    auto h_z           = Kokkos::create_mirror_view(z          );
    auto h_zh          = Kokkos::create_mirror_view(zh         );
    auto h_GPREC       = Kokkos::create_mirror_view(GPREC      );
    auto h_CBMFX       = Kokkos::create_mirror_view(CBMFX      );
    auto h_qd          = Kokkos::create_mirror_view(qd         );
    auto h_rceff       = Kokkos::create_mirror_view(rceff      );
    auto h_rceff_solid = Kokkos::create_mirror_view(rceff_solid);
    auto h_rceff_cld   = Kokkos::create_mirror_view(rceff_cld  );

    auto h_rctop  = Kokkos::create_mirror_view(rctop );
    auto h_rwtop  = Kokkos::create_mirror_view(rwtop );
    auto h_tctop  = Kokkos::create_mirror_view(tctop );

    auto h_GDCLW  = Kokkos::create_mirror_view(GDCLW );
    auto h_GDCFRC = Kokkos::create_mirror_view(GDCFRC);

    read_data_3d("data/rhog.dat",   h_rhog);
    read_data_3d("data/rhogvx.dat", h_rhogvx);
    read_data_3d("data/rhogvy.dat", h_rhogvy);
    read_data_3d("data/rhogvz.dat", h_rhogvz);
    read_data_3d("data/rhogw.dat",  h_rhogw);
    read_data_3d("data/rhoge.dat",  h_rhoge);
    read_data_3d("data/vx.dat",     h_vx);
    read_data_3d("data/vy.dat",     h_vy);
    read_data_3d("data/vz.dat",     h_vz);
    read_data_3d("data/w.dat",      h_w);
    read_data_3d("data/unccn.dat",  h_unccn);
    read_data_3d("data/rho.dat",    h_rho);
    read_data_3d("data/pre.dat",    h_pre);
    read_data_3d("data/tem.dat",    h_tem);

    read_data_4d("data/rhogq_Lswp.dat", h_rhogq_Lswp);
    read_data_4d("data/q_Lswp.dat",     h_q_Lswp);

    read_data_4d("data/precip_mp.dat",  h_precip_mp);
    read_data_4d("data/precip1_mp.dat", h_precip1_mp);
    read_data_4d("data/precip2_mp.dat", h_precip2_mp);

    read_data_3d("data/rhoein_precip_mp.dat", h_rhoein_precip_mp);
    read_data_3d("data/lh_precip_mp.dat",     h_lh_precip_mp);
    read_data_3d("data/rhophi_precip_mp.dat", h_rhophi_precip_mp);
    read_data_3d("data/rhokin_precip_mp.dat", h_rhokin_precip_mp);

    read_data_3d("data/frhoge_af.dat",  h_frhoge_af);
    read_data_3d("data/frhogqv_af.dat", h_frhogqv_af);
    read_data_3d("data/frhoge_rad.dat", h_frhoge_rad);
    read_data_3d("data/qke.dat",     h_qke);
    read_data_3d("data/gsgam2.dat",  h_gsgam2);
    read_data_3d("data/gsgam2h.dat", h_gsgam2h);
    read_data_3d("data/gam2.dat",    h_gam2);
    read_data_3d("data/gam2h.dat",   h_gam2h);

    read_data_2d("data/ix.dat", h_ix);
    read_data_2d("data/iy.dat", h_iy);
    read_data_2d("data/iz.dat", h_iz);
    read_data_2d("data/jx.dat", h_jx);
    read_data_2d("data/jy.dat", h_jy);
    read_data_2d("data/jz.dat", h_jz);

    read_data_3d("data/z.dat",  h_z);
    read_data_3d("data/zh.dat", h_zh);
    read_data_3d("data/GPREC.dat", h_GPREC);
    read_data_3d("data/CBMFX.dat", h_CBMFX);

    // copy data from host to device
    Kokkos::fence();
    Kokkos::deep_copy(rhog  , h_rhog   );
    Kokkos::deep_copy(rhogvx, h_rhogvx );
    Kokkos::deep_copy(rhogvy, h_rhogvy );
    Kokkos::deep_copy(rhogvz, h_rhogvz );
    Kokkos::deep_copy(rhogw , h_rhogw  );
    Kokkos::deep_copy(rhoge , h_rhoge  );
    Kokkos::deep_copy(vx    , h_vx     );
    Kokkos::deep_copy(vy    , h_vy     );
    Kokkos::deep_copy(vz    , h_vz     );
    Kokkos::deep_copy(w     , h_w      );
    Kokkos::deep_copy(unccn , h_unccn  );
    Kokkos::deep_copy(rho   , h_rho    );
    Kokkos::deep_copy(pre   , h_pre    );
    Kokkos::deep_copy(tem   , h_tem    );

    Kokkos::deep_copy(rhogq_Lswp, h_rhogq_Lswp);
    Kokkos::deep_copy(q_Lswp    , h_q_Lswp    );

    Kokkos::deep_copy(precip_mp , h_precip_mp );
    Kokkos::deep_copy(precip1_mp, h_precip1_mp);
    Kokkos::deep_copy(precip2_mp, h_precip2_mp);

    Kokkos::deep_copy(rhoein_precip_mp, h_rhoein_precip_mp );
    Kokkos::deep_copy(lh_precip_mp    , h_lh_precip_mp     );
    Kokkos::deep_copy(rhophi_precip_mp, h_rhophi_precip_mp );
    Kokkos::deep_copy(rhokin_precip_mp, h_rhokin_precip_mp );

    Kokkos::deep_copy(frhoge_af , h_frhoge_af );
    Kokkos::deep_copy(frhogqv_af, h_frhogqv_af);
    Kokkos::deep_copy(frhoge_rad, h_frhoge_rad);
    Kokkos::deep_copy(qke       , h_qke       );
    Kokkos::deep_copy(gsgam2    , h_gsgam2    );
    Kokkos::deep_copy(gsgam2h   , h_gsgam2h   );
    Kokkos::deep_copy(gam2      , h_gam2      );
    Kokkos::deep_copy(gam2h     , h_gam2h     );

    Kokkos::deep_copy(ix, h_ix);
    Kokkos::deep_copy(iy, h_iy);
    Kokkos::deep_copy(iz, h_iz);
    Kokkos::deep_copy(jx, h_jx);
    Kokkos::deep_copy(jy, h_jy);
    Kokkos::deep_copy(jz, h_jz);

    Kokkos::deep_copy(z          , h_z          );
    Kokkos::deep_copy(zh         , h_zh         );
    Kokkos::deep_copy(GPREC      , h_GPREC      );
    Kokkos::deep_copy(CBMFX      , h_CBMFX      );
    Kokkos::deep_copy(qd         , h_qd         );
    Kokkos::deep_copy(rceff      , h_rceff      );
    Kokkos::deep_copy(rceff_solid, h_rceff_solid);
    Kokkos::deep_copy(rceff_cld  , h_rceff_cld  );

    Kokkos::deep_copy(rctop, h_rctop);
    Kokkos::deep_copy(rwtop, h_rwtop);
    Kokkos::deep_copy(tctop, h_tctop);

    Kokkos::deep_copy(GDCLW , h_GDCLW );
    Kokkos::deep_copy(GDCFRC, h_GDCFRC);
    Kokkos::fence();
    /**
     * Vertical grid setup
     */
    GRD_Setup();
    CVW_CPW_Setup();
    /**
     * Saturation set_up
     */
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
    //  * create sub views:
    //  */
    // auto sub_rhog   = subview(rhog  , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rhogvx = subview(rhogvx, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rhogvy = subview(rhogvy, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rhogvz = subview(rhogvz, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rhogw  = subview(rhogw , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rhoge  = subview(rhoge , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_vx     = subview(vx    , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_vy     = subview(vy    , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_vz     = subview(vz    , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_w      = subview(w     , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_unccn  = subview(unccn , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rho    = subview(rho   , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_pre    = subview(pre   , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_tem    = subview(tem   , 0, Kokkos::ALL(), Kokkos::ALL());

    // auto sub_rhogq_Lswp = subview(rhogq_Lswp, 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    // auto sub_q_Lswp     = subview(q_Lswp    , 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    // auto sub_precip_mp  = subview(precip_mp , Kokkos::ALL(), 0, 0, Kokkos::ALL());
    // auto sub_precip1_mp = subview(precip1_mp, Kokkos::ALL(), 0, 0, Kokkos::ALL());
    // auto sub_precip2_mp = subview(precip2_mp, Kokkos::ALL(), 0, 0, Kokkos::ALL());

    // auto sub_rhoein_precip_mp = subview(rhoein_precip_mp, 0, 0, Kokkos::ALL());
    // auto sub_lh_precip_mp     = subview(lh_precip_mp    , 0, 0, Kokkos::ALL());
    // auto sub_rhophi_precip_mp = subview(rhophi_precip_mp, 0, 0, Kokkos::ALL());
    // auto sub_rhokin_precip_mp = subview(rhokin_precip_mp, 0, 0, Kokkos::ALL());

    // auto sub_frhoge_af  = subview(frhoge_af , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_frhogqv_af = subview(frhogqv_af, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_frhoge_rad = subview(frhoge_rad, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_qke        = subview(qke       , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_gsgam2     = subview(gsgam2    , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_gsgam2h    = subview(gsgam2h   , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_gam2       = subview(gam2      , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_gam2h      = subview(gam2h     , 0, Kokkos::ALL(), Kokkos::ALL());

    // auto sub_ix = subview(ix, 0, Kokkos::ALL());
    // auto sub_iy = subview(iy, 0, Kokkos::ALL());
    // auto sub_iz = subview(iz, 0, Kokkos::ALL());
    // auto sub_jx = subview(jx, 0, Kokkos::ALL());
    // auto sub_jy = subview(jy, 0, Kokkos::ALL());
    // auto sub_jz = subview(jz, 0, Kokkos::ALL());

    // auto sub_z           = subview(z           , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_zh          = subview(zh          , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_GPREC       = subview(GPREC       , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_CBMFX       = subview(CBMFX       , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_qd          = subview(qd          , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rceff       = subview(rceff       , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rceff_solid = subview(rceff_solid , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rceff_cld   = subview(rceff_cld   , 0, Kokkos::ALL(), Kokkos::ALL());

    // auto sub_rctop = subview(rctop, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_rwtop = subview(rwtop, 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_tctop = subview(tctop, 0, Kokkos::ALL(), Kokkos::ALL());

    // auto sub_GDCLW  = subview(GDCLW , 0, Kokkos::ALL(), Kokkos::ALL());
    // auto sub_GDCFRC = subview(GDCFRC, 0, Kokkos::ALL(), Kokkos::ALL());

    // for(int i = 0; i < SET_iteration; i++)
    // {
    //     // Call mp_driver
    //     mp_driver(
    //               l,
    //               sub_rhog             ,
    //               sub_rhogvx           ,
    //               sub_rhogvy           ,
    //               sub_rhogvz           ,
    //               sub_rhogw            ,
    //               sub_rhoge            ,
    //               sub_rhogq_Lswp       ,
    //               sub_vx               ,
    //               sub_vy               ,
    //               sub_vz               ,
    //               sub_w                ,
    //               sub_unccn            ,
    //               sub_rho              ,
    //               sub_tem              ,
    //               sub_pre              ,
    //               sub_q_Lswp           ,
    //               sub_qd               ,
    //               sub_precip_mp        ,  
    //               sub_precip1_mp       ,
    //               sub_precip2_mp       ,
    //               sub_rhoein_precip_mp ,
    //               sub_lh_precip_mp     ,
    //               sub_rhophi_precip_mp ,
    //               sub_rhokin_precip_mp ,
    //               sub_rceff            ,
    //               sub_rceff_solid      ,
    //               sub_rceff_cld        ,
    //               sub_rctop            ,
    //               sub_rwtop            ,
    //               sub_tctop            ,
    //               sub_frhoge_af        ,
    //               sub_frhogqv_af       ,
    //               sub_frhoge_rad       ,
    //               sub_qke              ,
    //               sub_gsgam2           ,
    //               sub_gsgam2h          ,
    //               sub_gam2             ,
    //               sub_gam2h            ,
    //               sub_ix               ,
    //               sub_iy               ,
    //               sub_iz               ,
    //               sub_jx               ,
    //               sub_jy               ,
    //               sub_jz               ,
    //               sub_z                ,
    //               sub_zh               ,
    //               TIME_DTL,
    //               TIME_CTIME,
    //               sub_GDCLW            ,
    //               sub_GDCFRC           ,
    //               sub_GPREC            ,
    //               sub_CBMFX           
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
    //     PROF_val_check("QV", rhogq_Lswp, I_QV, CHECK_QV1);
    //     PROF_val_check("QC", rhogq_Lswp, I_QC, CHECK_QC2);
    //     PROF_val_check("QR", rhogq_Lswp, I_QR, CHECK_QR3);
    //     PROF_val_check("QI", rhogq_Lswp, I_QI, CHECK_QI4);
    //     PROF_val_check("QS", rhogq_Lswp, I_QS, CHECK_QS5);
    //     PROF_val_check("QG", rhogq_Lswp, I_QG, CHECK_QG6);
    // }

    std::cout << "============= All process finished =============== \n";

    // try deallocate GRD views
    GRD_gz   = View1D<double, Kokkos::CudaSpace>();
    GRD_gzh  = View1D<double, Kokkos::CudaSpace>();
    GRD_dgz  = View1D<double, Kokkos::CudaSpace>();
    GRD_dgzh = View1D<double, Kokkos::CudaSpace>();
    GRD_rdgz = View1D<double, Kokkos::CudaSpace>();
    GRD_rdgzh= View1D<double, Kokkos::CudaSpace>();
    GRD_afact= View1D<double, Kokkos::CudaSpace>();
    GRD_bfact= View1D<double, Kokkos::CudaSpace>();
    GRD_cfact= View1D<double, Kokkos::CudaSpace>();
    GRD_dfact= View1D<double, Kokkos::CudaSpace>();

    CVW= View1D<double, Kokkos::CudaSpace>();
    CPW= View1D<double, Kokkos::CudaSpace>();

}
Kokkos::finalize();
}
