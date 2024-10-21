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

// declare all variables
double rhog  [ADM_lall][ADM_kall][ADM_gall_in];
double rhogvx[ADM_lall][ADM_kall][ADM_gall_in];
double rhogvy[ADM_lall][ADM_kall][ADM_gall_in];
double rhogvz[ADM_lall][ADM_kall][ADM_gall_in];
double rhogw [ADM_lall][ADM_kall][ADM_gall_in];
double rhoge [ADM_lall][ADM_kall][ADM_gall_in];
double vx    [ADM_lall][ADM_kall][ADM_gall_in];
double vy    [ADM_lall][ADM_kall][ADM_gall_in];
double vz    [ADM_lall][ADM_kall][ADM_gall_in];
double w     [ADM_lall][ADM_kall][ADM_gall_in];
double unccn [ADM_lall][ADM_kall][ADM_gall_in];
double rho   [ADM_lall][ADM_kall][ADM_gall_in];
double pre   [ADM_lall][ADM_kall][ADM_gall_in];
double tem   [ADM_lall][ADM_kall][ADM_gall_in];

double rhogq_Lswp [ADM_lall][TRC_VMAX][ADM_kall][ADM_gall_in];
double q_Lswp     [ADM_lall][TRC_VMAX][ADM_kall][ADM_gall_in];

double precip_mp  [2][ADM_lall][ADM_KNONE][ADM_gall_in];
double precip1_mp [2][ADM_lall][ADM_KNONE][ADM_gall_in];
double precip2_mp [2][ADM_lall][ADM_KNONE][ADM_gall_in];

double rhoein_precip_mp [ADM_lall][ADM_KNONE][ADM_gall_in];
double lh_precip_mp     [ADM_lall][ADM_KNONE][ADM_gall_in];
double rhophi_precip_mp [ADM_lall][ADM_KNONE][ADM_gall_in];
double rhokin_precip_mp [ADM_lall][ADM_KNONE][ADM_gall_in];

double frhoge_af [ADM_lall][ADM_kall][ADM_gall_in];
double frhogqv_af[ADM_lall][ADM_kall][ADM_gall_in];
double frhoge_rad[ADM_lall][ADM_kall][ADM_gall_in];
double qke       [ADM_lall][ADM_kall][ADM_gall_in];
double gsgam2    [ADM_lall][ADM_kall][ADM_gall_in];
double gsgam2h   [ADM_lall][ADM_kall][ADM_gall_in];
double gam2      [ADM_lall][ADM_kall][ADM_gall_in];
double gam2h     [ADM_lall][ADM_kall][ADM_gall_in];

double ix [ADM_lall][ADM_gall_in];
double iy [ADM_lall][ADM_gall_in];
double iz [ADM_lall][ADM_gall_in];
double jx [ADM_lall][ADM_gall_in];
double jy [ADM_lall][ADM_gall_in];
double jz [ADM_lall][ADM_gall_in];

double z           [ADM_lall][ADM_kall][ADM_gall_in];
double zh          [ADM_lall][ADM_kall][ADM_gall_in];
double GPREC       [ADM_lall][ADM_kall][ADM_gall_in];
double CBMFX       [ADM_lall][ADM_kall][ADM_gall_in];
double qd          [ADM_lall][ADM_kall][ADM_gall_in];
double rceff       [ADM_lall][ADM_kall][ADM_gall_in];
double rceff_solid [ADM_lall][ADM_kall][ADM_gall_in];
double rceff_cld   [ADM_lall][ADM_kall][ADM_gall_in];

double rctop [ADM_lall][ADM_KNONE][ADM_gall_in];
double rwtop [ADM_lall][ADM_KNONE][ADM_gall_in];
double tctop [ADM_lall][ADM_KNONE][ADM_gall_in];

double GDCLW  [ADM_lall][ADM_kall][ADM_gall_in];
double GDCFRC [ADM_lall][ADM_kall][ADM_gall_in];


int main(int argc, char* argv[])
{

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
    /**
     * Saturation set_up
     */
    SATURATION_Setup();
    /**
     * microphysics initialization
     */
    mp_init(MP_TYPE);

    int l = SET_l;

    std::cout << "============= Finish Initialize =============== \n";

    /**
     * Start Simulation
    */
    std::cout << "============= Start Kernel =============== \n";

    /**
     * create temp arrays to stroe :
     *      precip_mp(:,ADM_KNONE,l,:)
     *      precip1_mp(:,ADM_KNONE,l,:)
     *      precip2_mp(:,ADM_KNONE,l,:)
     */
    double precip_mp_tmp [2][ADM_gall_in];
    double precip1_mp_tmp[2][ADM_gall_in];
    double precip2_mp_tmp[2][ADM_gall_in];

    for(int k = 0; k < 2; k++)
    {
        for(int ij = 0; ij < ADM_gall_in; ij++)
        {
            precip_mp_tmp [k][ij] = precip_mp [k][0][0][ij];
            precip1_mp_tmp[k][ij] = precip1_mp[k][0][0][ij];
            precip2_mp_tmp[k][ij] = precip2_mp[k][0][0][ij];
        }
    }
    // Copy end

    for(int i = 0; i < SET_iteration; i++)
    {
        // Call mp_driver
        mp_driver(
                  l,
                  rhog            [0],
                  rhogvx          [0],
                  rhogvy          [0],
                  rhogvz          [0],
                  rhogw           [0],
                  rhoge           [0],
                  rhogq_Lswp      [0],
                  vx              [0],
                  vy              [0],
                  vz              [0],
                  w               [0],
                  unccn           [0],
                  rho             [0],
                  tem             [0],
                  pre             [0],
                  q_Lswp          [0],
                  qd              [0],
                  precip_mp_tmp,
                  precip1_mp_tmp,
                  precip2_mp_tmp,
                  rhoein_precip_mp[0][0],
                  lh_precip_mp    [0][0],
                  rhophi_precip_mp[0][0],
                  rhokin_precip_mp[0][0],
                  rceff           [0],
                  rceff_solid     [0],
                  rceff_cld       [0],
                  rctop           [0],
                  rwtop           [0],
                  tctop           [0],
                  frhoge_af       [0],
                  frhogqv_af      [0],
                  frhoge_rad      [0],
                  qke             [0],
                  gsgam2          [0],
                  gsgam2h         [0],
                  gam2            [0],
                  gam2h           [0],
                  ix              [0],
                  iy              [0],
                  iz              [0],
                  jx              [0],
                  jy              [0],
                  jz              [0],
                  z               [0],
                  zh              [0],
                  TIME_DTL,
                  TIME_CTIME,
                  GDCLW           [0],
                  GDCFRC          [0],
                  GPREC           [0],
                  CBMFX           [0]
                  );
    }
    std::cout << "============= Finish Kernel =============== \n";

    /**
     * Checking Results
    */

    return 0;
}
