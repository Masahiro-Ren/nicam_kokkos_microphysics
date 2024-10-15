#include "problemsize.h"
#include "data_io.h"

using namespace PROBLEM_SIZE;
using namespace DATA_IO;

// declare all variables
double GRD_gz   [ADM_kall];
double GRD_gzh  [ADM_kall];
double GRD_dgz  [ADM_kall];
double GRD_dgzh [ADM_kall];
double GRD_rdgz [ADM_kall];
double GRD_rdgzh[ADM_kall];
double GRD_afact[ADM_kall];
double GRD_bfact[ADM_kall];
double GRD_cfact[ADM_kall];
double GRD_dfact[ADM_kall];

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

double rhogeq_Lswp[ADM_lall][TRC_VMAX][ADM_kall][ADM_gall_in];
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

    std::cout << "============= Finish Initialize =============== \n";

    /**
     * Start Simulation
    */
    std::cout << "============= Start Kernel =============== \n";

    std::cout << "============= Finish Kernel =============== \n";

    /**
     * Checking Results
    */

    return 0;
}