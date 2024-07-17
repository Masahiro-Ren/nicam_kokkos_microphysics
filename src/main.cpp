#include "problemsize.h"
#include "data_io.h"
#include "mod_mp_driver.h"

using namespace PROBLEM_SIZE;

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