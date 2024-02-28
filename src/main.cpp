#include "problemsize.h"
#include "mod_mp_driver.h"

int main(int argc, char* argv[])
{
    // int kernel_num = 0;
    // /**
    //  * Read arguments
    // */
    // if(argc == 2)
    // {
    //     kernel_num = atoi(argv[1]);
    // }
    // else if (argc > 2)
    // {
    //     std::cerr << "Too many arguments \n";
    //     exit(1);
    // }

    /**
     * Display Simulation Configuerations
    */
    std::cout << "[KERNEL] physicskernel_microphysics \n";
    std::cout << "ADM_gall_in_orig = " << ADM_gall_in_orig << std::endl;
    std::cout << "ADM_kall = " << ADM_kall << std::endl;
    std::cout << "SET_l = " << SET_l << std::endl;
    std::cout << "MP_TYPE = " << MP_TYPE << std::endl;
    std::cout << "PI = " << PI << std::endl;
    std::cout << "EPS = " << EPS << std::endl;
    std::cout << "============= Start Initialize =============== \n";

    std::cout << "============= Finish Initialize =============== \n";

    /**
     * Start Simulation
    */
    std::cout << "============= Start Kernel =============== \n";

    mp_init(MP_TYPE);

    mp_driver();


    std::cout << "============= Finish Kernel =============== \n";

    /**
     * Checking Results
    */

    return 0;
}