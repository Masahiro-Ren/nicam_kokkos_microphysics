#include "problemsize.h"
#include "data_io.h"
#include "mod_mp_driver.h"


Vec2d<double> rhog;


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
    std::cout << "ADM_gall_in_orig = " <<  PROBLEM_SIZE::ADM_gall_in_orig << std::endl;
    std::cout << "ADM_kall = " << PROBLEM_SIZE::ADM_kall << std::endl;
    std::cout << "SET_l = " << PROBLEM_SIZE::SET_l << std::endl;
    std::cout << "MP_TYPE = " << PROBLEM_SIZE::MP_TYPE << std::endl;
    std::cout << "PI = " << PROBLEM_SIZE::PI << std::endl;
    std::cout << "EPS = " << PROBLEM_SIZE::EPS << std::endl;
    std::cout << "============= Start Initialize =============== \n";

    Data_IO dumpio;
    MOD_MP_Driver driver;
    dumpio.read(RHOG, rhog);


    std::cout << "============= Finish Initialize =============== \n";

    /**
     * Start Simulation
    */
    std::cout << "============= Start Kernel =============== \n";

    driver.mp_init(PROBLEM_SIZE::MP_TYPE);

    driver.mp_driver();


    std::cout << "============= Finish Kernel =============== \n";

    /**
     * Checking Results
    */

    return 0;
}