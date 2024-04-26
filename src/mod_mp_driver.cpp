#include "mod_mp_driver.h"

void MOD_MP_Driver::mp_init(std::string& MP_TYPE)
{
    if(MP_TYPE == "NSW6")
    {
        std::cout << "************ microphysics type = NSW6 *************** \n";
        // call mp_nsw6_init;
    }
}

void MOD_MP_Driver::mp_driver()
{

}
