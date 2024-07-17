#include "mod_mp_driver.h"

namespace MOD_MP_DIRVER{

void mp_init(std::string& MP_TYPE)
{
    if(MP_TYPE == "NSW6")
    {
        std::cout << "************ microphysics type = NSW6 *************** \n";
        // call mp_nsw6_init;
    }
}

void mp_driver()
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

};
