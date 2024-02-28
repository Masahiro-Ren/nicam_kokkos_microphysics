#include "mod_mp_driver.h"
#include "mod_mp_nsw6.h"

bool opt_radisu_explicit = false;
bool opt_volume_explicit = false;

double TSICE = 273.15;
double TWICE = 258.15;


void mp_init(std::string& MP_TYPE)
{
    if(MP_TYPE == "NSW6")
    {
        // call mp_nsw6_init
    }
}
