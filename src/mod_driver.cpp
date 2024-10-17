#include "mod_driver.h"

bool opt_radius_explicit = false;
bool opt_volume_explicit = false;

double TSICE = 273.15;
double TWICE = 258.15;

namespace DRIVER{

void mp_init(const std::string& MP_TYPE_in)
{
	if(MP_TYPE == "NSW6")
	{
		std::cout << "*** microphysics type: NSW6 *** \n";
		std::cout << "Call mp_nsw6_init \n";
	}
	else
	{
		std::cerr << __PRETTY_FUNCTION__ << " NOT appropriate type. Type:  " << MP_TYPE << std::endl;
		ADM_Proc_stop();
	}
}

}