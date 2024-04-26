#ifndef MOD_MP_DRIVER
#define MOD_MP_DRIVER

#include "problemsize.h"
#include "mod_debug.h"
#include "mod_mp_nsw6.h"

class MOD_MP_Driver {
private:
    bool opt_radisu_explicit;
    bool opt_volume_explicit;

    double TSICE;
    double TWICE;
public:
    MOD_MP_Driver()
    {
        opt_radisu_explicit = false;
        opt_volume_explicit = false;

        TSICE = 273.15;
        TWICE = 258.15;
    }

    void mp_init(std::string& MP_TYPE);
    void mp_driver();
};


#endif