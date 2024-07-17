#pragma once

#include "problemsize.h"
#include "mod_debug.h"
#include "mod_mp_nsw6.h"

namespace MOD_MP_DRIVER {
    
bool opt_radisu_explicit = false;
bool opt_volume_explicit = false;

double TSICE = 273.15;
double TWICE = 258.15;

void mp_init(std::string& MP_TYPE);
void mp_driver();

};
