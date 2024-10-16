#pragma once

#include "problemsize.h"
#include "data_io.h"

namespace DEBUG {

    int EX_STEP = 49;
    int EX_rgnid;

    int sat_ite_sum = 0;
    int sat_ite_count = 0;
    int sat_ite_max = 0;
    int sat_ite_min = -1;


    void ADM_Proc_stop();

    void PROF_val_check();

    void GRD_Setup();

    void GRD_Input_vgrid(const std::string& fname);

    void cnvvar_rhogkin_in();

    double MISC_gammafunc(double xx);

};