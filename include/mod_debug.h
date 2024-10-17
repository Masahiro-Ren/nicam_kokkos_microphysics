#pragma once

#include "problemsize.h"
#include "data_io.h"

namespace DEBUG {

    extern int EX_STEP;
    extern int EX_rgnid;

    extern int sat_ite_sum;
    extern int sat_ite_count;
    extern int sat_ite_max;
    extern int sat_ite_min;


    void ADM_Proc_stop();

    void PROF_val_check();

    void GRD_Setup();

    void GRD_Input_vgrid();

    void cnvvar_rhogkin_in();

    double MISC_gammafunc(double xx);

};