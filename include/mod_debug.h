#pragma once

#include "problemsize.h"
#include "data_io.h"

using namespace PROBLEM_SIZE;

namespace DEBUG {

extern int EX_STEP;
extern int EX_rgnid;

extern int sat_ite_sum;
extern int sat_ite_count;
extern int sat_ite_max;
extern int sat_ite_min;


// Host code
void ADM_Proc_stop();

void PROF_val_check(const std::string& val_name, double arr2d[ADM_kall][ADM_gall_in], double CHECK_arr2d[ADM_kall][ADM_gall_in]);

void PROF_val_check(const std::string& val_name, double arr3d[ADM_lall][ADM_kall][ADM_gall_in], double CHECK_arr3d[ADM_lall][ADM_kall][ADM_gall_in]);

void CVW_CPW_Setup();

void GRD_Setup();

void GRD_Input_vgrid();

void cnvvar_rhogkin_in(
    double rhog     [kdim][ijdim],
    double rhogvx   [kdim][ijdim],
    double rhogvy   [kdim][ijdim],
    double rhogvz   [kdim][ijdim],
    double rhogw    [kdim][ijdim],
    double C2Wfact  [2][kdim][ijdim],
    double W2Cfact  [2][kdim][ijdim],
    double rhogkin  [kdim][ijdim],
    double rhogkin_h[kdim][ijdim],
    double rhogkin_v[kdim][ijdim] );

/**
 * Kokkos ver.
 */
void cnvvar_rhogkin_in(
    const View2D<double, DEFAULT_MEM>&  rhog     ,
    const View2D<double, DEFAULT_MEM>&  rhogvx   ,
    const View2D<double, DEFAULT_MEM>&  rhogvy   ,
    const View2D<double, DEFAULT_MEM>&  rhogvz   ,
    const View2D<double, DEFAULT_MEM>&  rhogw    ,
    const View3D<double, DEFAULT_MEM>&  C2Wfact  ,
    const View3D<double, DEFAULT_MEM>&  W2Cfact  ,
    const View2D<double, DEFAULT_MEM>&  rhogkin  ,
    const View2D<double, DEFAULT_MEM>&  rhogkin_h,
    const View2D<double, DEFAULT_MEM>&  rhogkin_v );

double MISC_gammafunc(double xx);

void PROF_val_check(const std::string& val_name, const View4D<double, DEFAULT_MEM>& arr4d, const size_t idx_arr2d, const View3D<double, HOST_MEM>& CHECK_arr3d);

void PROF_val_check(const std::string& val_name, const View3D<double, DEFAULT_MEM>& arr3d, const View3D<double, HOST_MEM>& CHECK_arr3d);
}
