/**------------------------------------------------------------------------------
 
  Debug utility module
 
  @par Description
          This module is for debug.
 
  @author  H.Tomita
 
  @par History
  @li      2012-06-29 (H.Yashiro)  [NEW]
  
  C++ version implemented by Ren (2023-11-24)
*/
#ifndef MOD_DEBUG_H
#define MOD_DEBUG_H

#include "problem_size.h"
#include <cstdlib>


//++ Public parameters & variables
int EX_STEP = 49;
int EX_rgnid;
int EX_fid;
int EX_err;
char* EX_fname;
// change
int sat_ite_sum = 0;
int sat_ite_count = 0;
int sat_ite_max = 0;
int sat_ite_min = -1;

// void PROF_setup();
// void PROF_setprefx();
// void PROF_rapstart();
// void PROF_rapend();
// void PROF_rapreport();
// 
// void PROF_valcheck();
// void PROF_valcheck_DP_3D_ToText();
// void PROF_valcheck_DP_3D_ToText_max(); 
// void PROF_valcheck_DP_4D_ToText();

/**
 * Value Check functions
*/
// void PROF_valcheck_SP_1D();
// void PROF_valcheck_SP_2D();
// void PROF_valcheck_SP_3D();
// void PROF_valcheck_SP_4D();
// void PROF_valcheck_SP_5D();
// void PROF_valcheck_SP_6D();
// void PROF_valcheck_DP_1D();
// void PROF_valcheck_DP_2D();
// void PROF_valcheck_DP_3D();
// void PROF_valcheck_DP_4D();
// void PROF_valcheck_DP_5D();
// void PROF_valcheck_DP_6D();

// // make file name with a number
// void MISC_make_idstr();
// // get an available file ID
// void IO_get_available_fid();
// // Gamma function
// void MISC_gammafunc();
// void MISC_gammafunc_h();
// void MISC_gammafunc_s();
// void MISC_gammafunc_d();

void ADM_proc_stop();
void ADM_MPItime();
void GRD_Setup();
void GRD_input_vgrid();

void cnvvar_rhogkin_in();
void cnvvar_rhogkin_in__FROM__precip_transport_new_3578();
void cnvvar_rhogkin_in__FROM__precip_transport_new_3951();

/* ============== public parameters & variables ================  */
int EX_STEP = 49;
int EX_rgnid;
int EX_fid;
int EX_err;
int EX_fname;
int sat_ite_sum = 0;
int sat_ite_count = 0;
int sat_ite_max = 0;
int sat_ite_min = -1;

#endif