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

namespace mod_debug {

using namespace problem_size;

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

// ++ Private parameters & variables
// 
const int PROF_rapnlimit = 300;
int PROF_rapnmax   = 0;
int PROF_grpnmax   = 0;
int PROF_grpid   [PROF_rapnlimit];
int PROF_rapnstr [PROF_rapnlimit];
int PROF_rapnend [PROF_rapnlimit];
int PROF_raplevel[PROF_rapnlimit];
char* PROF_grpname  ;
char* PROF_rapname  ;
char* PROF_prefix   ;
double PROF_raptstr [PROF_rapnlimit];
double PROF_rapttot [PROF_rapnlimit];

int PROF_default_rap_level = 2;
int PROF_rap_level         = 2;
bool PROF_mpi_barrier       = false;

char* PROF_header;
char* PROF_item;
double PROF_max;
double PROF_min;
double PROF_sum;

const int min_fid = 7;
const int max_fid = 99;
const bool NSTR_ZERO_START = true;
const int NSTR_MAX_DIGIT  = 5;

// void PROF_setup();
// void PROF_setprefx();
// void PROF_rapstart();
// void PROF_rapend();
// void PROF_rapreport();
// 
// void PROF_valcheck();
// void PROF_valcheck_DP_3D_ToText();
// void PROF_valcheck_DP_3D_ToText_max(); // add yamanashi
// void PROF_valcheck_DP_4D_ToText();

void ADM_proc_stop();
void ADM_MPItime();
void GRD_Setup();
void GRD_input_vgrid();

void cnvvar_rhogkin_in();

}

#endif