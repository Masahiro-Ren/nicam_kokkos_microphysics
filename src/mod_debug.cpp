#include "mod_debug.h"

/* Some Private variables & parameters */
const int PROF_rapnlimit = 300;
int PROF_rapnmax = 0;
int PROF_grpnmax = 0;
int PROF_grpid[PROF_rapnlimit];
int PROF_rapnstr[PROF_rapnlimit];
int PROF_rapnend[PROF_rapnlimit];
int PROF_raplevel[PROF_rapnlimit];
char *PROF_grpname;
char *PROF_rapname;
char *PROF_prefix;
double PROF_raptstr[PROF_rapnlimit];
double PROF_rapttot[PROF_rapnlimit];

int PROF_default_rap_level = 2;
int PROF_rap_level = 2;
bool PROF_mpi_barrier = false;

char *PROF_header;
char *PROF_item;
double PROF_max;
double PROF_min;
double PROF_sum;

const int min_fid = 7;
const int max_fid = 99;
const bool NSTR_ZERO_START = true;
const int NSTR_MAX_DIGIT = 5;

void ADM_proc_stop()
{
    // Failure
    exit(1);
}

void ADM_MPItime()
{
    // to do
}

void GRD_Setup()
{

}

void get_rapid()
{

}