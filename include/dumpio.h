#ifndef DUMPIO_H
#define DUMPIO_H
/******************************************************************//**
 *functions
 *<utilities>
 *void           dumpio_ednchg                   : endian changer
 *void           dumpio_mk_fname                 : filename generator
 *void           dumpio_set_str, dumpio_trim_str    : string preprocess
 *int32_t        dumpio_syscheck                 : check endian
 *<file r/w>
 *int32_t        dumpio_fopen                    : open file IO stream
 *int32_t        dumpio_fclose                   : close file IO stream
 *int32_t        dumpio_write_data               : write data array
 *int32_t        dumpio_read_data                : read data array
 **********************************************************************/
#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>

#define H_SHORT  16
#define H_MID    64
#define H_LONG  256

/* action type */
#define H_FREAD   0
#define H_FWRITE  1
#define H_FAPPEND 2

/* return type */
#define ERROR_CODE   -1
#define SUCCESS_CODE  1

typedef float (real32_t);
typedef double(real64_t);

/* file structure */
typedef struct{
  char fname[H_LONG];
  int32_t opened;
  FILE *fp;
} fileinfo_t;

extern fileinfo_t *dumpio_finfo;

/******************************************************************************/
/* C functions                                                                */
/******************************************************************************/

/** endian change *****************************************************/
extern void dumpio_ednchg( void* avp_pointer,
                           const int32_t ai_size,
                           const int32_t ai_num   );

/** filename generator ************************************************/
extern void dumpio_mk_fname( char *fname,
                             const char *base,
                             const char *ext,
                             int32_t i,
                             int32_t y );

/** string preprocess *************************************************/
extern void dumpio_set_str( char *_str,
                            const char *str,
                            const int str_len );

/** string postprocess *************************************************/
extern void dumpio_trim_str( char *_str,
                             const char *str,
                             const int str_len );

/** check system & initialze ******************************************/
extern int32_t dumpio_syscheck( void );

/** open file IO stream ***********************************************/
extern int32_t dumpio_fopen( const char *vname_in, int32_t mode );

/** close file IO stream **********************************************/
extern int32_t dumpio_fclose( int32_t fid );

/** write data array **************************************************/
extern int32_t dumpio_write_data( int32_t fid,
                                  int32_t idxsize,
                                  void *data );

/** read data array (full size) ***************************************/
extern int32_t dumpio_read_data( int32_t fid,
                                 int32_t idxsize,
                                 void *data );

/* file ID counter */
int32_t dumpio_num_of_file = 0;

/* package+data+status container */
fileinfo_t *dumpio_finfo = NULL;

/* system information */
int32_t dumpio_system_ednchg = 0;

#endif