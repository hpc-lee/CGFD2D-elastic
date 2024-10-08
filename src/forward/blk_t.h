#ifndef BLK_T_H
#define BLK_T_H

#include "constants.h"
#include "fd_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
#include "io_funcs.h"

/*******************************************************************************
 * structure
 ******************************************************************************/

typedef struct
{
  // name for output file name
  char name[CONST_MAX_STRLEN];

  //// flag of medium
  //int medium_type;

  // fd
  fd_t    *fd;     // collocated grid fd
  
  // coordnate: x2d, z2d
  gd_t *gd;

  // grid metrics: jac, xi_x, etc
  gd_metric_t *gd_metric;

  // media: rho, lambda, mu etc
  md_t *md;

  // wavefield:
  wav_t *wav;
  
  // source term
  src_t *src;
  
  // bdry 
  bdry_t *bdry;
  
  // io
  iorecv_t  *iorecv;
  ioline_t  *ioline;
  iosnap_t  *iosnap;

  // fname and dir
  char output_fname_part[CONST_MAX_STRLEN];
  // wavefield output
  char output_dir[CONST_MAX_STRLEN];
  // seperate grid output to save grid for repeat simulation
  char grid_export_dir[CONST_MAX_STRLEN];
  // seperate medium output to save medium for repeat simulation
  char media_export_dir[CONST_MAX_STRLEN];

  // mem usage
  size_t number_of_float;
  size_t number_of_btye;
} blk_t;

/*******************************************************************************
 * function prototype
 ******************************************************************************/

int
blk_init(blk_t *blk, const int verbose);

// set str
int
blk_set_output(blk_t *blk,
               char *output_dir,
               char *grid_export_dir,
               char *media_export_dir,
               const int verbose);

int
blk_print(blk_t *blk);

int
blk_dt_esti_curv(gd_t *gd, md_t *md,
                 float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
                 int *dtmaxi, int *dtmaxk);

int
blk_dt_esti_cart(gd_t *gd, md_t *md,
                 float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
                 int *dtmaxi, int *dtmaxk);

float
blk_keep_three_digi(float dt);

#endif
