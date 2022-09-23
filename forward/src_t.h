#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

#include "constants.h"
#include "gd_info.h"
#include "gd_t.h"

// cal force_vec_stf/moment_ten_rate 1d index for icmp,it,istage
//  with respect to the start pointer of this source point
#define M_SRC_IND(icmp,it,istage,nt,num_stage) \
  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int total_number;
  int max_nt; // max nt of stf and mrf per src
  int max_stage; // max number of rk stages
  int max_ext; // max extened points

  // for getting value in calculation
  int it;
  int istage;

  // for output
  char evtnm[CONST_MAX_STRLEN];

  // time independent
  int *si; // local i index 
  int *sk; // local k index 
  int *it_begin; // start t index
  int *it_end;   // end   t index
  int   *ext_num; // valid extend points for this src
  int   *ext_indx; // max_ext * total_number
  float *ext_coef;

  // force and/or moment
  int force_actived;
  int moment_actived;

  // time dependent
  // force stf
  float *Fx; // max_stage * max_nt * total_number;
  float *Fz;
  // moment rate
  float *Mxx; // max_stage *max_nt * total_number;
  float *Mzz;
  float *Mxz;
} src_t;

/*************************************************
 * function prototype
 *************************************************/

int
src_glob_ext_ishere(int si, int sk, int half_ext, gdinfo_t *gdinfo);

int
src_set_by_par(gdinfo_t *gdinfo,
               gd_t *gdcurv,
               src_t    *src,
               float t0,
               float dt,
               int   max_stage,
               float *rk_stage_time,
               int   npoint_half_ext,
               char  *in_source_name,
               int   in_num_of_src,
               int   **source_index,
               float **source_inc,
               float **source_coords,
               float **force_vector, 
               int   *source_force_actived,
               float **moment_tensor,
               int   *source_moment_actived,
               char  **wavelet_name,
               float **wavelet_coefs,
               float *wavelet_tstart,
               float *wavelet_tend,
               int verbose);

float 
fun_ricker(float t, float fc, float t0);

float 
fun_ricker_deriv(float t, float fc, float t0);

float
fun_gauss(float t, float a, float t0);

float
fun_gauss_deriv(float t, float a, float t0);

void 
angle2moment(float strike, float dip, float rake, float* source_moment_tensor);

int
src_set_time(src_t *src, int it, int istage);

void
src_cal_norm_delt2d(float *delt, float x0, float z0,
                    float rx0, float rz0, int LenDelt);

#endif
