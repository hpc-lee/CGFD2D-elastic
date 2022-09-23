#ifndef SV_EQ1ST_CURV_COL_EL_ANISO_H
#define SV_EQ1ST_CURV_COL_EL_ANISO_H

#include "fd_t.h"
#include "gd_info.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_free.h"
#include "bdry_pml.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

void
sv_eq1st_curv_col_el_aniso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  gdcurv_metric_t  *metric,
  md_t *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  // include different order/stentil
  int num_of_fdx_op, fd_op_t *fdx_op,
  int num_of_fdz_op, fd_op_t *fdz_op,
  int fdz_max_len, 
  const int verbose);

void
sv_eq1st_curv_col_el_aniso_rhs_inner(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict c11d, float *restrict c13d,
    float *restrict c15d, float *restrict c33d,
    float *restrict c35d, float *restrict c55d,
    float *restrict slw3d,
    int ni1, int ni2, int nk1, int nk2,
    size_t siz_iz,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int verbose);

void
sv_eq1st_curv_col_el_aniso_rhs_timg_z2(
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    int ni1, int ni2, int nk1, int nk2,
    size_t siz_iz,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int verbose);

void
sv_eq1st_curv_col_el_aniso_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict c11d, float *restrict c13d,
    float *restrict c15d, float *restrict c33d,
    float *restrict c35d, float *restrict c55d,
    float *restrict slw3d,
    float *restrict vecVx2Vz,
    int ni1, int ni2, int nk1, int nk2,
    size_t siz_iz,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
    const int verbose);

void
sv_eq1st_curv_col_el_aniso_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict c11d, float *restrict c13d,
    float *restrict c15d, float *restrict c33d,
    float *restrict c35d, float *restrict c55d,
    float *restrict slw3d,
    int nk2, size_t siz_iz,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    bdrypml_t *bdrypml, bdryfree_t *bdryfree,
    const int verbose);

int
sv_eq1st_curv_col_el_aniso_dvh2dvz(gdinfo_t        *gdinfo,
                                   gdcurv_metric_t *metric,
                                   md_t       *md,
                                   bdryfree_t      *bdryfree,
                                   const int verbose);

#endif
