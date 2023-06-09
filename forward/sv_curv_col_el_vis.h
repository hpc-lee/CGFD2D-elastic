#ifndef SV_CURV_COL_EL_VIS_H
#define SV_CURV_COL_EL_VIS_H

#include "fd_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

int
sv_curv_col_el_vis_onestage(
             float *restrict w_cur,
             float *restrict rhs, 
             wav_t  *wav,
             gd_t   *gd,
             gdcurv_metric_t  *metric,
             md_t *md,
             bdry_t *bdry,
             src_t *src,
             // include different order/stentil
             int num_of_fdx_op, fd_op_t *fdx_op,
             int num_of_fdz_op, fd_op_t *fdz_op,
             int fdz_max_len, 
             const int verbose);

int
sv_curv_col_el_vis_rhs_vlow_z2(
             float *restrict  Vx , float *restrict  Vz ,
             float *restrict hTxx, float *restrict hTzz,
             float *restrict hTxz, 
             float **restrict  Jxx, float **restrict  Jzz,
             float **restrict  Jxz, 
             float **restrict hJxx, float **restrict hJzz,
             float **restrict hJxz, 
             float *restrict xi_x, float *restrict xi_z,
             float *restrict zt_x, float *restrict zt_z,
             float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
             float *restrict wl, float **restrict Ylam, float **restrict Ymu,  
             float *restrict vecVx2Vz,float *restrict vecA, 
             int ni1, int ni2, int nk1, int nk2, size_t siz_iz, int nmaxwell,
             int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
             int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
             const int verbose);

int
sv_curv_col_el_vis_dvh2dvz(gd_t        *gd,
                           gdcurv_metric_t *metric,
                           md_t       *md,
                           bdry_t      *bdry,
                           const int verbose);


int
sv_curv_col_el_vis_J(
                float *restrict hVx , float *restrict hVz ,
                float *restrict hTxx, float *restrict hTzz,
                float *restrict hTxz, 
                float **restrict  Jxx, float **restrict  Jzz,
                float **restrict  Jxz, 
                float **restrict hJxx, float **restrict hJzz,
                float **restrict hJxz, 
                float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
                float *restrict wl, float **restrict Ylam, float **restrict Ymu,  
                int ni1, int ni2, int nk1, int nk2, size_t siz_iz, int nmaxwell, 
                const int verbose);

#endif
