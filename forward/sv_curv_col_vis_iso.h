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
sv_curv_col_vis_iso_onestage(
             float *restrict w_cur,
             float *restrict rhs, 
             wav_t  *wav,
             gd_t   *gd,
             gd_metric_t  *metric,
             md_t *md,
             bdry_t *bdry,
             src_t *src,
             // include different order/stentil
             int num_of_fdx_op, fd_op_t *fdx_op,
             int num_of_fdz_op, fd_op_t *fdz_op,
             int fdz_max_len, 
             const int verbose);

int
sv_curv_col_vis_iso_atten(
           float *restrict hTxx, float *restrict hTzz, float *restrict hTxz,
           float **restrict  Jxx, float **restrict  Jzz, float **restrict  Jxz,
           float **restrict hJxx, float **restrict hJzz, float **restrict hJxz,
           float *restrict lam2d, float *restrict mu2d, float *restrict slw2d,
           float *restrict wl, float **restrict Ylam, float **restrict Ymu,
           int ni1, int ni2, int nk1, int nk2, size_t siz_iz, int nmaxwell,
           const int verbose);

int
sv_curv_col_vis_iso_dvh2dvz(gd_t        *gd,
                            gd_metric_t *metric,
                            md_t        *md,
                            bdry_t      *bdry,
                            const int verbose);

int
sv_curv_col_vis_iso_free(float *restrict w_end,
                         wav_t  *wav,
                         gd_t   *gdinfo,
                         gd_metric_t  *metric,
                         md_t *md,
                         bdry_t      *bdry,
                         const int verbose);

#endif
