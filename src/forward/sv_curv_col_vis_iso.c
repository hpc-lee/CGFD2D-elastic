/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_curv_col_el.h"
#include "sv_curv_col_el_iso.h"
#include "sv_curv_col_vis_iso.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

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
                  const int verbose)
{
  // local pointer get each vars
  float *restrict Vx    = w_cur + wav->Vx_pos ;
  float *restrict Vz    = w_cur + wav->Vz_pos ;
  float *restrict Txx   = w_cur + wav->Txx_pos;
  float *restrict Tzz   = w_cur + wav->Tzz_pos;
  float *restrict Txz   = w_cur + wav->Txz_pos;
  float *restrict hVx   = rhs   + wav->Vx_pos ; 
  float *restrict hVz   = rhs   + wav->Vz_pos ; 
  float *restrict hTxx  = rhs   + wav->Txx_pos; 
  float *restrict hTzz  = rhs   + wav->Tzz_pos; 
  float *restrict hTxz  = rhs   + wav->Txz_pos; 

  int nmaxwell = md->nmaxwell;

  float **Jxx  = (float **) malloc(nmaxwell * sizeof(float *)); 
  float **Jzz  = (float **) malloc(nmaxwell * sizeof(float *));
  float **Jxz  = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJxx = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJzz = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJxz = (float **) malloc(nmaxwell * sizeof(float *));

  for(int n=0; n < nmaxwell; n++)
  {
    Jxx[n]  = w_cur + wav->Jxx_pos[n];
    Jzz[n]  = w_cur + wav->Jzz_pos[n];
    Jxz[n]  = w_cur + wav->Jxz_pos[n];
    hJxx[n] = rhs   + wav->Jxx_pos[n];
    hJzz[n] = rhs   + wav->Jzz_pos[n];
    hJxz[n] = rhs   + wav->Jxz_pos[n];
  }

  float *restrict xi_x  = metric->xi_x;
  float *restrict xi_z  = metric->xi_z;
  float *restrict zt_x  = metric->zeta_x;
  float *restrict zt_z  = metric->zeta_z;
  float *restrict jac2d = metric->jac;

  float *restrict  lam2d = md->lambda;
  float *restrict  mu2d  = md->mu;
  float *restrict  slw2d = md->rho;
  float *restrict  wl    = md->wl;
  float **restrict Ylam  = md->Ylam;
  float **restrict Ymu   = md->Ymu;

  // grid size
  int ni1 = gd->ni1;
  int ni2 = gd->ni2;
  int nk1 = gd->nk1;
  int nk2 = gd->nk2;

  int ni  = gd->ni;
  int nk  = gd->nk;
  int nx  = gd->nx;
  int nz  = gd->nz;
  size_t siz_iz  = gd->siz_iz;

  float *vecVx2Vz = bdry->vecVx2Vz2;

  // local fd op
  int              fdx_inn_len;
  int    *restrict fdx_inn_indx;
  float  *restrict fdx_inn_coef;
  int              fdz_inn_len;
  int    *restrict fdz_inn_indx;
  float  *restrict fdz_inn_coef;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op
  fdx_inn_len  = fdx_op[num_of_fdx_op-1].total_len;
  fdx_inn_indx = fdx_op[num_of_fdx_op-1].indx;
  fdx_inn_coef = fdx_op[num_of_fdx_op-1].coef;

  fdz_inn_len  = fdz_op[num_of_fdz_op-1].total_len;
  fdz_inn_indx = fdz_op[num_of_fdz_op-1].indx;
  fdz_inn_coef = fdz_op[num_of_fdz_op-1].coef;

  // inner points
  sv_curv_col_el_iso_rhs_inner(Vx,Vz,Txx,Tzz,Txz,
                               hVx,hVz,hTxx,hTzz,hTxz,
                               xi_x, xi_z, zt_x, zt_z,
                               lam2d, mu2d, slw2d,
                               ni1,ni2,nk1,nk2,siz_iz,
                               fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                               fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                               verbose);

  // free, abs, source in turn
  // free surface at z2
  if (bdry->is_sides_free[CONST_NDIM-1][1] == 1)
  {
    // tractiong
    sv_curv_col_el_rhs_timg_z2(Txx,Tzz,Txz,hVx,hVz,
                               xi_x, xi_z, zt_x, zt_z,
                               jac2d, slw2d,
                               ni1,ni2,nk1,nk2,siz_iz,
                               fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                               fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                               verbose);

    // velocity: vlow
    sv_curv_col_el_iso_rhs_vlow_z2(Vx,Vz,hTxx,hTzz,hTxz,
                                   xi_x, xi_z, zt_x, zt_z,
                                   lam2d, mu2d, slw2d,
                                   vecVx2Vz,
                                   ni1,ni2,nk1,nk2,siz_iz,
                                   fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                   num_of_fdz_op,fdz_op,fdz_max_len,
                                   verbose);
  }

// cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_curv_col_el_iso_rhs_cfspml(Vx,Vz,Txx,Tzz,Txz,
                                  hVx,hVz,hTxx,hTzz,hTxz,
                                  xi_x, xi_z, zt_x, zt_z,
                                  lam2d, mu2d, slw2d,
                                  nk2, siz_iz,
                                  fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                  fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                  bdry,
                                  verbose);
    
  }


  sv_curv_col_vis_iso_atten(hTxx,hTzz,hTxz,
                            Jxx,Jzz,Jxz,
                            hJxx,hJzz,hJxz,
                            lam2d, mu2d, slw2d,
                            wl,Ylam,Ymu,
                            ni1,ni2,nk1,nk2,siz_iz,nmaxwell,
                            verbose);

  // add source term
  if (src->total_number > 0)
  {
    sv_curv_col_el_rhs_src(hVx,hVz,hTxx,hTzz,hTxz,
                           jac2d, slw2d, 
                           src,
                           verbose);
  }
  // end func

  return 0;
}

/*******************************************************************************
 * add the attenuation term
******************************************************************************/
int
sv_curv_col_vis_iso_atten(
           float *restrict hTxx, float *restrict hTzz, float *restrict hTxz,
           float **restrict  Jxx, float **restrict  Jzz, float **restrict  Jxz,
           float **restrict hJxx, float **restrict hJzz, float **restrict hJxz,
           float *restrict lam2d, float *restrict mu2d, float *restrict slw2d,
           float *restrict wl, float **restrict Ylam, float **restrict Ymu,
           int ni1, int ni2, int nk1, int nk2, size_t siz_iz, int nmaxwell,
           const int verbose)
{
  float lam,mu;
  float mem_Txx,mem_Tzz,mem_Txz;
  float sum_Jxxzz,sum_Jxx,sum_Jzz,sum_Jxz;
  float EVxx,EVzz,EVxz;
  float sum_hxz;
  size_t iptr;
  for (int k=nk1; k<=nk2; k++)
  {
    for (int i=ni1; i<=ni2; i++)
    {
      iptr = i + k*siz_iz;
      lam  = lam2d[iptr];
      mu   =  mu2d[iptr];

      sum_hxz = (hTxx[iptr]+hTzz[iptr])/(3*lam+2*mu);

      EVxx = ((2.0*hTxx[iptr]-hTzz[iptr])/(2*mu)+sum_hxz)/3;
      EVzz = ((2.0*hTzz[iptr]-hTxx[iptr])/(2*mu)+sum_hxz)/3;
      EVxz = hTxz[iptr]/mu*0.5;

      for (int n = 0; n < nmaxwell; n++)
      {
        hJxx[n][iptr] = wl[n]*(EVxx-Jxx[n][iptr]);
        hJzz[n][iptr] = wl[n]*(EVzz-Jzz[n][iptr]);
        hJxz[n][iptr] = wl[n]*(EVxz-Jxz[n][iptr]);
      }

      // sum of memory variable for attenuation
      sum_Jxx = 0.0;
      sum_Jzz = 0.0;
      sum_Jxz = 0.0;
      sum_Jxxzz = 0.0;
    
      for(int n=0; n<nmaxwell; n++)
      {
        sum_Jxxzz  += Ylam[n][iptr]  * (Jxx[n][iptr] + Jzz[n][iptr]);
        sum_Jxx  += Ymu[n][iptr]  * Jxx[n][iptr];
        sum_Jzz  += Ymu[n][iptr]  * Jzz[n][iptr];
        sum_Jxz  += Ymu[n][iptr]  * Jxz[n][iptr];
      }
    
      mem_Txx = lam*sum_Jxxzz + 2.0*mu*sum_Jxx;
      mem_Tzz = lam*sum_Jxxzz + 2.0*mu*sum_Jzz;
      mem_Txz = 2.0*mu*sum_Jxz;

      hTxx[iptr] -= mem_Txx;
      hTzz[iptr] -= mem_Tzz;
      hTxz[iptr] -= mem_Txz;
    }
  }
  return 0;
}


/*******************************************************************************
 * free surface coef
 ******************************************************************************/

int
sv_curv_col_vis_iso_dvh2dvz(gd_t        *gd,
                            gd_metric_t *metric,
                            md_t        *md,
                            bdry_t      *bdry,
                            const int verbose)
{
  int ierr = 0;

  int ni1 = gd->ni1;
  int ni2 = gd->ni2;
  int nk1 = gd->nk1;
  int nk2 = gd->nk2;
  int nx  = gd->nx;
  int nz  = gd->nz;
  size_t siz_iz  = gd->siz_iz;

  // point to each var
  float *restrict xi_x = metric->xi_x;
  float *restrict xi_z = metric->xi_z;
  float *restrict zt_x = metric->zeta_x;
  float *restrict zt_z = metric->zeta_z;

  float *restrict lam2d = md->lambda;
  float *restrict  mu2d = md->mu;

  float *vecVx2Vz = bdry->vecVx2Vz2;
  float *vecA = bdry->vecA;
  
  float A[2][2], B[2][2], C[2][2];
  float AB[2][2];

  for (int i = ni1; i <= ni2; i++)
  {
    size_t iptr = i + nk2 * siz_iz;

    float e11 = xi_x[iptr];
    float e12 = xi_z[iptr];
    float e21 = zt_x[iptr];
    float e22 = zt_z[iptr];

    float lam    = lam2d[iptr];
    float miu    =  mu2d[iptr];
    float lam2mu = lam + 2.0f * miu;
    
    // add some code






    // first dim: irow; sec dim: jcol, as Fortran code
    A[0][0]=lam2mu*e21*e21+miu*e22*e22;
    A[0][1]=lam*e21*e22+miu*e22*e21;
    A[1][0]=lam*e22*e21+miu*e21*e22;
    A[1][1]=lam2mu*e22*e22+miu*e21*e21;

    fdlib_math_invert2x2(A);

    B[0][0]=-lam2mu*e21*e11-miu*e22*e12;
    B[0][1]=-lam*e21*e12-miu*e22*e11;
    B[1][0]=-lam*e22*e11-miu*e21*e12;
    B[1][1]=-lam2mu*e22*e12-miu*e21*e11;

    fdlib_math_matmul2x2(A, B, AB);

    size_t ij = i * 4;

    // save into mat
    for(int irow = 0; irow < 2; irow++){
      for(int jcol = 0; jcol < 2; jcol++){
        vecVx2Vz[ij + irow*2 + jcol] = AB[irow][jcol];
        vecA[ij + irow*2 + jcol] = A[irow][jcol]; // A is DZ invert
      }
    }
  }

  return ierr;
}

int
sv_curv_col_vis_iso_free(float *restrict w_end,
                         wav_t  *wav,
                         gd_t   *gd,
                         gd_metric_t  *metric,
                         md_t *md,
                         bdry_t      *bdryfree, 
                         const int verbose)
{
  int ierr;

  return ierr;
}
