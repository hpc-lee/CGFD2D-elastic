/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_eq1st_curv_col_el_iso.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eq1st_curv_col_el_iso_onestage(
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

  float *restrict xi_x  = metric->xi_x;
  float *restrict xi_z  = metric->xi_z;
  float *restrict zt_x  = metric->zeta_x;
  float *restrict zt_z  = metric->zeta_z;
  float *restrict jac3d = metric->jac;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;
  float *restrict slw3d = md->rho;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;

  int ni  = gdinfo->ni;
  int nk  = gdinfo->nk;
  int nx  = gdinfo->nx;
  int nz  = gdinfo->nz;
  size_t siz_iz  = gdinfo->siz_iz;

  float *vecVx2Vz = bdryfree->vecVx2Vz2;

  // local fd op
  int              fdx_inn_len;
  int    *restrict fdx_inn_indx;
  float *restrict fdx_inn_coef;
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
  sv_eq1st_curv_col_el_iso_rhs_inner(Vx,Vz,Txx,Tzz,Txz,
                                    hVx,hVz,hTxx,hTzz,hTxz,
                                    xi_x, xi_z, zt_x, zt_z,
                                    lam3d, mu3d, slw3d,
                                    ni1,ni2,nk1,nk2,siz_iz,
                                    fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                    fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                    verbose);

  // free, abs, source in turn

  // free surface at z2
  if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1)
  {
    // tractiong
    sv_eq1st_curv_col_el_iso_rhs_timg_z2(Txx,Tzz,Txz,hVx,hVz,
                                        xi_x, xi_z, zt_x, zt_z,
                                        jac3d, slw3d,
                                        ni1,ni2,nk1,nk2,siz_iz,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        verbose);

    // velocity: vlow
    sv_eq1st_curv_col_el_iso_rhs_vlow_z2(Vx,Vz,hTxx,hTzz,hTxz,
                                        xi_x, xi_z, zt_x, zt_z,
                                        lam3d, mu3d, slw3d,
                                        vecVx2Vz,
                                        ni1,ni2,nk1,nk2,siz_iz,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        num_of_fdz_op,fdz_op,fdz_max_len,
                                        verbose);
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_curv_col_el_iso_rhs_cfspml(Vx,Vz,Txx,Tzz,Txz,
                                       hVx,hVz,hTxx,hTzz,hTxz,
                                       xi_x, xi_z, zt_x, zt_z,
                                       lam3d, mu3d, slw3d,
                                       nk2, siz_iz,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       bdrypml, bdryfree,
                                       verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_eq1st_curv_col_el_iso_rhs_src(hVx,hVz,hTxx,hTzz,hTxz,
                                    jac3d, slw3d, 
                                    src,
                                    verbose);
  }
  // end func

  return;
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

void
sv_eq1st_curv_col_el_iso_rhs_inner(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    int ni1, int ni2, int nk1, int nk2,
    size_t siz_iz, 
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTxx,DxTzz,DxTxz,DxVx,DxVz;
  float DzTxx,DzTzz,DzTxz,DzVx,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiz,ztx,ztz;

  float *restrict Vx_ptr;
  float *restrict Vz_ptr;
  float *restrict Txx_ptr;
  float *restrict Txz_ptr;
  float *restrict Tzz_ptr;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
        Vx_ptr = Vx + iptr;
        Vz_ptr = Vz + iptr;
        Txx_ptr = Txx + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;

        // Vx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxx, Txx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxx, Txx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tzz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTzz, Tzz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // metric
        xix = xi_x[iptr];
        xiz = xi_z[iptr];
        ztx = zt_x[iptr];
        ztz = zt_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

        // moment equation
        hVx[iptr] = slw*( xix*DxTxx + xiz*DxTxz  
                         +ztx*DzTxx + ztz*DzTxz );
        hVz[iptr] = slw*( xix*DxTxz + xiz*DxTzz 
                         +ztx*DzTxz + ztz*DzTzz );

        // Hooke's equatoin
        hTxx[iptr] =  lam2mu * ( xix*DxVx + ztx*DzVx)
                    + lam    * ( xiz*DxVz + ztz*DzVz);

        hTzz[iptr] = lam2mu * ( xiz*DxVz  + ztz*DzVz)
                    +lam    * ( xix*DxVx  + ztx*DzVx);

        hTxz[iptr] = mu *(
                     xiz*DxVx + xix*DxVz
                    +ztz*DzVx + ztx*DzVz
                    );

        iptr += 1;
      }
  }

  return;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_eq1st_curv_col_el_iso_rhs_timg_z2(
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
    const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTx,DzTz;
  float slwjac;
  float xix,xiz,ztx,ztz;

  // to save traction and other two dir force var
  float vecxi[fdx_len];
  float veczt[fdz_len];
  int n, iptr4vec;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // last indx, free surface force Tx/Ty/Tz to 0 in cal
  int k_min = nk2 - fdz_indx[fdz_len-1];

  // point affected by timg
  for (size_t k=k_min; k <= nk2; k++)
  {
    // k corresponding to 0 index of the fd op

    // index of free surface
    int n_free = nk2 - k - fdz_indx[0]; // first indx is negative

    // for 1d index
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
          // metric
          xix = xi_x[iptr];
          xiz = xi_z[iptr];
          ztx = zt_x[iptr];
          ztz = zt_z[iptr];

          // slowness and jac
          slwjac = slw3d[iptr] / jac3d[iptr];

          //
          // for hVx
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txx[iptr4vec]
                                          + xi_z[iptr4vec] * Txz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_iz;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                          + zt_z[iptr4vec] * Txz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++)
          {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            //int n_img = n_free-(n-n_free);
            iptr4vec = iptr + n_img * siz_iz;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                           + zt_z[iptr4vec] * Txz[iptr4vec] );
            //veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVx[iptr] = ( DxTx+DzTz ) * slwjac;

          //
          // for hVz
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txz[iptr4vec]
                                          + xi_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_iz;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                          + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            iptr4vec = iptr + n_img * siz_iz;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                           + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // for hVx 
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVz[iptr] = ( DxTx+DzTz ) * slwjac;

          // next
          iptr += 1;
      }
  }

  return;
}

/*
 * implement vlow boundary
 */

void
sv_eq1st_curv_col_el_iso_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict vecVx2Vz, 
    int ni1, int ni2, int nk1, int nk2,
    size_t siz_iz, 
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
    const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];

  // allocate max_len because fdz may have different lens
  float  lfdz_coef [fdz_max_len];
  int    lfdz_shift[fdz_max_len];

  // local var
  int i,k;
  int n_fd; // loop var for fd
  int fdz_len;

  // local var
  float DxVx,DxVz;
  float DzVx,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiz,ztx,ztz;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }

  // loop near surface layers
  //for (size_t n=0; n < 1; n++)
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // get pos and len for this point
    int  lfdz_len  = fdz_op[n].total_len;
    // point to indx/coef for this point
    int   *p_fdz_indx  = fdz_op[n].indx;
    float *p_fdz_coef  = fdz_op[n].coef;
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_shift[n_fd] = p_fdz_indx[n_fd] * siz_iz;
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd];
    }

    // for index
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (i=ni1; i<=ni2; i++)
      {
        // metric
        xix = xi_x[iptr];
        xiz = xi_z[iptr];
        ztx = zt_x[iptr];
        ztz = zt_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        if (k==nk2) // at surface, convert
        {
          size_t ij = (i )*4;
          DzVx = vecVx2Vz[ij+2*0+0] * DxVx
               + vecVx2Vz[ij+2*0+1] * DxVz;

          DzVz = vecVx2Vz[ij+2*1+0] * DxVx
               + vecVx2Vz[ij+2*1+1] * DxVz;
        }
        else // lower than surface, lower order
        {
          M_FD_SHIFT(DzVx, Vx, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVz, Vz, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
        }

        // Hooke's equatoin
        hTxx[iptr] =  lam2mu * ( xix*DxVx  + ztx*DzVx)
                    + lam    * ( xiz*DxVz  + ztz*DzVz);

        hTzz[iptr] = lam2mu * ( xiz*DxVz  + ztz*DzVz)
                    +lam    * ( xix*DxVx  + ztx*DzVx);

        hTxz[iptr] = mu *(
                     xiz*DxVx + xix*DxVz
                    +ztz*DzVx + ztx*DzVz
                    );

        iptr += 1;
      }
  }

  return;
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

void
sv_eq1st_curv_col_el_iso_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict xi_x, float *restrict xi_z,
    float *restrict zt_x, float *restrict zt_z,
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    int nk2, size_t siz_iz, 
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    bdrypml_t *bdrypml, bdryfree_t *bdryfree,
    const int verbose)
{

  float *vecVx2Vz = bdryfree->vecVx2Vz2;

  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // val on point
  float DxTxx,DxTzz,DxTxz,DxVx,DxVz;
  float DzTxx,DzTzz,DzTxz,DzVx,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiz,ztx,ztz;
  float hVx_rhs,hVz_rhs;
  float hTxx_rhs,hTzz_rhs,hTxz_rhs;
  // for free surface
  float Dx_DzVx,Dx_DzVz;

  // local
  int i,k;
  int iptr, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_B_minus_1;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      // skip to next face if not cfspml
      if (bdrypml->is_at_sides[idim][iside] == 0) continue;

      // get index into local var
      int abs_ni1 = bdrypml->ni1[idim][iside];
      int abs_ni2 = bdrypml->ni2[idim][iside];
      int abs_nk1 = bdrypml->nk1[idim][iside];
      int abs_nk2 = bdrypml->nk2[idim][iside];

      // get coef for this face
      float *restrict ptr_coef_A = bdrypml->A[idim][iside];
      float *restrict ptr_coef_B = bdrypml->B[idim][iside];
      float *restrict ptr_coef_D = bdrypml->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);

      // get pml vars
      float *restrict abs_vars_cur = auxvar->cur;
      float *restrict abs_vars_rhs = auxvar->rhs;

      float *restrict pml_Vx   = abs_vars_cur + auxvar->Vx_pos;
      float *restrict pml_Vz   = abs_vars_cur + auxvar->Vz_pos;
      float *restrict pml_Txx  = abs_vars_cur + auxvar->Txx_pos;
      float *restrict pml_Tzz  = abs_vars_cur + auxvar->Tzz_pos;
      float *restrict pml_Txz  = abs_vars_cur + auxvar->Txz_pos;

      float *restrict pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
      float *restrict pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
      float *restrict pml_hTxx = abs_vars_rhs + auxvar->Txx_pos;
      float *restrict pml_hTzz = abs_vars_rhs + auxvar->Tzz_pos;
      float *restrict pml_hTxz = abs_vars_rhs + auxvar->Txz_pos;

      // for each dim
      if (idim == 0 ) // x direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_iz;
            iptr = iptr_k + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // pml coefs
              int abs_i = i - abs_ni1;
              coef_D = ptr_coef_D[abs_i];
              coef_A = ptr_coef_A[abs_i];
              coef_B = ptr_coef_B[abs_i];
              coef_B_minus_1 = coef_B - 1.0;

              // metric
              xix = xi_x[iptr];
              xiz = xi_z[iptr];
              //xiz = fabs(xiz);

              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              slw = slw3d[iptr];
              lam2mu = lam + 2.0 * mu;

              // xi derivatives
              M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTzz, Tzz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( xix*DxTxx + xiz*DxTxz );
               hVz_rhs = slw * ( xix*DxTxz + xiz*DxTzz );
              hTxx_rhs = lam2mu*xix*DxVx + lam*xiz*DxVz;
              hTzz_rhs = lam*xix*DxVx + lam2mu*xiz*DxVz;
              hTxz_rhs = mu*( xiz*DxVx + xix*DxVz );

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
              hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
              pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
              pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];

              // add contributions from free surface condition
              //  not consider timg because conflict with main cfspml,
              //     need to revise in the future if required
              if (bdryfree->is_at_sides[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i )*4;
                Dx_DzVx = vecVx2Vz[ij+2*0+0] * DxVx
                        + vecVx2Vz[ij+2*0+1] * DxVz;

                Dx_DzVz = vecVx2Vz[ij+2*1+0] * DxVx
                        + vecVx2Vz[ij+2*1+1] * DxVz;

                // metric
                ztx = zt_x[iptr];
                ztz = zt_z[iptr];

                // keep xi derivative terms, including free surface convered
                hTxx_rhs =    lam2mu * (            ztx*Dx_DzVx)
                            + lam    * (            ztz*Dx_DzVz);

                hTzz_rhs =   lam2mu * (            ztz*Dx_DzVz)
                            +lam    * (            ztx*Dx_DzVx);

                hTxz_rhs = mu *(
                             ztz*Dx_DzVx + ztx*Dx_DzVz
                            );

                // make corr to Hooke's equatoin
                hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
                hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
                hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;

                // aux var
                //   a1 = alpha + d / beta, dealt in abs_set_cfspml
                pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
                pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
                pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
              }

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
        } // k
      }
      else // z direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_iz;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_B_minus_1 = coef_B - 1.0;

            iptr = iptr_k + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // metric
              ztx = zt_x[iptr];
              ztz = zt_z[iptr];
              //ztx = fabs(ztx);

              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              slw = slw3d[iptr];
              lam2mu = lam + 2.0 * mu;

              // zt derivatives
              M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxx, Txx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( ztx*DzTxx + ztz*DzTxz );
               hVz_rhs = slw * ( ztx*DzTxz + ztz*DzTzz );
              hTxx_rhs = lam2mu*ztx*DzVx + lam*ztz*DzVz;
              hTzz_rhs = lam*ztx*DzVx + lam2mu*ztz*DzVz;
              hTxz_rhs = mu*( ztz*DzVx + ztx*DzVz );

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
              hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
              pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
              pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
        } // k
      } // if which dim
    } // iside
  } // idim

  return;
}

/*******************************************************************************
 * free surface coef
 ******************************************************************************/

int
sv_eq1st_curv_col_el_iso_dvh2dvz(gdinfo_t        *gdinfo,
                                 gdcurv_metric_t *metric,
                                 md_t       *md,
                                 bdryfree_t      *bdryfree,
                                 const int verbose)
{
  int ierr = 0;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nx  = gdinfo->nx;
  int nz  = gdinfo->nz;
  size_t siz_iz  = gdinfo->siz_iz;

  // point to each var
  float *restrict xi_x = metric->xi_x;
  float *restrict xi_z = metric->xi_z;
  float *restrict zt_x = metric->zeta_x;
  float *restrict zt_z = metric->zeta_z;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *vecVx2Vz = bdryfree->vecVx2Vz2;
  
  float A[2][2], B[2][2];
  float AB[2][2];

  int k = nk2;

    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + k * siz_iz;

      float e11 = xi_x[iptr];
      float e12 = xi_z[iptr];
      float e21 = zt_x[iptr];
      float e22 = zt_z[iptr];

      float lam    = lam3d[iptr];
      float miu    =  mu3d[iptr];
      float lam2mu = lam + 2.0f * miu;

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

      size_t ij = (i) * 4;

      // save into mat
      for(int irow = 0; irow < 2; irow++)
        for(int jcol = 0; jcol < 2; jcol++){
          vecVx2Vz[ij + irow*2 + jcol] = AB[irow][jcol];
        }
    }

  return ierr;
}

/*******************************************************************************
 * add source terms
 ******************************************************************************/

int
sv_eq1st_curv_col_el_iso_rhs_src(
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict jac3d, float *restrict slw3d,
    src_t *src, // short nation for reference member
    const int verbose)
{
  int ierr = 0;

  // local var
  int si,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi / mij
  float fx, fz;
  float Mxx,Mzz,Mxz;

  int it     = src->it;
  int istage = src->istage;

  // add src; is is a commont iterater var
  for (int is=0; is < src->total_number; is++)
  {
    int   it_start = src->it_begin[is];
    int   it_end   = src->it_end  [is];

    if (it >= it_start && it <= it_end)
    {
      int   *ptr_ext_indx = src->ext_indx + is * max_ext;
      float *ptr_ext_coef = src->ext_coef + is * max_ext;
      int it_to_it_start = it - it_start;
      int iptr_cur_stage =   is * src->max_nt * src->max_stage // skip other src
                           + it_to_it_start * src->max_stage // skip other time step
                           + istage;
      if (src->force_actived == 1) {
        fx  = src->Fx [iptr_cur_stage];
        fz  = src->Fz [iptr_cur_stage];
      }
      if (src->moment_actived == 1) {
        Mxx = src->Mxx[iptr_cur_stage];
        Mzz = src->Mzz[iptr_cur_stage];
        Mxz = src->Mxz[iptr_cur_stage];
      }
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        //fprintf(stdout,"-----> src_iptr=%d,src_ceof=%g\n",iptr,coef);

        if (src->force_actived == 1) {
          float V = coef * slw3d[iptr] / jac3d[iptr];
          hVx[iptr] += fx * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / jac3d[iptr];
          hTxx[iptr] -= Mxx * rjac;
          hTzz[iptr] -= Mzz * rjac;
          hTxz[iptr] -= Mxz * rjac;
        }
      } // i_ext

    } // it
  } // is

  return ierr;
}

int
sv_eq1st_curv_col_el_iso_rhs_src_dbg(
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict jac3d, float *restrict slw3d,
    src_t *src, // short nation for reference member
    const int verbose)
{
  int ierr = 0;

  // local var
  int si,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi / mij
  float fx, fz;
  float Mxx,Mzz,Mxz;

  int it     = src->it;
  int istage = src->istage;

  // add src; is is a commont iterater var
  for (int is=0; is < src->total_number; is++)
  {
    int   it_start = src->it_begin[is];
    int   it_end   = src->it_end  [is];

    if (it >= it_start && it <= it_end)
    {
      int   *ptr_ext_indx = src->ext_indx + is * max_ext;
      float *ptr_ext_coef = src->ext_coef + is * max_ext;
      int it_to_it_start = it - it_start;
      int iptr_cur_stage =   is * src->max_nt * src->max_stage // skip other src
                           + it_to_it_start * src->max_stage // skip other time step
                           + istage;
      if (src->force_actived == 1) {
        fx  = src->Fx [iptr_cur_stage];
        fz  = src->Fz [iptr_cur_stage];
      }
      if (src->moment_actived == 1) {
        Mxx = src->Mxx[iptr_cur_stage];
        Mzz = src->Mzz[iptr_cur_stage];
        Mxz = src->Mxz[iptr_cur_stage];
      }
      
      // for extend points
        int   iptr = 39 + 3 + (49 + 3) * 106;
        float coef = 1.0;
        fprintf(stdout,"-----> src_iptr=%d,src_ceof=%g\n",iptr,coef);

        if (src->force_actived == 1) {
          float V = coef * slw3d[iptr] / jac3d[iptr];
          hVx[iptr] += fx * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / jac3d[iptr];
          hTxx[iptr] -= Mxx * rjac;
          hTzz[iptr] -= Mzz * rjac;
          hTxz[iptr] -= Mxz * rjac;
        }

    } // it
  } // is

  return ierr;
}

