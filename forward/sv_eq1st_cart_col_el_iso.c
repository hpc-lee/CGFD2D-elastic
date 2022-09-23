/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_eq1st_cart_col_el_iso.h"


/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eq1st_cart_col_el_iso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  gd_t    *gdcart,
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
  size_t siz_iz   = gdinfo->siz_iz;

  float *vecVx2Vz = bdryfree->vecVx2Vz2;

  float dx = gdcart->dx;
  float dz = gdcart->dz;

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

  // free surface at z2
  if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1)
  {
    // tractiong
    sv_eq1st_cart_col_el_iso_rhs_timg_z2(Tzz,Txz,
                                        ni1,ni2,nk1,nk2,nz,
                                        siz_iz,
                                        verbose);
  }

  // inner points
  sv_eq1st_cart_col_el_iso_rhs_inner(Vx,Vz,Txx,Tzz,Txz,
                                    hVx,hVz,hTxx,hTzz,hTxz,
                                    lam3d, mu3d, slw3d,
                                    dx, dz,
                                    ni1,ni2,nk1,nk2,siz_iz,
                                    fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                    fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                    verbose);

  // free, abs, source in turn

  // free surface at z2
  if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1)
  {
    // velocity: vlow
    sv_eq1st_cart_col_el_iso_rhs_vlow_z2(Vx,Vz,hTxx,hTzz,hTxz,
                                        lam3d, mu3d, slw3d,
                                        vecVx2Vz,
                                        dx, dz,
                                        ni1,ni2,nk1,nk2,siz_iz,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        num_of_fdz_op,fdz_op,fdz_max_len,
                                        verbose);
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_cart_col_el_iso_rhs_cfspml(Vx,Vz,Txx,Tzz,Txz,
                                       hVx,hVz,hTxx,hTzz,hTxz,
                                       lam3d, mu3d, slw3d,
                                       dx, dz,
                                       nk2, siz_iz,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       bdrypml, bdryfree,
                                       verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_eq1st_cart_col_el_iso_rhs_src(hVx,hVz,hTxx,hTzz,hTxz,
                                    slw3d, dx, dz,
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
sv_eq1st_cart_col_el_iso_rhs_inner(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float dx, float dz,
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
  float DxTxx,             DxTxz,DxVx,DxVz;
  float             DzTzz, DzTxz,DzVx,DzVz;
  float lam,mu,mu2,slw,Eii;

  float *restrict Vx_ptr;
  float *restrict Vz_ptr;
  float *restrict Txx_ptr;
  float *restrict Txz_ptr;
  float *restrict Tzz_ptr;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i] / dx;
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k] / dz;
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

        // Tzz derivatives
        M_FD_SHIFT_PTR_MACDRP(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // medium
        slw = slw3d[iptr];
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        // moment equation
        hVx[iptr] = slw*( DxTxx + DzTxz );
        hVz[iptr] = slw*( DxTxz + DzTzz );

        // Hooke's equatoin
        Eii = lam * (DxVx + DzVz);

        hTxx[iptr] = Eii + mu2 * DxVx;
        hTzz[iptr] = Eii + mu2 * DzVz;
        hTxz[iptr] = mu *( DxVz + DzVx );

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
sv_eq1st_cart_col_el_iso_rhs_timg_z2(
    float *restrict  Tzz, float *restrict  Txz, 
    int ni1, int ni2, int nk1, int nk2, int nz,
    size_t siz_iz,
    const int verbose)
{
  // nk2
  size_t iptr_k = nk2 * siz_iz;

    size_t iptr = iptr_k + ni1;
    for (size_t i=ni1; i<=ni2; i++)
    {
      Tzz[iptr] = 0.0;
      Txz[iptr] = 0.0;

      // next
      iptr += 1;
    }

  // mirror point
  for (size_t k=nk2+1; k<nz; k++)
  {
    int k_phy = nk2 - (k-nk2);
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr_gho = i + k     * siz_iz;
        size_t iptr_phy = i + k_phy * siz_iz;

        Tzz[iptr_gho] = -Tzz[iptr_phy];
        Txz[iptr_gho] = -Txz[iptr_phy];
      }
  }

  return;
}

/*
 * implement vlow boundary
 */

void
sv_eq1st_cart_col_el_iso_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict vecVx2Vz, 
    float dx, float dz,
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
  float lam,mu,mu2,slw,Eii;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i] / dx;
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
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    }

    // for index
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (i=ni1; i<=ni2; i++)
      {
        // medium
        slw = slw3d[iptr];
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        if (k==nk2) // at surface, convert
        {
          size_t ij = (i)*4;
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
        Eii = lam * (DxVx + DzVz);

        hTxx[iptr] = Eii + mu2 * DxVx;
        hTzz[iptr] = Eii + mu2 * DzVz;
        hTxz[iptr] = mu *( DxVz + DzVx );

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
sv_eq1st_cart_col_el_iso_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    float dx, float dz,
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
  float DxTxx,      DxTxz,DxVx,DxVz;
  float       DzTzz,DzTxz,DzVx,DzVz;
  float lam,mu,lam2mu,slw;
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
    lfdx_coef [i] = fdx_coef[i] / dx;
    lfdx_shift[i] = fdx_indx[i];
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k] / dz;
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

              // medium
              slw = slw3d[iptr];
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              lam2mu = lam + 2.0 * mu;

              // xi derivatives
              M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( DxTxx  );
               hVz_rhs = slw * ( DxTxz  );
              hTxx_rhs = lam2mu * DxVx;
              hTzz_rhs = lam * DxVx;
              hTxz_rhs = mu*( DxVz );

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

                // keep xi derivative terms, including free surface convered
                hTxx_rhs = lam    * Dx_DzVz;

                hTzz_rhs = lam2mu * Dx_DzVz;

                hTxz_rhs = mu     * Dx_DzVx ;

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
              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              slw = slw3d[iptr];
              lam2mu = lam + 2.0 * mu;

              // zt derivatives
              M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * DzTxz;
               hVz_rhs = slw * DzTzz;
              hTxx_rhs = lam    * DzVz;
              hTzz_rhs = lam2mu * DzVz;
              hTxz_rhs = mu     * DzVx; 

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
sv_eq1st_cart_col_el_iso_dvh2dvz(gdinfo_t   *gdinfo,
                                 md_t       *md,
                                 bdryfree_t      *bdryfree,
                                 const int verbose)
{
  int ierr = 0;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk2 = gdinfo->nk2;
  size_t siz_iz   = gdinfo->siz_iz;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *vecVx2Vz = bdryfree->vecVx2Vz2;
  
  float lam2mu, lam, mu;
 
  int k = nk2;

    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + k * siz_iz;

      lam    = lam3d[iptr];
      mu     =  mu3d[iptr];
      lam2mu = lam + 2.0f * mu;

      size_t ij = ( i) * 4;

      // save into mat
      for(int irow = 0; irow < 2; irow++)
        for(int jcol = 0; jcol < 2; jcol++){
          vecVx2Vz[ij + irow*2 + jcol] = 0.0;
        }

      // DzVx = -DxVz
      int DzVx_j = 0;
      int DxVz_i = 1;
      float coef = -1.0;
      vecVx2Vz[ij + DzVx_j * CONST_NDIM + DxVz_i] = coef;

      // DzVz = - (lam/lam2mu) * DxVx - (lam/lam2mu) * DyVy
      int DzVz_j = 1;
      int DxVx_i = 0;
      coef = - lam / lam2mu;
      vecVx2Vz[ij + DzVz_j * CONST_NDIM + DxVx_i] = coef;
    }

  return ierr;
}

/*******************************************************************************
 * add source terms
 ******************************************************************************/

int
sv_eq1st_cart_col_el_iso_rhs_src(
    float *restrict hVx , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTzz,
    float *restrict hTxz, 
    float *restrict slw3d,
    float dx, float dz,
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

  float vol = dx * dz;

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

        if (src->force_actived == 1) {
          float V = coef * slw3d[iptr] / vol;
          hVx[iptr] += fx * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / vol;
          hTxx[iptr] -= Mxx * rjac;
          hTzz[iptr] -= Mzz * rjac;
          hTxz[iptr] -= Mxz * rjac;
        }
      } // i_ext

    } // it
  } // is

  return ierr;
}

