/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using cart grid and staggerd scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"
#include "sv_eq1st_cart_stg_el_iso.h"

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 *  simple MPI exchange without computing-communication overlapping
 ******************************************************************************/

void
sv_eq1st_cart_stg_el_iso_allstep(
  fdstg_t    *fd,
  gdinfo_t   *gdinfo,
  gd_t       *gdcart,
  md_t       *md,
  src_t      *src,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  wav_t      *wav,
  iorecv_t   *iorecv,
  ioline_t   *ioline,
  iosnap_t   *iosnap,
  // time
  float dt, int nt_total, float t0,
  char *output_dir,
  int qc_check_nan_num_of_step,
  const int output_all, // qc all var
  const int verbose)
{
  // retrieve from struct
  float dx = gdcart->dx;
  float dz = gdcart->dz;

  float t_cur;
  float t_end; // time after this loop for nc output

  // create snapshot nc output files
  if (verbose>0) fprintf(stdout,"prepare snap nc output ...\n"); 
  iosnap_nc_t  iosnap_nc;
  io_snap_nc_create(iosnap, &iosnap_nc);

  // alloc working space for slice and snap ouotput
  size_t siz_max_wrk = iosnap->siz_max_wrk;
  float *restrict wrk = (float *) fdlib_mem_calloc_1d_float(siz_max_wrk,
                        0.0, "wrk in stg el iso");

  // calculate conversion matrix for free surface
  if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1)
  {
     sv_eq1st_cart_stg_el_iso_dvh2dvz(gdinfo,md,bdryfree,verbose);
  }

  //--------------------------------------------------------
  // time loop
  //--------------------------------------------------------

  if (verbose>0) fprintf(stdout,"start time loop ...\n"); 

  for (int it=0; it<nt_total; it++)
  {
    //
    // Update Tij first, considering force Green function 
    //

    t_cur  = it * dt + t0;
    t_end = t_cur + 0.5 * dt;

    if (verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // set src_t time
    src_set_time(src, it, 0);

    sv_eq1st_cart_stg_el_iso_hook(fd, wav,
        gdinfo, md, bdryfree, bdrypml, src,
        dt, dx, dz,
        verbose);

    //
    // Update Vi
    //

    t_cur  = it * dt + t0 + 0.5 * dt;
    t_end = t_cur + 0.5 * dt;
    if (verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // set src_t time
    src_set_time(src, it, 0);

    // from Tij to Vi
    sv_eq1st_cart_stg_el_iso_moment(fd,wav,
        gdinfo, md, bdryfree, bdrypml, src,
        dt, dx, dz,
        verbose);

    //--------------------------------------------
    // QC
    //--------------------------------------------

    if (qc_check_nan_num_of_step >0  && (it % qc_check_nan_num_of_step) == 0) {
      if (verbose>10) fprintf(stdout,"-> check value nan\n");
        //wav_check_value(wav->v5d);
    }

    //--------------------------------------------
    // save results
    //--------------------------------------------

    //-- recv by interp
    io_recv_keep(iorecv, wav->v5d, it, wav->ncmp, wav->siz_icmp);

    //-- line values
    io_line_keep(ioline, wav->v5d, it, wav->ncmp, wav->siz_icmp);

    // snapshot
    io_snap_nc_put(iosnap, &iosnap_nc, gdinfo, wav, 
                   wav->v5d, wrk, nt_total, it, t_end, 0, 1, 0);

    io_snap_nc_put(iosnap, &iosnap_nc, gdinfo, wav, 
                   wav->v5d, wrk, nt_total, it, t_end, 1, 0, 1);

    // debug output
    if (output_all==1)
    {
      char ou_file[CONST_MAX_STRLEN];

        io_build_fname_time(output_dir,"wav",".nc",it,ou_file);
        io_var2d_export_nc(ou_file,
                           wav->v5d,
                           wav->cmp_pos,
                           wav->cmp_name,
                           wav->ncmp,
                           gdinfo->index_name,
                           gdinfo->nx,
                           gdinfo->nz);
    }

  } // time loop

  // postproc

  // close nc
  io_snap_nc_close(&iosnap_nc);

  if (wrk) free(wrk);

  return;
}

/*******************************************************************************
 * perform one stage calculation
 ******************************************************************************/

void
sv_eq1st_cart_stg_el_iso_hook(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  float dt, float dx, float dz,
  const int verbose)
{
  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  float *restrict Txx   = wav->v5d + wav->Txx_pos;
  float *restrict Tzz   = wav->v5d + wav->Tzz_pos;
  float *restrict Txz   = wav->v5d + wav->Txz_pos;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *vecVx2Vz = bdryfree->vecVx2Vz2;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nz  = gdinfo->nz ;

  size_t siz_iz  = gdinfo->siz_iz;

  // local var
  float lam, mu, mu2;
  float Eii;
  float DxVx,DzVx;
  float DxVz,DzVz;

  // allocate max_len because fdz may have different lens
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  int n_fd;
  float *restrict Vx_ptr;
  float *restrict Vz_ptr;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op

  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  // free surface at z2
  if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1)
  {
    sv_eq1st_cart_stg_el_iso_free_vzero(Vx,Vz,
                                  ni1,ni2,nk1,nk2,nz,
                                  siz_iz,
                                  verbose);
  }

  // at surface layers
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    int k = nk2 - n;

    // use normal op if not free surface
    int lfdz_op_n = num_of_fdz_op - 1; 
    // use lower order free surface at z2
    if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1) {
      lfdz_op_n = n;
    }

    // get pos and len for this point
    int  lfdz_len      = fd->lay_fdz_op[lfdz_op_n].total_len;
    int   *p_fdz_indx  = fd->lay_fdz_op[lfdz_op_n].indx;
    float *p_fdz_coef  = fd->lay_fdz_op[lfdz_op_n].coef;
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
      lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_iz;
      lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_iz;
    }

    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        Vx_ptr = Vx + iptr;
        Vz_ptr = Vz + iptr;

        // Derivatives for Tii
        M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // derivatives for Txz
        M_FD_SHIFT_PTR(DxVz, Vz_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzVx, Vx_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // if surface
        if (k==nk2 && bdryfree->is_at_sides[CONST_NDIM-1][1]) // at surface, convert
        {
          size_t ij = (i )*4;
          DzVz = vecVx2Vz[ij+1*CONST_NDIM+0] * DxVx;
        }

        // Hooke's equatoin
        Eii = lam * (DxVx + DzVz);

        Txx[iptr] += dt * (Eii + mu2 * DxVx);
        Tzz[iptr] += dt * (Eii + mu2 * DzVz);

        Txz[iptr] += dt * mu *( DxVz + DzVx );

        iptr += 1;
      }
  }

  // for inner layers
  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]    = p_fdz_coef[n_fd] / dz;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_iz;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_iz;
  }

  //for (size_t k=nk1; k <= nk2; k++)
  for (size_t k=nk1; k <= nk2 - num_of_fdz_op + 1; k++)
  {
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        Vx_ptr = Vx + iptr;
        Vz_ptr = Vz + iptr;

        // Derivatives for Tii
        M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // derivatives for Txz
        M_FD_SHIFT_PTR(DxVz, Vz_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzVx, Vx_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // Hooke's equatoin
        Eii = lam * (DxVx + DzVz);

        Txx[iptr] += dt * (Eii + mu2 * DxVx);
        Tzz[iptr] += dt * (Eii + mu2 * DzVz);

        Txz[iptr] += dt * mu *( DxVz + DzVx );

        iptr += 1;
      }
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_cart_stg_el_iso_hook_cfspml(fd,Vx,Vz,Txx,Tzz,Txz,
                                       lam3d, mu3d, 
                                       dt, dx, dz,
                                       nk2, siz_iz,
                                       bdrypml, bdryfree,
                                       verbose);
  }

  // add source term
  if (src->moment_actived == 1)
  {
    sv_eq1st_cart_stg_el_iso_hook_src(Txx,Tzz,Txz,
                                     dt, dx, dz, src,
                                     verbose);
  }
  // end func

  return;
}

void
sv_eq1st_cart_stg_el_iso_moment(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  float dt, float dx, float dz,
  const int verbose)
{
  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  float *restrict Txx   = wav->v5d + wav->Txx_pos;
  float *restrict Tzz   = wav->v5d + wav->Tzz_pos;
  float *restrict Txz   = wav->v5d + wav->Txz_pos;

  float *restrict slw3d = md->rho;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nz  = gdinfo->nz ;

  size_t siz_iz  = gdinfo->siz_iz;

  // local var
  float slwdt;
  float DxTxx,            DxTxz;
  float             DzTzz,DzTxz;

  // allocate max_len because fdz may have different lens
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  int n_fd;
  float *restrict Txx_ptr;
  float *restrict Tzz_ptr;
  float *restrict Txz_ptr;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op

  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  // for inner layers
  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_iz;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_iz;
  }

  // free surface at z2
  if (bdryfree->is_at_sides[CONST_NDIM-1][1] == 1)
  {
    sv_eq1st_cart_stg_el_iso_free_simg(Tzz,Txz,
                                  ni1,ni2,nk1,nk2,nz,
                                  siz_iz,
                                  verbose);
  }

  for (size_t k=nk1; k <= nk2; k++)
  {
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        slwdt = slw3d[iptr] * dt;

        Txx_ptr = Txx + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;

        // Tii deriv
        M_FD_SHIFT_PTR(DxTxx, Txx_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzTzz, Tzz_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // Txz deriv
        M_FD_SHIFT_PTR(DxTxz, Txz_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzTxz, Txz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // Moment equatoin
        Vx[iptr] += slwdt * (DxTxx + DzTxz);
        Vz[iptr] += slwdt * (DxTxz + DzTzz);

        iptr += 1;
      }
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_cart_stg_el_iso_moment_cfspml(fd,Vx,Vz,Txx,Tzz,Txz,
                                       slw3d, dt, dx, dz,
                                       nk2, siz_iz,
                                       bdrypml,
                                       verbose);
  }

  // add source term
  if (src->force_actived == 1)
  {
    sv_eq1st_cart_stg_el_iso_moment_src(Vx, Vz,
                                     slw3d, dt, dx, dz, src,
                                    verbose);
  }
  // end func

  return;
}


/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * updating Tij
 */

void
sv_eq1st_cart_stg_el_iso_hook_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict lam3d, float *restrict  mu3d,
    float dt, float dx, float dz,
    int nk2, size_t siz_iz,
    bdrypml_t *bdrypml, bdryfree_t *bdryfree,
    const int verbose)
{
  float *vecVx2Vz = bdryfree->vecVx2Vz2;

  // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  int n_fd;
  float *restrict Vx_ptr;
  float *restrict Vz_ptr;

  // val on point
  float lam,mu,lam2mu;

  // local
  int i,k;
  int iptr, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_rB_minus_1;
  float coef_Am, coef_Bm, coef_Dm, coef_rBm_minus_1;

  // put fd op into local array
  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_iz;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_iz;
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
      float *restrict ptr_coef_Am = bdrypml->Am[idim][iside];
      float *restrict ptr_coef_Bm = bdrypml->Bm[idim][iside];
      float *restrict ptr_coef_Dm = bdrypml->Dm[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);

      // for each dim
      if (idim == 0 ) // x direction
      {
        float DxVx,DxVz;
        // get pml vars, some vars are not used
        float *restrict pml_DxVx   = auxvar->var + auxvar->Vx_pos;
        float *restrict pml_DxVz   = auxvar->var + auxvar->Vz_pos;

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
              coef_Dm = ptr_coef_Dm[abs_i];
              coef_Am = ptr_coef_Am[abs_i];
              coef_Bm = ptr_coef_Bm[abs_i];

              coef_rB_minus_1 = coef_B - 1.0;
              coef_rBm_minus_1 = coef_Bm - 1.0;

              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              lam2mu = lam + 2.0 * mu;

              Vx_ptr = Vx + iptr;
              Vz_ptr = Vz + iptr;

              // for x deriv
              M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
              M_FD_SHIFT_PTR(DxVz, Vz_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
              float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

              // correct Txx,Tyy,Tzz
              float Vx1 = (2.0 * pml_DxVx[iptr_a] + dt* coef_D *DxVx ) * r_dt_AD;
              float Vx1_corr = coef_rB_minus_1 * DxVx - Vx1 * coef_B;
              Txx[iptr] += dt * lam2mu * Vx1_corr;
              Tzz[iptr] += dt * lam    * Vx1_corr;

              // correct Tij
              float Vz1 =  (2.0 * pml_DxVz[iptr_a] + dt * coef_Dm * DxVz ) * r_dt_ADm; 
              float Vz1_corr = coef_rBm_minus_1 * DxVz - Vz1 * coef_Bm;
              Txz[iptr] += dt * mu * Vz1_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DxVx[iptr_a]  = 2.0 * Vx1  - pml_DxVx[iptr_a];
              pml_DxVz[iptr_a]  = 2.0 * Vz1  - pml_DxVz[iptr_a];

              // add contributions from free surface condition
              //     only consider Tii on surface
              if (bdryfree->is_at_sides[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i )*4;
                float Dx_DzVz = vecVx2Vz[ij+1*CONST_NDIM+0] * DxVx;

                // make corr to Hooke's equatoin
                float lamsq_lam2u = lam*lam/lam2mu;
                Txx[iptr] +=   dt * lam * Dx_DzVz * coef_rB_minus_1
                             + dt * lamsq_lam2u * Vx1 * coef_B;
              }

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
        } // k
      }
      else // z direction
      {
        float DzVx,DzVz;

        // get pml vars, some vars are not used
        float *restrict pml_DzVx   = auxvar->var + auxvar->Vx_pos;
        float *restrict pml_DzVz   = auxvar->var + auxvar->Vz_pos;

        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_iz;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_rB_minus_1 = coef_B - 1.0;

          coef_Dm = ptr_coef_Dm[abs_k];
          coef_Am = ptr_coef_Am[abs_k];
          coef_Bm = ptr_coef_Bm[abs_k];
          coef_rBm_minus_1 = coef_Bm - 1.0;

          float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
          float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

            iptr = iptr_k + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              lam2mu = lam + 2.0 * mu;

              Vx_ptr = Vx + iptr;
              Vz_ptr = Vz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DzVx, Vx_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);
              M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

              // correct Txx,Tyy,Tzz
              float Vz3 =  (2.0 * pml_DzVz[iptr_a] + dt * coef_D * DzVz ) * r_dt_AD;
              float Vz3_corr = coef_rB_minus_1 * DzVz - Vz3 * coef_B;

              Txx[iptr] += dt *lam   * Vz3_corr;
              Tzz[iptr] += dt *lam2mu* Vz3_corr;

              pml_DzVz[iptr_a] = 2.0 * Vz3 - pml_DzVz[iptr_a];

              // correct Txz
              float Vx3 =  (2.0 * pml_DzVx[iptr_a] + dt * coef_Dm * DzVx ) * r_dt_ADm;
              float Vx3_corr = coef_rBm_minus_1 * DzVx - Vx3 * coef_Bm;

              Txz[iptr] += dt * mu * Vx3_corr;

              pml_DzVx[iptr_a] = 2.0 * Vx3 - pml_DzVx[iptr_a];

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

void
sv_eq1st_cart_stg_el_iso_moment_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tzz,
    float *restrict  Txz, 
    float *restrict slw3d,
    float dt, float dx, float dz,
    int nk2, size_t siz_iz,
    bdrypml_t *bdrypml,
    const int verbose)
{
  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  // local
  int i,k;
  int iptr, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_rB_minus_1;
  float coef_Am, coef_Bm, coef_Dm, coef_rBm_minus_1;
  float slw;

  float *restrict Txx_ptr;
  float *restrict Tzz_ptr;
  float *restrict Txz_ptr;

  // put fd op into local array

  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_iz;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_iz;
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
      float *restrict ptr_coef_Am = bdrypml->Am[idim][iside];
      float *restrict ptr_coef_Bm = bdrypml->Bm[idim][iside];
      float *restrict ptr_coef_Dm = bdrypml->Dm[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);

      // for each dim
      if (idim == 0 ) // x direction
      {
        float DxTxx, DxTxz;

        // get pml vars, some vars are not used
        float *restrict pml_DxTxx   = auxvar->var + auxvar->Txx_pos;
        float *restrict pml_DxTxz   = auxvar->var + auxvar->Txz_pos;

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
              coef_rB_minus_1 = coef_B - 1.0;

              coef_Dm = ptr_coef_Dm[abs_i];
              coef_Am = ptr_coef_Am[abs_i];
              coef_Bm = ptr_coef_Bm[abs_i];
              coef_rBm_minus_1 = coef_Bm - 1.0;

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
              float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

              // medium
              slw = slw3d[iptr];

              Txx_ptr = Txx + iptr;
              Txz_ptr = Txz + iptr;

              //  x deriv
              M_FD_SHIFT_PTR(DxTxx, Txx_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
              M_FD_SHIFT_PTR(DxTxz, Txz_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);

              // correct Vx
              float Txx1 = (2.0 * pml_DxTxx[iptr_a] + dt* coef_Dm *DxTxx ) * r_dt_ADm;
              float Txx1_corr = coef_rBm_minus_1 * DxTxx - Txx1 * coef_Bm;
              Vx[iptr] += dt * slw * Txx1_corr;

              // correct Vz
              float Txz1 = (2.0 * pml_DxTxz[iptr_a] + dt* coef_D *DxTxz ) * r_dt_AD;
              float Txz1_corr = coef_rB_minus_1 * DxTxz - Txz1 * coef_B;
              Vz[iptr] += dt * slw * Txz1_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DxTxx[iptr_a]  = 2.0 * Txx1  - pml_DxTxx[iptr_a];
              pml_DxTxz[iptr_a]  = 2.0 * Txz1  - pml_DxTxz[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
        } // k
      }
      else // z direction
      {
        float DzTxz, DzTzz;

        // get pml vars, some vars are not used
        float *restrict pml_DzTxz   = auxvar->var + auxvar->Txz_pos;
        float *restrict pml_DzTzz   = auxvar->var + auxvar->Tzz_pos;

        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_iz;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_rB_minus_1 = coef_B - 1.0;

          coef_Dm = ptr_coef_Dm[abs_k];
          coef_Am = ptr_coef_Am[abs_k];
          coef_Bm = ptr_coef_Bm[abs_k];
          coef_rBm_minus_1 = coef_Bm - 1.0;

          float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
          float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

            iptr = iptr_k + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              slw = slw3d[iptr];

              Txz_ptr = Txz + iptr;
              Tzz_ptr = Tzz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DzTxz, Txz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);
              M_FD_SHIFT_PTR(DzTzz, Tzz_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

              // correct Vx
              float Txz3 = (2.0 * pml_DzTxz[iptr_a] + dt* coef_D *DzTxz ) * r_dt_AD;
              float Txz3_corr = coef_rB_minus_1 * DzTxz - Txz3 * coef_B;
              Vx[iptr] += dt * slw * Txz3_corr;

              // correct Vz
              float Tzz3 = (2.0 * pml_DzTzz[iptr_a] + dt* coef_Dm *DzTzz ) * r_dt_ADm;
              float Tzz3_corr = coef_rBm_minus_1 * DzTzz - Tzz3 * coef_Bm;
              Vz[iptr] += dt * slw * Tzz3_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DzTxz[iptr_a]  = 2.0 * Txz3  - pml_DzTxz[iptr_a];
              pml_DzTzz[iptr_a]  = 2.0 * Tzz3  - pml_DzTzz[iptr_a];

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
 * add source terms
 ******************************************************************************/

int
sv_eq1st_cart_stg_el_iso_hook_src(
    float *restrict Txx, float *restrict Tzz,
    float *restrict Txz, 
    float dt, float dx, float dz,
    src_t *src, // short nation for reference member
    const int verbose)
{
  int ierr = 0;

  // local var
  int si,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi / mij
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

      Mxx = src->Mxx[iptr_cur_stage];
      Mzz = src->Mzz[iptr_cur_stage];
      Mxz = src->Mxz[iptr_cur_stage];
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        float rjac = dt * coef / vol;

        Txx[iptr] -= Mxx * rjac;
        Tzz[iptr] -= Mzz * rjac;
        Txz[iptr] -= Mxz * rjac;
      } // i_ext

    } // it
  } // is

  return ierr;
}

int
sv_eq1st_cart_stg_el_iso_moment_src(
    float *restrict Vx, float *restrict Vz,
    float *restrict slw3d, float dt, float dx, float dz,
    src_t *src, // short nation for reference member
    const int verbose)
{
  int ierr = 0;

  // local var
  int si,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi 
  float fx, fz;

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

      fx  = src->Fx [iptr_cur_stage];
      fz  = src->Fz [iptr_cur_stage];
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        float V = dt * coef * slw3d[iptr] / vol;

        Vx[iptr] += fx * V;
        Vz[iptr] += fz * V;

      } // i_ext

    } // it
  } // is

  return ierr;
}

/*******************************************************************************
 * free surface coef
 ******************************************************************************/

int
sv_eq1st_cart_stg_el_iso_dvh2dvz(gdinfo_t   *gdinfo,
                                 md_t       *md,
                                 bdryfree_t      *bdryfree,
                                 const int verbose)
{
  int ierr = 0;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nk2 = gdinfo->nk2;
  size_t siz_iz  = gdinfo->siz_iz;

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

      size_t ij = (i) * 4;

      // save into mat
      for(int irow = 0; irow < 2; irow++)
        for(int jcol = 0; jcol < 2; jcol++){
          vecVx2Vz[ij + irow*2 + jcol] = 0.0;
        }

      // DzVx = -DxVz from Txz=0, not used
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
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_eq1st_cart_stg_el_iso_free_simg(
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

      Txz[iptr] = -Txz[iptr - siz_iz];

      // next
      iptr += 1;
    }

  // mirror point
  for (size_t k=nk2+1; k<nz; k++)
  {
    int k_phy = nk2 - (k-nk2);
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr_gho = i +  k     * siz_iz;
        size_t iptr_phy = i +  k_phy * siz_iz;

        Tzz[iptr_gho] = -Tzz[iptr_phy];
        Txz[iptr_gho] = -Txz[iptr_phy - siz_iz];
      }
  }

  return;
}

/*
 * set zero velocity above free surface
 */

void
sv_eq1st_cart_stg_el_iso_free_vzero(
    float *restrict  Vx, float *restrict  Vz,
    int ni1, int ni2, int nk1, int nk2, int nz,
    size_t siz_iz, 
    const int verbose)
{
  // nk2
  size_t iptr_k = nk2 * siz_iz;

    size_t iptr = iptr_k + ni1;
    for (size_t i=ni1; i<=ni2; i++)
    {
      Vz[iptr] = 0.0;

      // next
      iptr += 1;
    }

  // above surface
  for (size_t k=nk2+1; k<nz; k++)
  {
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr = i + k * siz_iz;

        Vz[iptr] = 0.0;
      }
  }

  return;
}
