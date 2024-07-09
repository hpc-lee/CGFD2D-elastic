/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"

#include "constants.h"
#include "fdlib_mem.h"
#include "md_t.h"

int
md_init(gd_t *gd, md_t *md, int media_type, int visco_type, int nmaxwell)
{
  int ierr = 0;

  md->nx = gd->nx;
  md->nz = gd->nz;

  md->siz_iz   = md->nx;
  md->siz_icmp = md->nx * md->nz;

  // media type
  md->medium_type = media_type;
  if (media_type == CONST_MEDIUM_ACOUSTIC_ISO)
  {
    md->ncmp = 2; // rho + kappa
  }
  else if (media_type == CONST_MEDIUM_ELASTIC_ISO)
  {
    md->ncmp = 3; // rho + lambda + mu
  } else if (media_type == CONST_MEDIUM_ELASTIC_VTI) {
    md->ncmp = 5; // c11 13 33 55 + rho
  } else if (media_type == CONST_MEDIUM_ELASTIC_ANISO){
    md->ncmp = 7; // 11, 13, 15, 33, 35, 55, rho
  } else if (media_type == CONST_MEDIUM_VISCOELASTIC_ISO) {
    // visco
    md->visco_type = visco_type;
    if (visco_type == CONST_VISCO_GRAVES) {
     md->ncmp = 3 + 1;
    } else if(visco_type == CONST_VISCO_GMB) {
      md->nmaxwell = nmaxwell;
      md->ncmp = 3 + 2*md->nmaxwell + 2;
    }
  } else{
    fprintf(stderr,"ERROR: media_type=%d is not implemented\n",media_type);
    exit(1);
  }

  /*
   * 0: rho
   * 1: lambda
   * 2: mu
   */
  
  // vars
  md->v3d = (float *) fdlib_mem_calloc_1d_float(
                          md->siz_icmp * md->ncmp,
                          0.0, "md_init");
  if (md->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc medium_el_iso\n");
      fflush(stderr);
  }

  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(md->ncmp,
                                                         0,
                                                         "medium_init");

  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(md->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "medium_init");

 if (visco_type == CONST_VISCO_GMB) 
 {
    md->Ylam = (float **) malloc(md->nmaxwell*sizeof(float *));
    md->Ymu  = (float **) malloc(md->nmaxwell*sizeof(float *));
    md->wl   = (float *)  malloc(md->nmaxwell*sizeof(float ));
 }


  // set pos
  for (int icmp=0; icmp < md->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * md->siz_icmp;
  }

  // init
  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","rho");
  md->rho = md->v3d + cmp_pos[icmp];

  // acoustic iso
  if (media_type == CONST_MEDIUM_ACOUSTIC_ISO) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","kappa");
    md->kappa = md->v3d + cmp_pos[icmp];
  }

  // iso
  if (media_type == CONST_MEDIUM_ELASTIC_ISO) 
  {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","lambda");
    md->lambda = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","mu");
    md->mu = md->v3d + cmp_pos[icmp];
  }

  // vti
  if (media_type == CONST_MEDIUM_ELASTIC_VTI)
  {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c11");
    md->c11 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c13");
    md->c13 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c33");
    md->c33 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c55");
    md->c55 = md->v3d + cmp_pos[icmp];
  }

  // aniso
  if (media_type == CONST_MEDIUM_ELASTIC_ANISO)
  {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c11");
    md->c11 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c13");
    md->c13 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c15");
    md->c15 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c33");
    md->c33 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c35");
    md->c35 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c55");
    md->c55 = md->v3d + cmp_pos[icmp];
  }

  // vis_iso
  if (media_type == CONST_MEDIUM_VISCOELASTIC_ISO)
  {
    if (visco_type == CONST_VISCO_GRAVES) {
      icmp += 1;
      sprintf(cmp_name[icmp],"%s","Qs");
      md->Qs = md->v3d + cmp_pos[icmp];
    } else if(visco_type == CONST_VISCO_GMB) {
      icmp += 1;
      sprintf(cmp_name[icmp],"%s","lambda");
      md->lambda = md->v3d + cmp_pos[icmp];

      icmp += 1;
      sprintf(cmp_name[icmp],"%s","mu");
      md->mu = md->v3d + cmp_pos[icmp];

      icmp += 1;
      sprintf(cmp_name[icmp],"%s","Qp");
      md->Qp = md->v3d + cmp_pos[icmp];

      icmp += 1;
      sprintf(cmp_name[icmp],"%s","Qs");
      md->Qs = md->v3d + cmp_pos[icmp];

      for(int i=0; i < md->nmaxwell; i++)
      { 
        icmp += 1;
        sprintf(cmp_name[icmp],"%s%d","Ylam",i+1);
        md->Ylam[i] = md->v3d + cmp_pos[icmp];

        icmp += 1;
        sprintf(cmp_name[icmp],"%s%d","Ymu",i+1);
        md->Ymu[i] = md->v3d + cmp_pos[icmp];
      }
    }
  }

  
  // set pointer
  md->cmp_pos  = cmp_pos;
  md->cmp_name = cmp_name;

  return ierr;
}

//
//
//

int
md_import(gd_t *gd, md_t *md, char *in_dir)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/media_px0_pz0.nc", in_dir);

  int ni1 = gd->ni1;
  int nk1 = gd->nk1;
  int ni2 = gd->ni2;
  int nk2 = gd->nk2;
  int ni  = gd->ni;
  int nk  = gd->nk;
  size_t  siz_iz = gd->siz_iz;
  
  size_t iptr, iptr1;
  
  float *var_in = (float *) malloc(sizeof(float)*ni*nk);
  size_t start[] = {0, 0};
  size_t count[] = {nk, ni};
  
  // read in nc
  int ncid;
  int varid;

  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  handle_nc_err(ierr);
  
  for (int icmp=0; icmp < md->ncmp; icmp++) 
  {
    ierr = nc_inq_varid(ncid, md->cmp_name[icmp], &varid);
    handle_nc_err(ierr);
  
    ierr = nc_get_var(ncid,varid,var_in);
    handle_nc_err(ierr);
    float *ptr = md->v3d + md->cmp_pos[icmp];
    for(int k=nk1; k<=nk2; k++) {
      for(int i=ni1; i<=ni2; i++)
      {
        iptr = i + k*siz_iz; 
        iptr1 = (i-3) + (k-3)*ni; 
        ptr[iptr] = var_in[iptr1];
      }
    }
  }

  geometric_symmetry(gd,md->v3d,md->ncmp);
  
  // close file
  ierr = nc_close(ncid);  
  handle_nc_err(ierr);

  free(var_in);

  return 0;
}

int
md_export(gd_t  *gd,
          md_t  *md,
          char *output_dir)
{
  int  number_of_vars = md->ncmp;
  int  nx = gd->nx;
  int  nz = gd->nz;
  int  ni1 = gd->ni1;
  int  nk1 = gd->nk1;
  int  ni2 = gd->ni2;
  int  nk2 = gd->nk2;
  int  ni  = gd->ni;
  int  nk  = gd->nk;
  int  gni1 = gd->gni1;
  int  gnk1 = gd->gnk1;
  size_t  siz_iz = gd->siz_iz;
  size_t iptr, iptr1;

  float *var_out = (float *) malloc(sizeof(float)*ni*nk);


  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/media_px0_pz0.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", ni, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nk, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, md->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
    handle_nc_err(ierr);
  }

  // attribute:
  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = md->v3d + md->cmp_pos[ivar];
    for(int k=nk1; k<=nk2; k++) {
      for(int i=ni1; i<=ni2; i++)
      {
        iptr = i + k*siz_iz; 
        iptr1 = (i-3) + (k-3)*ni; 
        var_out[iptr1] = ptr[iptr];
      }
    }
    ierr = nc_put_var_float(ncid, varid[ivar], var_out); 
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  free(var_out);

  return 0;
}



/*
 * test
 */

int
md_gen_test_ac_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *kappa2d = md->kappa;
  float *rho2d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
    for (size_t i=0; i<nx; i++)
    {
      size_t iptr = i + k * siz_iz;
      float Vp=3000.0;
      float rho=1500.0;
      float kappa = Vp*Vp*rho;
      kappa2d[iptr] = kappa;
      rho2d[iptr] = rho;
    }
  }

  return ierr;
}

int
md_gen_test_el_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *lam2d = md->lambda;
  float  *mu2d = md->mu;
  float *rho2d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        float Vp=3000.0;
        float Vs=2000.0;
        float rho=1500.0;
        float mu = Vs*Vs*rho;
        float lam = Vp*Vp*rho - 2.0*mu;
        lam2d[iptr] = lam;
         mu2d[iptr] = mu;
        rho2d[iptr] = rho;
      }
  }

  return ierr;
}

int
md_gen_test_el_vti(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float rho=1500.0;

        md->rho[iptr] = rho;

	      md->c11[iptr] = 25.2*1e9;//lam + 2.0f*mu;
	      md->c13[iptr] = 10.9620*1e9;//lam;
	      md->c33[iptr] = 18.0*1e9;//lam + 2.0f*mu;
	      md->c55[iptr] = 5.12*1e9;//mu;
        //-- Vp ~ sqrt(c11/rho) = 4098
      }
  }

  return ierr;
}

int
md_gen_test_el_aniso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float rho=1500.0;

        md->rho[iptr] = rho;

	      md->c11[iptr] = 25.2*1e9;//lam + 2.0f*mu;
	      md->c13[iptr] = 10.9620*1e9;//lam;
	      md->c15[iptr] = 0.0;
	      md->c33[iptr] = 18.0*1e9;//lam + 2.0f*mu;
	      md->c35[iptr] = 0.0;
	      md->c55[iptr] = 5.12*1e9;//mu;

        //-- Vp ~ sqrt(c11/rho) = 4098
      }
  }

  return ierr;
}

int
md_gen_test_Qs(md_t *md, float Qs_freq)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  md->visco_Qs_freq = Qs_freq;

  float *Qs = md->Qs;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        Qs[iptr] = 20;
      }
  }

  return ierr;
}

/*
 * convert rho to slowness to reduce number of arithmetic cal
 */

int
md_rho_to_slow(float *restrict rho, size_t siz_icmp)
{
  int ierr = 0;

  for (size_t iptr=0; iptr<siz_icmp; iptr++) {
    if (rho[iptr] > 1e-10) {
      rho[iptr] = 1.0 / rho[iptr];
    } else {
      rho[iptr] = 0.0;
    }
  }

  return ierr;
}

int 
md_vis_GMB_cal_Y(md_t *md, float freq, float fmin, float fmax)
{
  int ierr = 0;

  md->visco_GMB_freq = freq;
  md->visco_GMB_fmin = fmin;
  md->visco_GMB_fmax = fmax;

  int kmax = 2*md->nmaxwell-1;
  float wr = 2*PI*md->visco_GMB_freq;
  float wmin = 2*PI*md->visco_GMB_fmin;
  float wmax = 2*PI*md->visco_GMB_fmax;
  float wratio = wmax/wmin;

  int nmaxwell = md->nmaxwell;
  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *wk = (float *) fdlib_mem_calloc_1d_float(kmax,0,
                                                  "medium_visco_iso_cal");
  float *YP = (float *) fdlib_mem_calloc_1d_float(nmaxwell,0,
                                                  "medium_visco_iso_cal");
  float *YS = (float *) fdlib_mem_calloc_1d_float(nmaxwell,0,
                                                  "medium_visco_iso_cal");
  float **GP = (float **) fdlib_mem_calloc_2l_float(kmax, nmaxwell, 0,
                                                   "medium_visco_iso_cal");
  float **GS = (float **) fdlib_mem_calloc_2l_float(kmax, nmaxwell, 0,
                                                   "medium_visco_iso_cal");
  for(int k=0; k<kmax; k++)
  {
    wk[k] = wmin*pow(wratio,k/(float)(kmax-1));
  }

  float *wl = md->wl;
  float *Qp = md->Qp;
  float *Qs = md->Qs; 
  float *lambda = md->lambda;
  float *mu = md->mu;
  float QP1,QS1,theta1,theta2,thetatmp,lam,muu;
  float R, kappa,muunrelax;

  for(int n=0; n<nmaxwell; n++)
  {
    wl[n] = wk[2*n];
  }

  for (int k=0; k<nz; k++)
  {
    for(int i=0; i<nx; i++)
    {
      size_t iptr = i + k * siz_iz;
      QP1 = 1.0/Qp[iptr];
      QS1 = 1.0/Qs[iptr];
      lam = lambda[iptr];
      muu = mu[iptr];
      for(int m=0; m<kmax; m++)
      {
        for(int n=0; n<nmaxwell; n++) 
        {
          GP[m][n] = (wl[n]*wk[m]+pow(wl[n],2)*QP1)/(pow(wl[n],2)+pow(wk[m],2));
          GS[m][n] = (wl[n]*wk[m]+pow(wl[n],2)*QS1)/(pow(wl[n],2)+pow(wk[m],2));
        }
      }

      md_visco_LS(GP,YP,QP1,kmax,nmaxwell);
      md_visco_LS(GS,YS,QS1,kmax,nmaxwell);

      //P
      theta1 = 0.0;
      theta2 = 0.0;
      thetatmp = 0.0;

      for(int n=0; n<nmaxwell; n++)
      {
        thetatmp = thetatmp+YP[n]/(1+pow(wr/wl[n],2));
      }
      theta1 = 1-thetatmp;
  
      thetatmp = 0.0;
      for(int n=0; n<nmaxwell; n++)
      {
        thetatmp = thetatmp+YP[n]*wr/wl[n]/(1+pow(wr/wl[n],2));
      }
      theta2 = thetatmp;

      R = sqrt(pow(theta1,2)+pow(theta2,2));
      kappa = (lam+2*muu)*(R+theta1)/(2*pow(R,2));

      //S 
      theta1 = 0.0;
      theta2 = 0.0;
      thetatmp = 0.0;
    
      for(int n=0; n<nmaxwell; n++)
      {
        thetatmp = thetatmp+YS[n]/(1+pow(wr/wl[n],2));
      }
      theta1 = 1-thetatmp;
    
      thetatmp = 0.0;
      for(int n=0; n<nmaxwell; n++)
      {
        thetatmp = thetatmp+YS[n]*wr/wl[n]/(1+pow(wr/wl[n],2));
      }
      theta2 = thetatmp;
    
      R = sqrt(pow(theta1,2)+pow(theta2,2));
      muunrelax = muu*(R+theta1)/(2*pow(R,2));
    
      md->mu[iptr] = muunrelax;
    
      md->lambda[iptr] = kappa-2*muunrelax;
    
      for(int n=0; n<nmaxwell; n++)
      {
        md->Ylam[n][iptr] = (1+2*md->mu[iptr]/md->lambda[iptr])*YP[n]-2*md->mu[iptr]/md->lambda[iptr]*YS[n];
        md->Ymu[n][iptr]  = YS[n];
        if(md->Ylam[n][iptr] > 1 || md->Ymu[n][iptr] > 1)
        {
          fprintf(stdout,"attention, the coef over the normal range!");
          fprintf(stdout,"Ylam[%d][%d][%d]=%f,Ymu[%d][%d][%d]=%f\n",
                          n,i,k,md->Ylam[n][iptr],n,i,k,md->Ymu[n][iptr]);
        }
      }
    }
  }
  
  free(YP);
  free(YS);
  free(wk);

  fdlib_mem_free_2l_float(GP, kmax, "visco_GP");
  fdlib_mem_free_2l_float(GS, kmax, "visco_GS");
  
  return ierr;

}

int 
md_visco_LS(float **restrict input, float *restrict output, float d, int m, int n)
{
  // G=md,m=(G^TG)^(-1)G^Td

  int ierr = 0;

  float trans[VISCO_LS_MAXSIZE][VISCO_LS_MAXSIZE];
  float multi[VISCO_LS_MAXSIZE][VISCO_LS_MAXSIZE];
  float inver[VISCO_LS_MAXSIZE][VISCO_LS_MAXSIZE];

  for(int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++)
    {
      if(i==j)
          inver[i][j] = 1;
      else
          inver[i][j] = 0;
    }
  }
  
  // transposition
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<m; j++)
    {
      trans[i][j] = input[j][i];
    }
  }

  for(int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++)
    {
      multi[i][j] = 0;
      for(int k=0; k<m; k++)
      {
        multi[i][j] = multi[i][j]+trans[i][k]*input[k][j];
      }
    }
  }

  md_visco_LS_mat_inv(multi,inver,n);

  for(int i=0; i<n; i++)
  {
    for(int j=0; j<m; j++)
    {
      multi[i][j] = 0;
      for(int k=0; k<n; k++)
      {
        multi[i][j] = multi[i][j]+inver[i][k]*trans[k][j];
      }
    }
  }

  float sum;
  for(int i=0; i<n; i++)
  {
    sum = 0.0;
    for(int j=0; j<m; j++)
    {
      sum = sum+multi[i][j]*d;
    }
    output[i] = sum;
  }

  return ierr;
}

int 
md_visco_LS_mat_inv(float matrix[][VISCO_LS_MAXSIZE], float inverse[][VISCO_LS_MAXSIZE], int n)
{
  int ierr = 0;
  
  float tmp;
  int j=0;

  for(int k=0; k<n; k++)
  {
    if(matrix[k][k] == 0)
    {
      for(int jj=k+1; jj<n; jj++)
      {
        j = jj;
        if(matrix[j][k]!=0) break;
      }
      if(j == n)
      {
        fprintf(stderr,"Error: Matrix is not inversible(Cal of visco coef)\n");
        fflush(stderr);
        exit(-1);
      }
      for(int i=0; i<n; i++)
      {
        tmp = matrix[k][i];
        matrix[k][i] = matrix[j][i];
        matrix[j][i] = tmp;
        tmp = inverse[k][i];
        inverse[k][i] = inverse[j][i];
        inverse[j][i] = tmp;
      }
    }
            
    tmp = matrix[k][k];
    for(int j=0; j<n; j++)
    {
      matrix[k][j] = matrix[k][j]/tmp;
      inverse[k][j] = inverse[k][j]/tmp;
    }
    
    for(int i=0; i<n; i++)
    {
      tmp = matrix[i][k];
      for(j=0;j<n;j++)
      {
        if(i==k) break;
        matrix[i][j] = matrix[i][j]-matrix[k][j]*tmp;
        inverse[i][j] = inverse[i][j]-inverse[k][j]*tmp;
      }
    }
  }

  return ierr;
}

int
md_gen_test_vis_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  size_t siz_iz  = md->siz_iz;

  float *lam2d = md->lambda;
  float  *mu2d = md->mu;
  float *rho2d = md->rho;
  float *Qs = md->Qs;
  float *Qp = md->Qp;

  float Vp=3000.0;
  float Vs=2000.0;
  float rho=1500.0;
  float mu = Vs*Vs*rho;
  float lam = Vp*Vp*rho - 2.0*mu;

  for (int k=0; k<nz; k++)
  {
    for (int i=0; i<nx; i++)
    {
      size_t iptr = i + k * siz_iz;
      lam2d[iptr] = lam;
       mu2d[iptr] = mu;
      rho2d[iptr] = rho;
      Qs[iptr] = 40;
      Qp[iptr] = 80;
    }
  }

  return ierr;
}
