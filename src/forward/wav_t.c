/*********************************************************************
 * wavefield for 2d elastic 1st-order equations
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "fdlib_mem.h"
#include "wav_t.h"

int 
wav_init(gd_t *gd,
         wav_t *V,
         int number_of_levels,
         int visco_type,
         int nmaxwell)
{
  int ierr = 0;

  // Vx,Vz,Txx,Tzz,Txz
  V->ncmp = 5;

  V->nx   = gd->nx;
  V->nz   = gd->nz;
  V->nlevel = number_of_levels;

  V->siz_iz    = V->nx;
  V->siz_icmp  = V->nx * V->nz;

  V->nmaxwell = nmaxwell;
  V->visco_type = visco_type ;

  if(visco_type == CONST_VISCO_GMB)
  {
    V->ncmp += 3*nmaxwell;  // Jxx Jzz Jxz
  }

  V->siz_ilevel = V->siz_icmp * V->ncmp;

  // vars
  // 2 Vi, 3 Tij, 4 rk stages
  V->v4d = (float *) fdlib_mem_calloc_1d_float(V->siz_ilevel * V->nlevel,
                        0.0, "v4d, wf_el2d_1st");
  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                      V->ncmp, 0, "w2d_pos, wf_el2d_1st");
  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(
                      V->ncmp, CONST_MAX_STRLEN, "w2d_name, wf_el2d_1st");
  if (visco_type == CONST_VISCO_GMB)
  {
    V->Jxx_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                nmaxwell, 0, "Jxx_pos, wf_el2d_1st");
    V->Jzz_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                nmaxwell, 0, "Jzz_pos, wf_el2d_1st");
    V->Jxz_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                nmaxwell, 0, "Jxz_pos, wf_el2d_1st");
    V->Jxx_seq = (size_t *) fdlib_mem_calloc_1d_sizet(
                nmaxwell, 0, "Jxx_seq, wf_el2d_1st");
    V->Jzz_seq = (size_t *) fdlib_mem_calloc_1d_sizet(
                nmaxwell, 0, "Jzz_seq, wf_el2d_1st");
    V->Jxz_seq = (size_t *) fdlib_mem_calloc_1d_sizet(
                nmaxwell, 0, "Jxz_seq, wf_el2d_1st");
  }

  // set value
  for (int icmp=0; icmp < V->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * V->siz_icmp;
  }

  // set values
  int icmp = 0;

  /*
   * 0-1: Vx,Vz
   * 2-4: Txx,Tzz,Txz
   */

  sprintf(cmp_name[icmp],"%s","Vx");
  V->Vx_pos = cmp_pos[icmp];
  V->Vx_seq = 0;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vz");
  V->Vz_pos = cmp_pos[icmp];
  V->Vz_seq = 1;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txx");
  V->Txx_pos = cmp_pos[icmp];
  V->Txx_seq = 2;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Tzz");
  V->Tzz_pos = cmp_pos[icmp];
  V->Tzz_seq = 3;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txz");
  V->Txz_pos = cmp_pos[icmp];
  V->Txz_seq = 4;
  icmp += 1;

  if (visco_type == CONST_VISCO_GMB) 
  {
    for(int i=0; i < nmaxwell; i++)
    {
      sprintf(cmp_name[icmp],"%s%d","Jxx",i+1);
      V->Jxx_pos[i] = cmp_pos[icmp];
      V->Jxx_seq[i] = icmp;
      icmp += 1;
  
      sprintf(cmp_name[icmp],"%s%d","Jzz",i+1);
      V->Jzz_pos[i] = cmp_pos[icmp];
      V->Jzz_seq[i] = icmp;
      icmp += 1;
  
      sprintf(cmp_name[icmp],"%s%d","Jxz",i+1);
      V->Jxz_pos[i] = cmp_pos[icmp];
      V->Jxz_seq[i] = icmp;
      icmp += 1;
    }
  }
  // set pointer
  V->cmp_pos  = cmp_pos;
  V->cmp_name = cmp_name;

  return ierr;
}

int 
wav_ac_init(gd_t *gd,
            wav_t *V,
            int number_of_levels)
{
  int ierr = 0;

  // Vx,Vz,P
  V->ncmp = 3;

  V->nx   = gd->nx;
  V->nz   = gd->nz;
  V->nlevel = number_of_levels;

  V->siz_iz   = V->nx;
  V->siz_icmp  = V->nx * V->nz;
  V->siz_ilevel = V->siz_icmp * V->ncmp;

  // vars
  // 2 Vi, 1 P, 4 rk stages
  V->v4d = (float *) fdlib_mem_calloc_1d_float(V->siz_ilevel * V->nlevel,
                        0.0, "v5d, wf_ac2d_1st");
  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                      V->ncmp, 0, "w2d_pos, wf_ac2d_1st");
  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(
                      V->ncmp, CONST_MAX_STRLEN, "w2d_name, wf_ac2d_1st");
  
  // set value
  for (int icmp=0; icmp < V->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * V->siz_icmp;
  }

  // set values
  int icmp = 0;

  /*
   * 0-1: Vx,Vz
   * 2: P
   */

  sprintf(cmp_name[icmp],"%s","Vx");
  V->Vx_pos = cmp_pos[icmp];
  V->Vx_seq = 0;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vz");
  V->Vz_pos = cmp_pos[icmp];
  V->Vz_seq = 1;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","P");
  V->Txx_pos = cmp_pos[icmp];
  V->Txx_seq = 2;
  icmp += 1;

  // set pointer
  V->cmp_pos  = cmp_pos;
  V->cmp_name = cmp_name;

  return ierr;
}

int
wav_check_value(float *restrict w, wav_t *wav)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *ptr = w + icmp * wav->siz_icmp;
    for (size_t iptr=0; iptr < wav->siz_icmp; iptr++)
    {
      if (ptr[iptr] != ptr[iptr])
      {
        fprintf(stderr, "ERROR: NaN occurs at iptr=%d icmp=%d\n", iptr, icmp);
        fflush(stderr);
        exit(-1);
      }
    }
  }

  return ierr;
}

int
wav_zero_edge(gd_t *gd, wav_t *wav, float *restrict w4d)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *restrict var = w4d + wav->cmp_pos[icmp];

    // z1
    for (int k=0; k < gd->nk1; k++)
    {
      size_t iptr_k = k * gd->siz_iz;
        for (int i=0; i < gd->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    }

    // z2
    for (int k=gd->nk2+1; k < gd->nz; k++)
    {
      size_t iptr_k = k * gd->siz_iz;
        for (int i=0; i < gd->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    }

    // x1
    for (int k = gd->nk1; k <= gd->nk2; k++)
    {
      size_t iptr_k = k * gd->siz_iz;
        for (int i=0; i < gd->ni1; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    } 

    // x2
    for (int k = gd->nk1; k <= gd->nk2; k++)
    {
      size_t iptr_k = k * gd->siz_iz;
        for (int i = gd->ni2+1; i < gd->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    } 

  } // icmp

  return ierr;
}
