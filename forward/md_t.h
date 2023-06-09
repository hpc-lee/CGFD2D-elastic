#ifndef MD_EL_ISO_H
#define MD_EL_ISO_H

#include "gd_t.h"
#include "par_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int n1, n2, n3, n4;
  int nx, nz, ncmp;
  float *v3d; // allocated var

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

  // flag to determine medium type
  int medium_type;
  int visco_type;

  // rho for all media
  float *rho;

  // for acustic
  float *kappa; // pointer to var

  // for isotropic media
  float *lambda; // pointer to var
  float *mu;

  // for visco attenuation
  int nmaxwell;
  float visco_GMB_freq;
  float visco_GMB_fmin;
  float visco_GMB_fmax;
  float *Qs;
  float *Qp;
  float **Ylam; // float *Ylam[nmaxwell];
  float **Ymu;  // float *Ymu[nmaxwell];
  float *wl;

  // for anisotropic media
  float *c11;
  float *c13;
  float *c15;
  float *c33;
  float *c35;
  float *c55;

  float visco_Qs_freq;

} md_t;

/*************************************************
 * function prototype
 *************************************************/

int
md_init(gd_t *gd, md_t *md, int media_type, int visco_type, int nmaxwell);

int
md_import(md_t *md, char *in_dir);

int
md_export(gd_t *gd,
          md_t *md,
          char *output_dir);

int
md_gen_test_ac_iso(md_t *md);

int
md_gen_test_el_iso(md_t *md);

int
md_gen_test_el_vti(md_t *md);

int
md_gen_test_el_aniso(md_t *md);

int
md_gen_test_Qs(md_t *md, float Qs_freq);

int
md_gen_test_GMB(md_t *md);

int
md_rho_to_slow(float *restrict rho, size_t siz_volume);

int 
md_vis_GMB_cal_Y(md_t *md, float freq, float fmin, float fmax);

int 
md_visco_LS(float **restrict input, float *restrict output, float d, int m, int n);

int 
md_visco_LS_mat_inv(float matrix[][VISCO_LS_MAXSIZE], float inverse[][VISCO_LS_MAXSIZE], int n);

#endif
