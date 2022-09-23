/*******************************************************************************
 * Staggered Grid Finite Difference Wave Propagation Simulation 
 *
 * Copyright (c) 2021 ZHANG Wei. All rights reserved.
 *
 * Author(s): ZHANG Wei <zhangwei@sustech.edu.cn>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "constants.h"
#include "par_t.h"
#include "blk_t.h"

#include "media_discrete_model.h"
#include "sv_eq1st_cart_stg_ac_iso.h"
#include "sv_eq1st_cart_stg_el_iso.h"

int main(int argc, char** argv)
{
  int verbose = 1; // default fprint
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

  // get commond-line argument

  // argc checking
  if (argc < 2) {
    fprintf(stdout,"usage: main_cart_stg_2d <par_file> <opt: verbose>\n");
    exit(1);
  }

  //strncpy(par_fname, argv[1], sizeof(argv[1]));
  par_fname = argv[1];

  if (argc >= 3) {
    verbose = atoi(argv[2]); // verbose number
    fprintf(stdout,"verbose=%d\n", verbose); fflush(stdout);
  }

  fprintf(stdout,"par file =  %s\n", par_fname); 

  // read par

  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_read_from_file(par_fname, par, verbose);

  if (verbose>0) par_print(par);

//-------------------------------------------------------------------------------
// init blk_t
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"create blk ...\n"); 

  // malloc blk
  blk_t *blk = (blk_t *) malloc(sizeof(blk_t));

  // malloc inner vars
  blk_init(blk, verbose);

  fdstg_t     *fd        = blk->fdstg;
  gdinfo_t    *gdinfo    = blk->gdinfo;
  gd_t        *gdcart    = blk->gd;
  md_t        *md        = blk->md;
  wav_t       *wav       = blk->wav;
  src_t       *src       = blk->src;
  bdryfree_t  *bdryfree  = blk->bdryfree;
  bdrypml_t   *bdrypml   = blk->bdrypml;
  iorecv_t    *iorecv    = blk->iorecv;
  ioline_t    *ioline    = blk->ioline;
  iosnap_t    *iosnap    = blk->iosnap;

  // set up fd_t
  //    not support selection scheme by par file yet
  if (verbose>0) fprintf(stdout,"set scheme ...\n"); 
  fd_set_stg4(fd);

  // set gdinfo
  gd_info_set(gdinfo, 
              par->number_of_total_grid_points_x,
              par->number_of_total_grid_points_z,
              par->abs_num_of_layers,
              fd->fdx_nghosts,
              fd->fdz_nghosts,
              verbose);

  // set str in blk
  blk_set_output(blk, 
                 par->output_dir,
                 par->grid_export_dir,
                 par->media_export_dir,
                 verbose);

//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"allocate grid vars ...\n"); 

  // generate grid coord
  if (par->grid_generation_itype != PAR_GRID_CARTESIAN) {
    fprintf(stderr, "ERROR: grid type wrong\n");
    exit(1);
  }

  // malloc var in gdcart
  float dx = par->cartesian_grid_stepsize[0];
  float dz = par->cartesian_grid_stepsize[1];

  float x0 = par->cartesian_grid_origin[0];
  float z0 = par->cartesian_grid_origin[1];

  gd_cart_init_set(gdinfo,gdcart,dx,x0,dz,z0);

  fprintf(stdout, " --> done\n"); fflush(stdout);

  // output
  if (par->is_export_grid==1)
  {
    if (verbose>0) fprintf(stdout,"export coord to file ...\n"); 
    gd_cart_coord_export(gdinfo, gdcart,
                         blk->grid_export_dir);
  }
  fprintf(stdout, " --> done\n"); fflush(stdout);

//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

  // allocate media vars
  if (verbose>0) fprintf(stdout,"allocate media vars ...\n"); 
  md_init(gdinfo, md, par->media_itype, par->aniso_itype, par->visco_itype);

  // read or discrete velocity model
  switch (par->media_input_itype)
  {
    case PAR_MEDIA_CODE :
        if (verbose>0) fprintf(stdout,"generate simple medium in code ...\n"); 

        if (md_is_el_iso(md)==1) {
            md_gen_test_el_iso(md);
        }

        if (md_is_ac_iso(md)==1) {
          md_gen_test_ac_iso(md);
        }

        break;

    case PAR_MEDIA_IMPORT :
        if (verbose>0) fprintf(stdout,"import discrete medium file ...\n"); 
        md_import(md, par->grid_import_dir);

        break;

    case PAR_MEDIA_3LAY : {
        if (verbose>0) fprintf(stdout,"read and discretize 3D layer medium file ...\n"); 

        float *lam3d = md->lambda;
        float  *mu3d = md->mu;
        float *rho3d = md->rho;
        //media_el_iso_layer2model(lam3d, mu3d, rho3d,
        //                         gdcurv->x3d,
        //                         gdcurv->y3d,
        //                         gdcurv->z3d,
        //                         gdcurv->nx,
        //                         gdcurv->ny,
        //                         gdcurv->nz,
        //                         par->media_input_file,
        //                         par->equivalent_medium_method);
        break;
    }

    case PAR_MEDIA_3GRD : {
        if (verbose>0) fprintf(stdout,"read and descretize 3D grid medium file ...\n"); 

        float *lam3d = md->lambda;
        float  *mu3d = md->mu;
        float *rho3d = md->rho;

        //media_el_iso_grid2model(lam3d, mu3d, rho3d,
        //                        gdcurv->x3d,
        //                        gdcurv->y3d,
        //                        gdcurv->z3d,
        //                        gdcurv->nx,
        //                        gdcurv->ny,
        //                        gdcurv->nz,
        //                        gdcurv->xmin, gdcurv->xmax,   //float Xmin, float Xmax,
        //                        gdcurv->ymin, gdcurv->ymax,   //float Ymin, float Ymax, 
        //                        par->media_input_file,
        //                        par->equivalent_medium_method); 
        break;
    }
  }

  // export grid media
  if (par->is_export_media==1)
  {
    if (verbose>0) fprintf(stdout,"export discrete medium to file ...\n"); 

    md_export(gdinfo, md, blk->media_export_dir);
  } else {
    if (verbose>0) fprintf(stdout,"do not export medium\n"); 
  }

//-------------------------------------------------------------------------------
//-- estimate/check/set time step
//-------------------------------------------------------------------------------

  float   t0 = par->time_start;
  float   dt = par->size_of_time_step;
  int     nt_total = par->number_of_time_steps+1;

  if (par->time_check_stability==1)
  {
    float dtmax, dtmaxVp, dtmaxL;
    int   dtmaxi, dtmaxk;

    //-- estimate time step
    fprintf(stdout,"   estimate time step ...\n"); 
    blk_dt_esti_cart(gdinfo, gdcart,md,fd->CFL,
            &dtmax, &dtmaxVp, &dtmaxL, &dtmaxi, &dtmaxk);
    
    //-- print for QC
    fprintf(stdout, "-> dtmax=%f, Vp=%f, L=%f, i=%d, k=%d\n",
            dtmax, dtmaxVp, dtmaxL, dtmaxi, dtmaxk);
    
    // check valid
    if (dtmax <= 0.0) {
       fprintf(stderr,"ERROR: maximum dt <= 0, stop running\n");
       exit(1);
    }

    //-- auto set stept
    if (dt < 0.0) {
       dt       = blk_keep_two_digi(dtmax);
       nt_total = (int) (par->time_window_length / dt + 0.5);

       fprintf(stdout, "-> Set dt       = %f according to maximum allowed value\n", dt);
       fprintf(stdout, "-> Set nt_total = %d\n", nt_total);
    }

    //-- if input dt, check value
    if (dtmax < dt) {
       fprintf(stdout, "Serious Error: dt=%f > dtmax=%f, stop!\n", dt, dtmax);
       exit(1);
    }
  }

//-------------------------------------------------------------------------------
//-- source import or locate on fly
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"set source using info from par file ...\n"); 
  
  // temp solution for src_set
  float rk_rhs_time[] = { 0.0 };

  src_set_by_par(gdinfo, gdcart, src,
                 t0, dt,
                 1, rk_rhs_time,
                 fd->fdx_max_half_len,
                 par->source_name,
                 par->source_number,
                 par->source_index,
                 par->source_inc,
                 par->source_coords,
                 par->source_force_vector,
                 par->source_force_actived,
                 par->source_moment_tensor,
                 par->source_moment_actived,
                 par->wavelet_name,
                 par->wavelet_coefs,
                 par->wavelet_tstart,
                 par->wavelet_tend,
                 verbose);

  /*
  if (par->is_export_source==1)
  {
      ierr = src_export();
  }
  */

//-------------------------------------------------------------------------------
//-- allocate main var
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"allocate solver vars ...\n"); 
  if (md->medium_type == CONST_MEDIUM_ACOUSTIC)
  {
    wav_ac_init(gdinfo, wav, 1);
  } else
  {
    wav_init(gdinfo, wav, 1);
  }

//-------------------------------------------------------------------------------
//-- setup output, may require coord info
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"setup output info ...\n"); 

  // receiver: need to do
  io_recv_read_locate(gdinfo, gdcart, iorecv,
                      nt_total, wav->ncmp, par->in_station_file);

  // line
  io_line_locate(gdinfo, gdcart, ioline,
                 wav->ncmp,
                 nt_total,
                 par->number_of_receiver_line,
                 par->receiver_line_index_start,
                 par->receiver_line_index_incre,
                 par->receiver_line_count,
                 par->receiver_line_name);
  
  // snapshot
  io_snapshot_locate(gdinfo, iosnap,
                     par->number_of_snapshot,
                     par->snapshot_name,
                     par->snapshot_index_start,
                     par->snapshot_index_count,
                     par->snapshot_index_incre,
                     par->snapshot_time_start,
                     par->snapshot_time_incre,
                     par->snapshot_save_velocity,
                     par->snapshot_save_stress,
                     par->snapshot_save_strain,
                     blk->output_dir);

//-------------------------------------------------------------------------------
//-- absorbing boundary etc auxiliary variables
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"setup absorbingg boundary ...\n"); 
  
  if (par->bdry_has_cfspml == 1)
  {
    bdry_pml_set_stg(gdinfo, gdcart, wav, bdrypml,
                 par->cfspml_is_sides,
                 par->abs_num_of_layers,
                 par->cfspml_alpha_max,
                 par->cfspml_beta_max,
                 par->cfspml_velocity,
                 verbose);
  }

//-------------------------------------------------------------------------------
//-- free surface preproc
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"cal free surface matrix ...\n"); 

  if (par->bdry_has_free == 1)
  {
    bdry_free_set(gdinfo,bdryfree, par->free_is_sides, verbose);
  }

//-------------------------------------------------------------------------------
//-- qc
//-------------------------------------------------------------------------------
  
  //fd_print(fd);

  blk_print(blk);

  gd_info_print(gdinfo);

  iosnap_print(iosnap);

//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------
  
  // convert rho to 1 / rho to reduce number of arithmetic cal
  md_rho_to_slow(md->rho, md->siz_icmp);

  if (verbose>0) fprintf(stdout,"start solver ...\n"); 
  
  time_t t_start = time(NULL);
  
  // compute rhs
  if (md_is_el_iso(md)==1)
  {
    sv_eq1st_cart_stg_el_iso_allstep(fd,gdinfo,gdcart,md,
                            src,bdryfree,bdrypml,
                            wav, 
                            iorecv,ioline,iosnap,
                            dt,nt_total,t0,
                            blk->output_dir,
                            par->check_nan_every_nummber_of_steps,
                            par->output_all,
                            verbose);
  }
  else if (md_is_ac_iso(md)==1)
  {
    sv_eq1st_cart_stg_ac_iso_allstep(fd,gdinfo,gdcart,md,
                            src,bdryfree,bdrypml,
                            wav, 
                            iorecv,ioline,iosnap,
                            dt,nt_total,t0,
                            blk->output_dir,
                            par->check_nan_every_nummber_of_steps,
                            par->output_all,
                            verbose);
  }
  else
  {
      fprintf(stderr,"ERROR: stg solver execpt el_iso and ac_iso not implemented\n");
      exit(1);
  }
  
  time_t t_end = time(NULL);
  
  fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));

//-------------------------------------------------------------------------------
//-- save station and line seismo to sac
//-------------------------------------------------------------------------------

  io_recv_output_sac(iorecv,dt,wav->ncmp,wav->cmp_name,
                      src->evtnm,blk->output_dir,err_message);

  io_line_output_sac(ioline,dt,wav->cmp_name,src->evtnm,blk->output_dir);

//-------------------------------------------------------------------------------
//-- postprocess
//-------------------------------------------------------------------------------

  return 0;
}
