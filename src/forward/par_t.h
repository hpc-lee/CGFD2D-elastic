#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "constants.h"
#include "par_t.h"

#define PAR_MAX_STRLEN 1000
#define PAR_TYPE_STRLEN 50

#define PAR_GRID_IMPORT       1
#define PAR_GRID_CARTESIAN    2

#define PAR_METRIC_CALCULATE 1
#define PAR_METRIC_IMPORT    2

#define PAR_MEDIA_IMPORT 1
#define PAR_MEDIA_CODE   2
#define PAR_MEDIA_2LAY   3
#define PAR_MEDIA_2GRD   4
#define PAR_MEDIA_2BIN   5

#define PAR_MEDIA_CMP_VELOCITY 1
#define PAR_MEDIA_CMP_THOMSEN  2
#define PAR_MEDIA_CMP_CIJ      3

#define PAR_SOURCE_JSON  1
#define PAR_SOURCE_FILE  3

typedef struct{

  //-- dirs and file name
  //char project_dir  [PAR_MAX_STRLEN];
  char output_dir   [PAR_MAX_STRLEN];
  char out_grid_dir     [PAR_MAX_STRLEN];
  char media_dir    [PAR_MAX_STRLEN];
  //char source_dir   [PAR_MAX_STRLEN];
  //char station_dir  [PAR_MAX_STRLEN];
  //char log_file_name[PAR_MAX_STRLEN];
  char tmp_dir[PAR_MAX_STRLEN];

  // time step
  int   number_of_time_steps;
  float size_of_time_step ;
  int   time_start_index;
  int   time_end_index;
  float time_start;
  //float time_end  ;
  int time_check_stability;
  float time_window_length;

  // for each block
  //char grid_name[PAR_MAX_STRLEN];
  
  // grid size
  int  number_of_total_grid_points_x;
  int  number_of_total_grid_points_z;

  // dis grid
  int disg_num_level;
  int *disg_at_zindx;
  int *disg_factor;

  // boundary, CONST_NDIM_2
  char **boundary_type_name;
  
  // abs layer-based, for pml or exp
  int   abs_num_of_layers[CONST_NDIM][2];

  // pml
  int   cfspml_is_sides[CONST_NDIM][2];
  float cfspml_alpha_max[CONST_NDIM][2];
  float cfspml_beta_max[CONST_NDIM][2];
  float cfspml_velocity[CONST_NDIM][2];
  int   bdry_has_cfspml;

  // exp
  int   ablexp_is_sides[CONST_NDIM][2];
  float ablexp_velocity[CONST_NDIM][2];
  int   bdry_has_ablexp;

  // free
  int   free_is_sides[CONST_NDIM][2];
  int   bdry_has_free;

  // grid
  int grid_generation_itype;
  int is_export_grid;

  char grid_export_dir[PAR_MAX_STRLEN];
  char grid_import_dir[PAR_MAX_STRLEN];

  float cartesian_grid_origin[CONST_NDIM];
  float cartesian_grid_stepsize[CONST_NDIM];

  // metric
  int metric_method_itype;
  int is_export_metric;
  char metric_export_dir[PAR_MAX_STRLEN];
  char metric_import_dir[PAR_MAX_STRLEN];

  // medium
  char media_type[PAR_MAX_STRLEN]; // iso, vti, or aniso
  int  media_itype; // iso, vti, or aniso
  char media_input_way[PAR_MAX_STRLEN]; // in_code, import, file
  int  media_input_itype;
  char media_input_cmptype[PAR_MAX_STRLEN]; // cij, thomson
  int  media_input_icmptype;

  int is_export_media;
  char equivalent_medium_method[PAR_MAX_STRLEN]; // For layer2model
  char media_export_dir[PAR_MAX_STRLEN];
  char media_import_dir[PAR_MAX_STRLEN];
  char media_input_file[PAR_MAX_STRLEN];

  // visco
  char visco_type[PAR_MAX_STRLEN]; // graves_Qs
  int  visco_itype; // graves_Qs
  float visco_Qs_freq;
  int nmaxwell;
  float fmax;
  float fmin;
  float fr;

  // source
  //int source_input_itype;
  char source_input_file[PAR_MAX_STRLEN];
  int is_export_source;
  char source_export_dir[PAR_MAX_STRLEN];

  char source_dd_input_file[PAR_MAX_STRLEN];
  int  source_dd_add_at_point;
  int  source_dd_nt_per_read;

  // output
  // receiver
  char in_station_file[PAR_MAX_STRLEN];
  // line
  int number_of_receiver_line;
  int *receiver_line_index_start;
  int *receiver_line_index_incre;
  int *receiver_line_count;
  //int *receiver_line_time_interval;
  char **receiver_line_name;
  // snapshot
  int number_of_snapshot;
  char **snapshot_name;
  int *snapshot_index_start;
  int *snapshot_index_count;
  int *snapshot_index_incre;
  int *snapshot_time_start;
  //int *snapshot_time_count; // should output to end 
  int *snapshot_time_incre;
  int *snapshot_save_velocity;
  int *snapshot_save_stress;
  int *snapshot_save_strain;

  // misc
  int check_nan_every_nummber_of_steps;
  int output_all;
} par_t;

int
par_read_from_file(char *par_fname, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int 
par_read_json_cfspml(cJSON *item,
      int *nlay, float *amax, float *bmax, float *vel);

int 
par_read_json_ablexp(cJSON *item, int *nlay, float *vel);

int
par_print(par_t *par);

#endif

