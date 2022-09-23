#!/bin/bash

#set -x
set -e

date

#-- source intel lib
#source /share/apps/intel-2019.3/intel/bin/compilervars.sh intel64
#export FI_PROVIDER=tcp
#MPIDIR=/share/apps/intel-2019.3/mpich-3.3

#-- program related dir
#EXEC_DIR=/home/zhangw/code/zwlab/CGFD2D-wave
#EXEC_WAVE=$EXEC_DIR/main_curv_col_2d
#EXEC_WAVE=$EXEC_DIR/main_cart_stg_2d
#EXEC_WAVE=$EXEC_DIR/main_cart_col_el_2d
EXEC_WAVE=`pwd`/../main_curv_col_2d

#-- project dir for output
OUTDIR=`pwd`
PROJDIR=${OUTDIR}/00
PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/output
MEDIA_DIR=${PROJDIR}/output
SOURCE_DIR=${PROJDIR}/output
OUTPUT_DIR=${PROJDIR}/output

#-- input dir for input files
INPUTDIR=`pwd`
IN_STATION_LIST_FILE=${INPUTDIR}/recv.list
#IN_GRID_LAY_FILE=/home/zhangw/code/zwlab/CGFD2D-wave/example/grid_hill2d.gdlay
#IN_MEDIA_LAY_FILE=/home/zhangw/code/zwlab/CGFD2D-wave/example/media_hill2d.mdlay
#IN_MEDIA_GRID_FILE=/home/zhangw/code/zwlab/CGFD2D-wave/example/media_hill2d.mdgrd
#IN_SOURCE_FILE=/home/zhangw/code/zwlab/CGFD2D-wave/example/source.anasrc

#-- create dir
mkdir -p $PROJDIR
mkdir -p $OUTPUT_DIR
mkdir -p $GRID_DIR
mkdir -p $MEDIA_DIR

#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $PAR_FILE
{
  "number_of_total_grid_points_x" : 300,
  "number_of_total_grid_points_z" : 201,

  "grid_generation_method" : {
      "#import" : "$GRID_DIR",
      "cartesian" : {
        "origin"  : [ -1000.0  , -4000.0],
        "inteval" : [ 20.0,  20.0 ]
      },
      "#layer_interp" : {
        "in_grid_layer_file" : "$IN_GRID_LAY_FILE",
        "refine_factor" : [ 1, 1, 1 ],
        "horizontal_start_index" : [ 3, 3 ],
        "vertical_ToFreeSurf_resample_index" : 0
      }
  },
  "is_export_grid" : 1,
  "grid_export_dir"   : "$GRID_DIR",

  "metric_calculation_method" : {
      "#import" : "$GRID_DIR",
      "calculate" : 1
  },
  "is_export_metric" : 1,


  "#size_of_time_step" : 0.008,
  "#size_of_time_step" : 0.02,
  "#number_of_time_steps" : 500,
  "time_window_length" : 15,
  "check_stability" : 1,

  "boundary_x_left" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_x_right" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_z_bottom" : {
      "cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          }
      },
  "boundary_z_top" : {
      "#cfspml" : {
          "number_of_layers" : 10,
          "alpha_max" : 3.14,
          "beta_max" : 2.0,
          "ref_vel"  : 5000.0
          },
      "free" : "timg"
      },


  "medium" : {
      "#type" : "elastic_iso",
      "#type" : "elastic_vti",
      "#type" : "acoustic_iso",
      "type" : "elastic_aniso",
      "#import" : "$MEDIA_DIR",
      "#code_generate" : 1,
      "infile_layer" : "$INPUTDIR/prep_medium/test.md2lay",
      "#infile_grid" : "$INPUTDIR/prep_medium/test.md2grd",
      "equivalent_medium_method" : "tti",
      "#in_2grd_file" : "$IN_MEDIA_GRID_FILE"
  },
  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "#visco_config" : {
      "type" : "graves_Qs",
      "Qs_freq" : 1.0
  },

  "source_input" : {
      "in_par" : {
         "name" : "evt_by_par",
         "source" : [
            {
                "index" : [ 100, 97 ],
                "#coord" : [ 4000, 200 ],
                "wavelet_name" : "ricker",
                "ricker_center_frequency" : 10.0,
                "ricker_peak_time" : 1.0,
                "#wavelet_name" : "gaussian",
                "#gaussian_rms_width" : 2.0,
                "#gaussian_peak_time" : 0.5,
                "start_time" : 0.0,
                "end_time"   : 3.0,
                "#force_vector" : [ 0.0, 1e16],
                "moment_tensor" : [ 1e16, 1e16, 0.0]
            }
         ]
      },
      "#in_source_file" : "$IN_SOURCE_FILE"
  },
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "output_dir" : "$OUTPUT_DIR",

  "in_station_file" : "$IN_STATION_LIST_FILE",

  "receiver_line" : [
    {
      "name" : "line_x_1",
      "grid_index_start"    : [   0, 100 ],
      "grid_index_incre"    : [  10,  0 ],
      "grid_index_count"    : 30
    }
  ],

  "snapshot" : [
    {
      "name" : "volume_vel",
      "grid_index_start" : [ 0, 0 ],
      "grid_index_count" : [ 300,101 ],
      "grid_index_incre" : [  1, 1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 0,
      "save_strain"   : 0
    }
  ],

  "check_nan_every_nummber_of_steps" : 0,
  "output_all" : 0 
}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#

#-- gen run script
cat << ieof > ${PROJDIR}/cgfd2d_wave_sim.sh
#!/bin/bash

#-- simulation
printf "\nStart simualtion ...\n";

time $EXEC_WAVE $PAR_FILE 100
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

ieof

#-------------------------------------------------------------------------------
#-- start run
#-------------------------------------------------------------------------------

chmod 755 ${PROJDIR}/cgfd2d_wave_sim.sh
${PROJDIR}/cgfd2d_wave_sim.sh
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

date

#
#-------------------------------------------------------------------------------
#-- plot results
#-------------------------------------------------------------------------------
#

#EXEC_MATLAB=/share/apps/Matlab/R2015b/bin/matlab2015
#EXEC_Plot_Script=./hill3d_plot.sh

#if [[ -n "${Plot_Result}" && ${Plot_Result} -eq 1 ]]; then
#   printf "\nPlot result using %s\n", ${EXEC_Plot_Script}
#   time ${EXEC_Plot_Script}
#else
#   printf "\nDo not plot result\n"
#fi
#
#date

# vim:ft=conf:ts=4:sw=4:nu:et:ai:
