################################################################################
#  Makefile for CGFD2D-wave package
#
################################################################################

# known issues:
#  .h dependency is not included, use make cleanall

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CC     := $(GNUHOME)/bin/gcc
CXX    := $(GNUHOME)/bin/g++

#-- 
CFLAGS := -I$(NETCDF)/include -I./src/lib/ -I./src/forward/ -I./src/media/  $(CFLAGS)

#- O3
CFLAGS   := -O3 -std=c99 $(CFLAGS)
CPPFLAGS := -O3 -std=c++11 $(CPPFLAGS)

#- dynamic
LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

skeldirs := obj
DIR_OBJ    := ./obj
#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

OBJS := cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry2d.o \
		media_read_file.o \
		fd_t.o par_t.o \
		gd_t.o md_t.o wav_t.o \
		bdry_t.o src_t.o io_funcs.o \
		blk_t.o interp.o\
		drv_rk_curv_col.o \
		sv_curv_col_el.o \
		sv_curv_col_ac_iso.o \
		sv_curv_col_el_iso.o \
		sv_curv_col_el_vti.o \
		sv_curv_col_el_aniso.o \
		sv_curv_col_vis_iso.o \
		main_curv_col_2d.o

OBJS := $(addprefix $(DIR_OBJ)/,$(OBJS))

vpath  %.cpp ./
vpath  %.c ./

all: skel main
skel:
	@mkdir -p $(skeldirs)

main: $(OBJS)
	$(CXX) $^ $(LDFLAGS) -o $@ 

$(DIR_OBJ)/%.o: src/media/%.cpp
	${CXX} $(CPPFLAGS) -c $^ -o $@ 
$(DIR_OBJ)/%.o: src/lib/%.c
	${CC} $(CFLAGS) -c $^ -o $@
$(DIR_OBJ)/%.o: src/forward/%.c
	${CC} $(CFLAGS) -c $^ -o $@


cleanexe:
	rm -f main

cleanobj:
	rm -f $(DIR_OBJ)/*.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
