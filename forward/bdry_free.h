#ifndef BDRY_FREE_H
#define BDRY_FREE_H

#include "constants.h"
#include "gd_info.h"
#include "gd_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct
{
  int is_enable; //
  int is_at_sides[CONST_NDIM][2];

  int nx;
  int nz;

  // top
  float *vecVx2Vz2; // [j,i, dzVi, dxVi]

  // bottom
  float *vecVx2Vz1;

  // left
  float *vecVy2Vx1;

  // right
  float *vecVy2Vx2;

  // front
  float *vecVx2Vy1;

  // back
  float *vecVx2Vy2;

} bdryfree_t;

/*************************************************
 * function prototype
 *************************************************/

int
bdry_free_set(gdinfo_t        *gdinfo,
              bdryfree_t      *bdryfree,
              int   in_is_sides[][2],
              const int verbose);

#endif
