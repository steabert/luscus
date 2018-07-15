/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck surface.h $Revision: 7.7 $ */
/*****
Project: Iso-surface drawing

Multi-surface support

*******/


/**** Header file ****/

#ifndef _surface_h_
#define _surface_h_

/* types from ogl.h */
#include<GL/gl.h>
/*#include "gveps.h"*/

typedef GLdouble point_t[3];

typedef struct triangle_t_s
{
  point_t vertex[3];
  point_t normal[3];
  double value[3];/*electrostatic potential value -> maps into a color!*/
  GLdouble distance;
  unsigned char clr_idx; /* color index: 0, 1 */
  struct triangle_t_s* lnk_prev; /* link field for linked list */
} triangle_t;

struct multisrf_t_s;

typedef struct
{
  triangle_t* head; /* list head */
  int n_triangles;

  point_t max_coor,min_coor; /* bounding box: for max. diameter calc. */

  struct multisrf_t_s* parent; /* parent msrf structure */
  int offset; /* index in parent */
} surface_t;

typedef struct multisrf_t_s
{
  triangle_t** index; /* array of pointers to triangles */
  int is_valid;
  int is_sorted;
  int n_triangles; /* total number of triangles */
  
  surface_t** child; /* array of child surfaces */
  int n_surfaces; /* max. number of surfaces */
} multisrf_t;


#define msrf_Set_Unsorted(self) ((self)->is_sorted=0)

multisrf_t* msrf_Init(multisrf_t*, int n_srf);
void msrf_Clean(multisrf_t*);

surface_t* msrf_New_Surface(multisrf_t*);
void msrf_Delete_Surface(surface_t*);

/* add triangle with normals */
triangle_t* srf_Add_Triangle(surface_t*,
                             double txyz[3][3],
                             double tgrd[3][3],
                             double val[3],
                             int clridx);

/* draw all surfaces in current OpenGL context */
void msrf_Draw(multisrf_t*, color_t[],int, FILE *, int, FILE *, gvepsdef *ctl);

double grayscale(double,double,double);

#endif /* _surface_h_ */
