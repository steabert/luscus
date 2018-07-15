/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck gveps.h $Revision: 7.7 $ */
#ifndef _GVEPS_H_
#define _GVEPS_H_
/**************************************************************************/
/*                                                                        */
/* Header file for gveps utilities.                                       */
/*                                                                        */
/**************************************************************************/
/*========================================================================*/
/* define's                                                               */
/*========================================================================*/
/* internal */
#define KIND_ATOM    1
#define KIND_BOND    2
#define KIND_TRI     3
#define KIND_TEXT    4
/* for gveps_text */
#define TEXT_TOP     1
#define TEXT_MIDDLE  2
#define TEXT_BOTTOM  3
#define TEXT_LEFT    4
#define TEXT_CENTER  8
#define TEXT_RIGHT  12
/* for gveps_init */
#define SIZE_SMALL   1
#define SIZE_LARGE   2
#define SIZE_HUGE    3
#define DO_BACKG     4
/* for gveps_bond */
#define BOND_DASH    1
/*========================================================================*/
/* typedef's                                                              */
/*========================================================================*/
typedef struct {
   double x,y,z;
   double x1,y1,z1;
   double x2,y2,z2;
   double nx1,ny1,nz1;
   double nx2,ny2,nz2;
   double nx3,ny3,nz3;
   double radius;
   double red,green,blue;
   double red1,green1,blue1;
   int    kind;
   int    dash;
   int    vert,hori;
   char   str[8];
} obj;
typedef struct {
   obj    *X;
   double  red_bg;
   double  green_bg;
   double  blue_bg;
   double  xmax,xmin;
   double  ymax,ymin;
   double  zmax,zmin;
   double  U[3][3];
   int     do_bg;
   int     do_view;
   int     maxobj;
   int     nobj;
   int     size;
} gvepsdef;

/* static gvepsdef ctl; */
/*========================================================================*/
/* external references                                                    */
/*========================================================================*/
/* Basic routines */
extern void gveps_init(FILE *f, int opt, double red, double green, double blue, gvepsdef *ctl);
extern void gveps_atom(FILE *f, int opt, double radius, double x, double y, double z,
   double red, double green, double blue, gvepsdef *ctl);
extern void gveps_bond(FILE *f, int opt,
   double x1, double y1, double z1, double red1, double green1, double blue1,
   double x2, double y2, double z2, double red2, double green2, double blue2, gvepsdef *ctl);
extern void gveps_triangle(FILE *f, int opt,
   double x1, double y1, double z1,
   double x2, double y2, double z2,
   double x3, double y3, double z3,
   double nx1, double ny1, double nz1,
   double nx2, double ny2, double nz2,
   double nx3, double ny3, double nz3,
   double red, double green, double blue,  gvepsdef *ctl);
extern void gveps_text(FILE *f, int opt, char *str, double x, double y, double z, gvepsdef *ctl);
extern void gveps_view(FILE *f, int opt, double U[3][3], gvepsdef *ctl);
extern void gveps_paint(FILE *f, int opt, gvepsdef *ctl);
/* Advanced routines */
extern void gveps_move(FILE *f, int opt, double x, double y, double z, gvepsdef *ctl);
extern void gveps_euler(FILE *f, int opt, double alpha, double beta, double gamma, gvepsdef *ctl);
/* Internal routines */
extern void gveps_allo(gvepsdef *ctl);
extern void gveps_sort(obj X[], int n);
extern void gveps_trim(obj X[], int n);
#endif
