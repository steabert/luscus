/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*fordeck findsym.h $Revision: 7.7 Patch(7.7): 279_gv $ */
#ifndef __FIND_SYM_H__
#define __FIND_SYM_H__
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
#include <stdio.h>
#define VERSION "0.11"
#ifdef _IN_GATE_
#ifdef _CAPITALS_
#define findsym     FINDSYM 
#else
#ifndef ADD_
#define findsym findsym_
#endif
#endif
#else
#define INT int
#endif

/*--------------------------------------------------------------------------*/
/* Type definitions.                                                        */
/*--------------------------------------------------------------------------*/
typedef struct {
   int     n_atoms;
   int     *Z;
   char  **label;
   char  *point_group;
   double *r[3];
   double *mass;
   double MoI[4];
} nuclear_frame;
/*--------------------------------------------------------------------------*/
/* Global variables.                                                        */
/*--------------------------------------------------------------------------*/

int	read_input_findsym(int, int, char [][4],  FILE*, double [][3], nuclear_frame *, double uthr);
void	print_output_findsym(int, FILE *, nuclear_frame *);
void	scramble_findsym(nuclear_frame *);
void	CM_adjust_findsym(nuclear_frame *);
void	inertia_adjust_findsym(nuclear_frame *);
void	jacobi_findsym(double *, double *, int, int);
int	find_point_group(int, FILE*, nuclear_frame *,char*, char*, int[],int);
void	Cn_axis_findsym(double *, double [3][3], int);
void	mirror_plane_findsym(double *, double [3][3], int);
void	Unit_findsym(double [3][3]);
void	Invert_findsym(double [3][3]);
void	transform_findsym(double [3][3], nuclear_frame *);
double  match_findsym(double [3][3], nuclear_frame *);
int	keep(nuclear_frame *, nuclear_frame *);
void	recover(nuclear_frame *, nuclear_frame *);
int	freeframe(nuclear_frame *);
int find_point_group_unscrambled(int IN_GATE, FILE *ou, nuclear_frame *x, char *what, char *outgroup, int outsym[],int skip_print);

#endif
