/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*fordeck findsym.c $Revision: 7.7 Patch(7.7): 279_gv 571_sync $ */
/****************************************************************************/
/*                                                                          */
/*                              F i n d S y m                               */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* This program determines the point group of a molecule. The full group or */
/* a subgroup of D2h is deternined. The molecule is aligned along principal */
/* axes.                                                                    */
/*                                                                          */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* Author:  Per-Olof Widmark                                                */
/*          Dept. of Theoretical Chemistry                                  */
/*          Lund University                                                 */
/*          Sweden                                                          */
/* Modified: Valera Veryazov                                                */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* History: 0.1 VV: check also symmetry for non scrambled structure         */
/*          0.11 VV: remove MAXATOM limits                                  */
/*                                                                          */
/****************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "findsym.h"
#include "vstring.h"
#ifdef _IN_GATE_
#include "molcastype.h"
#endif

#define MY_PI 3.14159265358979323846e0

double THR_FINDSYM=0.01e0;
typedef struct {
   int   n_gen;
   int   gen[3];
   char group[8];
} point_group;

nuclear_frame new_coords;

/*#ifdef _IN_GATE_*/
void printlog(char * str)
{ fprintf(stderr,"%s",str); return; }
/*#endif*/
/*==========================================================================*/
/*  findsym works different in gv and Gateway.
in gv input is given via parameters, in gate - read from external file.

list of parameters:
natoms: number of atoms (gv)
atomNames: atomic labels (gv)
xyz: coordinates (gv)
FileIn: input file (gate)
lenFileIn: length of FileIn (gate)
what: string to set up task
lenwhat: length of what (gate)
utht: threshold
scramble: flag to scramble the structure (gv)
outgroup: group (gv)
outsym: (gv)
ierr: (gate)



*/
/*==========================================================================*/

INT  findsym(INT IN_GATE, INT natoms, char atomNames[][4], double xyz[][3],
             char *FileIn, INT *lenFileIn,
             char *what, INT *lenwhat,
             double *uthr, INT scramble,
	     char *outgroup, INT outsym[], INT *ierr ) {
   nuclear_frame x;
   nuclear_frame xinit;
   nuclear_frame xxx;
   FILE *iu, *ou;
   char *FileIn0,*what0;
   int i, j;
   if(*uthr==0.0) *uthr=THR_FINDSYM;
if(IN_GATE)
{
   FileIn0=(char*) malloc( (*lenFileIn+1)*sizeof(char));
   strncpy(FileIn0,FileIn, *lenFileIn);
   FileIn0[*lenFileIn]=0;

   what0=(char*) malloc( (*lenwhat+1)*sizeof(char));
   strncpy(what0,what, *lenwhat);
   what0[*lenwhat]=0;


   iu=fopen(FileIn0,"r");
   if(iu==NULL) {fprintf (stderr,"Can't open file %s\n",FileIn0); *ierr=1; return 0;}
   free(FileIn0);

   ou=fopen("findsym.out","w");
}
else
{
  what0=0;
  iu=0;
  ou=0;
}
/*--------------------------------------------------------------------------*/
/* Read input structure.                                                    */
/*--------------------------------------------------------------------------*/
   if(read_input_findsym(IN_GATE,natoms,atomNames,iu,xyz,&x,*uthr)) {
      printlog("Aborting after reading input\n");
      *ierr=1; return 0;
   }
#ifdef DEBUG
   fprintf(stdout,"Input structure\n");
   print_output_findsym(IN_GATE,stdout,&x);
   fprintf(stdout,"\n");
#endif
   if(keep(&x,&xinit)) {
      printlog("malloc problems..\n");
      *ierr=1; return 0;
   }
/*--------------------------------------------------------------------------*/
/* Optionally scramble the input structure.                                 */
/*--------------------------------------------------------------------------*/
/*   if(scramble) scramble_findsym(&x); */
#ifdef DEBUG
   fprintf(stdout,"Scrambled structure\n");
   print_output_findsym(IN_GATE,stdout,&x);
   fprintf(stdout,"\n");
#endif
/*--------------------------------------------------------------------------*/
/* Adjust center of mass.                                                   */
/*--------------------------------------------------------------------------*/
   CM_adjust_findsym(&x);
#ifdef DEBUG
   fprintf(stdout,"CM adjusted structure\n");
   print_output_findsym(IN_GATE,stdout,&x);
   fprintf(stdout,"\n");
#endif
/*--------------------------------------------------------------------------*/
/* Adjust moment of inertia.                                                */
/*--------------------------------------------------------------------------*/
   inertia_adjust_findsym(&x);
#ifdef DEBUG
   fprintf(stdout,"Moment of inertia adjusted structure\n");
   print_output_findsym(IN_GATE,stdout,&x);
   fprintf(stdout,"\n");
#endif
/*--------------------------------------------------------------------------*/
/* Align an atom with a coordinate axis if possible/desirable               */
/* Move this code to a separate file.                                       */
/*--------------------------------------------------------------------------*/
   {
      double U[3][3];
      double r[3];
      double trMoI;
      double r2,r2max,r2near;
      double x0,y0;
      double x1,y1;
      double p,q;
      double q1,q2;
      double g;
      double xthr;
      int    debug;
      int    adjust;
      int    nequiv;
      int    kfar;
      int    knear;
      int    ix,iy,iz;
      int    n[3],try[3];
      int    k,m;

      r2near=0.0;
      debug=0;
      xthr=1.0e-3;
#ifdef DEBUG
      debug=1;
#endif

      for(k=0; k<3; k++) n[k]=0;
      for(k=0; k<3; k++) try[k]=0;
      trMoI=x.MoI[3];
      if(trMoI<1.0) trMoI=1.0;
      if(fabs(x.MoI[0]-x.MoI[1])/trMoI<xthr && fabs(x.MoI[1]-x.MoI[2])/trMoI<xthr) {
         if(debug) printf("This is a qubic point group!\n");
         try[0]=1;
         try[1]=1;
         try[2]=1;
      } else if(fabs(x.MoI[0]-x.MoI[1])/trMoI<xthr) {
         if(debug) printf("z-axis is highest axis!\n");
         try[2]=1;
      } else if(fabs(x.MoI[0]-x.MoI[2])/trMoI<xthr) {
         if(debug) printf("y-axis is highest axis!\n");
         try[1]=1;
      } else if(fabs(x.MoI[1]-x.MoI[2])/trMoI<xthr) {
         if(debug) printf("x-axis is highest axis!\n");
         try[0]=1;
      } else {
         if(debug) printf("All three axes have different eigenvalues!\n");
      }
      for(k=0; k<3; k++) if(try[k]) {
#ifdef DEBUG
         if(k==0)      printf("Checking x-axis!\n");
         else if(k==1) printf("Checking y-axis!\n");
         else if(k==2) printf("Checking z-axis!\n");
#endif
         r[0]=0.0e0; r[1]=0.0e0; r[2]=0.0e0;
         r[k]=1.0e0;
         Unit_findsym(U);
         n[k]=0;
         for(m=2; m<9; m++) {
            Cn_axis_findsym(r,U,m);
            q=match_findsym(U,&x);
            if(fabs(q)<*uthr) n[k]=m;
#ifdef DEBUG
            if(fabs(q)<*uthr) printf("Is a %i-fold axis, q=%.3f\n",m,q);
            else              printf("Is not a %i-fold axis, q=%.3f\n",m,q);
#endif
         }
#ifdef DEBUG
         if(n[k]>0) printf("Highest axis is %i\n",n[k]);
#endif
      }
      adjust=0;
      if(n[0]>0 && n[1]==0 && n[2]==0) {
         if(debug) printf("Adjust around x-axis!\n");
         ix=1;
         iy=2;
         iz=0;
         adjust=1;
      } else if(n[1]>0 && n[2]==0 && n[0]==0) {
         if(debug) printf("Adjust around y-axis!\n");
         ix=2;
         iy=0;
         iz=1;
         adjust=1;
      } else if(n[2]>0 && n[0]==0 && n[1]==0) {
         if(debug) printf("Adjust around z-axis!\n");
         ix=0;
         iy=1;
         iz=2;
         adjust=1;
      }
      if(adjust) {
         if(debug) printf("Adjustment indices: %i, %i\n",ix,iy);
         if(debug) printf("Axis order: %i\n",n[iz]);
         r2max=0.0;
         kfar=0;
         for(k=0; k<x.n_atoms; k++) {
            r2=x.r[ix][k]*x.r[ix][k]+x.r[iy][k]*x.r[iy][k];
            if(r2>r2max) { r2max=r2; kfar=k; }
         }
         if(debug) printf("Atom %i is %.3f away from axis!\n",kfar,sqrt(r2max));
         if(n[iz]==2 && r2max>0.1) { /* 2-fold axis mean accidental degeneracy! */
            if(debug) printf("Axis order is 2, do something else\n");
            nequiv=0;
            for(k=0; k<x.n_atoms; k++) {
               r2=x.r[ix][k]*x.r[ix][k]+x.r[iy][k]*x.r[iy][k];
               if(fabs(r2-r2max)<*uthr) nequiv=nequiv+1;
            }
            if(debug) printf("%i equivalent atoms\n",nequiv);
            if(nequiv==2) {
               if(debug) printf("2 equivalent atoms, must be on one axis!\n");
            } else {
               if(debug) printf("4 equivalent atoms, can be on one axis, try both!\n");
               knear=kfar;
               r2near=1.0e12;
               for(k=0; k<x.n_atoms; k++) {
                  r2=x.r[ix][k]*x.r[ix][k]+x.r[iy][k]*x.r[iy][k];
                  if(fabs(r2-r2max)<*uthr && k!=kfar) {
                     r2=(x.r[0][kfar]-x.r[0][k])*(x.r[0][kfar]-x.r[0][k])
                       +(x.r[1][kfar]-x.r[1][k])*(x.r[1][kfar]-x.r[1][k])
                       +(x.r[2][kfar]-x.r[2][k])*(x.r[2][kfar]-x.r[2][k]);
                     if(debug) printf("[%i] %.6f\n",k,r2);
                     if(r2<r2near) { r2near=r2; knear=k; }
                  }
               }
               if(debug) printf("kfar=%i, knear=%i\n",kfar,knear);
               if(debug) printf("Try aligning midpoint with axis and see if mirrors are ok\n");
               x0=0.5*(x.r[ix][kfar]+x.r[ix][knear]);
               y0=0.5*(x.r[iy][kfar]+x.r[iy][knear]);
               if(fabs(x0)<fabs(y0)) { g=x0/y0; p=g/sqrt(1.0+g*g); q=1.0/sqrt(1.0+g*g); }
               else                  { g=y0/x0; p=1.0/sqrt(1.0+g*g); q=g/sqrt(1.0+g*g); }
               for(k=0; k<x.n_atoms; k++) {
                  x1 = p*x.r[ix][k]+q*x.r[iy][k];
                  y1 =-q*x.r[ix][k]+p*x.r[iy][k];
                  x.r[ix][k]=x1;
                  x.r[iy][k]=y1;
               }
               r[0]=0.0; r[1]=0.0; r[2]=0.0;
               r[ix]=1.0;
               Unit_findsym(U); mirror_plane_findsym(r,U,2); q1=match_findsym(U,&x);
               Unit_findsym(U); Cn_axis_findsym(r,U,2); q2=match_findsym(U,&x);
               if(debug) printf("q1=%.3f, q2=%.3f\n",q1,q2);
               if(q1<*uthr) {
                  if(debug) printf("Did find mirror plane, q1=%.3f\n",q1);
               } else if(q2<*uthr) {
                  if(debug) printf("Did find C2 axis, q2=%.3f\n",q2);
               } else {
                  if(debug) printf("Did not find mirror plane, q1=%.3f, q2=%.3f\n",q1,q2);
                  x0=x.r[ix][kfar];
                  y0=x.r[iy][kfar];
                  if(fabs(x0)<fabs(y0)) { g=x0/y0; p=g/sqrt(1.0+g*g); q=1.0/sqrt(1.0+g*g); }
                  else                  { g=y0/x0; p=1.0/sqrt(1.0+g*g); q=g/sqrt(1.0+g*g); }
                  for(k=0; k<x.n_atoms; k++) {
                     x1 = p*x.r[ix][k]+q*x.r[iy][k];
                     y1 =-q*x.r[ix][k]+p*x.r[iy][k];
                     x.r[ix][k]=x1;
                     x.r[iy][k]=y1;
                  }
                  r[0]=0.0; r[1]=0.0; r[2]=0.0;
                  r[ix]=1.0;
                  Unit_findsym(U); mirror_plane_findsym(r,U,2); q1=match_findsym(U,&x);
                  Unit_findsym(U); Cn_axis_findsym(r,U,2); q2=match_findsym(U,&x);
                  if(debug) printf("q1=%.3f, q2=%.3f\n",q1,q2);
                  if(q1<*uthr) {
                     if(debug) printf("Did find mirror plane second time, q1=%.3f\n",q1);
                  } else if(q2<*uthr) {
                     if(debug) printf("Did find C2 axis second time, q2=%.3f\n",q2);
                  } else {
                     if(debug) printf("Did not find mirror plane second time, q1=%.3f, q2=%.3f\n",q1,q2);
                  }
               }
            }
         } else if(r2max>0.1) { /* Not 2-fold axis */
            x0=x.r[ix][kfar];
            y0=x.r[iy][kfar];
            if(fabs(x0)<fabs(y0)) {
               g=x0/y0;
               p=g/sqrt(1.0+g*g);
               q=1.0/sqrt(1.0+g*g);
#ifdef DEBUG
               printf("p/q=x0/y0=%.3f, p=%.3f, q=%.3f\n",g,p,q);
#endif
            } else {
               g=y0/x0;
               p=1.0/sqrt(1.0+g*g);
               q=g/sqrt(1.0+g*g);
#ifdef DEBUG
               printf("q/p=y0/x0=%.3f, p=%.3f, q=%.3f\n",g,p,q);
#endif
            }
            for(k=0; k<x.n_atoms; k++) {
               x1 = p*x.r[ix][k]+q*x.r[iy][k];
               y1 =-q*x.r[ix][k]+p*x.r[iy][k];
               x.r[ix][k]=x1;
               x.r[iy][k]=y1;
            }
         }
      }
   }
#ifdef DEBUG
   fprintf(stdout,"Align on axis adjusted structure\n");
   print_output_findsym(IN_GATE,stdout,&x);
   fprintf(stdout,"\n");
#endif
/*--------------------------------------------------------------------------*/
/* Find pointgroup.                                                         */
/*--------------------------------------------------------------------------*/
  if(IN_GATE)
  {
   find_point_group(IN_GATE,ou,&x,what0,outgroup,outsym,0);
  }
  else
  {
keep(&x, &xxx);
   if(find_point_group(IN_GATE,ou,&x,what,outgroup,outsym,1))
   {
     return 1;
   }


  }
#ifdef DEBUG
   fprintf(stdout,"Point group adjusted structure\n");
   print_output_findsym(IN_GATE,stdout,&x);
   fprintf(stdout,"\n");
#endif
/*--------------------------------------------------------------------------*/
/* Find unscrambled pointgroup.                                             */
/*--------------------------------------------------------------------------*/
   recover(&x,&xinit);
if(IN_GATE)
{
   find_point_group_unscrambled(IN_GATE,ou,&x,what0,outgroup,outsym,0);
}
else
{ int rank_m=0, rank_f=0;
/* VV:: here we have to remember our selection, and if needed use unscrambled structure
      FIXME: this is very durty solution, by recalling symmetry analyzer twice.

*/
   rank_m=mycount(outgroup,' ');
   find_point_group_unscrambled(IN_GATE,stderr,&x,what,outgroup,outsym,1);
   rank_f=mycount(outgroup,' ');
   if(rank_m>rank_f)
   {
     fprintf(stderr,"Symmetry can be gained by transformation\n");
     recover(&x,&xxx);
     find_point_group(IN_GATE,ou,&x,what,outgroup,outsym,0);
     freeframe(&xxx);
   }
   else
   {
   recover(&x,&xinit);
   find_point_group_unscrambled(IN_GATE,stderr,&x,what,outgroup,outsym,0);
   }
}
/*--------------------------------------------------------------------------*/
/* Clean up!      .                                                         */
/*--------------------------------------------------------------------------*/
 freeframe(&x);
 freeframe(&xinit);

   /*copying coordinates back -Goran*/
     for(i = 0; i < x.n_atoms; i++)
       for(j = 0; j < 3; j++)
         xyz[i][j] = new_coords.r[j][i];
      freeframe(&new_coords);
if(IN_GATE)
{
   free(what0);
   ierr=0;
   fclose(ou);
   return 0;
}
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   return 0;
}
int  read_input_findsym(int IN_GATE, int natoms, char atomNames[][4], FILE *io, double xyz[][3],
nuclear_frame *x, double uthr) {

   char   line[256],*ptr;
   double  r1,r2,r3;
   double Bohr=1.0e0;

   int     i,n, i3;
   int  iv, im, nMendeleev, ll;
   char Mendeleev[][3]={"X","H", "He",
            "Li","Be","B", "C", "N", "O", "F", "Ne",
            "Na","Mg","Al","Si","P", "S", "Cl","Ar",
            "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co","Ni",
            "Cu","Zn","Ga","Ge","As","Se","Br","Kr",
            "Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd",
            "Ag","Cd","In","Sn","Sb","Te","I", "Xe",
            "Cs","Ba","La","Ce","Pr","Nb","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                    "Hf","Ta","W", "Re","Os","Ir","Pt",
            "Au","Hg","Tl","Pb","Bi","Po","At","Rn",
            "Fr","Ra","Ac","Th","Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"
             };

   nMendeleev=sizeof(Mendeleev)/(sizeof(char)*3);
if(IN_GATE)
{
   if(fgets(line,sizeof(line),io)==NULL) {
      printlog("Premature end of input file\n");
      return 1;
   }

   if( sscanf(line,"%d",&natoms)!= 1 )
        { printlog("Can't read Number of atom \n"); return 1;}
   if(fgets(line,sizeof(line),io)==NULL)
     {
       printlog("Error during reading comment line in XYZ file\n"); return 1;
     }
     mystrlwr(line);
     if(strstr(line,"a.u.")!=NULL) Bohr=0.52917721067e0;
     if(strstr(line,"bohr")!=NULL) Bohr=0.52917721067e0;

}

   if(uthr<1.0e-6) uthr=1.0e-6;
   if(uthr>1.0e0)  uthr=1.0e0;
   THR_FINDSYM=uthr;
/*==========================================================================*/
/*   Allocate memory                                                        */
/*==========================================================================*/
   n=natoms;

#ifdef TRACE_MALLOC
 puts("allocate x->Z");
#endif
   if( (x->Z=(int  *) malloc(sizeof(int) * n))==NULL)
            return 1;

   for (i3=0; i3<3; i3++)
    {

   if( (x->r[i3]=(double *) malloc(sizeof(double) * n))==NULL)
            return 1;
    }
   if( (x->label=(char **) malloc(sizeof(char*) * n))==NULL)
            return 1;
   if( (x->mass=(double *) malloc(sizeof(double) * n))==NULL)
            return 1;

   if( (x->point_group=(char *) malloc(sizeof(char) * 64))==NULL)
            return 1;

   x->n_atoms=0;

if(IN_GATE)
{
  for (n=0; n<natoms; n++)
  {
   if(fgets(line,sizeof(line),io)==NULL)
      {
      printlog("Premature end of xyz file\n");
      return 1;
      }
   ptr=strchr(line,' ');
   if(ptr!=NULL)
   {
         ll=ptr-line;
         x->label[n]=(char *) malloc((ll+1)*sizeof(char));
         strncpy(x->label[n],line,ll);
         x->label[n][ll]=0;

      if( sscanf(ptr,"%lf %lf %lf", &r1, &r2, &r3) != 3)
           return 1;

           iv=strlen(x->label[n]);
           for (i=0; i<iv; i++)
           {
             if(isalpha(x->label[n][i])==0) {iv=i; break; }
           }
             if(iv==0) {printlog("Input file is not complete\n"); return 1; }
              x->label[n][0]=toupper(x->label[n][0]);
           if(iv>1)
              x->label[n][1]=tolower(x->label[n][1]);
           for(im=0; im<nMendeleev; im++)
           {
           if(strncmp(x->label[n], Mendeleev[im], iv)==0) { x->Z[n]=im; break;}
           }

                   x->r[0][n]=r1*Bohr;
                   x->r[1][n]=r2*Bohr;
                   x->r[2][n]=r3*Bohr;
      x->mass[n]=(double)(x->Z[n]);

   }
   }


}
else
{
   for(i=0; i<natoms; i++)
   {

         x->label[i]=(char *) malloc(sizeof(char)*4);
         strcpy(x->label[i],atomNames[i]);
         iv=(int)strlen(x->label[i]);
           for(im=0; im<nMendeleev; im++)
           {
           if(strncmp(x->label[i], Mendeleev[im],(unsigned int)iv)==0) { x->Z[i]=im; break;}
           }


                   x->r[0][i]=xyz[i][0];
                   x->r[1][i]=xyz[i][1];
                   x->r[2][i]=xyz[i][2];
      x->mass[i]=(double)(x->Z[i]);

   }
}
   for(i=0; i<4; i++) x->MoI[i]=-1.0e0;
   x->point_group=NULL;
   x->n_atoms=natoms;


   return 0;
}

/*==========================================================================*/
void CM_adjust_findsym(nuclear_frame *x) {
   double T[3],M;
   int     i,k;
/*--------------------------------------------------------------------------*/
/* Find center of mass.                                                     */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) T[i]=0.0e0;
   for(k=0; k<x->n_atoms; k++) {
      for(i=0; i<3; i++) T[i]=T[i]+x->mass[k]*x->r[i][k];
   }
   M=0.0e0;
   for(k=0; k<x->n_atoms; k++) M=M+x->mass[k];
   for(i=0; i<3; i++) T[i]=T[i]/M;
#ifdef DEBUG
   fprintf(stdout,"Center of mass\n");
   for(i=0; i<3; i++) fprintf(stdout,"%12.6f ",T[i]); fprintf(stdout,"\n");
   fprintf(stdout,"\n");
#endif
/*--------------------------------------------------------------------------*/
/* Perform translation.                                                     */
/*--------------------------------------------------------------------------*/
   for(k=0; k<x->n_atoms; k++) {
      for(i=0; i<3; i++) x->r[i][k]=x->r[i][k]-T[i];
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
void Cn_axis_findsym(double r[3], double U[3][3], int  n) {
   double u[3][3],R[3][3],C[3][3];
   double t,c,s;
   int  i,j,k,m;
#ifdef DEBUG
   printf("Cn_axis: Rotation axis\n");
   for(i=0; i<3; i++) printf("%12.6f ",r[i]); printf("\n\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct orthonormal directions.                                        */
/*--------------------------------------------------------------------------*/
/* First direction = rotation axis */
   for(i=0; i<3; i++) u[i][0]=r[i];
/* Second direction */
   k=0;
   for(i=0; i<3; i++) if(fabs(u[i][0])<fabs(u[k][0])) k=i;
   if(k==0) {
      u[0][1]= 0.0e0;
      u[1][1]=-u[2][0];
      u[2][1]= u[1][0];
   } else if(k==1) {
      u[0][1]=-u[2][0];
      u[1][1]= 0.0e0;
      u[2][1]= u[0][0];
   } else if(k==2) {
      u[0][1]=-u[1][0];
      u[1][1]= u[0][0];
      u[2][1]= 0.0e0;
   }
/* Third direction */
   u[0][2]= u[1][0]*u[2][1]-u[2][0]*u[1][1];
   u[1][2]=-u[0][0]*u[2][1]+u[2][0]*u[0][1];
   u[2][2]= u[0][0]*u[1][1]-u[1][0]*u[0][1];
/* Normalize */
   for(k=0; k<3; k++) {
      t=0.0e0;
      for(i=0; i<3; i++) t=t+u[i][k]*u[i][k];
      t=1.0e0/sqrt(t);
      for(i=0; i<3; i++) u[i][k]=t*u[i][k];
   }
#ifdef DEBUG
   printf("Cn_axis: Orthonormal directions\n");
   for(i=0; i<3; i++) {
      for(k=0; k<3; k++) printf("%12.6f ",u[i][k]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct core transformation matrix.                                    */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         C[i][j]=0.0e0;
      }
   }
   c=cos(2.0e0*MY_PI/(double)(n));
   s=sin(2.0e0*MY_PI/(double)(n));
   C[0][0]=1.0e0;
   C[1][1]= c;
   C[1][2]=-s;
   C[2][1]= s;
   C[2][2]= c;
#ifdef DEBUG
   printf("Cn_axis: Core transformation matrix\n");
   for(i=0; i<3; i++) {
      for(k=0; k<3; k++) printf("%12.6f ",C[i][k]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct transformation matrix.                                         */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         t=0.0e0;
         for(k=0; k<3; k++) {
            for(m=0; m<3; m++) {
               t=t+u[i][k]*C[k][m]*u[j][m];
            }
         }
         R[i][j]=t;
      }
   }
#ifdef DEBUG
   printf("Cn_axis: Transformation matrix\n");
   for(i=0; i<3; i++) {
      for(k=0; k<3; k++) printf("%12.6f ",R[i][k]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Copy local transformation matrix to parameter.                           */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         U[i][j]=R[i][j];
      }
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
/*==========================================================================*/
/* Symmetrize                                                               */
/*==========================================================================*/
int symmetrize_d2h(point_group *p, nuclear_frame *full, nuclear_frame *reduced,
                   int outsym[], int IN_GATE) {
   double  sx,sy,sz;
   double  ax,ay,az;
   double  ri[3],rj[3];
   double  r2;
   double  thr;
   int     n_op,op[8];
   int     i,j,k,m;
   double *rr0, *rr1, *rr2;
   int isfound;
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   thr=THR_FINDSYM;
/*--------------------------------------------------------------------------*/
/* Create operator list.                                                    */
/*--------------------------------------------------------------------------*/
   op[0]=0; n_op=1;
   for(k=0; k<p->n_gen; k++) {
      for(i=0; i<n_op; i++) op[n_op+i]=op[i]^p->gen[k];
      n_op=2*n_op;
   }
/*--------------------------------------------------------------------------*/
/* Symmetrize coordinate list.                                              */
/*--------------------------------------------------------------------------*/
   for(k=0; k<n_op; k++) {
      if((op[k]&1)==1) sx=-1.0e0;
      else             sx= 1.0e0;
      if((op[k]&2)==2) sy=-1.0e0;
      else             sy= 1.0e0;
      if((op[k]&4)==4) sz=-1.0e0;
      else             sz= 1.0e0;
      for(i=0; i<full->n_atoms; i++) {
         ri[0]=sx*full->r[0][i];
         ri[1]=sy*full->r[1][i];
         ri[2]=sz*full->r[2][i];
         for(j=0; j<full->n_atoms; j++) {
            rj[0]=full->r[0][j];
            rj[1]=full->r[1][j];
            rj[2]=full->r[2][j];
            r2=0.0e0; for(m=0; m<3; m++) r2=r2+(ri[m]-rj[m])*(ri[m]-rj[m]);
            if(r2<thr*thr) {
               if(full->Z[i]!=full->Z[j]) {
                  printlog("Internal error in symmetrize_d2h\n");
                  fprintf(stderr,"Atoms %i and %i are not equivalent\n",i,j);
                  return 1;
               }
               ax=0.5e0*(ri[0]+rj[0]);
               ay=0.5e0*(ri[1]+rj[1]);
               az=0.5e0*(ri[2]+rj[2]);
               full->r[0][i]=sx*ax;
               full->r[1][i]=sy*ay;
               full->r[2][i]=sz*az;
               full->r[0][j]=   ax;
               full->r[1][j]=   ay;
               full->r[2][j]=   az;
            }
         }
      }
   }
/*--------------------------------------------------------------------------*/
/* Set all small values to zero.                                            */
/*--------------------------------------------------------------------------*/
   for(i=0; i<full->n_atoms; i++) {
      for(m=0; m<3; m++) {
         if(fabs(full->r[m][i])<1.0e-8) full->r[m][i]=0.0e0;
      }
   }
/*--------------------------------------------------------------------------*/
/* Compute reduced coordinate list.                                         */
/*--------------------------------------------------------------------------*/
   *reduced=*full;

   rr0=(double *) malloc(sizeof(double) * full->n_atoms);
   rr1=(double *) malloc(sizeof(double) * full->n_atoms);
   rr2=(double *) malloc(sizeof(double) * full->n_atoms);
   for(i=0; i<full->n_atoms; i++) {
    rr0[i]=full->r[0][i];
    rr1[i]=full->r[1][i];
    rr2[i]=full->r[2][i];
   }


   for(k=0; k<n_op; k++) {
      if((op[k]&1)==1) sx=-1.0e0;
      else             sx= 1.0e0;
      if((op[k]&2)==2) sy=-1.0e0;
      else             sy= 1.0e0;
      if((op[k]&4)==4) sz=-1.0e0;
      else             sz= 1.0e0;
      for(i=0; i<reduced->n_atoms; i++) {
         ri[0]=sx*reduced->r[0][i];
         ri[1]=sy*reduced->r[1][i];
         ri[2]=sz*reduced->r[2][i];
         for(j=i+1; j<reduced->n_atoms; j++) {
            rj[0]=reduced->r[0][j];
            rj[1]=reduced->r[1][j];
            rj[2]=reduced->r[2][j];
            r2=0.0e0; for(m=0; m<3; m++) r2=r2+(ri[m]-rj[m])*(ri[m]-rj[m]);
            if(r2<1.0e-6) {
               if(reduced->Z[i]!=reduced->Z[j]) {
                  printlog("Internal error in symmetrize_d2h\n");
                  fprintf(stderr,"Atoms %i and %i are not equivalent\n",i,j);
                  return 1;
               }
               for(m=j; m<reduced->n_atoms-1; m++) {
                  reduced->r[0][m]=reduced->r[0][m+1];
                  reduced->r[1][m]=reduced->r[1][m+1];
                  reduced->r[2][m]=reduced->r[2][m+1];
                  reduced->Z[m]=reduced->Z[m+1];
                  reduced->mass[m]=reduced->mass[m+1];
                  reduced->label[m]=reduced->label[m+1];
               }
               reduced->n_atoms--;
            }
         }
      }
   }

/* here we have to compare full and reduce lists */
   for(i=0;i<full->n_atoms; i++)
   {
     isfound=0;
      for(j=0; j<reduced->n_atoms; j++)
      {
        if(fabs(reduced->r[0][j]-rr0[i])<1.0e-6 &&
           fabs(reduced->r[1][j]-rr1[i])<1.0e-6 &&
           fabs(reduced->r[2][j]-rr2[i])<1.0e-6)
           isfound=1;
      }
      if(IN_GATE==0 && isfound==0) outsym[i]=0;
   }
  free(rr0);
  free(rr1);
  free(rr2);

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
return 0;
}


/*==========================================================================*/
/*                                                                          */
/* Find subgroups of D2h for molecules with noncubic pointgroups.           */
/*                                                                          */
/*==========================================================================*/
int find_point_group_2(int IN_GATE, FILE *ou, nuclear_frame *x, unsigned int TM, char *what, char *outgroup, int outsym[],
int skip_print) {
   char   *oplab[]={"E","x","y","xy","z","xz","yz","xyz"};
   char   *optyp[]={"e","s","s","c" ,"s","c", "c", "i"};
   point_group   p[16];
   nuclear_frame full,reduced;
   double  U[3][3];
   double  r[3];
   double  thr,err,q;
   int      op[8];
   int      op1;
   int      gen[3];
   int      n_op,n_gen,n_group,kase;
   int      i,j,k,m;
   int      requested_n_gen=0;
   char     r_gen[3][4];
   int      found=0,im;
   int      user_input=0;
   int      foundsome=0;
   FILE *fsym;
/*--------------------------------------------------------------------------*/
/* Initialize                                                               */
/*--------------------------------------------------------------------------*/
   n_op=0; op[n_op]=0; n_op++;
   n_gen=0;
   thr=THR_FINDSYM;
#ifdef DEBUG
   printf("# thr=%.6f\n",thr);
#endif
/*--------------------------------------------------------------------------*/
/* Make a decision                                                          */
/*--------------------------------------------------------------------------*/
/*   FIXME!!!!!!!!!!
    if(strchr(what,'*')==NULL && TM=='*') return 0;
*/
if(IN_GATE)
{
    mydiet(what);
    mystrlwr(what);
}
    if(strcmp(what,"full")!=0 && strcmp(what,"nosym")!=0 && strcmp(what,"list")!=0)
    {
    requested_n_gen=mycount(what,' ')+1;
/*    if(requested_n_gen==1) --  easy case only one element to check
      if(requested_n_gen==3) --  easy case - assuming D2h
      if(requested_n_gen==2) --  lets calculate missing generator
*/
      user_input=1;
      if(requested_n_gen==2)
        {
          mytoken(what,' ',r_gen[0],r_gen[1]);
          j=0;
          if(mycount(what,'x')==1) { r_gen[2][j]='x'; j++;}
          if(mycount(what,'y')==1) { r_gen[2][j]='y'; j++;}
          if(mycount(what,'z')==1) { r_gen[2][j]='z'; j++;}
          r_gen[2][j]=0;
        }
/*   printf(">%s< >%s< >%s< \n",r_gen[0],r_gen[1],r_gen[2]);	 
*/
   }
/*--------------------------------------------------------------------------*/
/* Find generators and group elements.                                      */
/*--------------------------------------------------------------------------*/
   for(k=0; k<=6; k++) {
      if(k==0) {
         r[0]=0.0e0; r[1]=0.0e0; r[2]=1.0e0;
         op1=3;
         kase=0;
      } else if(k==1) {
         r[0]=0.0e0; r[1]=1.0e0; r[2]=0.0e0;
         op1=5;
         kase=0;
      } else if(k==2) {
         r[0]=1.0e0; r[1]=0.0e0; r[2]=0.0e0;
         op1=6;
         kase=0;
      } else if(k==3) {
         r[0]=1.0e0; r[1]=0.0e0; r[2]=0.0e0;
         op1=7;
         kase=2;
      } else if(k==4) {
         r[0]=0.0e0; r[1]=0.0e0; r[2]=1.0e0;
         op1=4;
         kase=1;
      } else if(k==5) {
         r[0]=0.0e0; r[1]=1.0e0; r[2]=0.0e0;
         op1=2;
         kase=1;
      } else if(k==6) {
         r[0]=1.0e0; r[1]=0.0e0; r[2]=0.0e0;
         op1=1;
         kase=1;
      } else {
         printlog("Internal error in find_point_group_2\n");
         fprintf(stderr,"k=%i\n",k);
         return 1;
      }
      Unit_findsym(U);
      if(kase==0)      Cn_axis_findsym(r,U,2);
      else if(kase==1) mirror_plane_findsym(r,U,2);
      else if(kase==2) Invert_findsym(U);
      else {
         printlog("Internal error in find_point_group_2\n");
         fprintf(stderr,"kase=%i\n",kase);
         return 1;
      }

      if(match_findsym(U,x)<thr) {
         m=-1; for(i=0; i<n_op; i++) if(op[i]==op1) m=i;
         if(m==-1) {
            for(i=0; i<n_op; i++) op[i+n_op]=op[i]^op1; n_op=2*n_op;
            gen[n_gen]=op1; n_gen++;
         }
      }
   }
#ifdef DEBUG
   printf("# n_op=%i ",n_op);
   for(k=0; k<n_op; k++) printf("[%s] ",oplab[op[k]]); printf("\n");
   printf("# n_gen=%i ",n_gen);
   for(k=0; k<n_gen; k++) printf("[%s] ",oplab[gen[k]]); printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Do we have a point group?                                                */
/*--------------------------------------------------------------------------*/
   for(i=0; i<n_op; i++) {
      for(j=0; j<n_op; j++) {
         k=-1;
         for(m=0; m<n_op; m++) if((op[i]^op[j])==op[m]) k=m;
         if(k==-1) {
            printlog("Internal error in find_point_group_2\nNot a point group\n");
            return 1;
         }
      }
   }
/*--------------------------------------------------------------------------*/
/* Generate all subgroups.                                                  */
/*--------------------------------------------------------------------------*/
   n_group=0;
   if(n_gen==3) {
      p[n_group].n_gen=3;
      p[n_group].gen[0]=gen[0];
      p[n_group].gen[1]=gen[1];
      p[n_group].gen[2]=gen[2];
      n_group++;

      p[n_group].n_gen=2; p[n_group].gen[0]=op[1]; p[n_group].gen[1]=op[2]; n_group++;
      p[n_group].n_gen=2; p[n_group].gen[0]=op[1]; p[n_group].gen[1]=op[4]; n_group++;
      p[n_group].n_gen=2; p[n_group].gen[0]=op[1]; p[n_group].gen[1]=op[6]; n_group++;
      p[n_group].n_gen=2; p[n_group].gen[0]=op[2]; p[n_group].gen[1]=op[4]; n_group++;
      p[n_group].n_gen=2; p[n_group].gen[0]=op[2]; p[n_group].gen[1]=op[5]; n_group++;
      p[n_group].n_gen=2; p[n_group].gen[0]=op[3]; p[n_group].gen[1]=op[4]; n_group++;
      p[n_group].n_gen=2; p[n_group].gen[0]=op[3]; p[n_group].gen[1]=op[5]; n_group++;

      p[n_group].n_gen=1; p[n_group].gen[0]=op[1]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[2]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[3]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[4]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[5]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[6]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[7]; n_group++;

      p[n_group].n_gen=0; n_group++;
   }
   if(n_gen==2) {
      p[n_group].n_gen=2;
      p[n_group].gen[0]=gen[0];
      p[n_group].gen[1]=gen[1];
      n_group++;

      p[n_group].n_gen=1; p[n_group].gen[0]=op[1]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[2]; n_group++;
      p[n_group].n_gen=1; p[n_group].gen[0]=op[3]; n_group++;

      p[n_group].n_gen=0; n_group++;
   }
   if(n_gen==1) {
      p[n_group].n_gen=1;
      p[n_group].gen[0]=gen[0];
      n_group++;

      p[n_group].n_gen=0; n_group++;
   }
   if(n_gen==0) {
      p[n_group].n_gen=0; n_group++;
   }
/*--------------------------------------------------------------------------*/
/* Characterize all groups.                                                 */
/*--------------------------------------------------------------------------*/
   for(k=0; k<n_group; k++) {
      strcpy(p[k].group,"<group>");
   }
   for(k=0; k<n_group; k++) {
      if(p[k].n_gen==3) {
         strcpy(p[k].group,"D2h");
      } else if(p[k].n_gen==2) {
         if(strcmp(optyp[p[k].gen[0]],"c")==0) {
            if(strcmp(optyp[p[k].gen[1]],"c")==0) {
               strcpy(p[k].group,"D2");
            } else if(strcmp(optyp[p[k].gen[1]],"s")==0) {
               strcpy(p[k].group,"C2v");
            } else if(strcmp(optyp[p[k].gen[1]],"i")==0) {
               strcpy(p[k].group,"C2h");
            }
         }
      } else if(p[k].n_gen==1) {
         if(strcmp(optyp[p[k].gen[0]],"c")==0) {
            strcpy(p[k].group,"C2");
         } else if(strcmp(optyp[p[k].gen[0]],"s")==0) {
            strcpy(p[k].group,"Cs");
         } else if(strcmp(optyp[p[k].gen[0]],"i")==0) {
            strcpy(p[k].group,"Ci");
         }
      } else if(p[k].n_gen==0) {
         strcpy(p[k].group,"C1");
      }
   }
/*--------------------------------------------------------------------------*/
/* Print all passing point groups.                                          */
/*--------------------------------------------------------------------------*/
   for(k=0; k<n_group; k++) {
      op[0]=0; n_op=1;
      for(m=0; m<p[k].n_gen; m++) {
         for(i=0; i<n_op; i++) op[n_op+i]=op[i]^p[k].gen[m]; n_op=2*n_op;
      }
      err=0.0;
      for(m=1; m<n_op; m++) {
         if(strcmp(oplab[op[m]],"xy")==0)  { r[0]=0.0e0; r[1]=0.0e0; r[2]=1.0e0; Unit_findsym(U); Cn_axis_findsym(r,U,2); }
         if(strcmp(oplab[op[m]],"xz")==0)  { r[0]=0.0e0; r[1]=1.0e0; r[2]=0.0e0; Unit_findsym(U); Cn_axis_findsym(r,U,2); }
         if(strcmp(oplab[op[m]],"yz")==0)  { r[0]=1.0e0; r[1]=0.0e0; r[2]=0.0e0; Unit_findsym(U); Cn_axis_findsym(r,U,2); }
         if(strcmp(oplab[op[m]],"z")==0)   { r[0]=0.0e0; r[1]=0.0e0; r[2]=1.0e0; Unit_findsym(U); mirror_plane_findsym(r,U,2); }
         if(strcmp(oplab[op[m]],"y")==0)   { r[0]=0.0e0; r[1]=1.0e0; r[2]=0.0e0; Unit_findsym(U); mirror_plane_findsym(r,U,2); }
         if(strcmp(oplab[op[m]],"x")==0)   { r[0]=1.0e0; r[1]=0.0e0; r[2]=0.0e0; Unit_findsym(U); mirror_plane_findsym(r,U,2); }
         if(strcmp(oplab[op[m]],"xyz")==0) { r[0]=0.0e0; r[1]=0.0e0; r[2]=0.0e0; Unit_findsym(U); Invert_findsym(U); }
         q=match_findsym(U,x);
#ifdef DEBUG
         printf("# err[%s]=%.6f\n",oplab[op[m]],q);
#endif
         if(q>err) err=q;
      }
      if(err<thr) {
         /* if NOSYM - skip it until C1 is found */
         if(strcmp(what,"nosym")==0 && p[k].n_gen!=0) continue;

              if(strcmp(what,"full")!=0 && strcmp(what,"nosym")!=0)
              {
                   /* Waiting for D2h                     */

                   if(requested_n_gen==3 && p[k].n_gen!=3) continue;

                   if(requested_n_gen==1 && p[k].n_gen!=1) continue;
                   if(requested_n_gen==1 && strcmp(what,oplab[p[k].gen[0]])!=0 ) continue;

                   /* most complicated case: 2 generators must match */
                   if(requested_n_gen==2 && p[k].n_gen!=2) continue;
                   if(requested_n_gen==2)
                   {
                     im=0;
                     if(strcmp(oplab[p[k].gen[0]],r_gen[0])==0) im++;
                     if(strcmp(oplab[p[k].gen[0]],r_gen[1])==0) im++;
                     if(strcmp(oplab[p[k].gen[0]],r_gen[2])==0) im++;
                     if(strcmp(oplab[p[k].gen[1]],r_gen[0])==0) im++;
                     if(strcmp(oplab[p[k].gen[1]],r_gen[1])==0) im++;
                     if(strcmp(oplab[p[k].gen[1]],r_gen[2])==0) im++;

                     if(im!=2) continue;
                   }
               }
if(IN_GATE)
{
         fprintf(ou,"%c %s ",TM,p[k].group);
         fprintf(ou,"%.6f\n",err);
         fprintf(ou,"%i ",p[k].n_gen);
}
         if(strcmp(what,"full")==0 && p[k].n_gen==3)
          {
if(IN_GATE)
{
          fprintf(ou,"x y z ");
}
          strcpy(outgroup,"x y z");
          foundsome=1;
          }
         else
          {
            if(user_input)
              {
if(IN_GATE)
{
                fprintf(ou,"%s ",what);
}
else
{
                fprintf(stderr,"%s ",what);
}
              }
            else
              {
              outgroup[0]=0;
                for(i=0; i<p[k].n_gen; i++)
                     {
if(IN_GATE)
{
                      fprintf(ou,"%s ",oplab[p[k].gen[i]]);
}
                     strcat(outgroup,oplab[p[k].gen[i]]);
                     mystrcatl(outgroup,' ');
                     foundsome=1;
                     }
              }
          }
if(IN_GATE)
{
         fprintf(ou,"\n");
}
         found=1;
         keep(x,&full);
         keep(x,&reduced);
         if(symmetrize_d2h(&p[k],&full,&reduced,outsym,IN_GATE)) return 1;

if(IN_GATE)
{
         fprintf(ou,"%i\n",reduced.n_atoms);
         print_output_findsym(IN_GATE,ou,&reduced);
}
if(!IN_GATE)
{

         if(foundsome==0)
           {
             fprintf(stderr,"Sorry, no symmetry element has been found\n");
	    /*copy coordinates here! -Goran*/
	    keep(x, &new_coords);
             return 0;
           }
  if(skip_print==0)
  {
         fprintf(stderr,"Reduced coordinates\n");
         fprintf(stderr,"%i\nGroup generators: %s\n",reduced.n_atoms, outgroup);
         print_output_findsym(IN_GATE,stderr,&reduced);
  }
         if(foundsome>0)
           {
            keep(x,&full);
#ifdef DEBUG
            fsym= fopen("symmetry.xyz","w");
            fprintf(fsym,"%i\nGroup generators: %s\n",full.n_atoms, outgroup);
            print_output_findsym(IN_GATE,fsym,&full);
            fclose(fsym);
#endif
	    /*copy coordinates here! -Goran*/
	    keep(x, &new_coords);
           }
}
         if(strcmp(what,"list")!=0) return 0;

/*
         freeframe(&reduced);
         freeframe(&full);
*/
      }
   }
   return 1;
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
/* End of find_point_group_2 */
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
int find_point_group(int IN_GATE, FILE *ou, nuclear_frame *x,char *what, char *outgroup, int outsym[], int skip_print) {
   return find_point_group_2(IN_GATE, ou, x,'+',what,outgroup,outsym,skip_print);
}

/*==========================================================================*/
void inertia_adjust_findsym(nuclear_frame *x) {
   double I[6],u[9],U[3][3];
   double ri,rj,r2,g;
   int     i,j,k,ind;
/*--------------------------------------------------------------------------*/
/* Construct moment of inertia.                                             */
/*--------------------------------------------------------------------------*/
   for(i=0; i<6; i++) I[i]=0.0e0;
   for(k=0; k<x->n_atoms; k++) {
      r2=0.0e0;
      for(i=0; i<3; i++) {
         ri=x->r[i][k];
         r2=r2+ri*ri;
      }
      for(i=0; i<3; i++) {
         ri=x->r[i][k];
         for(j=0; j<=i; j++) {
            rj=x->r[j][k];
            if(i==j) g=ri*rj-r2;
            else     g=ri*rj;
            I[(i*(i+1))/2+j]=I[(i*(i+1))/2+j]-x->mass[k]*g;
         }
      }
   }
#ifdef DEBUG
   printf("Moment of inertia\n");
   for(i=0; i<6; i++) if(fabs(I[i])<1.0e-8) I[i]=0.0e0;
   for(i=0; i<3; i++) {
      for(j=0; j<=i; j++) printf("%12.6f ",I[(i*(i+1))/2+j]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Diagonalize moment of inertia.                                           */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         if(i==j) u[3*i+j]=1.0e0;
         else     u[3*i+j]=0.0e0;
      }
   }
   jacobi_findsym(I,u,3,3);
#ifdef DEBUG
   printf("Diagonalized moment of inertia\n");
   for(i=0; i<6; i++) if(fabs(I[i])<1.0e-8) I[i]=0.0e0;
   for(i=0; i<3; i++) {
      for(j=0; j<=i; j++) printf("%12.6f ",I[(i*(i+1))/2+j]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Copy eigenvectors to transformation matrix.                              */
/*--------------------------------------------------------------------------*/
   ind=0;
   for(j=0; j<3; j++) {
      for(i=0; i<3; i++) {
         U[j][i]=u[ind];
         ind++;
      }
   }
#ifdef DEBUG
   printf("Rotation matrix\n");
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) printf("%12.6f ",U[i][j]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Transform.                                                               */
/*--------------------------------------------------------------------------*/
   transform_findsym(U,x);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
   g=0.0e0;
   for(i=0; i<3; i++) {
      j=i;
      x->MoI[i]=I[(i*(i+1))/2+j];
      g=g+I[(i*(i+1))/2+j];
   }
   x->MoI[3]=g;
#ifdef DEBUG
   printf("Moments of inertia\n");
   for(i=0; i<4; i++) printf("%12.6f ",x->MoI[i]); printf("\n");
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
/* this is currently only in GATE we need to record outgroup2???? */
int find_point_group_unscrambled(int IN_GATE, FILE *ou, nuclear_frame *x, char *what, char *outgroup, int outsym[],int skip_print) {

   return find_point_group_2(IN_GATE, ou,x,'*',what,outgroup,outsym,skip_print);

}

/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void Invert_findsym(double U[3][3]) {
   int     i,j;
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) U[i][j]=-U[i][j];
   }
}
void jacobi_findsym(double f[], double u[], int  n, int  nr) {
   int  i,j,k;
   int  ii,jj,ij;
   int  ki,kj;
   double eps,err;
   double fii,fjj,fij;
   double diff,sign;
   double d,t,c,s,tmp;
   eps=1.0e-20;
   err=1.0e0;
   while (err>eps) {
      err=0.0e0;
      for (i=1; i<n; i++) {
         for (j=0; j<i; j++) {
            ii=i+i*(i+1)/2;
            jj=j+j*(j+1)/2;
            ij=j+i*(i+1)/2;
            fii=f[ii];
            fjj=f[jj];
            fij=f[ij];
            tmp=fij;
            if(fij<0) tmp=-fij;
            if(tmp>eps) {
               if(tmp>err) err=tmp;
               diff=fii-fjj;
               sign=1.0e0;
               if (diff<0) {
                  diff=-diff;
                  sign=-sign;
               }
               d=diff+sqrt(diff*diff+4*fij*fij);
               t=2*sign*fij/d;
               c=1/sqrt(1+t*t);
               s=c*t;
/* rotate u */
               for (k=0; k<nr; k++) {
                  d=-s*u[k+i*nr]+c*u[k+j*nr];
                  u[k+i*nr]=c*u[k+i*nr]+s*u[k+j*nr];
                  u[k+j*nr]=d;
               }
/* rotate f with on index in i or j */
               for (k=0; k<n; k++) {
                  if ( (k!=i) || (k!=j) ) {
                     ki=k+i*(i+1)/2;
                     kj=k+j*(j+1)/2;
                     if (k>i) ki=i+k*(k+1)/2;
                     if (k>j) kj=j+k*(k+1)/2;
                     d=c*f[kj]-s*f[ki];
                     f[ki]=s*f[kj]+c*f[ki];
                     f[kj]=d;
                  }
               }
               f[ii]=c*c*fii+s*s*fjj+2*c*s*fij;
               f[jj]=s*s*fii+c*c*fjj-2*c*s*fij;
               f[ij]=0;
            }
         }
      }
   }
}
/*==========================================================================*/
double match_findsym(double U[3][3], nuclear_frame *x) {
   double        dr,r2,r2min,r2max;
   int            i,ka,kb,n;
   double *newr[3];
   double V[3];
   int     j,k;
   n=x->n_atoms;
/*
    for(ka=0; ka<n; ka++)
            printf(">>%d %f %f %f\n",ka, x->r[0][ka],x->r[1][ka],x->r[2][ka]);
*/

   if( (newr[0]=(double *) malloc(sizeof(double) * n))==NULL)
            return -1;
   if( (newr[1]=(double *) malloc(sizeof(double) * n))==NULL)
            return -1;
   if( (newr[2]=(double *) malloc(sizeof(double) * n))==NULL)
            return -1;

   for(k=0; k<n; k++) {
      for(i=0; i<3; i++) {
         V[i]=0.0e0; for(j=0; j<3; j++) V[i]=V[i]+U[i][j]*x->r[j][k];
      }
      for(i=0; i<3; i++) newr[i][k]=V[i];
   }

   r2max=0.0e0;
   for(ka=0; ka<n; ka++) {
      r2min=1.0e6;
      for(kb=0; kb<n; kb++) {
         if(x->Z[ka]==x->Z[kb]) {
            r2=0.0e0;
            for(i=0; i<3; i++) {
               dr=x->r[i][ka]-newr[i][kb];
               r2=r2+dr*dr;
            }
            if(r2<r2min)r2min=r2;
         }
      }
      if(r2min>r2max) r2max=r2min;
   }



   free(newr[0]);
   free(newr[1]);
   free(newr[2]);
   return sqrt(r2max);
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
void mirror_plane_findsym(double r[3], double U[3][3], int  n) {
   double u[3][3],R[3][3],C[3][3];
   double t;
   int  i,j,k,m;
#ifdef DEBUG
   printf("mirror_plane: Rotation axis\n");
   for(i=0; i<3; i++) printf("%12.6f ",r[i]); printf("\n\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct orthonormal directions.                                        */
/*--------------------------------------------------------------------------*/
/* First direction = rotation axis */
   for(i=0; i<3; i++) u[i][0]=r[i];
/* Second direction */
   k=0;
   for(i=0; i<3; i++) if(fabs(u[i][0])<fabs(u[k][0])) k=i;
   if(k==0) {
      u[0][1]= 0.0e0;
      u[1][1]=-u[2][0];
      u[2][1]= u[1][0];
   } else if(k==1) {
      u[0][1]=-u[2][0];
      u[1][1]= 0.0e0;
      u[2][1]= u[0][0];
   } else if(k==2) {
      u[0][1]=-u[1][0];
      u[1][1]= u[0][0];
      u[2][1]= 0.0e0;
   }
/* Third direction */
   u[0][2]= u[1][0]*u[2][1]-u[2][0]*u[1][1];
   u[1][2]=-u[0][0]*u[2][1]+u[2][0]*u[0][1];
   u[2][2]= u[0][0]*u[1][1]-u[1][0]*u[0][1];
/* Normalize */
   for(k=0; k<3; k++) {
      t=0.0e0;
      for(i=0; i<3; i++) t=t+u[i][k]*u[i][k];
      t=1.0e0/sqrt(t);
      for(i=0; i<3; i++) u[i][k]=t*u[i][k];
   }
#ifdef DEBUG
   printf("mirror_plane: Orthonormal directions\n");
   for(i=0; i<3; i++) {
      for(k=0; k<3; k++) printf("%12.6f ",u[i][k]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct core transformation matrix.                                    */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         C[i][j]=0.0e0;
      }
   }
   C[0][0]=-1.0e0;
   C[1][1]= 1.0e0;
   C[2][2]= 1.0e0;
#ifdef DEBUG
   printf("mirror_plane: Core transformation matrix\n");
   for(i=0; i<3; i++) {
      for(k=0; k<3; k++) printf("%12.6f ",C[i][k]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct transformation matrix.                                         */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         t=0.0e0;
         for(k=0; k<3; k++) {
            for(m=0; m<3; m++) {
               t=t+u[i][k]*C[k][m]*u[j][m];
            }
         }
         R[i][j]=t;
      }
   }
#ifdef DEBUG
   printf("mirror_plane: Transformation matrix\n");
   for(i=0; i<3; i++) {
      for(k=0; k<3; k++) printf("%12.6f ",R[i][k]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Copy local transformation matrix to parameter.                           */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         U[i][j]=R[i][j];
      }
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}


int  keep(nuclear_frame *x, nuclear_frame *xinit) {
int  i3,i,n,k;

    n=x->n_atoms;

  if( (xinit->Z=(int  *) malloc(sizeof(int) * n))==NULL)
            return 1;
  if( (xinit->mass=(double *) malloc(sizeof(double) * n))==NULL)
            return 1 ;

   if( (xinit->point_group=(char *) malloc(sizeof(char) * 64))==NULL)
            return 1;

   for (i3=0; i3<3; i3++)
    {

   if( (xinit->r[i3]=(double *) malloc(sizeof(double) * n))==NULL)
            return 1;
    }
   if( (xinit->label=(char **) malloc(sizeof(char*) * n))==NULL)
            return 1;

   for (i=0; i<n; i++)
   if( (xinit->label[i]=(char *) malloc( sizeof(char) * (strlen(x->label[i])+1)))==NULL)
            return 1;


   for(i=0; i<n; i++)
   {
   for(k=0; k<3; k++)
     {
       xinit->r[k][i]=x->r[k][i];
     }
     strcpy(xinit->label[i],x->label[i]);
   xinit->Z[i]=x->Z[i];
   xinit->mass[i]=x->mass[i];

   }
   xinit->n_atoms=n;
   for (i=0; i<4; i++)
    xinit->MoI[i]=x->MoI[i];


 return 0;

}
void recover(nuclear_frame *x, nuclear_frame *xinit) {

int  n,i,k;
    n=xinit->n_atoms;
   for(i=0; i<n; i++)
   {
   for(k=0; k<3; k++)
     {
       x->r[k][i]=xinit->r[k][i];
     }
     strcpy(x->label[i],xinit->label[i]);
   x->Z[i]=xinit->Z[i];
   x->mass[i]=xinit->mass[i];

   }
   x->n_atoms=n;

 return;
}

int  freeframe(nuclear_frame *x) {
int  i3,n,i;
   for (i3=0; i3<3; i3++)
      free(x->r[i3]);

   n=x->n_atoms;
   for (i=0; i<n; i++)
    free(x->label[i]);
   free(x->Z);
   free(x->mass);
   free(x->point_group);

 return 0;
}

#define DISORDER
/*==========================================================================*/
/*                                                                          */
/*==========================================================================*/
void scramble_findsym(nuclear_frame *x) {
   double U[3][3],T[3],V[3];
   double s;
   int     i,j,k;
/*--------------------------------------------------------------------------*/
/* Construct rotation matrix.                                               */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         U[i][j]=(double)rand();
      }
   }
   for(i=0; i<3; i++) {
      s=0.0e0; for(k=0; k<3; k++) s=s+U[k][i]*U[k][i];
      s=1.0e0/sqrt(s);
      for(k=0; k<3; k++) U[k][i]=s*U[k][i];
      for(j=i+1; j<3; j++) {
         s=0.0e0; for(k=0; k<3; k++) s=s+U[k][j]*U[k][i];
         for(k=0; k<3; k++) U[k][j]=U[k][j]-s*U[k][i];
      }
   }
#ifdef DEBUG
   printf("Rotation matrix\n");
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) printf("%12.6f ",U[i][j]); printf("\n");
   }
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Construct translation vector.                                            */
/*--------------------------------------------------------------------------*/
   for(i=0; i<3; i++) {
      T[i]=(double)rand();
   }
   s=0.0e0; for(k=0; k<3; k++) s=s+T[k]*T[k];
   s=1.0e0/sqrt(s);
   for(k=0; k<3; k++) T[k]=s*T[k];
   for(k=0; k<3; k++) T[k]=1.0e0*T[k];
#ifdef DEBUG
   printf("Translation vector\n");
   for(i=0; i<3; i++) printf("%12.6f ",T[i]); printf("\n");
   printf("\n");
#endif
/*--------------------------------------------------------------------------*/
/* Perform transformation.                                                  */
/*--------------------------------------------------------------------------*/
   for(k=0; k<x->n_atoms; k++) {
      for(i=0; i<3; i++) {
         V[i]=T[i]; for(j=0; j<3; j++) V[i]=V[i]+U[i][j]*x->r[j][k];
      }
      for(i=0; i<3; i++) {
         x->r[i][k]=V[i];
      }
   }
/*--------------------------------------------------------------------------*/
/* Add random disorder.                                                     */
/*--------------------------------------------------------------------------*/
#ifdef DISORDER
   for(k=0; k<x->n_atoms; k++) {
      for(i=0; i<3; i++) {
         s=0.001e0*((double)(rand())/(double)(RAND_MAX)-0.5);
         x->r[i][k]=s+x->r[i][k];
      }
   }
#endif
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
void transform_findsym(double U[3][3], nuclear_frame *x) {
   double V[3];
   int     i,j,k;
   for(k=0; k<x->n_atoms; k++) {
      for(i=0; i<3; i++) {
         V[i]=0.0e0; for(j=0; j<3; j++) V[i]=V[i]+U[i][j]*x->r[j][k];
      }
      for(i=0; i<3; i++) x->r[i][k]=V[i];
   }
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
}
void Unit_findsym(double U[3][3]) {
   int  i,j;
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         if(i==j) U[i][j]=1.0e0;
         else     U[i][j]=0.0e0;
      }
   }
#ifdef DEBUG
   printf("Unit matrix\n");
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) printf("%12.6f ",U[i][j]); printf("\n");
   }
   printf("\n");
#endif
}
/*==========================================================================*/
void print_output_findsym(int IN_GATE, FILE *io, nuclear_frame *x) {
   int   k,n;
   n=x->n_atoms;
   for(k=0; k<n; k++) {
      fprintf(io,"%-8s ",x->label[k]);
if(IN_GATE)
{
      fprintf(io,"%3i ",x->Z[k]);
}
      fprintf(io,"%15.9f ",x->r[0][k]);
      fprintf(io,"%15.9f ",x->r[1][k]);
      fprintf(io,"%15.9f ",x->r[2][k]);
if(IN_GATE)
{
      fprintf(io,"%9.3f ",x->mass[k]);
}
      fprintf(io,"\n");
   }
/*
   if(x->point_group!=NULL) {
      fprintf(io,"Point group: %s\n",x->point_group);

   }
*/
}
