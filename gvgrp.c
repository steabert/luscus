/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck gvgrp.c $Revision: 7.7 $ */

/**************************************************************************/
/*                                                                        */
/* This routine takes group generators and produces objects to be         */
/* painted by gv. Only D2h and subgroups are handeled.                    */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/* Written: May 2007                                                      */
/*                                                                        */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gvgrp.h"
/*========================================================================*/
/*                                                                        */
/*========================================================================*/
int gvgrp(char *gen, SymmetryElement X[]) {
   char *line;
   char *tok;
   char *list[]={"E","x","y","xy","z","xz","yz","xyz"};
   int   val[8];
   int   indx=0;
   int   i,k,m,n;

#ifdef DEBUG
   fprintf(stderr,"Got generators '%s'\n",gen);
#endif

   for(i=0; i<8; i++) val[i]=-1;

   n=0;
   line=(char *)malloc(strlen(gen)+1);
   strcpy(line,gen);
   tok=strtok(line," \t\n");
   while(tok!=NULL) {
      k=-1;
      for(i=0; i<8; i++) if(strcmp(list[i],tok)==0) k=i;
      if(k==-1) return 0;
      if(n==0) indx=1;
      else if(n==1) indx=2;
      else if(n==2) indx=4;
      val[indx]=k;
#ifdef DEBUG
      fprintf(stderr,"Generator '%s', index %i\n",tok,k);
#endif
      tok=strtok(NULL," \t\n");
      n++;
   }
   free(line);
   val[0]=0;

   if(n==1) {
      m=2;
   } else if(n==2) {
      m=4;
      val[3]=val[1]^val[2];
   } else if(n==3) {
      m=8;
      val[3]=val[1]^val[2];
      for(i=0; i<3; i++) val[5+i]=val[4]^val[i+1];
   } else {
      return 0;
   }

   for(i=0; i<m; i++) for(k=0; k<i; k++) if(val[i]==val[k]) return 0;

   for(i=0; i<m; i++) {
      if(val[i]==0) {        /* Unit operation         */
         X[i].type=unit;
         strcpy(X[i].label,"E");
      } else if(val[i]==1) { /* yz mirror plane, gen=x */
         X[i].type=plane;
         strcpy(X[i].label,"x");
         X[i].x[0]= 0.0; X[i].y[0]= 1.0; X[i].z[0]= 1.0;
         X[i].x[1]= 0.0; X[i].y[1]=-1.0; X[i].z[1]= 1.0;
         X[i].x[2]= 0.0; X[i].y[2]=-1.0; X[i].z[2]=-1.0;
         X[i].x[3]= 0.0; X[i].y[3]= 1.0; X[i].z[3]=-1.0;
      } else if(val[i]==2) { /* xz mirror plane, gen=y */
         X[i].type=plane;
         strcpy(X[i].label,"y");
         X[i].y[0]= 0.0; X[i].z[0]= 1.0; X[i].x[0]= 1.0;
         X[i].y[1]= 0.0; X[i].z[1]=-1.0; X[i].x[1]= 1.0;
         X[i].y[2]= 0.0; X[i].z[2]=-1.0; X[i].x[2]=-1.0;
         X[i].y[3]= 0.0; X[i].z[3]= 1.0; X[i].x[3]=-1.0;
      } else if(val[i]==3) { /* z C2 axiz, gen=xy      */
         X[i].type=axis;
         strcpy(X[i].label,"xy");
         X[i].x[0]= 0.0; X[i].y[0]= 0.0; X[i].z[0]=-1.0;
         X[i].x[1]= 0.0; X[i].y[1]= 0.0; X[i].z[1]= 1.0;
      } else if(val[i]==4) { /* xy mirror plane, gen=z */
         X[i].type=plane;
         strcpy(X[i].label,"z");
         X[i].z[0]= 0.0; X[i].x[0]= 1.0; X[i].y[0]= 1.0;
         X[i].z[1]= 0.0; X[i].x[1]=-1.0; X[i].y[1]= 1.0;
         X[i].z[2]= 0.0; X[i].x[2]=-1.0; X[i].y[2]=-1.0;
         X[i].z[3]= 0.0; X[i].x[3]= 1.0; X[i].y[3]=-1.0;
      } else if(val[i]==5) { /* y C2 axiz, gen=xz      */
         X[i].type=axis;
         strcpy(X[i].label,"xz");
         X[i].z[0]= 0.0; X[i].x[0]= 0.0; X[i].y[0]=-1.0;
         X[i].z[1]= 0.0; X[i].x[1]= 0.0; X[i].y[1]= 1.0;
      } else if(val[i]==6) { /* x C2 axiz, gen=yz      */
         X[i].type=axis;
         strcpy(X[i].label,"yz");
         X[i].y[0]= 0.0; X[i].z[0]= 0.0; X[i].x[0]=-1.0;
         X[i].y[1]= 0.0; X[i].z[1]= 0.0; X[i].x[1]= 1.0;
      } else if(val[i]==7) { /* inversion              */
         X[i].type=inversion;
         strcpy(X[i].label,"xyz");
      } else {
         return 0;
      }
   }

#ifdef DEBUG
   fprintf(stderr,"values");
   for(i=0; i<8; i++) fprintf(stderr," %2i",val[i]);
   fprintf(stderr,"\n");
#endif

   return m;

}
