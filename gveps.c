/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck gveps.c $Revision: 7.7 $ */

/**************************************************************************/
/*                                                                        */
/* Libriry to paint  molecules and orbitals in eps format                 */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/* Written: April 2007                                                    */
/*                                                                        */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gveps.h"
/**************************************************************************/
/*                                                                        */
/* This routine initialize the painting of molecules and orbitals in eps  */
/* format.                                                                */
/*                                                                        */
/**************************************************************************/

void gveps_init(FILE *f, int opt, double red, double green, double blue, gvepsdef *ctl) {

   ctl->X=NULL;
   ctl->red_bg=red;
   ctl->green_bg=green;
   ctl->blue_bg=blue;
   ctl->xmax=-1.0e6;
   ctl->xmin= 1.0e6;
   ctl->ymax=-1.0e6;
   ctl->ymin= 1.0e6;
   ctl->zmax=-1.0e6;
   ctl->zmin= 1.0e6;
   ctl->do_bg=0;
   ctl->maxobj=0;
   ctl->nobj=0;
   ctl->size=SIZE_SMALL;
   ctl->do_view=0;

   if((opt&3)!=0) ctl->size=(opt&3);
   if((opt&DO_BACKG)==DO_BACKG) ctl->do_bg=1;

   gveps_allo(ctl);

}
/**************************************************************************/
/*                                                                        */
/* This routine defines atoms to be drawn in eps format.                  */
/*                                                                        */
/**************************************************************************/
void gveps_atom(FILE *f, int opt, double radius, double x, double y, double z,
   double red, double green, double blue, gvepsdef *ctl) {
   double llx,lly,llz,urx,ury,urz;
   double tx,ty,tz;
   int n;

#ifdef DEBUG
   fprintf(stderr,"atom %6.3f %6.3f %6.3f\n",x,y,z);
#endif
   if(ctl->nobj==ctl->maxobj) gveps_allo(ctl);

   if(ctl->do_view) {
      tx=ctl->U[0][0]*x+ctl->U[0][1]*y+ctl->U[0][2]*z;
      ty=ctl->U[1][0]*x+ctl->U[1][1]*y+ctl->U[1][2]*z;
      tz=ctl->U[2][0]*x+ctl->U[2][1]*y+ctl->U[2][2]*z;
      x=tx;
      y=ty;
      z=tz;
   }

   llx=x-radius;
   lly=y-radius;
   llz=z-radius;
   urx=x+radius;
   ury=y+radius;
   urz=z+radius;
   if(llx<ctl->xmin) ctl->xmin=llx;
   if(lly<ctl->ymin) ctl->ymin=lly;
   if(llz<ctl->zmin) ctl->zmin=llz;
   if(urx>ctl->xmax) ctl->xmax=urx;
   if(ury>ctl->ymax) ctl->ymax=ury;
   if(urz>ctl->zmax) ctl->zmax=urz;

   n=ctl->nobj;
   ctl->X[n].x=x;
   ctl->X[n].y=y;
   ctl->X[n].z=z;
   ctl->X[n].radius=radius;
   ctl->X[n].red=red;
   ctl->X[n].green=green;
   ctl->X[n].blue=blue;

   ctl->X[n].kind=KIND_ATOM;
   ctl->nobj=n+1;
}
/**************************************************************************/
/*                                                                        */
/* This routine defines bonds to be drawn in eps format.                  */
/*                                                                        */
/**************************************************************************/
void gveps_bond(FILE *f, int opt,
   double x1, double y1, double z1, double red1, double green1, double blue1,
   double x2, double y2, double z2, double red2, double green2, double blue2, gvepsdef *ctl) {
   double tx,ty,tz;
   double t;
   int    n;

#ifdef DEBUG
   fprintf(stderr,"bond %6.3f %6.3f %6.3f\n",x1,y1,z1);
   fprintf(stderr,"     %6.3f %6.3f %6.3f\n",x2,y2,z2);
#endif
   if(ctl->nobj==ctl->maxobj) gveps_allo(ctl);

   if(ctl->do_view) {
      tx=ctl->U[0][0]*x1+ctl->U[0][1]*y1+ctl->U[0][2]*z1;
      ty=ctl->U[1][0]*x1+ctl->U[1][1]*y1+ctl->U[1][2]*z1;
      tz=ctl->U[2][0]*x1+ctl->U[2][1]*y1+ctl->U[2][2]*z1;
      x1=tx;
      y1=ty;
      z1=tz;
      tx=ctl->U[0][0]*x2+ctl->U[0][1]*y2+ctl->U[0][2]*z2;
      ty=ctl->U[1][0]*x2+ctl->U[1][1]*y2+ctl->U[1][2]*z2;
      tz=ctl->U[2][0]*x2+ctl->U[2][1]*y2+ctl->U[2][2]*z2;
      x2=tx;
      y2=ty;
      z2=tz;
   }


   if(z2<z1) {
      t=x1; x1=x2; x2=t;
      t=y1; y1=y2; y2=t;
      t=z1; z1=z2; z2=t;
      t=red1; red1=red2; red2=t;
      t=green1; green1=green2; green2=t;
      t=blue1; blue1=blue2; blue2=t;
   }

   n=ctl->nobj;
   ctl->X[n].x=x1;
   ctl->X[n].y=y1;
   ctl->X[n].z=z1;
   ctl->X[n].x1=x2;
   ctl->X[n].y1=y2;
   ctl->X[n].z1=z2;
   ctl->X[n].red=red1;
   ctl->X[n].green=green1;
   ctl->X[n].blue=blue1;
   ctl->X[n].red1=red2;
   ctl->X[n].green1=green2;
   ctl->X[n].blue1=blue2;
   if((opt&BOND_DASH)==BOND_DASH) ctl->X[n].dash=1;
   else                           ctl->X[n].dash=0;

   ctl->X[n].kind=KIND_BOND;
   ctl->nobj=n+1;

   if(x1<ctl->xmin) ctl->xmin=x1;
   if(y1<ctl->ymin) ctl->ymin=y1;
   if(z1<ctl->zmin) ctl->zmin=z1;
   if(x1>ctl->xmax) ctl->xmax=x1;
   if(y1>ctl->ymax) ctl->ymax=y1;
   if(z1>ctl->zmax) ctl->zmax=z1;

   if(x2<ctl->xmin) ctl->xmin=x2;
   if(y2<ctl->ymin) ctl->ymin=y2;
   if(z2<ctl->zmin) ctl->zmin=z2;
   if(x2>ctl->xmax) ctl->xmax=x2;
   if(y2>ctl->ymax) ctl->ymax=y2;
   if(z2>ctl->zmax) ctl->zmax=z2;
}
/**************************************************************************/
/*                                                                        */
/* This routine defines coloured triangled to be drawn in eps format.     */
/*                                                                        */
/**************************************************************************/
void gveps_triangle(FILE *f, int opt,
   double x1, double y1, double z1,
   double x2, double y2, double z2,
   double x3, double y3, double z3,
   double nx1, double ny1, double nz1,
   double nx2, double ny2, double nz2,
   double nx3, double ny3, double nz3,
   double red, double green, double blue,  gvepsdef *ctl) {

   double tx,ty,tz;
   double scale;
   double t;
   int    n;

#ifdef DEBUG
   fprintf(stderr,"tri  %6.3f %6.3f %6.3f\n",x1,y1,z1);
   fprintf(stderr,"     %6.3f %6.3f %6.3f\n",x2,y2,z2);
   fprintf(stderr,"     %6.3f %6.3f %6.3f\n",x3,y3,z3);
   fprintf(stderr,"     %6.3f %6.3f %6.3f\n",nx1,ny1,nz1);
   fprintf(stderr,"     %6.3f %6.3f %6.3f\n",nx2,ny2,nz2);
   fprintf(stderr,"     %6.3f %6.3f %6.3f\n",nx3,ny3,nz3);
#endif
   if(ctl->nobj==ctl->maxobj) gveps_allo(ctl);

   if(ctl->do_view) {
      tx=ctl->U[0][0]*x1+ctl->U[0][1]*y1+ctl->U[0][2]*z1;
      ty=ctl->U[1][0]*x1+ctl->U[1][1]*y1+ctl->U[1][2]*z1;
      tz=ctl->U[2][0]*x1+ctl->U[2][1]*y1+ctl->U[2][2]*z1;
      x1=tx;
      y1=ty;
      z1=tz;
      tx=ctl->U[0][0]*x2+ctl->U[0][1]*y2+ctl->U[0][2]*z2;
      ty=ctl->U[1][0]*x2+ctl->U[1][1]*y2+ctl->U[1][2]*z2;
      tz=ctl->U[2][0]*x2+ctl->U[2][1]*y2+ctl->U[2][2]*z2;
      x2=tx;
      y2=ty;
      z2=tz;
      tx=ctl->U[0][0]*x3+ctl->U[0][1]*y3+ctl->U[0][2]*z3;
      ty=ctl->U[1][0]*x3+ctl->U[1][1]*y3+ctl->U[1][2]*z3;
      tz=ctl->U[2][0]*x3+ctl->U[2][1]*y3+ctl->U[2][2]*z3;
      x3=tx;
      y3=ty;
      z3=tz;
      tx=ctl->U[0][0]*nx1+ctl->U[0][1]*ny1+ctl->U[0][2]*nz1;
      ty=ctl->U[1][0]*nx1+ctl->U[1][1]*ny1+ctl->U[1][2]*nz1;
      tz=ctl->U[2][0]*nx1+ctl->U[2][1]*ny1+ctl->U[2][2]*nz1;
      nx1=tx;
      ny1=ty;
      nz1=tz;
      tx=ctl->U[0][0]*nx2+ctl->U[0][1]*ny2+ctl->U[0][2]*nz2;
      ty=ctl->U[1][0]*nx2+ctl->U[1][1]*ny2+ctl->U[1][2]*nz2;
      tz=ctl->U[2][0]*nx2+ctl->U[2][1]*ny2+ctl->U[2][2]*nz2;
      nx2=tx;
      ny2=ty;
      nz2=tz;
      tx=ctl->U[0][0]*nx3+ctl->U[0][1]*ny3+ctl->U[0][2]*nz3;
      ty=ctl->U[1][0]*nx3+ctl->U[1][1]*ny3+ctl->U[1][2]*nz3;
      tz=ctl->U[2][0]*nx3+ctl->U[2][1]*ny3+ctl->U[2][2]*nz3;
      nx3=tx;
      ny3=ty;
      nz3=tz;
   }

   scale=1.0/sqrt(nx1*nx1+ny1*ny1+nz1*nz1);
   nx1=scale*nx1; ny1=scale*ny1; nz1=scale*nz1;
   scale=1.0/sqrt(nx2*nx2+ny2*ny2+nz2*nz2);
   nx2=scale*nx2; ny2=scale*ny2; nz2=scale*nz2;
   scale=1.0/sqrt(nx3*nx3+ny3*ny3+nz3*nz3);
   nx3=scale*nx3; ny3=scale*ny3; nz3=scale*nz3;

   if(z2<z1) {
      if(z3<z2) {
         t=x1;  x1=x3;   x3=t;
         t=y1;  y1=y3;   y3=t;
         t=z1;  z1=z3;   z3=t;
         t=nx1; nx1=nx3; nx3=t;
         t=ny1; ny1=ny3; ny3=t;
         t=nz1; nz1=nz3; nz3=t;
      } else {
         t=x1;  x1=x2;   x2=t;
         t=y1;  y1=y2;   y2=t;
         t=z1;  z1=z2;   z2=t;
         t=nx1; nx1=nx2; nx2=t;
         t=ny1; ny1=ny2; ny2=t;
         t=nz1; nz1=nz2; nz2=t;
      }
   } else if(z3<z1) {
      t=x1;  x1=x3;   x3=t;
      t=y1;  y1=y3;   y3=t;
      t=z1;  z1=z3;   z3=t;
      t=nx1; nx1=nx3; nx3=t;
      t=ny1; ny1=ny3; ny3=t;
      t=nz1; nz1=nz3; nz3=t;
   }

   n=ctl->nobj;
   ctl->X[n].x=x1;
   ctl->X[n].y=y1;
   ctl->X[n].z=z1;
   ctl->X[n].x1=x2;
   ctl->X[n].y1=y2;
   ctl->X[n].z1=z2;
   ctl->X[n].x2=x3;
   ctl->X[n].y2=y3;
   ctl->X[n].z2=z3;
   ctl->X[n].nx1=nx1;
   ctl->X[n].ny1=ny1;
   ctl->X[n].nz1=nz1;
   ctl->X[n].nx2=nx2;
   ctl->X[n].ny2=ny2;
   ctl->X[n].nz2=nz2;
   ctl->X[n].nx3=nx3;
   ctl->X[n].ny3=ny3;
   ctl->X[n].nz3=nz3;
   ctl->X[n].red=red;
   ctl->X[n].green=green;
   ctl->X[n].blue=blue;

   ctl->X[n].kind=KIND_TRI;
   ctl->nobj=n+1;

   if(x1<ctl->xmin) ctl->xmin=x1;
   if(y1<ctl->ymin) ctl->ymin=y1;
   if(z1<ctl->zmin) ctl->zmin=z1;
   if(x1>ctl->xmax) ctl->xmax=x1;
   if(y1>ctl->ymax) ctl->ymax=y1;
   if(z1>ctl->zmax) ctl->zmax=z1;

   if(x2<ctl->xmin) ctl->xmin=x2;
   if(y2<ctl->ymin) ctl->ymin=y2;
   if(z2<ctl->zmin) ctl->zmin=z2;
   if(x2>ctl->xmax) ctl->xmax=x2;
   if(y2>ctl->ymax) ctl->ymax=y2;
   if(z2>ctl->zmax) ctl->zmax=z2;

   if(x3<ctl->xmin) ctl->xmin=x3;
   if(y3<ctl->ymin) ctl->ymin=y3;
   if(z3<ctl->zmin) ctl->zmin=z3;
   if(x3>ctl->xmax) ctl->xmax=x3;
   if(y3>ctl->ymax) ctl->ymax=y3;
   if(z3>ctl->zmax) ctl->zmax=z3;
}
/**************************************************************************/
/*                                                                        */
/* This routine defines text to be drawn in eps format. Note that only up */
/* to 7 characters is supported.                                          */
/*                                                                        */
/**************************************************************************/
void gveps_text(FILE *f, int opt, char *str, double x, double y, double z, gvepsdef *ctl) {
   int n;

#ifdef DEBUG
   fprintf(stderr,"text %6.3f %6.3f %6.3f '%s'\n",x,y,z,str);
#endif
   if(ctl->nobj==ctl->maxobj) gveps_allo(ctl);

   if(x<ctl->xmin) ctl->xmin=x;
   if(y<ctl->ymin) ctl->ymin=y;
   if(z<ctl->zmin) ctl->zmin=z;

   n=ctl->nobj;
   ctl->X[n].x=x;
   ctl->X[n].y=y;
   ctl->X[n].z=z;
   strncpy(ctl->X[n].str,str,7);
   ctl->X[n].vert=TEXT_MIDDLE;
   ctl->X[n].hori=TEXT_CENTER;
   if((opt&3)!=0) {
      ctl->X[n].vert=(opt&3);
#ifdef DEBUG
      fprintf(stderr,"Vertical for '%s' is %i\n",str,(opt&3));
#endif
   }
   if((opt&12)!=0) {
      ctl->X[n].hori=(opt&12);
#ifdef DEBUG
      fprintf(stderr,"Horizontal for '%s' is %i\n",str,(opt&12));
#endif
   }

   ctl->X[n].kind=KIND_TEXT;
   ctl->nobj=n+1;
}
/**************************************************************************/
/*                                                                        */
/* This routine defines the viewpoint of the viewer.                      */
/*                                                                        */
/**************************************************************************/
void gveps_view(FILE *f, int opt, double U[3][3], gvepsdef *ctl) {
   int i,j;

   ctl->do_view=1;
   for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
         ctl->U[i][j]=U[i][j];
      }
   }

}
/**************************************************************************/
/*                                                                        */
/* This routine finally paint molecules and orbitals in eps format.       */
/*                                                                        */
/* opt=1 - use eps level2                                                 */
/*                                                                        */
/**************************************************************************/
void gveps_paint(FILE *f, int opt, gvepsdef *ctl) {
   double llx,lly,urx,ury;
   double scale,tx,ty;
   double xoff,yoff;
   double dx,dy,dz;
   double slant;
   double nz1,nz2,nz3;
   double margin;
   double x1,y1;
   double x2,y2;
   double x3,y3;
   double bondwidth;
   double notzero;
   double fs;
   int    i,n;

   bondwidth=0.15;
   notzero=0.0001;
   fs=0.2;

   n=ctl->nobj;
   gveps_sort(ctl->X,n);
   gveps_trim(ctl->X,n);
#ifdef DEBUG
/*   dump_def(stderr,*ctl); */
#endif

/*------------------------------------------------------------------------*/
/* Figure out translate and scale                                         */
/*------------------------------------------------------------------------*/
   margin=0.1;
   llx=ctl->xmin-margin;
   lly=ctl->ymin-margin;
   urx=ctl->xmax+margin;
   ury=ctl->ymax+margin;
   xoff=-llx;
   yoff=-lly;
   tx=451.0/(urx-llx);
   ty=681.0/(ury-lly);
   if(tx>ty) scale=ty;
   else      scale=tx;
   ty=1.0e-2;
   for(i=-3; i<4; i++) {
      tx=exp((double)(i)*log(10.0));
      if(scale>1.0*tx) ty=1.0*tx;
      if(scale>2.0*tx) ty=2.0*tx;
      if(scale>5.0*tx) ty=5.0*tx;
   }
   scale=ty;
   llx=72.0;
   lly=72.0;
   urx=llx+scale*(ctl->xmax-ctl->xmin+2*margin);
   ury=lly+scale*(ctl->ymax-ctl->ymin+2*margin);
/*------------------------------------------------------------------------*/
/* Print prologue                                                         */
/*------------------------------------------------------------------------*/
   fprintf(f,"%s\n","%!PS-Adobe-2.0 EPSF");
   fprintf(f,"%s","%%BoundingBox: ");
   fprintf(f,"%i %i %i %i\n",(int)llx,(int)lly,(int)urx,(int)ury);
   fprintf(f,"50 dict begin\n");
   fprintf(f,"/languagelevel where {\n");
   fprintf(f,"  pop languagelevel 3 ge {\n");
   fprintf(f,"    /islevel3 true def\n");
   fprintf(f,"  } {\n");
   fprintf(f,"    /islevel3 false def\n");
   fprintf(f,"  } ifelse\n");
   fprintf(f,"} {\n");
   fprintf(f,"  /islevel3 false def\n");
   fprintf(f,"} ifelse\n");
   fprintf(f,"%s\n","%--- start configure");
   if(opt==0)
   {
    fprintf(f,"%s\n","%/islevel3 false def");
   }
   else
   {
    fprintf(f,"%s\n","/islevel3 false def");  
   }
   fprintf(f,"/bondwidth %.3f def\n",bondwidth);
   fprintf(f,"/notzero %.6f def\n",notzero);
   fprintf(f,"/fs %.3f def\n",fs);
   fprintf(f,"/RedColor   { 1.0 0.0 0.0 } def\n");
   fprintf(f,"/GreenColor { 0.0 1.0 0.0 } def\n");
   fprintf(f,"/BlueColor  { 0.0 0.0 1.0 } def\n");
   fprintf(f,"%s\n","%--- end configure");
   fprintf(f,"/Helvetica findfont fs scalefont setfont\n");
   fprintf(f,"/Bond {\n");
   fprintf(f,"  /dodash exch def\n");
   fprintf(f,"  /slant exch def\n");
   fprintf(f,"  /y2 exch def\n");
   fprintf(f,"  /x2 exch def\n");
   fprintf(f,"  /y1 exch def\n");
   fprintf(f,"  /x1 exch def\n");
   fprintf(f,"  gsave x1 y1 translate\n");
   fprintf(f,"    dodash { [0.07 0.10] 0 setdash } if\n");
   fprintf(f,"    x2 x1 neg add\n");
   fprintf(f,"    y2 y1 neg add\n");
   fprintf(f,"    atan neg rotate\n");
   fprintf(f,"    /dist x2 x1 neg add dup mul y2 y1 neg add dup mul add sqrt def\n");
   fprintf(f,"    90 -1 1 {\n");
   fprintf(f,"      /angle exch def\n");
   fprintf(f,"      /w angle sin bondwidth mul def\n");
   fprintf(f,"      /s angle cos def\n");
   fprintf(f,"      /y bondwidth notzero add dup mul w dup mul neg add sqrt neg 0.15 add slant mul 0.5 mul def\n");
   fprintf(f,"      s setgray\n");
   fprintf(f,"      w setlinewidth\n");
   fprintf(f,"      0 y moveto 0 dist lineto stroke\n");
   fprintf(f,"    } for\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/Atom {\n");
   fprintf(f,"  /blue exch def\n");
   fprintf(f,"  /green exch def\n");
   fprintf(f,"  /red exch def\n");
   fprintf(f,"  /radius exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    90 -1 1 {\n");
   fprintf(f,"      /angle exch def\n");
   fprintf(f,"      /r angle sin radius mul def\n");
   fprintf(f,"      /s angle cos def\n");
   fprintf(f,"      red s mul green s mul blue s mul setrgbcolor\n");
   fprintf(f,"      newpath\n");
   fprintf(f,"      0 0 r 0 360 arc\n");
   fprintf(f,"      closepath fill\n");
   fprintf(f,"    } for\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"islevel3 {\n");
   fprintf(f,"  /Triangle {\n");
   fprintf(f,"    /blue exch def\n");
   fprintf(f,"    /green exch def\n");
   fprintf(f,"    /red exch def\n");
   fprintf(f,"    /s3 exch def\n");
   fprintf(f,"    /s2 exch def\n");
   fprintf(f,"    /s1 exch def\n");
   fprintf(f,"    /y3 exch def\n");
   fprintf(f,"    /x3 exch def\n");
   fprintf(f,"    /y2 exch def\n");
   fprintf(f,"    /x2 exch def\n");
   fprintf(f,"    /y1 exch def\n");
   fprintf(f,"    /x1 exch def\n");
   fprintf(f,"    /r1 red s1 mul def\n");
   fprintf(f,"    /r2 red s2 mul def\n");
   fprintf(f,"    /r3 red s3 mul def\n");
   fprintf(f,"    /g1 green s1 mul def\n");
   fprintf(f,"    /g2 green s2 mul def\n");
   fprintf(f,"    /g3 green s3 mul def\n");
   fprintf(f,"    /b1 blue s1 mul def\n");
   fprintf(f,"    /b2 blue s2 mul def\n");
   fprintf(f,"    /b3 blue s3 mul def\n");
   fprintf(f,"    gsave\n");
   fprintf(f,"      /DeviceRGB setcolorspace\n");
   fprintf(f,"      <<\n");
   fprintf(f,"        /ShadingType 4\n");
   fprintf(f,"        /ColorSpace [/DeviceRGB]\n");
   fprintf(f,"        /DataSource [ 0 x1 y1 r1 g1 b1\n");
   fprintf(f,"                      0 x2 y2 r2 g2 b2\n");
   fprintf(f,"                      0 x3 y3 r3 g3 b3 ]\n");
   fprintf(f,"      >>\n");
   fprintf(f,"      shfill\n");
   fprintf(f,"    grestore\n");
   fprintf(f,"  } def\n");
   fprintf(f,"} {\n");
   fprintf(f,"  /Triangle {\n");
   fprintf(f,"    /blue exch def\n");
   fprintf(f,"    /green exch def\n");
   fprintf(f,"    /red exch def\n");
   fprintf(f,"    add add 0.333 mul\n");
   fprintf(f,"    /s exch def\n");
   fprintf(f,"    red s mul green s mul blue s mul setrgbcolor\n");
   fprintf(f,"    moveto lineto lineto\n");
   fprintf(f,"    closepath fill\n");
   fprintf(f,"  } def\n");
   fprintf(f,"} ifelse\n");
   fprintf(f,"/TextCM {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str stringwidth pop -0.5 mul fs -0.35 mul rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextLM {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto 0 fs -0.35 mul rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextRM {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str stringwidth pop neg fs -0.35 mul rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextCT {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str stringwidth pop -0.5 mul fs -0.70 mul rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextLT {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto 0 fs -0.70 mul rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextRT {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str stringwidth pop neg fs -0.70 mul rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextCB {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str stringwidth pop -0.5 mul 0 rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextLB {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   fprintf(f,"/TextRB {\n");
   fprintf(f,"  /str exch def\n");
   fprintf(f,"  gsave translate\n");
   fprintf(f,"    0 setgray\n");
   fprintf(f,"    0 0 moveto str stringwidth pop neg 0 rmoveto str show\n");
   fprintf(f,"  grestore\n");
   fprintf(f,"} def\n");
   if(ctl->do_bg) {
      fprintf(f,"%.3f %.3f %.3f setrgbcolor\n",ctl->red_bg,ctl->green_bg,ctl->blue_bg);
      fprintf(f,"%i %i moveto\n",(int)llx,(int)lly);
      fprintf(f,"%i %i lineto\n",(int)urx,(int)lly);
      fprintf(f,"%i %i lineto\n",(int)urx,(int)ury);
      fprintf(f,"%i %i lineto\n",(int)llx,(int)ury);
      fprintf(f,"closepath fill\n");
   }
   fprintf(f,"72 72 translate %.3f %.3f scale\n",scale,scale);
/*------------------------------------------------------------------------*/
/* Paint objects                                                          */
/*------------------------------------------------------------------------*/
   for (i=0; i<n; i++) {
      if(ctl->X[i].kind==KIND_ATOM) {
         x1=ctl->X[i].x+xoff;
         y1=ctl->X[i].y+yoff;
         fprintf(f,"%.3f ",x1);
         fprintf(f,"%.3f ",y1);
         fprintf(f,"%.3f ",ctl->X[i].radius);
         fprintf(f,"%.3f ",ctl->X[i].red);
         fprintf(f,"%.3f ",ctl->X[i].green);
         fprintf(f,"%.3f ",ctl->X[i].blue);
         fprintf(f,"Atom\n");
      } else if(ctl->X[i].kind==KIND_BOND) {
         dx=ctl->X[i].x1-ctl->X[i].x;
         dy=ctl->X[i].y1-ctl->X[i].y;
         dz=ctl->X[i].z1-ctl->X[i].z;
         slant=dz/sqrt(dx*dx+dy*dy+dz*dz);
         x1=ctl->X[i].x+xoff;
         y1=ctl->X[i].y+yoff;
         x2=ctl->X[i].x1+xoff;
         y2=ctl->X[i].y1+yoff;
         if(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) < 0.001) {
            x2=x2+0.001;
            y2=y2+0.001;
         }
         fprintf(f,"%.3f ",x1);
         fprintf(f,"%.3f ",y1);
         fprintf(f,"%.3f ",x2);
         fprintf(f,"%.3f ",y2);
         fprintf(f,"%.3f ",slant);
         if(ctl->X[i].dash) fprintf(f,"true ");
         else               fprintf(f,"false ");
         fprintf(f,"Bond\n");
      } else if(ctl->X[i].kind==KIND_TRI) {
         x1=ctl->X[i].x+xoff;
         y1=ctl->X[i].y+yoff;
         x2=ctl->X[i].x1+xoff;
         y2=ctl->X[i].y1+yoff;
         x3=ctl->X[i].x2+xoff;
         y3=ctl->X[i].y2+yoff;
         if(ctl->X[i].nz1<0.0) nz1=0.0;
         else                  nz1=ctl->X[i].nz1;
         if(ctl->X[i].nz2<0.0) nz2=0.0;
         else                  nz2=ctl->X[i].nz2;
         if(ctl->X[i].nz3<0.0) nz3=0.0;
         else                  nz3=ctl->X[i].nz3;
         fprintf(f,"%.3f ",x1);
         fprintf(f,"%.3f ",y1);
         fprintf(f,"%.3f ",x2);
         fprintf(f,"%.3f ",y2);
         fprintf(f,"%.3f ",x3);
         fprintf(f,"%.3f ",y3);
         fprintf(f,"%.3f ",nz1);
         fprintf(f,"%.3f ",nz2);
         fprintf(f,"%.3f ",nz3);
         if(fabs(ctl->X[i].red-1.0)+fabs(ctl->X[i].green)+fabs(ctl->X[i].blue)<0.01) {
            fprintf(f,"RedColor ");
         } else if(fabs(ctl->X[i].red)+fabs(ctl->X[i].green-1.0)+fabs(ctl->X[i].blue)<0.01) {
            fprintf(f,"GreenColor ");
         } else if(fabs(ctl->X[i].red)+fabs(ctl->X[i].green)+fabs(ctl->X[i].blue-1.0)<0.01) {
            fprintf(f,"BlueColor ");
         } else {
            fprintf(f,"%.3f ",ctl->X[i].red);
            fprintf(f,"%.3f ",ctl->X[i].green);
            fprintf(f,"%.3f ",ctl->X[i].blue);
         }
         fprintf(f,"Triangle\n");
      } else if(ctl->X[i].kind==KIND_TEXT) {
         x1=ctl->X[i].x+xoff;
         y1=ctl->X[i].y+yoff;
         fprintf(f,"%.3f ",x1);
         fprintf(f,"%.3f ",y1);
         fprintf(f,"(%s) ",ctl->X[i].str);
         if(ctl->X[i].vert==TEXT_TOP) {
            if(ctl->X[i].hori==TEXT_LEFT) {
               fprintf(f,"TextLT\n");
            } else if(ctl->X[i].hori==TEXT_RIGHT) {
               fprintf(f,"TextRT\n");
            } else {
               fprintf(f,"TextCT\n");
            }
         } else if(ctl->X[i].vert==TEXT_BOTTOM) {
            if(ctl->X[i].hori==TEXT_LEFT) {
               fprintf(f,"TextLB\n");
            } else if(ctl->X[i].hori==TEXT_RIGHT) {
               fprintf(f,"TextRB\n");
            } else {
               fprintf(f,"TextCB\n");
            }
         } else {
            if(ctl->X[i].hori==TEXT_LEFT) {
               fprintf(f,"TextLM\n");
            } else if(ctl->X[i].hori==TEXT_RIGHT) {
               fprintf(f,"TextRM\n");
            } else {
               fprintf(f,"TextCM\n");
            }
         }
      } else {
         fprintf(f,"%s\n","%What");
      }
   }
/*------------------------------------------------------------------------*/
/* Pprint epilogue                                                        */
/*------------------------------------------------------------------------*/
   fprintf(f,"end\n");
   fprintf(f,"showpage\n");
}
/**************************************************************************/
/*                                                                        */
/* This routine rotates all objects defined so far. Transformation is     */
/* done using euler angles.                                               */
/*                                                                        */
/**************************************************************************/
void gveps_euler(FILE *f, int opt, double alpha, double beta, double gamma, gvepsdef *ctl) {
   double U[3][3];
   double r1[3],r2[3],r3[3],r4[3],r5[3],r6[3];
   double v1[3],v2[3],v3[3],v4[3],v5[3],v6[3];
   double llx,lly,llz,urx,ury,urz;
   double t;
   int    i,k,m,n;

   U[0][0]= cos(gamma)*cos(beta)*cos(alpha)-sin(gamma)*sin(alpha);
   U[0][1]= cos(gamma)*cos(beta)*sin(alpha)+sin(gamma)*cos(alpha);
   U[0][2]=-cos(gamma)*sin(beta);
   U[1][0]=-sin(gamma)*cos(beta)*cos(alpha)-cos(gamma)*sin(alpha);
   U[1][1]=-sin(gamma)*cos(beta)*sin(alpha)+cos(gamma)*cos(alpha);
   U[1][2]= sin(gamma)*sin(beta);
   U[2][0]= sin(beta)*cos(alpha);
   U[2][1]= sin(beta)*sin(alpha);
   U[2][2]= cos(beta);

   n=ctl->nobj;
   ctl->xmax=ctl->X[0].x;
   ctl->ymax=ctl->X[0].y;
   ctl->zmax=ctl->X[0].z;
   ctl->xmin=ctl->X[0].x;
   ctl->ymin=ctl->X[0].y;
   ctl->zmin=ctl->X[0].z;

   for (i=0; i<n; i++) {
      if(ctl->X[i].kind==KIND_ATOM) {
         r1[0]=ctl->X[i].x;
         r1[1]=ctl->X[i].y;
         r1[2]=ctl->X[i].z;
         for(k=0; k<3; k++) {
            v1[k]=0; for(m=0; m<3; m++) v1[k]=v1[k]+U[k][m]*r1[m];
         }
         ctl->X[i].x=v1[0];
         ctl->X[i].y=v1[1];
         ctl->X[i].z=v1[2];
         llx=v1[0]-ctl->X[i].radius;
         lly=v1[1]-ctl->X[i].radius;
         llz=v1[2]-ctl->X[i].radius;
         urx=v1[0]+ctl->X[i].radius;
         ury=v1[1]+ctl->X[i].radius;
         urz=v1[2]+ctl->X[i].radius;
         if(llx<ctl->xmin) ctl->xmin=llx;
         if(lly<ctl->ymin) ctl->ymin=lly;
         if(llz<ctl->zmin) ctl->zmin=llz;
         if(urx>ctl->xmax) ctl->xmax=urx;
         if(ury>ctl->ymax) ctl->ymax=ury;
         if(urz>ctl->zmax) ctl->zmax=urz;
      } else if(ctl->X[i].kind==KIND_BOND) {
         r1[0]=ctl->X[i].x;
         r1[1]=ctl->X[i].y;
         r1[2]=ctl->X[i].z;
         r2[0]=ctl->X[i].x1;
         r2[1]=ctl->X[i].y1;
         r2[2]=ctl->X[i].z1;
         for(k=0; k<3; k++) {
            v1[k]=0; for(m=0; m<3; m++) v1[k]=v1[k]+U[k][m]*r1[m];
            v2[k]=0; for(m=0; m<3; m++) v2[k]=v2[k]+U[k][m]*r2[m];
         }
         if(v2[2]<v1[2]) {
            for(m=0; m<3; m++) { t=v1[m]; v1[m]=v2[m]; v2[m]=t; }
            t=ctl->X[i].red; ctl->X[i].red=ctl->X[i].red1; ctl->X[i].red1=t;
            t=ctl->X[i].green; ctl->X[i].green=ctl->X[i].green1; ctl->X[i].green1=t;
            t=ctl->X[i].blue; ctl->X[i].blue=ctl->X[i].blue1; ctl->X[i].blue1=t;
         }
         if(v1[0]<ctl->xmin) ctl->xmin=v1[0];
         if(v1[1]<ctl->ymin) ctl->ymin=v1[1];
         if(v1[2]<ctl->zmin) ctl->zmin=v1[2];
         if(v2[0]<ctl->xmin) ctl->xmin=v2[0];
         if(v2[1]<ctl->ymin) ctl->ymin=v2[1];
         if(v2[2]<ctl->zmin) ctl->zmin=v2[2];
         if(v1[0]>ctl->xmax) ctl->xmax=v1[0];
         if(v1[1]>ctl->ymax) ctl->ymax=v1[1];
         if(v1[2]>ctl->zmax) ctl->zmax=v1[2];
         if(v2[0]>ctl->xmax) ctl->xmax=v2[0];
         if(v2[1]>ctl->ymax) ctl->ymax=v2[1];
         if(v2[2]>ctl->zmax) ctl->zmax=v2[2];
         ctl->X[i].x=v1[0];
         ctl->X[i].y=v1[1];
         ctl->X[i].z=v1[2];
         ctl->X[i].x1=v2[0];
         ctl->X[i].y1=v2[1];
         ctl->X[i].z1=v2[2];
      } else if(ctl->X[i].kind==KIND_TRI) {
         r1[0]=ctl->X[i].x;
         r1[1]=ctl->X[i].y;
         r1[2]=ctl->X[i].z;
         r2[0]=ctl->X[i].x1;
         r2[1]=ctl->X[i].y1;
         r2[2]=ctl->X[i].z1;
         r3[0]=ctl->X[i].x2;
         r3[1]=ctl->X[i].y2;
         r3[2]=ctl->X[i].z2;
         r4[0]=ctl->X[i].nx1;
         r4[1]=ctl->X[i].ny1;
         r4[2]=ctl->X[i].nz1;
         r5[0]=ctl->X[i].nx2;
         r5[1]=ctl->X[i].ny2;
         r5[2]=ctl->X[i].nz2;
         r6[0]=ctl->X[i].nx3;
         r6[1]=ctl->X[i].ny3;
         r6[2]=ctl->X[i].nz3;
         for(k=0; k<3; k++) {
            v1[k]=0; for(m=0; m<3; m++) v1[k]=v1[k]+U[k][m]*r1[m];
            v2[k]=0; for(m=0; m<3; m++) v2[k]=v2[k]+U[k][m]*r2[m];
            v3[k]=0; for(m=0; m<3; m++) v3[k]=v3[k]+U[k][m]*r3[m];
            v4[k]=0; for(m=0; m<3; m++) v4[k]=v4[k]+U[k][m]*r4[m];
            v5[k]=0; for(m=0; m<3; m++) v5[k]=v5[k]+U[k][m]*r5[m];
            v6[k]=0; for(m=0; m<3; m++) v6[k]=v6[k]+U[k][m]*r6[m];
         }
         if(v2[2]<v1[2]) {
            if(v3[2]<v2[2]) {
               for(m=0; m<3; m++) { t=v1[m]; v1[m]=v3[m]; v3[m]=t; }
               for(m=0; m<3; m++) { t=v4[m]; v4[m]=v6[m]; v6[m]=t; }
            } else {
               for(m=0; m<3; m++) { t=v1[m]; v1[m]=v2[m]; v2[m]=t; }
               for(m=0; m<3; m++) { t=v4[m]; v4[m]=v5[m]; v5[m]=t; }
            }
         } else if(v3[2]<v1[2]) {
            for(m=0; m<3; m++) { t=v1[m]; v1[m]=v3[m]; v3[m]=t; }
            for(m=0; m<3; m++) { t=v4[m]; v4[m]=v6[m]; v6[m]=t; }
         }
         if(v1[0]<ctl->xmin) ctl->xmin=v1[0];
         if(v1[1]<ctl->ymin) ctl->ymin=v1[1];
         if(v1[2]<ctl->zmin) ctl->zmin=v1[2];
         if(v2[0]<ctl->xmin) ctl->xmin=v2[0];
         if(v2[1]<ctl->ymin) ctl->ymin=v2[1];
         if(v2[2]<ctl->zmin) ctl->zmin=v2[2];
         if(v3[0]<ctl->xmin) ctl->xmin=v3[0];
         if(v3[1]<ctl->ymin) ctl->ymin=v3[1];
         if(v3[2]<ctl->zmin) ctl->zmin=v3[2];
         if(v1[0]>ctl->xmax) ctl->xmax=v1[0];
         if(v1[1]>ctl->ymax) ctl->ymax=v1[1];
         if(v1[2]>ctl->zmax) ctl->zmax=v1[2];
         if(v2[0]>ctl->xmax) ctl->xmax=v2[0];
         if(v2[1]>ctl->ymax) ctl->ymax=v2[1];
         if(v2[2]>ctl->zmax) ctl->zmax=v2[2];
         if(v3[0]>ctl->xmax) ctl->xmax=v3[0];
         if(v3[1]>ctl->ymax) ctl->ymax=v3[1];
         if(v3[2]>ctl->zmax) ctl->zmax=v3[2];
         ctl->X[i].x=v1[0];
         ctl->X[i].y=v1[1];
         ctl->X[i].z=v1[2];
         ctl->X[i].x1=v2[0];
         ctl->X[i].y1=v2[1];
         ctl->X[i].z1=v2[2];
         ctl->X[i].x2=v3[0];
         ctl->X[i].y2=v3[1];
         ctl->X[i].z2=v3[2];
         ctl->X[i].nx1=v4[0];
         ctl->X[i].ny1=v4[1];
         ctl->X[i].nz1=v4[2];
         ctl->X[i].nx2=v5[0];
         ctl->X[i].ny2=v5[1];
         ctl->X[i].nz2=v5[2];
         ctl->X[i].nx3=v6[0];
         ctl->X[i].ny3=v6[1];
         ctl->X[i].nz3=v6[2];
      } else if(ctl->X[i].kind==KIND_TEXT) {
         r1[0]=ctl->X[i].x;
         r1[1]=ctl->X[i].y;
         r1[2]=ctl->X[i].z;
         for(k=0; k<3; k++) {
            v1[k]=0;
            for(m=0; m<3; m++) v1[k]=v1[k]+U[k][m]*r1[m];
         }
         if(v1[0]<ctl->xmin) ctl->xmin=v1[0];
         if(v1[1]<ctl->ymin) ctl->ymin=v1[1];
         if(v1[2]<ctl->zmin) ctl->zmin=v1[2];
         if(v1[0]>ctl->xmax) ctl->xmax=v1[0];
         if(v1[1]>ctl->ymax) ctl->ymax=v1[1];
         if(v1[2]>ctl->zmax) ctl->zmax=v1[2];
         ctl->X[i].x=v1[0];
         ctl->X[i].y=v1[1];
         ctl->X[i].z=v1[2];
      }
   }

}
/**************************************************************************/
/*                                                                        */
/* This routine translates all objects defined so far.                    */
/*                                                                        */
/**************************************************************************/
void gveps_move(FILE *f, int opt, double x, double y, double z, gvepsdef *ctl) {
   int i,n;

   n=ctl->nobj;
   for (i=0; i<n; i++) {
      if(ctl->X[i].kind==KIND_ATOM) {
         ctl->X[i].x=ctl->X[i].x+x;
         ctl->X[i].y=ctl->X[i].y+y;
         ctl->X[i].z=ctl->X[i].z+z;
      } else if(ctl->X[i].kind==KIND_BOND) {
         ctl->X[i].x=ctl->X[i].x+x;
         ctl->X[i].y=ctl->X[i].y+y;
         ctl->X[i].z=ctl->X[i].z+z;
         ctl->X[i].x1=ctl->X[i].x1+x;
         ctl->X[i].y1=ctl->X[i].y1+y;
         ctl->X[i].z1=ctl->X[i].z1+z;
      } else if(ctl->X[i].kind==KIND_TRI) {
         ctl->X[i].x=ctl->X[i].x+x;
         ctl->X[i].y=ctl->X[i].y+y;
         ctl->X[i].z=ctl->X[i].z+z;
         ctl->X[i].x1=ctl->X[i].x1+x;
         ctl->X[i].y1=ctl->X[i].y1+y;
         ctl->X[i].z1=ctl->X[i].z1+z;
         ctl->X[i].x2=ctl->X[i].x2+x;
         ctl->X[i].y2=ctl->X[i].y2+y;
         ctl->X[i].z2=ctl->X[i].z2+z;
      } else if(ctl->X[i].kind==KIND_TEXT) {
         ctl->X[i].x=ctl->X[i].x+x;
         ctl->X[i].y=ctl->X[i].y+y;
         ctl->X[i].z=ctl->X[i].z+z;
      }
   }

   ctl->xmax=ctl->xmax+x;
   ctl->ymax=ctl->ymax+y;
   ctl->zmax=ctl->zmax+z;
   ctl->xmin=ctl->xmin+x;
   ctl->ymin=ctl->ymin+y;
   ctl->zmin=ctl->zmin+z;

}
/**************************************************************************/
/*                                                                        */
/* This routine allocate space for list of postscript objects. If space   */
/* is already allocated, a new list is allocated, the contents copied and */
/* the old list is deallocated.                                           */
/*                                                                        */
/**************************************************************************/
void gveps_allo(gvepsdef *ctl) {
   obj *X;
   obj *Y;
   int  i,m,n;

   if(ctl->X==NULL) {
      if(ctl->size==SIZE_HUGE) {
         n=1000000;
      } else if(ctl->size==SIZE_LARGE) {
         n=10000;
      } else {
         n=100;
      }
#ifdef DEBUG
      fprintf(stderr,"Allocating list for %i objects\n",n);
#endif
      ctl->X=(obj *)malloc(sizeof(obj)*n);
      ctl->maxobj=n;
   } else {
      n=ctl->maxobj;
      m=2*n;
      ctl->maxobj=m;
#ifdef DEBUG
      fprintf(stderr,"Reallocating list for %i objects\n",m);
#endif
      X=(obj *)malloc(sizeof(obj)*m);
      for(i=0; i<n; i++) X[i]=ctl->X[i];
      Y=ctl->X;
      ctl->X=X;
      free(Y);
   }

}
/**************************************************************************/
/*                                                                        */
/* This routine sorts objects with respect to z coordinate.               */
/*                                                                        */
/**************************************************************************/
void gveps_sort(obj X[], int n) {
   obj Z;
   int i,j,k;

   for(i=0; i<n; i++) {
      k=i;
      for(j=i; j<n; j++) if(X[j].z<X[k].z) k=j;
      Z=X[i];
      X[i]=X[k];
      X[k]=Z;
   }

}
/**************************************************************************/
/*                                                                        */
/* This routine trims the to the surface of the atoms.                    */
/*                                                                        */
/**************************************************************************/
void gveps_trim(obj X[], int n) {
   double dx,dy,dz;
   double ux,uy,uz;
   double r;
   int    k1,k2;
   int    i,j;

   for(i=0; i<n; i++) {
      if(X[i].kind==KIND_BOND) {
#ifdef DEBUG
         fprintf(stderr,"Object %i is bond\n",i);
#endif
         k1=-1;
         for(j=0; j<n; j++) {
            if(X[j].kind==KIND_ATOM) {
               dx=X[i].x-X[j].x;
               dy=X[i].y-X[j].y;
               dz=X[i].z-X[j].z;
               r=sqrt(dx*dx+dy*dy+dz*dz);
               if(r<0.5) k1=j;
#ifdef DEBUG
               fprintf(stderr,"Distance for (%i,%i) is %.3f\n",i,j,r);
#endif
            }
         }
#ifdef DEBUG
         fprintf(stderr,"Object k1=%i is atom\n",k1);
#endif
         k2=-1;
         for(j=0; j<n; j++) {
            if(X[j].kind==KIND_ATOM) {
               dx=X[i].x1-X[j].x;
               dy=X[i].y1-X[j].y;
               dz=X[i].z1-X[j].z;
               r=sqrt(dx*dx+dy*dy+dz*dz);
               if(r<0.5) k2=j;
#ifdef DEBUG
               fprintf(stderr,"Distance for (%i,%i) is %.3f\n",i,j,r);
#endif
            }
         }
#ifdef DEBUG
         fprintf(stderr,"Object k2=%i is atom\n",k2);
#endif
         if(k1!=-1 && k2!=-1) {
#ifdef DEBUG
            fprintf(stderr,"Bond %i connects atoms %i and %i\n",i,k1,k2);
#endif
            ux=X[k2].x-X[k1].x;
            uy=X[k2].y-X[k1].y;
            uz=X[k2].z-X[k1].z;
            r=sqrt(ux*ux+uy*uy+uz*uz);
            ux=ux/r;
            uy=uy/r;
            uz=uz/r;
            X[i].x=X[k1].x+(X[k1].radius-0.015)*ux;
            X[i].y=X[k1].y+(X[k1].radius-0.015)*uy;
            X[i].z=X[k1].z+(X[k1].radius-0.015)*uz;
            X[i].x1=X[k2].x-(X[k2].radius-0.015)*ux;
            X[i].y1=X[k2].y-(X[k2].radius-0.015)*uy;
            X[i].z1=X[k2].z-(X[k2].radius-0.015)*uz;
         }
      }
   }
}
