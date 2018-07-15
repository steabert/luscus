/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<gtk/gtk.h>
#ifdef GTK_GLEXT
#include<gtk/gtkgl.h>
#else
#include<GL/glu.h>
#include<math.h>
#endif
#include<GL/gl.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

void luscus_gtk_draw_sphere(char solid, double radius, int slices, int stacks)
{
#ifdef GTK_GLEXT
  gdk_gl_draw_sphere((gboolean) solid, radius, slices, stacks);
#else
  GLUquadricObj *qobj;
  qobj = gluNewQuadric();
  if (solid) gluQuadricDrawStyle(qobj, GLU_FILL);
  else gluQuadricDrawStyle(qobj, GLU_LINE);
  
  gluSphere(qobj, radius, slices, stacks);

  gluDeleteQuadric(qobj);
#endif
}

#ifdef GTK_GLEXT
void luscus_gtk_draw_torus(char solid, double in_rad, double out_rad, int nsides, int nrings)
{

  gdk_gl_draw_torus((gboolean) solid, in_rad, out_rad, nsides, nrings);
/*  int i, j;
  double phi1, chi1;
  double cphi0, sphi0, cchi0, schi0;
  double cphi1, sphi1, cchi1, schi1;
  double dphi, dchi;

  cphi1 = 0.0;
  sphi1 = 1.0;
  cchi1 = 0.0;
  schi1 = 1.0;
  for (i = 0; i < nrings; i++, chi1 = M_PI * (double) (i+1))
  {
    cchi0 = cchi1;
    schi0 = schi1;
    cchi1 = cos(chi1);
    schi1 = sin(chi1);
    for(j = 0; j < nsides; j++, phi1 = M_PI * (double) (j+1))
    {
      cphi0 = cphi1;
      sphi0 = sphi1;
      cphi1 = cos(phi1);
      sphi1 = sin(phi1);

      if (solid) gl_Begin(GL_QUADS);
      else glBegin(GL_LINE_LOOP);

      glNormal3f();
      glVertex3f();

      glEnd();
    }
  }*/
}

void luscus_gtk_draw_teapot(char solid, double scale)
{
  gdk_gl_draw_teapot((gboolean) solid, scale);
}

void luscus_gtk_draw_cube(char solid, double size)
{
  gdk_gl_draw_cube((gboolean) solid, size);
}

#endif

