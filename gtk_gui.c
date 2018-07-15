/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdlib.h>
#include<string.h>
#include<gtk/gtk.h>
#include<math.h>
#ifdef GTK_GLEXT
#include<gtk/gtkgl.h>
#else
#ifdef WINDOWS
#include<gdk/gdkwin32.h>
#else
#include<GL/glx.h>
#include<gdk/gdkx.h>
#endif
#endif
#include<GL/gl.h>
#include<GL/glu.h>
#include"luscus.h"
#include"gv.h"
#include"gv_gtk.h"
#include"gv_glx.h"
#include"gv_functions.h"
#include"molcas024001.xpm"
#include "gveps.h"
#include"surface.h"

#define MAXCONNECTED 16

/*#ifdef GTK_GLEXT
#define BEGIN_OGL_STUFF \
        GdkGLContext *glcontext = gtk_widget_get_gl_context(da); \
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(da); \
        if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) return FALSE

#define END_OGL_STUFF gdk_gl_drawable_gl_end(gldrawable);
#define SWAP_OGL_BUFFERS gdk_gl_drawable_swap_buffers(gldrawable);
#else

#ifdef LINUX
typedef struct glx_conf
{
  Display *display;
  GLXContext context;
  Colormap xcolormap;
  GdkVisual* visual;
} GLXVIS;*/

#ifndef GTK_GLEXT
#ifdef LINUX

GLXVIS *tmpglxvis;
GLXVIS glxvis;

#endif
#endif

/*#ifdef GTK2
#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(da); \
        Display *display = gdk_x11_display_get_xdisplay(gdk_window_get_display(window)); \
        int id = gdk_x11_drawable_get_xid(window); \
        if (!glXMakeCurrent(display, id, glxvis.context)) return FALSE
#endif
#ifdef GTK3
#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(da); \
        Display *display = gdk_x11_display_get_xdisplay(gdk_window_get_display(window)); \
        int id = gdk_x11_window_get_xid(window); \
        if (!glXMakeCurrent(display, id, glxvis.context)) return FALSE
#endif
#define END_OGL_STUFF

#define SWAP_OGL_BUFFERS glXSwapBuffers(display, id);

#else
typedef struct glx_conf
{
  Display *display;
  int id;
  this_can_not_be_compiled context;
} GLXVIS;
#endif

#endif*/


static gboolean luscus_draw_init(GtkWidget*, gpointer);
#ifdef GTK3
static gboolean luscus_draw_reshape(GtkWidget*, GdkEvent*, gpointer);
static gboolean luscus_draw_display(GtkWidget*, cairo_t*, gpointer);
#endif
#ifdef GTK2
static gboolean luscus_draw_reshape(GtkWidget*, GdkEventConfigure*, gpointer);
static gboolean luscus_draw_display(GtkWidget*, GdkEventExpose*, gpointer);
#endif
static gboolean luscus_gtk_mouse_button_press(GtkWidget*, GdkEventButton*, gpointer);
static gboolean luscus_gtk_mouse_button_release(GtkWidget*, GdkEventButton*, gpointer);
static gboolean luscus_gtk_mouse_move(GtkWidget*, GdkEventMotion*, gpointer);
static gboolean luscus_gtk_mouse_move(GtkWidget*, GdkEventMotion*, gpointer);
void luscus_gtk_mouse_scroll(GtkWidget*, GdkEventScroll*, gpointer);

int mcubes(surface_t*, int, int, int, double[], double[], double[3], double, double[3], int);
void msrf_Draw1(multisrf_t*, color_t[], int, FILE*, int, FILE*, gvepsdef*);


GtkWidget *drawingarea;
GtkWidget *statusbar1;
GtkWidget *statusbar2;
GdkCursor *write_cursor;

static int vibr_count = 0;
static int is_is;
static double rvibr = 0.0;

static int slices=15;
static int stacks=15;
static int shift_move=0;
static int shift_move_x_b;
static int shift_move_y_b;
static int shift_move_x_e;
static int shift_move_y_e;

static int last_fragment = 0;

static int vec_rot = 0;
static int connected[MAXCONNECTED];
static int selected_go_num = 0;
static int selected_go_type = 0;
static int fullscreen = 0;
static double molecule_diameter;
static double x_size, y_size;
static double scale = 1.0;
static double expand_or_retract;
static double screen_width, screen_height;
static double camera_x = 0.0, camera_y = 0.0;
static double x_prev, y_prev;
static double tot_ang = 0.;
static double keepmvm[16]; /*modelview matrix for the orbital display*/
static double ang_rot=0.0;
GLdouble mvm[16];
multisrf_t* All_Surfaces;
surface_t** surf; /* dynamically allocated array */
static GLuint list_3d;
static char redraw_3d;
FILE *fl_povray;
FILE *fl_ps;
int textbox_state;
static int itextbox = 0;

typedef char atn[4];

gvepsdef ctl;

/*communication with other functions*/

void get_sizes_data(double *ret_x_size, double *ret_y_size, double *ret_scale, double *ret_camera_x, double *ret_camera_y, double *ret_mol_diameter)
{
  *ret_x_size = x_size;
  *ret_y_size = y_size;
  *ret_scale = scale;
  *ret_camera_x = camera_x;
  *ret_camera_y = camera_y;
  *ret_mol_diameter = molecule_diameter;
}

#ifndef GTK_GLEXT
GLXVIS *get_glx_conf(void)
{
  return tmpglxvis;
}
#endif

void rerender_3d(void)
{
  redraw_3d = 1;
}

unsigned int get_list_3d(void)
{
  return list_3d;
}

int find_connected(int iat)
{
  int i;
  int iconnected = 0;

  for(i = 0; i < m->nbond && iconnected < MAXCONNECTED; i++)
  {
    if (iat == m->bond[i].iat1)
      connected[iconnected++] = m->bond[i].iat2;
    if (iat == m->bond[i].iat2)
      connected[iconnected++] = m->bond[i].iat1;
  }
  return iconnected;
}

int get_vec_rot(void)
{
  return vec_rot;
}

GLdouble *get_main_mvm(void)
{
  return mvm;
}

void luscus_gtk_set_vec_rot(int new_rot)
{
  vec_rot = new_rot;
}

int get_number_of_slices(void)
{
  return slices;
}

int get_number_of_stacks(void)
{
  return stacks;
}

void set_number_of_slices(int newslices)
{
  slices = newslices;
}

void set_number_of_stacks(int newstacks)
{
  stacks = newstacks;
}

void gv_gtk_get_screen_size(gint *width, gint *height)
{
  GtkAllocation allocation;
  gtk_widget_get_allocation(drawingarea, &allocation);

  *width = allocation.width;
  *height = allocation.height;
  return;
}

void change_scale(int isc)
{
  if (isc == 1) scale*=1.1;
  else if (isc == -1) scale/=1.1;
  else if (isc == 0) scale = 1.0;
  set_scale();
}

void add_surface(MOL *m, XYZ p0, XYZ p1, XYZ p2)
{
  int i;
  double v1len, v1v2;

/*  rr_cross(&norm[0], &norm[1], &norm[2],
           p0[0], p1[0], p2[0],
           p0[1], p1[1], p2[1],
           p0[2], p1[2], p2[2]);

  rr_cross(&p2[0], &p2[1], &p2[2],
           p0[0], p1[0], p2[0],
           p0[1], p1[1], p2[1],
          norm[0], norm[1], norm[2]);*/

  allocate_surfaces(m, m->nsurf + 1);
/*  v1len = 0.0;
  v1v2 = 0.0;
  for(i = 0; i < 3; i++)
  {
    v1v2 += (p2[i] - p0[i])*(p1[i] - p0[i]);
    v1len += (p1[i] - p0[i])*(p1[i] - p0[i]);
  }*/
  for(i = 0; i < 3; i++)
  {
    m->surf1[m->nsurf-1][i] = p0[i];
/*    m->surf2[m->nsurf-1][i] = p0[i] + (p1[i] - p0[i]) * 2.0 * molecule_diameter / v1len;*/
/*    m->surf3[m->nsurf-1][i] = p0[i] + (p2[i] - p0[i] - (p1[i] - p0[i]) * v1v2 / v1len);*/ /*Gram-Schmidt 23.01.2014*/
    m->surf2[m->nsurf-1][i] = p1[i];/*Gram-Schmidt is moved to pl3to4*//*2016.09.15*/
    m->surf3[m->nsurf-1][i] = p2[i];/*Gram-Schmidt is moved to pl3to4*//*2016.09.15*/
  }
  surface_set_default_values(m, m->nsurf-1);
  luscus_gtk_update_3Dobject_info();
}

void add_vector(MOL *m, XYZ p0, XYZ p1)
{
  int i;
  allocate_vectors(m, m->nvector + 1);
  for(i = 0; i < 3; i++)
  {
    m->vector1[m->nvector-1][i] = p0[i];
    m->vector2[m->nvector-1][i] = p1[i];
  }
  vector_set_default_values(m, m->nvector-1);
  m->vector_color[m->nvector-1][3] = 1.0;
  luscus_gtk_update_3Dobject_info();
}

void add_sphere(MOL *m, XYZ p0, double r)
{
  int i;
  allocate_spheres(m, m->nsphere+1);
  for(i = 0; i < 3; i++)
    m->sphere_center[m->nsphere-1][i] = p0[i];
  sphere_set_defeault_values(m, m->nsphere-1);
  m->sphere_radius[m->nsphere-1] = r;
  luscus_gtk_update_3Dobject_info();
}

void add_triangle(MOL *m, XYZ p0, XYZ p1, XYZ p2)
{
  int i;
  allocate_triangles(m, m->ntriangle+1);
  for(i = 0; i < 3; i++)
  {
    m->triangle1[m->ntriangle-1][i] = p0[i];
    m->triangle2[m->ntriangle-1][i] = p1[i];
    m->triangle3[m->ntriangle-1][i] = p2[i];
  }
  triangle_set_defeault_values(m, m->ntriangle-1);
  luscus_gtk_update_3Dobject_info();
}

void add_cell(MOL *m, XYZ p0,  XYZ p1,  XYZ p2,  XYZ p3)
{
  int i;
  allocate_cells(m, m->ncells+1);
  for(i = 0; i < 3; i++)
  {
    m->cell1[m->ncells-1][i] = p0[i];
    m->cell2[m->ncells-1][i] = p1[i];
    m->cell3[m->ncells-1][i] = p2[i];
    m->cell4[m->ncells-1][i] = p3[i];
  }
  cell_set_default_values(m, m->ncells-1);
  luscus_gtk_update_3Dobject_info(); 
}

void add_fragment(int ifrag)
{
  MOL f;
  int i, j, k;
  int i0;
  double xm, rm;
  double d, dist;
  int iconnected;
  int shift = 0; /*in the case the shift button is pressed, this should be 1*/
  double A[3][3], xi3[3];
  XYZ cm;
  XYZ tmp;
  XYZ vec, vecf;
  XYZ xyz0;

  last_fragment = ifrag;
/*  int natom_old = m->natom;
  if (m->n_selected == 0)
  {
  }*/
  f = load_fragment(ifrag);
/*  printf("bonds:\n");
  for(i = 0; i < f.nbond; i++)
    printf("bond#%d: %d - %d   %d\n", i, f.bond[i].iat1, f.bond[i].iat2, f.bond[i].bond_type);*/

  if (m->n_selected == 0)
  {
    i0 = 0;
    xm = 0.F;

    if (m->natom)
    {
      for (i = 0; i < m->natom; i++)
      {
        if (m->xyz[i][0] < xm)
        {
          xm = m->xyz[i][0];
          i0 = i;
        }
      }
      get_center(&cm[0], &cm[1], &cm[2]);
      for(i = 0; i < 3; i++) cm[i] -= m->xyz[i0][i];
      rm = sqrt(cm[0]*cm[0] + cm[1]*cm[1] + cm[2]*cm[2]);
      if(rm < 0.1)
      {
        xm=0.5;
        rm=1;
      }

      for(k = 1; k < 100; k++)
      {
        dist = 10000.F;
        for(i = 1; i < f.natom; i++)
          for(j = 0; j < 3; j++) tmp[j] = f.xyz[i][j] + (double) k;
          for(j = 0; j < m->natom; j++)
          {
            d = rr_dist(m->xyz[j], tmp);
            if (d < dist) dist = d;
          }
          if (dist > 2.5F) break;
      }

      allocate_atoms(m, m->natom + f.natom - 1);

      for(i = 1; i < f.natom; i++)
      {
        for(j = 0; j < 3; j++)
          m->xyz[m->natom - f.natom + i][j] = (double) k + f.xyz[i][j] + m->xyz[i0][j];
        m->elem[m->natom - f.natom + i].name = strdup(f.elem[i].name); /*--------------------*/
        for(j = 0; j < 4; j++)
          m->elem[m->natom - f.natom + i].color[j] = f.elem[i].color[j];
        m->elem[m->natom - f.natom + i].vdw_rad = f.elem[i].vdw_rad;
        m->elem[m->natom - f.natom + i].bond_rad = f.elem[i].bond_rad;
        m->elem[m->natom - f.natom + i].valency = f.elem[i].valency;

        /*insert atom into atom_list*/
        insert_atom_into_list(m->natom - f.natom + i);
      }

      /*add bonds*/
      for(i = 0; i < f.nbond; i++)
        add_bond(m, m->natom - f.natom + f.bond[i].iat1, m->natom - f.natom + f.bond[i].iat2, f.bond[i].bond_type);
    }
    else
    {
      allocate_atoms(m, f.natom - 1);

      for(i = 1; i < f.natom; i++)
      {
        for(j = 0; j < 3; j++)
          m->xyz[i - 1][j] = f.xyz[i][j];
        m->elem[i - 1].name = strdup(f.elem[i].name);  /*-----------------------------------*/
        for(j = 0; j < 4; j++)
          m->elem[i - 1].color[j] = f.elem[i].color[j];
        m->elem[i - 1].vdw_rad = f.elem[i].vdw_rad;
        m->elem[i - 1].bond_rad = f.elem[i].bond_rad;
        m->elem[i - 1].valency = f.elem[i].valency;

        /*insert atom into atom_list*/
        insert_atom_into_list(i-1);
      }

      for(i = 0; i < f.nbond; i++)
        add_bond(m, m->natom - f.natom + f.bond[i].iat1, m->natom - f.natom + f.bond[i].iat2, f.bond[i].bond_type);

      /*add bonds*/
    }

    free_fragment_data(f);
  }
  if(m->n_selected > 0)
  {
    i0=m->selected[0];

    /*I. define direction vector of the existing molecule*/
    if (m->n_selected != 2)
    {
/*
   1. find connected to this atom.
   2. if one is connected - use this bond as reference
   3. if several connected - compute 'vector'
   4. orient fragment accordingly
*/
      vec[0]=0; vec[1]=0; vec[2]=0;
      iconnected=find_connected(i0);
      for(i=0; i < iconnected; i++)
        for(j = 0; j < 3; j++)
          vec[j] += m->xyz[connected[i]][j] - m->xyz[i0][j];

      if (iconnected == 0)
      {
        vec[0]=1.0; vec[1]=0.0; vec[2]=0.0;
      }
    }
    else
    {
       /*find direction of the bond*/
      for(i = 0; i < 3; i++)
       	vec[i] = m->xyz[i0][i] - m->xyz[m->selected[1]][i];
       /*delete second atom*/
       if (i0 > m->selected[1]) i0--;
       delete_atom(m->selected[1]);
    }
    unselect_all();

    /*II. find direction vector in fragment */
    for(i = 0; i < 3; i++)
      vecf[i] = f.xyz[0][i] - f.xyz[2][i];

    /*III. define translation vector*/
    for(i = 0; i < 3; i++)
      xyz0[i] = m->xyz[i0][i] - f.xyz[0][i];

    /*IV. define rotation matrix*/
/*    for(i = 0; i < 3; i++)
    {
      xyz[0][i] = vecf[i];
      xyz[1][i] = 0.0;
      xyz[2][i] = vec[i];
    }*/

    build_rotation_matrix(vec, vecf, A);

    /*V. do rotation of the fragment*/
    for(j = 2; j < f.natom; j++)
    {
      for(i = 0; i < 3; i++)
        xi3[i] = f.xyz[j][i] + xyz0[i];

      translate_to_origin(xi3,m->xyz[i0],A);

      for(i = 0; i < 3; i++)
        f.xyz[j][i] = xi3[i];
    }

    /*VI. add atoms from the fragment + translate for the translation vector*/
    allocate_atoms(m, m->natom + f.natom - 2);

    for(i = 2; i < f.natom; i++)
    {
      for(j = 0; j < 3; j++)
        m->xyz[m->natom - f.natom + i][j] = f.xyz[i][j];
      m->elem[m->natom - f.natom + i].name = strdup(f.elem[i].name); /*-------------------*/
      for(j = 0; j < 4; j++)
        m->elem[m->natom - f.natom + i].color[j] = f.elem[i].color[j];
      m->elem[m->natom - f.natom + i].vdw_rad = f.elem[i].vdw_rad;
      m->elem[m->natom - f.natom + i].bond_rad = f.elem[i].bond_rad;
      m->elem[m->natom - f.natom + i].valency = f.elem[i].valency;

      /*insert atom into atom_list*/
      insert_atom_into_list(m->natom - f.natom + i);
    }

    /*add bonds from fragment*/
    for(i = 0; i < f.nbond; i++)
      if (f.bond[i].iat1 != 1 && f.bond[i].iat2 != 1)
        add_bond(m, m->natom - f.natom + f.bond[i].iat1,
                    m->natom - f.natom + f.bond[i].iat2, f.bond[i].bond_type);
    /*add bond between new fragment and selected atom*/
    add_bond(m, i0, m->natom - f.natom + 2, 1);

    free_fragment_data(f); 
  }
  draw_all_pixdata();
  append_backup();
  set_scale();
  if (Input_Data.automatic_rebonding) rebond();
}

int get_last_fragment(void)
{
  return last_fragment;
}

void textbox_delete_last_char(MOL *m)
{
  int ichar;
  if (!m) return;
  if (!m->ntextboxes) return;
  if (!m->textboxes[m->ntextboxes-1].message) return;
  if (m->textboxes[m->ntextboxes-1].message[0] == 0) return;

  ichar = strlen(m->textboxes[m->ntextboxes-1].message);
  m->textboxes[m->ntextboxes-1].message[ichar-1] = 0;
  m->textboxes[m->ntextboxes-1].message = (char*) realloc(m->textboxes[m->ntextboxes-1].message, ichar * sizeof(char));
}

void textbox_delete_insert_char(MOL *m, char c)
{
  int ichar;
  if (!m) return;
  if (!m->ntextboxes) return;
  if (!m->textboxes[m->ntextboxes-1].message)
  {
    m->textboxes[m->ntextboxes-1].message = (char*) malloc(2 * sizeof(char));
    m->textboxes[m->ntextboxes-1].message[0] = c;
    m->textboxes[m->ntextboxes-1].message[1] = 0;
  }
  else
  {
    ichar = strlen(m->textboxes[m->ntextboxes-1].message) + 2;
    m->textboxes[m->ntextboxes-1].message = (char*) realloc(m->textboxes[m->ntextboxes-1].message, ichar * sizeof(char));
    m->textboxes[m->ntextboxes-1].message[ichar-2] = c;
    m->textboxes[m->ntextboxes-1].message[ichar-1] = 0;
  }
}

void luscus_set_textbox_number(int niboxnum)
{
  if (niboxnum >= m->ntextboxes)
  {
    textbox_state = 0;
    return;
  }
  itextbox = niboxnum;
}

void set_go_selected(int go_type, int go_num)
{
  selected_go_num = go_num;
  selected_go_type = go_type;
}

void unset_go_selected(void)
{
  selected_go_num = 0;
  selected_go_type = 0;
}

void luscus_set_cursor_text(void)
{
  GdkWindow *window = gtk_widget_get_window(drawingarea);
  GdkDisplay *display = gtk_widget_get_display(drawingarea);
  GdkCursor *cursor = gdk_cursor_new_for_display(display, GDK_XTERM);

  gdk_window_set_cursor(window, cursor);
}

void luscus_set_cursor_default(void)
{
  GdkWindow *window = gtk_widget_get_window(drawingarea);
  gdk_window_set_cursor(window, NULL);
}

void animate(int flag)
{
  float ang=0.5;
  if (n_geometries > 1)
  {
    if (!m) return;
    if (igeo == n_geometries-1) callback_geo_play(NULL, NULL);
    igeo++;
    if (igeo >= n_geometries) igeo = n_geometries - 1;
    get_new_section();
    luscus_gtk_update_geo_info();
    set_current_graph_data();
    set_scale();
    rerender_3d();
  }
  else
  {
    ang_rot=ang;
    glRotated(ang_rot, 0.5, 0.5, 1.0);
  }
  redraw();
}

void do_symmetry(void)
{
  int i;
  int rc;
  char what[5]="full";
  double fthr = 0.05;
  int *outsym;
  static char outgroup[16];
  atn *atnam;
  int n_sym_elements;
  SymmetryElement X[8];
  XYZ p0, p1, p2;

  if (m->n_selected == 0)
  {
/*    printf("is_findsym = %d\n", is_findsym); fflush(stdout);
    is_findsym=!is_findsym;
    if (is_findsym)
    {*/
      /*place atom names in temporary array of correct format*/
      atnam = (atn*) malloc(sizeof(atn) * m->natom);
      for(i = 0; i < m->natom; i++)
      {
        strncpy(atnam[i], m->elem[i].name, 3);
        atnam[i][3] = 0;
      }

      outsym = (int*) malloc(sizeof(int) * m->natom);

      if (Input_Data.force_symmetry) 
      {
        fthr = 0.15;
        rc=findsym(0, m->natom, atnam, m->xyz, " ", 0,
                   what, 4, &fthr, 0, outgroup, outsym, 0);
        fthr = 0.05;
        Input_Data.force_symmetry = 0;
      }
      else
        rc=findsym(0, m->natom, atnam, m->xyz, " ", 0,
                   what, 4, &fthr, 0, outgroup, outsym, 0);

      if (outgroup[0] == 0)
      {
        luscus_gtk_pop_message_from_statusbar2();
       	luscus_gtk_push_message_to_statusbar2("No symmetry element has been found");
      }
      else luscus_gtk_pop_message_from_statusbar2();
      free(atnam);
      free(outsym);

      if (rc) return;
      if (outgroup[0])
      {
        luscus_gtk_pop_message_from_statusbar2();
       	luscus_gtk_push_message_to_statusbar2(outgroup);
      }

      /*destroy all geometry elements that allready exist!*/
      deallocate_vectors(m);
      deallocate_triangles(m);
      deallocate_spheres(m);
      deallocate_surfaces(m);
      deallocate_cells(m);
      luscus_gtk_update_upon_select_or_mark();
      luscus_gtk_update_3Dobject_info();
      /*add new geometry elements*/
      n_sym_elements = gvgrp(outgroup,X);

      for(i = 0; i < n_sym_elements; i++)
      {
        if (X[i].type==plane)
        {
          p0[0] = X[i].x[0];
          p0[1] = X[i].y[0];
          p0[2] = X[i].z[0];
          p1[0] = X[i].x[1];
          p1[1] = X[i].y[1];
          p1[2] = X[i].z[1];
          p2[0] = X[i].x[2];
          p2[1] = X[i].y[2];
          p2[2] = X[i].z[2];

          add_surface(m, p1, p0, p2);
        }
        if (X[i].type==axis)
        {
          p0[0] = X[i].x[0] * 3.0;
          p0[1] = X[i].y[0] * 3.0;
          p0[2] = X[i].z[0] * 3.0;
          p1[0] = X[i].x[1] * 3.0;
          p1[1] = X[i].y[1] * 3.0;
          p1[2] = X[i].z[1] * 3.0;
          add_vector(m, p0, p1);
          m->sharpness[m->nvector-1] = 0.0;
          m->radius[m->nvector-1] = 0.025;
        }
      }
      luscus_gtk_update_upon_select_or_mark();
/*    }*/
  }
  else if (m->n_selected == 1) do_inversion();
  else if (m->n_selected == 2 && vec_rot != 0) do_rotation();
  else if (m->n_selected == 2 && vec_rot == 0) do_translation();
  else if (m->n_selected == 3) do_mirroring();
  append_backup();
}

void do_rotation(void)
{
  int i, j, k, l;
  int natom = m->natom;
  int nbond = m->nbond;
  int i0 = m->selected[0];
  int i1 = m->selected[1];
  XYZ newatom;
  double dist2;
  int noaccepted;
  int iatom = 0;
  int *copied_atoms1 = NULL;
  int *copied_atoms2 = NULL;

  if (m->n_selected != 2) return;

  for(i = 1; i < vec_rot; i++)
  {
    if (m->n_marked)
    {
      for(j = 0, iatom = 0; j < m->n_marked; j++)
      {
        if (m->marked[j] != i0 && m->marked[j] != i1)
        {
/*          allocate_atoms(m, m->natom + 1);
          copy_atom_data(m->marked[j], m->natom - 1);*/
          for(k = 0; k < 3; k++) newatom[k] = m->xyz[m->marked[j]][k];
          transform_n_fold(vec_rot, i, m->xyz[i0], m->xyz[i1], newatom/* m->xyz[m->natom-1]*/);
          /*check for collisions*/
          noaccepted = 0;
          for(k = 0; k < m->natom; k++)
          {
            dist2 = (m->xyz[k][0] - newatom[0])*(m->xyz[k][0] - newatom[0]) +
                    (m->xyz[k][1] - newatom[1])*(m->xyz[k][1] - newatom[1]) +
                    (m->xyz[k][2] - newatom[2])*(m->xyz[k][2] - newatom[2]);
            if (dist2 < (m->elem[k].vdw_rad + m->elem[m->marked[j]].vdw_rad) * (m->elem[k].vdw_rad + m->elem[m->marked[j]].vdw_rad)) noaccepted = 1;
          }
          /*No collisions; accept new atom*/
          if (!noaccepted)
          {
            copied_atoms1 = (int*) realloc(copied_atoms1, (iatom+1) * sizeof(int));
            copied_atoms2 = (int*) realloc(copied_atoms2, (iatom+1) * sizeof(int));
            copied_atoms1[iatom] = m->marked[j];
            copied_atoms2[iatom] = m->natom;

            allocate_atoms(m, m->natom + 1);
            copy_atom_data(m->marked[j], m->natom - 1);
            for (k = 0; k < 3; k++) m->xyz[m->natom-1][k] = newatom[k];
            insert_atom_into_list(m->natom-1);
            iatom++;
          }
        }
      }

      /*copy bonds with marked atoms*/

      if (!Input_Data.automatic_rebonding)
        for(l = 0; l < nbond; l++)
          for(k = 0; k < iatom; k++)
            if (m->bond[l].iat1 == copied_atoms1[k])
              for(j = 0; j < iatom; j++)
                if (m->bond[l].iat2 == copied_atoms1[j])
                  add_bond(m, copied_atoms2[k], copied_atoms2[j], m->bond[l].bond_type);

      if (copied_atoms1) free(copied_atoms1);
      if (copied_atoms2) free(copied_atoms2);
      copied_atoms1 = NULL;
      copied_atoms2 = NULL;
    }
    else/*rotate all atoms*/
    {
      for(j = 0, iatom = 0; j < natom; j++)
      {
        if (j != i0 && j != i1)
        {
          for(k = 0; k < 3; k++) newatom[k] = m->xyz[j][k];
          transform_n_fold(vec_rot, i, m->xyz[i0], m->xyz[i1], newatom/* m->xyz[m->natom-1]*/);
          /*check for collisions*/
          noaccepted = 0;
          for(k = 0; k < m->natom; k++)
          {
            dist2 = (m->xyz[k][0] - newatom[0])*(m->xyz[k][0] - newatom[0]) +
                    (m->xyz[k][1] - newatom[1])*(m->xyz[k][1] - newatom[1]) +
                    (m->xyz[k][2] - newatom[2])*(m->xyz[k][2] - newatom[2]);
            if (dist2 < (m->elem[k].vdw_rad + m->elem[j].vdw_rad) * (m->elem[k].vdw_rad + m->elem[j].vdw_rad)) noaccepted = 1;
          }
          /*No collisions; accept new atom*/
          if (!noaccepted)
          {
            copied_atoms1 = (int*) realloc(copied_atoms1, (iatom+1) * sizeof(int));
            copied_atoms2 = (int*) realloc(copied_atoms2, (iatom+1) * sizeof(int));
            copied_atoms1[iatom] = j;
            copied_atoms2[iatom] = m->natom;

            allocate_atoms(m, m->natom + 1);
            copy_atom_data(j, m->natom - 1);
            for (k = 0; k < 3; k++) m->xyz[m->natom-1][k] = newatom[k];
            insert_atom_into_list(m->natom-1);
            iatom++;
          }
        }
      }

      /*copy all bonds*/
      if (!Input_Data.automatic_rebonding)
        for(l = 0; l < nbond; l++)
          for(k = 0; k < iatom; k++)
            if (m->bond[l].iat1 == copied_atoms1[k])
              for(j = 0; j < iatom; j++)
                if (m->bond[l].iat2 == copied_atoms1[j])
                  add_bond(m, copied_atoms2[k], copied_atoms2[j], m->bond[l].bond_type);

      if (copied_atoms1) free(copied_atoms1);
      if (copied_atoms2) free(copied_atoms2);
      copied_atoms1 = NULL;
      copied_atoms2 = NULL;
    }
  }

  if (Input_Data.automatic_rebonding) rebond();
}

void do_translation(void)
{
  int natom = m->natom;
  int nbond = m->nbond;
  int i, j, k;
  double len;
  XYZ tv; /*translation vector*/
  XYZ newatom;
  double dist2;
  int iatom = 0;
  int *copied_atoms1 = NULL;
  int *copied_atoms2 = NULL;
  int noaccepted;

  if (m->n_selected != 2) return;

  for(i = 0; i < 3; i++)
    tv[i] = m->xyz[m->selected[0]][i] - m->xyz[m->selected[1]][i];

  len = sqrt(tv[0]*tv[0] + tv[1]*tv[1] + tv[2]*tv[2]);

  for(i = 0; i < 3; i++)
    tv[i] *= Input_Data.symmetry_translation / len;

  if (m->n_marked)
  {
    /*copy atoms*/
    for(i = 0; i < m->n_marked; i++)
    {
      noaccepted = 0;
      for(j = 0; j < 3; j++) newatom[j] = m->xyz[m->marked[i]][j] - tv[j];
      for(j = 0; j < m->natom; j++)
      {
        dist2 = (m->xyz[j][0] - newatom[0]) * (m->xyz[j][0] - newatom[0]) + 
                (m->xyz[j][1] - newatom[1]) * (m->xyz[j][1] - newatom[1]) + 
                (m->xyz[j][2] - newatom[2]) * (m->xyz[j][2] - newatom[2]);
        if (dist2 < (m->elem[j].vdw_rad + m->elem[m->marked[i]].vdw_rad) * 
                    (m->elem[j].vdw_rad + m->elem[m->marked[i]].vdw_rad)) noaccepted=1; 
      }
      if (!noaccepted)
      {
        copied_atoms1 = (int*) realloc(copied_atoms1, (iatom+1) * sizeof(int));
        copied_atoms2 = (int*) realloc(copied_atoms2, (iatom+1) * sizeof(int));
        copied_atoms1[iatom] = m->marked[i];
        copied_atoms2[iatom] = m->natom;

        allocate_atoms(m, m->natom + 1);
        copy_atom_data(m->marked[i], m->natom - 1);
        for(k = 0; k < 3; k++) m->xyz[m->natom-1][k] = newatom[k];
        insert_atom_into_list(m->natom-1);
        iatom++;
      }
    }
    /*copy bonds*/
    if (!Input_Data.automatic_rebonding)
      for(i = 0; i < nbond; i++)
        for(j = 0; j < iatom; j++)
          if (m->bond[i].iat1 == copied_atoms1[j])
            for(k = 0; k < iatom; k++)
              if (m->bond[i].iat2 == copied_atoms1[k])
                add_bond(m, copied_atoms2[j], copied_atoms2[k], m->bond[i].bond_type);
  }
  else
  {
    /*copy atoms*/
    for(i = 0; i < natom; i++)
    {
      noaccepted = 0;
      for(j = 0; j < 3; j++) newatom[j] = m->xyz[i][j] - tv[j];
      for(j = 0; j < m->natom; j++)
      {
        dist2 = (m->xyz[j][0] - newatom[0]) * (m->xyz[j][0] - newatom[0]) + 
                (m->xyz[j][1] - newatom[1]) * (m->xyz[j][1] - newatom[1]) + 
                (m->xyz[j][2] - newatom[2]) * (m->xyz[j][2] - newatom[2]);
        if (dist2 < (m->elem[j].vdw_rad + m->elem[i].vdw_rad) * 
                    (m->elem[j].vdw_rad + m->elem[i].vdw_rad)) noaccepted=1; 
      }
      if (!noaccepted)
      {
        copied_atoms1 = (int*) realloc(copied_atoms1, (iatom+1) * sizeof(int));
        copied_atoms2 = (int*) realloc(copied_atoms2, (iatom+1) * sizeof(int));
        copied_atoms1[iatom] = i;
        copied_atoms2[iatom] = m->natom;

        allocate_atoms(m, m->natom + 1);
        copy_atom_data(i, m->natom - 1);
        for(k = 0; k < 3; k++) m->xyz[m->natom-1][k] = newatom[k];
        insert_atom_into_list(m->natom-1);
        iatom++;
      }
    }
    /*copy bonds*/
    if (!Input_Data.automatic_rebonding)
      for(i = 0; i < nbond; i++)
        for(j = 0; j < iatom; j++)
          if (m->bond[i].iat1 == copied_atoms1[j])
            for(k = 0; k < iatom; k++)
              if (m->bond[i].iat2 == copied_atoms1[k])
                add_bond(m, copied_atoms2[j], copied_atoms2[k], m->bond[i].bond_type);
  }
  if (copied_atoms1) free(copied_atoms1);
  if (copied_atoms2) free(copied_atoms2);
  if (Input_Data.automatic_rebonding) rebond();
}

void do_expand_or_retract(void)
{
  int i, j;
  double center[3];
  /*calculate center*/
  for(j = 0; j < 3; j++) center[j] = 0.0;
  for(i = 0; i < m->natom; i++)
  {
    for(j = 0; j < 3; j++)
      center[j] += m->xyz[i][j];
  }
  /*increase/decrease the distance to the center*/
  for(i = 0; i < m->natom; i++)
    for(j=0; j < 3; j++)
      m->xyz[i][j] += (m->xyz[i][j]-center[j]) * expand_or_retract;

  accumulate_motion(expand_or_retract);
  if (Input_Data.automatic_rebonding) rebond();
}

void transform_mirror(XYZ p0, XYZ p1, XYZ p2, XYZ px, XYZ pu)
/*  transform px by mirroring around plane, specified by p0, p1 and p2 */
{
  double al, be, r10, r20, q, t;
  r10=(p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
  r20=(p2[0]-p0[0])*(p2[0]-p0[0]) + (p2[1]-p0[1])*(p2[1]-p0[1]) + (p2[2]-p0[2])*(p2[2]-p0[2]);
  q  =(p1[0]-p0[0])*(p2[0]-p0[0]) + (p1[1]-p0[1])*(p2[1]-p0[1]) + (p1[2]-p0[2])*(p2[2]-p0[2]);
  
  be=q*( (p0[0]-px[0])*(p1[0]-p0[0]) + (p0[1]-px[1])*(p1[1]-p0[1]) + (p0[2]-px[2])*(p1[2]-p0[2]) ) -
   r10*( (p0[0]-px[0])*(p2[0]-p0[0]) + (p0[1]-px[1])*(p2[1]-p0[1]) + (p0[2]-px[2])*(p2[2]-p0[2]) );
  t=r20*r10 - q*q;
  if(fabs(t) < 0.00001) be=0;  
  else be=be/t; 
  
  al=-( (p0[0]-px[0])*(p1[0]-p0[0]) + (p0[1]-px[1])*(p1[1]-p0[1]) + (p0[2]-px[2])*(p1[2]-p0[2]) ) -
    be*((p2[0]-p0[0])*(p1[0]-p0[0]) + (p2[1]-p0[1])*(p1[1]-p0[1]) + (p2[2]-p0[2])*(p1[2]-p0[2]) );

  if(fabs(r10) < 0.00001) al=0;
  else   al=al/r10;
  
  pu[0]=2*(p0[0]+al*(p1[0]-p0[0]) + be*(p2[0]-p0[0]) )-px[0];
  pu[1]=2*(p0[1]+al*(p1[1]-p0[1]) + be*(p2[1]-p0[1]) )-px[1];
  pu[2]=2*(p0[2]+al*(p1[2]-p0[2]) + be*(p2[2]-p0[2]) )-px[2];
  return;  
}

void do_mirroring(void)
{
  int i, j, k, l;
  int natom = m->natom;
  int nbond = m->nbond;
/*  int nmatch;*/
  int collision;
  XYZ tmpxyz, tmpatom;
  double d2;
  int *copied_atoms1 = NULL;
  int *copied_atoms2 = NULL;

  if (m->n_selected != 3) return;

  if (m->n_marked)
  {
    /*copy atoms*/
/*    allocate_atoms(m, natom + m->n_marked);*/
    j = 0;
    for(i = 0; i < m->n_marked; i++)
    {
      if (m->marked[i] != m->selected[0] && m->marked[i] != m->selected[1] && m->marked[i] != m->selected[2])
      {
        for(k = 0; k < 3; k++) tmpxyz[k] = m->xyz[m->marked[i]][k];
        transform_mirror(m->xyz[m->selected[0]], m->xyz[m->selected[1]], m->xyz[m->selected[2]], tmpxyz, tmpatom);
        collision = 0;
        for(k = 0; k < natom; k++)
        {
          d2 = (m->xyz[k][0] - tmpatom[0]) * (m->xyz[k][0] - tmpatom[0]) + (m->xyz[k][1] - tmpatom[1]) * (m->xyz[k][1] - tmpatom[1]) + (m->xyz[k][2] - tmpatom[2]) * (m->xyz[k][2] - tmpatom[2]);
          if (d2 < (m->elem[k].vdw_rad + m->elem[m->marked[i]].vdw_rad) * (m->elem[k].vdw_rad + m->elem[m->marked[i]].vdw_rad)) collision = 1;
        }
        if (!collision)
        {
          copied_atoms1 = (int*) realloc(copied_atoms1, (j+1) * sizeof(int));
          copied_atoms2 = (int*) realloc(copied_atoms2, (j+1) * sizeof(int));
          copied_atoms1[j] = m->marked[i];
          copied_atoms2[j] = natom+j;
          allocate_atoms(m, natom+j+1);
          copy_atom_data(m->marked[i], natom+j);
          for(k = 0; k < 3; k++) m->xyz[natom+j][k] = tmpatom[k];
          insert_atom_into_list(m->natom-1);
          j++;
        }
      }
    }
    /*copy bonds*/
    if (!Input_Data.automatic_rebonding)
      for(i = 0; i < nbond; i++)
        for(k = 0; k < j; k++)
          if (m->bond[i].iat1 == copied_atoms1[k])
            for(l = 0; l < j; l++)
              if (m->bond[i].iat2 == copied_atoms1[l])
                add_bond(m, copied_atoms2[k], copied_atoms2[l], m->bond[i].bond_type);
  }
  else
  {
    /*copy atoms*/
    for(i = 0, j = 0; i < natom; i++)
    {
      if (i != m->selected[0] && i != m->selected[1] && i != m->selected[2])
      {

        for(k = 0; k < 3; k++) tmpxyz[k] = m->xyz[i][k];
        transform_mirror(m->xyz[m->selected[0]], m->xyz[m->selected[1]], m->xyz[m->selected[2]], tmpxyz, tmpatom);
        collision = 0;
        for(k = 0; k < natom; k++) /*check for collisions*/
        {
          d2 = (m->xyz[k][0] - tmpatom[0]) * (m->xyz[k][0] - tmpatom[0]) + (m->xyz[k][1] - tmpatom[1]) * (m->xyz[k][1] - tmpatom[1]) + (m->xyz[k][2] - tmpatom[2]) * (m->xyz[k][2] - tmpatom[2]);
          if (d2 < (m->elem[k].vdw_rad + m->elem[i].vdw_rad) * (m->elem[k].vdw_rad + m->elem[i].vdw_rad)) collision = 1;
        }
        if (! collision)
        {
          copied_atoms1 = (int*) realloc(copied_atoms1, (j+1) * sizeof(int));
          copied_atoms2 = (int*) realloc(copied_atoms2, (j+1) * sizeof(int));
          copied_atoms1[j] = i;
          copied_atoms2[j] = natom+j;
          allocate_atoms(m, natom+j+1);
          copy_atom_data(i, natom+j);
          for(k = 0; k < 3; k++) m->xyz[natom+j][k] = tmpatom[k];
          insert_atom_into_list(m->natom-1);

          j++;
        }
      }
    }
    /*copy bonds*/
    if (!Input_Data.automatic_rebonding)
      for(i = 0; i < nbond; i++)
        for(k = 0; k < j; k++)
          if (m->bond[i].iat1 == copied_atoms1[k])
            for(l = 0; l < j; l++)
              if (m->bond[i].iat2 == copied_atoms1[l])
                add_bond(m, copied_atoms2[k], copied_atoms2[l], m->bond[i].bond_type);
  }
  if (copied_atoms1) free(copied_atoms1);
  if (copied_atoms2) free(copied_atoms2);

  if (Input_Data.automatic_rebonding) rebond();
}

void change_element_by_one(int ichange)
{
  int i;
  int ielem;
  int npos;
  gchar *tmp;
  if (!m) return;
  if (m->n_selected != 1) return;

  for (ielem = 0; ielem < number_of_elements && g_strcmp0(e[ielem].name, m->elem[m->selected[0]].name) != 0; ielem++); /*searching for correct element*/
  for (i = 0; i < def_elem.n_default_elem && ielem != def_elem.default_elem[i]; i++);
  if (i == def_elem.n_default_elem) i = 0;
/*  ielem += ichange;*/
  npos = (i + ichange) % def_elem.n_default_elem;
  if (npos < 0) npos += def_elem.n_default_elem;  /*for some implementations of modulo function negative values might appear!*/
  ielem = def_elem.default_elem[npos];
  if (ielem >= number_of_elements) ielem = 0;
  if (ielem < 0) ielem = number_of_elements - 1;
/**/
  m->elem[m->selected[0]].color[0] = e[ielem].color[0];
  m->elem[m->selected[0]].color[1] = e[ielem].color[1];
  m->elem[m->selected[0]].color[2] = e[ielem].color[2];
  m->elem[m->selected[0]].color[3] = e[ielem].color[3];
  m->elem[m->selected[0]].vdw_rad = e[ielem].vdw_rad;
  m->elem[m->selected[0]].bond_rad = e[ielem].bond_rad;
  m->elem[m->selected[0]].valency = e[ielem].valency;
  free(m->elem[m->selected[0]].name);
  m->elem[m->selected[0]].name = strdup(e[ielem].name);
  change_atom_parameters_in_list(m->selected[0]);

  luscus_gtk_pop_message_from_statusbar2();
  tmp = g_strdup_printf("%s %d (%f %f %f)", m->elem[m->selected[0]].name, m->selected[0]+1, m->xyz[m->selected[0]][0],  m->xyz[m->selected[0]][1], m->xyz[m->selected[0]][2]);
  luscus_gtk_push_message_to_statusbar2(tmp);
  g_free(tmp);
  append_backup();
}

void ViewPoint(int i)
{
  int j;
  FILE *status;
  char tmpstr[64];
  char err[]="Error during reading ViewPoint file\n";
  double z_depth;
  if (i == 1)
  {
    if((status=fopen(".ViewPoint","r"))==NULL)
    {
      make_warning("no saved ViewPoint\n");
      return;
    }
    else
    {
/*      double xx, yy;
      double ss, xsz, ysz;*/
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err); return;}
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr,"%lf %lf",&screen_width, &screen_height);
  
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr, "%lf %lf", &x_size, &y_size);

      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr,"%lf",&scale);

      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr, "%lf", &z_depth);

      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr, "%lf %lf", &camera_x, &camera_y);
  
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr,"%lf %lf %lf %lf",&keepmvm[0],&keepmvm[1],&keepmvm[2],&keepmvm[3]);
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr,"%lf %lf %lf %lf",&keepmvm[4],&keepmvm[5],&keepmvm[6],&keepmvm[7]);
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);}
      sscanf(tmpstr,"%lf %lf %lf %lf",&keepmvm[8],&keepmvm[9],&keepmvm[10],&keepmvm[11]);
      if(fgets(tmpstr,sizeof(tmpstr),status)==NULL) {printf ("%s",err);} 
      sscanf(tmpstr,"%lf %lf %lf %lf",&keepmvm[12],&keepmvm[13],&keepmvm[14],&keepmvm[15]);

      molecule_diameter=z_depth/2.0;
      /*load keepmvm into a mvm*/
      for(j = 0; j < 16; j++) mvm[j] = keepmvm[j];

     /*screen set size*/
      gtk_widget_set_size_request(GTK_WIDGET(drawingarea), screen_width, screen_height);


      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();

      glOrtho(-x_size - camera_x, x_size - camera_x,
              -y_size - camera_y, y_size - camera_y,
              -z_depth, z_depth);  
      glMatrixMode(GL_MODELVIEW);


      glMatrixMode(GL_MODELVIEW);
      glLoadMatrixd (keepmvm);
  
      if (m->ngrids) luscus_gtk_resort_surfaces();
      redraw();
    }
  }
  else
  {
    status=fopen(".ViewPoint","w");
  
    fprintf(status,"screen sizes\n %lf %lf\n", screen_width, screen_height);
    fprintf(status,"field sizes\n %lf %lf\n", x_size, y_size);
    fprintf(status,"scale\n %lf\n", scale);
    fprintf(status,"z_depth\n %lf\n", 2.0 * molecule_diameter);
    fprintf(status,"camera position\n %lf %lf\n", camera_x, camera_y);
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX,keepmvm);
    fprintf(status,"matrix \n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n",
            keepmvm[0], keepmvm[1], keepmvm[2], keepmvm[3],
            keepmvm[4], keepmvm[5], keepmvm[6], keepmvm[7],
            keepmvm[8], keepmvm[9], keepmvm[10],keepmvm[11],
            keepmvm[12],keepmvm[13],keepmvm[14],keepmvm[15]);
    fclose(status);
    printf("current ViewPoint saved in .ViewPoint file\n");
  }
  return;
}

/*select*/

void unselect_all(void)
{
  int i;
  for(i = 0; i < 4; i++) m->selected[i] = -1;
  m->n_selected = 0;
  rerender_3d();
  luscus_gtk_update_upon_select_or_mark();
}

/*mark*/

void unmark_all(void)
{
  m->n_marked = 0;
  if (m->marked) free(m->marked);
  m->marked = NULL;
  luscus_gtk_update_upon_select_or_mark();
  rerender_3d();
}

int is_atom_marked(int iatom)
{
  int i;
  if (m->n_marked <= 0) return 0;
  for(i = 0; i < m->n_marked; i++) if (m->marked[i] == iatom) return 1;
  return 0;
}

void mark_atom(int float_select)
{
  m->n_marked++;
  m->marked = (int*) realloc(m->marked, sizeof(int) * m->n_marked);
  m->marked[m->n_marked-1] = float_select;
  luscus_gtk_update_upon_select_or_mark();
  rerender_3d();
}

void mark_H(void)
{
  int i;
  for(i = 0; i < m->natom; i++)
    if (m->elem[i].name[0] == 'H' && m->elem[i].name[1] == 0)
      if (!is_atom_marked(i))
       	mark_atom(i);
  rerender_3d();
}

void mark_element_as_selected(void)
{
  int i;
  if (m->n_selected <= 0) return;
  for(i = 0; i < m->natom; i++)
    if (strcmp(m->elem[i].name, m->elem[m->selected[0]].name) == 0)
      if (!is_atom_marked(i))
       	mark_atom(i);
  rerender_3d();
}

void mark_neighbor(void)
{
  int i;
  if (m->n_selected <= 0) return;
  if (m->n_selected == 1)
  {
    for(i = 0; i < m->nbond; i++)
      if (m->bond[i].iat1 == m->selected[0])
      {
        if (!is_atom_marked(m->bond[i].iat2))
          mark_atom(m->bond[i].iat2);
      }
      else if (m->bond[i].iat2 == m->selected[0])
      {
        if (!is_atom_marked(m->bond[i].iat1))
          mark_atom(m->bond[i].iat1);
      }
    rerender_3d();
  }
  else if (m->n_selected == 2)
  {
    mark_atom(m->selected[0]);
    for(i = 0; i < m->nbond; i++)
      if (m->bond[i].iat1 == m->selected[0] && m->bond[i].iat2 != m->selected[1])
      {
        mark_atom(m->bond[i].iat2);
        mark_nei(m->bond[i].iat2);
      }
      else if (m->bond[i].iat2 == m->selected[0] && m->bond[i].iat1 != m->selected[1])
      {
        mark_atom(m->bond[i].iat1);
        mark_nei(m->bond[i].iat1);
      }
    rerender_3d();
  }
}

void mark_nei(int iatom)
{
  int i;
  if (iatom > m->natom) return;
  for(i = 0; i < m->nbond; i++)
  {
    if(m->bond[i].iat1 == iatom && !is_atom_marked(m->bond[i].iat2))
    {
      mark_atom(m->bond[i].iat2);
      mark_nei(m->bond[i].iat2);
    }
    else if(m->bond[i].iat2 == iatom && !is_atom_marked(m->bond[i].iat1))
    {
      mark_atom(m->bond[i].iat1);
      mark_nei(m->bond[i].iat1);
    }
    rerender_3d();
  }
}

void mark_one_side(void)
{
  int i, j;
  double center[3];
  double vector[3];
  double tvec1[3], tvec2[3];
  double tvec[3];
  double cphi;
  double vlen;
  int marked;

  if (m->n_selected == 2)
    for(i = 0, vlen = 0.0; i < 3; i++)
    {
      center[i] = 0.5 * (m->xyz[m->selected[0]][i] + m->xyz[m->selected[1]][i]);
      vector[i] = m->xyz[m->selected[1]][i] - m->xyz[m->selected[0]][i];
      vlen += vector[i] * vector[i];
    }
  else if (m->n_selected == 3)
  {
    for(i = 0, vlen = 0.0; i < 3; i++)
    {
      center[i] = (m->xyz[m->selected[0]][i] + 
                   m->xyz[m->selected[1]][i] + 
                   m->xyz[m->selected[2]][i])/3.0;
      tvec1[i] = m->xyz[m->selected[1]][i] - m->xyz[m->selected[0]][i];
      tvec2[i] = m->xyz[m->selected[2]][i] - m->xyz[m->selected[1]][i];
    }
    /*vector product*/
    vector[0] = tvec1[1] * tvec2[2] - tvec1[2] * tvec2[1];
    vector[1] = tvec1[2] * tvec2[0] - tvec1[0] * tvec2[2];
    vector[2] = tvec1[0] * tvec2[1] - tvec1[1] * tvec2[0];
    vlen += vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
  }

  vlen = 1.0/sqrt(vlen);
  for(i = 0; i < 3; i++) vector[i] *= vlen;

  for(i = 0; i < m->natom; i++)
  {
    cphi = 0.0;
    for(j = 0, vlen = 0.0; j < 3; j++)
    {
      tvec[j] = m->xyz[i][j] - center[j];
      vlen += tvec[j]*tvec[j];
    }
    vlen = 1.0/sqrt(vlen);
    for(j = 0; j < 3; j++) cphi += vector[j] * tvec[j] * vlen;
    if(cphi > 0.0) mark_atom(i);
  }
  if (m->n_selected == 3)
    marked = 0;
    for(i = 0; i < 3; i++)
    {
      for(j = 0; j < m->n_marked; j++)
        if (m->marked[j] == m->selected[i]) marked = 1;
      if (!marked)
        mark_atom(m->selected[i]);
    }
  rerender_3d();
}

void unmark_atom(int float_select)
{
  int i;
  for (i = float_select; i < m->n_marked-1; i++)
    m->marked[i] = m->marked[i + 1];
  m->n_marked--;
  m->marked = (int*) realloc(m->marked, sizeof(int) * m->n_marked);
  luscus_gtk_update_upon_select_or_mark();
  rerender_3d();
}

void renumber_marked(void)
{
  int i, j;
  int im, inext;
  XYZ *tmpxyz;
  ELEM_DATA *tmpe;
  double *tmp_charg_m, *tmp_charg_l;
  int *tmp_add_numer, *tmp_symm;
  char **tmp_name;
  int *index;

  tmpxyz = (XYZ*) malloc(m->natom * sizeof(XYZ));
  tmpe = (ELEM_DATA*) malloc(m->natom * sizeof(ELEM_DATA));
  tmp_charg_m = (double*) malloc(m->natom * sizeof(double));
  tmp_charg_l = (double*) malloc(m->natom * sizeof(double));
  tmp_add_numer = (int*) malloc(m->natom * sizeof(int));
  tmp_symm = (int*) malloc(m->natom * sizeof(int));
  tmp_name = (char**) malloc(m->natom * sizeof(char*));
  index = (int*) malloc(m->natom * sizeof(int));

  /*copy marked atoms into a temprary arrays*/
  for(i = 0; i < m->n_marked; i++)
  {
    index[i] = m->marked[i];
    for(j = 0; j < 3; j++) tmpxyz[i][j] = m->xyz[m->marked[i]][j];
    tmpe[i].name = m->elem[m->marked[i]].name;
    tmpe[i].vdw_rad = m->elem[m->marked[i]].vdw_rad;
    tmpe[i].bond_rad = m->elem[m->marked[i]].bond_rad;
    tmpe[i].valency = m->elem[m->marked[i]].valency;
    for(j = 0; j < 4; j++) tmpe[i].color[j] = m->elem[m->marked[i]].color[j];
    tmp_charg_m[i] = m->charge_m[m->marked[i]];
    tmp_charg_l[i] = m->charge_l[m->marked[i]];
    tmp_add_numer[i] = m->additional_numeration[m->marked[i]];
    tmp_symm[i] = m->symmetry[m->marked[i]];
    tmp_name[i] = m->name[m->marked[i]];
  }
  /*copy unmarked atoms into a temprary arrays*/
  im = m->n_marked;
  for(i = 0; i < m->natom; i++)
  {
    inext = 1;
    for(j = 0; j < m->n_marked; j++)
      if (m->marked[j] == i) inext = 0;
    if (inext)
    {
      index[im] = i;
      for(j = 0; j < 3; j++) tmpxyz[im][j] = m->xyz[i][j];
      tmpe[im].name = m->elem[i].name;
      tmpe[im].vdw_rad = m->elem[i].vdw_rad;
      tmpe[im].bond_rad = m->elem[i].bond_rad;
      tmpe[im].valency = m->elem[i].valency;
      for(j = 0; j < 4; j++) tmpe[im].color[j] = m->elem[i].color[j];
      tmp_charg_m[im] = m->charge_m[i];
      tmp_charg_l[im] = m->charge_l[i];
      tmp_add_numer[im] = m->additional_numeration[i];
      tmp_symm[im] = m->symmetry[i];
      tmp_name[im] = m->name[i];
      im++;
    }
  }

  /*update bonds*/
  for(i = 0; i < m->nbond; i++)
  {
    for(j = 0; j < m->natom; j++)
      if(index[j] == m->bond[i].iat1)
      {
        m->bond[i].iat1 = j;
        break;
      }
    for(j = 0; j < m->natom; j++)
      if(index[j] == m->bond[i].iat2)
      {
        m->bond[i].iat2 = j; 
        break;
      }
  }

  /*copy data back to the main arrays*/
  for(i = 0; i < m->natom; i++)
  {
    for(j = 0; j < 3; j++) m->xyz[i][j] = tmpxyz[i][j];
    m->elem[i].name = tmpe[i].name;
    m->elem[i].vdw_rad = tmpe[i].vdw_rad;
    m->elem[i].bond_rad = tmpe[i].bond_rad;
    m->elem[i].valency = tmpe[i].valency;
    for(j = 0; j < 4; j++) m->elem[i].color[j] = tmpe[i].color[j];
    m->charge_m[i] = tmp_charg_m[i];
    m->charge_l[i] = tmp_charg_l[i];
    m->additional_numeration[i] = tmp_add_numer[i];
    m->symmetry[i] = tmp_symm[i];
    m->name[i] = tmp_name[i];
  }

/*    swap_atoms(m->marked[i], j++);*/
  unselect_all();
  unmark_all();
  append_backup();
  deallocate_atom_list();
  insert_all_atoms_into_list();

  free(tmpxyz);
  free(tmpe);
  free(tmp_charg_m);
  free(tmp_charg_l);
  free(tmp_add_numer);
  free(tmp_symm);
  free(tmp_name);
  free(index);
}

void reverse_marked(void)
{
  int i, j;
  int new_n_marked;
  int *new_marked;

  new_n_marked = m->natom - m->n_marked;
  new_marked = (int*) malloc(sizeof(int) * new_n_marked);

  for (i = 0, j = 0; i < m->natom; i++)
    if (!is_atom_marked(i))
      new_marked[j++] = i;

  unmark_all();
  m->marked = new_marked;
  m->n_marked = new_n_marked;
  luscus_gtk_update_upon_select_or_mark();
  rerender_3d();
}

/*----- EPS -----*/

void initPS(int opt, char *fname)
{
  double U[3][3];
  double V[16];
#ifdef TRACE
  puts("Trace: initPS");
#endif

  Input_Data.ps=1;

/* nextfname(Input_Data.basename,"eps",fname);*/

  fl_ps = fopen(fname, "w");
  if(fl_ps==NULL)
  {
    make_warning("Error: Can't create file");
/*      sprintf(LOG,"Can't create file %s \n",fname); printlog(LOG);*/
    return;
  }

  gveps_init(fl_ps,0,1.0,0.0,0.0,&ctl);
  glGetDoublev(GL_MODELVIEW_MATRIX,V);

  U[0][0]=V[0];
  U[1][0]=V[1];
  U[2][0]=V[2];

  U[0][1]=V[4];
  U[1][1]=V[5];
  U[2][1]=V[6];

  U[0][2]=V[8];
  U[1][2]=V[9];
  U[2][2]=V[10];

  gveps_view(fl_ps,0,U,&ctl);

  draw_bonds();
  draw_atoms();
  draw_surfaces();

  gveps_paint(fl_ps,opt,&ctl);

  fclose(fl_ps);
  Input_Data.ps=0;
}

void print_bond_eps(int ia1, int ia2, double mid[3])
{
  double povcsize=0.08;

  /* povray starts here */
  if(Input_Data.povray)
  {
    fprintf(fl_povray, "\ncone {\n <%10.7f,%10.7f,%10.7f>, %10.7f\n  <%10.7f,%10.7f,%10.7f>, %10.7f\n pigment { color rgb <%7.4f %7.4f %7.4f>\n  }\n}\n",
    m->xyz[ia1][0], m->xyz[ia1][1],  m->xyz[ia1][2],
    povcsize,
    mid[0], mid[1], mid[2],
    povcsize,
    m->elem[ia1].color[0],
    m->elem[ia1].color[1],
    m->elem[ia1].color[2]);
    fprintf(fl_povray, "\ncone {\n <%10.7f,%10.7f,%10.7f>, %10.7f\n  <%10.7f,%10.7f,%10.7f>, %10.7f\n pigment { color rgb <%7.4f %7.4f %7.4f>\n  }\n}\n",
    m->xyz[ia2][0], m->xyz[ia2][1], m->xyz[ia2][2],
    povcsize,
    mid[0], mid[1], mid[2],
    povcsize,
    m->elem[ia2].color[0],
    m->elem[ia2].color[1],
    m->elem[ia2].color[2]);
  }
  /* povray ends here */

  /* gveps starts here */
  if(Input_Data.ps)
  {
    gveps_bond(fl_ps, 0,
              m->xyz[ia1][0], m->xyz[ia1][1], m->xyz[ia1][2],
              m->elem[ia1].color[0],
              m->elem[ia1].color[1],
              m->elem[ia1].color[2],
              m->xyz[ia2][0], m->xyz[ia2][1], m->xyz[ia2][2],
              m->elem[ia2].color[0],
              m->elem[ia2].color[1],
              m->elem[ia2].color[2],
              &ctl);
  }
}

void initPovray(char *fname)
{
  FILE *fl_header;
  char str[256];
  float camera_org[3];
#ifdef TRACE
  puts("Trace: initPovray");
#endif
  Input_Data.povray=1;

/* nextfname(Input_Data.basename,"pov",fname);*/ /*this is done in kbd_press function -goran*/

  fl_povray = fopen(fname, "w");
  if(fl_povray==NULL)
  {
/*     sprintf(LOG,"Can't create file %s\n",fname); printlog(LOG);*/
    make_warning("Can't create file ");
    return;
  }
 /* print header */
  fl_header=fopen("header.pov","r");
  if(fl_header != NULL)
  {
    while( fgets(str, 256, fl_header))
    {
      fputs(str,fl_povray);
    }
    fclose(fl_header);
  }
  else
  {
    camera_org[0] = -(mvm[0] * mvm[12] + mvm[1] * mvm[13] + mvm[2] * mvm[14]);
    camera_org[1] = -(mvm[4] * mvm[12] + mvm[5] * mvm[13] + mvm[6] * mvm[14]);
    camera_org[2] = -(mvm[8] * mvm[12] + mvm[9] * mvm[13] + mvm[10] * mvm[14]);

/*    fprintf(fl_povray, "#include \"colors.inc\"\n#include \"textures.inc\"\n#include \"skies.inc\"\n#include \"stones.inc\"\n#include \"woods.inc\"\n#include \"metals.inc\"\n#include \"ash.map\"\n#include \"bubinga.map\"\n#include \"cedar.map\"\n#include \"orngwood.map\"\n#include \"teak.map\"\n#include \"whiteash.map\"\n#include \"benediti.map\"\n#include \"marbteal.map\"\n#include \"pinkmarb.map\"\n#include \"rdgranit.map\"\n\n");*/
    fprintf(fl_povray, "global_settings {\n  assumed_gamma 2.2\n  ambient_light color rgb <1, 1, 1>\n}\n\n");
    fprintf(fl_povray, "camera {\n \torthographic\n \tlocation <%f, %f, %f>\n \tlook_at <%f, %f, %f>\n \tangle 35\n}\n\n",
    
    2*molecule_diameter * mvm[2], 2*molecule_diameter*mvm[6], 2*molecule_diameter*mvm[10], camera_org[0], camera_org[1], camera_org[2]);
/*    camera_org[0], camera_org[1], camera_org[2], mvm[2], mvm[6], mvm[10]);*/
    fprintf(fl_povray, "sky_sphere {pigment {rgb <0.859,0.910,0.831>}}\n\n");
    fprintf(fl_povray, "light_source\n{  // #1\n  <-0.5,-0.5,1.0>\n  color rgb <1,1,1>\n  parallel\n  point_at <0,0,0>\n  fade_distance 1\n  fade_power 0\n}\n");
  }

  /* update picture*/

  if (m->natom && !Input_Data.hide_atoms) draw_atoms();
  if (m->nbond && !Input_Data.hide_bonds) draw_bonds();

  if (Input_Data.show_axis) draw_axes();
  if (m->ishow & HAS_DIPOLE) draw_dipole();
  if (m->nvector) draw_vectors();
  if (m->ntriangle) draw_triangles();
  if (m->nsphere) draw_spheres();
  if (m->nsurf) draw_surfaces();
  if (m->ncells) draw_cells();
  if (m->ngrids) draw_grid();

  fclose(fl_povray);
  Input_Data.povray=0;
}


/*----- surface -----*/

void luscus_gtk_resort_surfaces(void)
{
  if (All_Surfaces) All_Surfaces->is_sorted = 0;
  rerender_3d();
}

void deallocate_msrf(void)
{
  if (All_Surfaces) msrf_Clean(All_Surfaces);
  All_Surfaces = NULL;
}

void allocate_msrf(int nngrids)
{
  int i;
  surf = (surface_t**) malloc(sizeof(surface_t*) * nngrids);
  for(i = 0; i < nngrids; i++) surf[i] = NULL;
  if (All_Surfaces) msrf_Clean(All_Surfaces);

  All_Surfaces = msrf_Init(All_Surfaces, nngrids);
}

void make_surfaces(void)
{
  surface_t* srf;
  double dummy,min,max;
  double *data_ptr;
  double rad[3];

  do_remove_grid();

  if (iorb < 0 || iorb > m->ngrids) return;

  current_grid_limits(&min, &max, &dummy);
  if (min < 0.0 && max < 0.0)
  {
    Input_Data.lev_min = -max;
    Input_Data.lev_max = -min;
  }
  else if (min < 0.0 && max > 0.0)
  {
    Input_Data.lev_min = 0.0;
    Input_Data.lev_max = max;
  }
  else
  {
    Input_Data.lev_min = min;
    Input_Data.lev_max = max;
  }

 /* calculate isolevel step, if necessary */
  if (Input_Data.n_lev_steps)
    Input_Data.lev_step =(Input_Data.lev_max - Input_Data.lev_min) / Input_Data.n_lev_steps;

  if (Input_Data.auto_lev)
    Input_Data.lev = 0.75 * Input_Data.lev_min + 0.25 * Input_Data.lev_max;

  current_grid_data(&data_ptr);
  rad[0]= m->grid.axisvec1[0];
  rad[1]= m->grid.axisvec2[1];
  rad[2]= m->grid.axisvec3[2];

  srf= surf[iorb] = msrf_New_Surface(All_Surfaces);

  mcubes(srf,
         m->grid.npt[0], m->grid.npt[1], m->grid.npt[2],
         data_ptr, m->epot, m->grid.origin, Input_Data.lev, rad, 1);
  Input_Data.max_diam_ok=0; /* diameter is changed */
}

void do_remove_grid(void)
{
  int i;

  for(i = 0; i < m->ngrids; i++)
    if (surf[i])
    {
      msrf_Delete_Surface(surf[i]);
      surf[i] = NULL;
    }
#ifdef _DEBUG_
  Debug("Surface %p deleted",srf);
#endif

  Input_Data.max_diam_ok=0; /* diameter is changed */
}

/*statusbars*/

void make_warning(gchar *warning_text)
{
  GtkWidget *dialog;

  dialog = gtk_message_dialog_new(GTK_WINDOW(window), GTK_DIALOG_DESTROY_WITH_PARENT,
                                  GTK_MESSAGE_ERROR, GTK_BUTTONS_CLOSE, "%s", warning_text);
  gtk_dialog_run(GTK_DIALOG (dialog));
  gtk_widget_destroy (dialog);

  return;
}

void print_reference(void)
{
  luscus_gtk_push_message_to_statusbar2("luscus reference: G. Kovačević, V. Veryazov, J. Cheminformatics, 7 (2015) 1-10");
/*  print_plugin_warnings(); */
}

void luscus_gtk_push_message_to_statusbar1(char* message)
{
  guint context_id;
  guint message_id;

  if (GTK_IS_STATUSBAR(statusbar1))
  {
    context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusbar1), "some message");
    message_id = gtk_statusbar_push(GTK_STATUSBAR(statusbar1), context_id, message);
  }
}

void luscus_gtk_push_message_to_statusbar2(gchar* message)
{
  guint context_id;
  guint message_id;

  if (GTK_IS_STATUSBAR(statusbar2))
  {
    context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusbar2), "some message");
    message_id = gtk_statusbar_push(GTK_STATUSBAR(statusbar2), context_id, message);
  }
}

void luscus_gtk_pop_message_from_statusbar1(void)
{
  guint context_id;
  context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusbar1), "some message");
  gtk_statusbar_pop(GTK_STATUSBAR(statusbar1), context_id);
}

void luscus_gtk_pop_message_from_statusbar2(void)
{
  guint context_id;
  context_id = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusbar2), "some message");
  gtk_statusbar_pop(GTK_STATUSBAR(statusbar2), context_id);
}

/*----- orientation/size -----*/

void Do_center(void)
{
  XYZ origin;

  if (move_camera)
  {
    get_center(&origin[0], &origin[1], &origin[2]);
    camera_x = -origin[0];
    camera_y = -origin[1];
    set_scale();
  }
  else
    set_origin_molecule();
}

void set_sizes(void)
{
  molecule_diameter = Calc_Diameter();
}

void set_scale(void)
{
/*  set_sizes();*/
  if (screen_width < screen_height)
  {
    x_size = molecule_diameter / (scale * 2.F);
    y_size = molecule_diameter * screen_height / (screen_width * scale * 2.F);
  }
  else
  {
    x_size = molecule_diameter * screen_width / (screen_height * scale * 2.F);
    y_size = molecule_diameter / (scale * 2.F);
  }
/*  if (m)
    if (m->ngrids)
    {
      x_size *= 2.0;
      y_size *= 2.0;
    }*/
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glOrtho(-x_size - camera_x, x_size - camera_x,
          -y_size - camera_y, y_size - camera_y,
          -2.0*molecule_diameter, 2.0*molecule_diameter); 
  glMatrixMode(GL_MODELVIEW);
  glViewport(0, 0, (int) screen_width, (int) screen_height);
}

void Do_key_ud(int iopt, int is_shift) /* 0 for down 1 for up */
{
  double sh=-5.0;
  GLint isd;
  XYZ p1, p2;
  XYZ trans;
  GLdouble modelMatrix[16], projMatrix[16];
  GLint viewport[4];

  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);

  if(iopt==1) sh=5.0;
  if(is_shift)
  {
    if (move_camera)
    {
      camera_y += sh * 0.2;
      set_scale();
    }
    else
    {
/*      move_molecule[1] += sh*5.0;*/
      isd = gluUnProject((viewport[2]-viewport[0])/2, (viewport[3]-viewport[1])/2 + sh*5.0, 0.0, modelMatrix, projMatrix, viewport, &p1[0], &p1[1], &p1[2]);
      isd = gluUnProject((viewport[2]-viewport[0])/2, (viewport[3]-viewport[1])/2, 0.0, modelMatrix, projMatrix, viewport, &p2[0], &p2[1], &p2[2]);
      trans[0] = p1[0] - p2[0];
      trans[1] = p1[1] - p2[1];
      trans[2] = p1[2] - p2[2];

      if (isd == GL_TRUE) translate(trans);
      else fprintf(stderr,"ERROR in GluUnProject\n");
/*      glMatrixMode(GL_PROJECTION);
      glTranslated(0.0,sh*5.0,0.0);
      glMatrixMode(GL_MODELVIEW);*/
    }
  }
  else
  {
    if (move_camera)
    {
      camera_y += sh * 0.02;
      set_scale();
    }
    else
    {
/*      move_molecule[1] += sh;*/
      isd = gluUnProject((viewport[2]-viewport[0])/2, (viewport[3]-viewport[1])/2 + sh, 0.0, modelMatrix, projMatrix, viewport, &p1[0], &p1[1], &p1[2]);
      isd = gluUnProject((viewport[2]-viewport[0])/2, (viewport[3]-viewport[1])/2, 0.0, modelMatrix, projMatrix, viewport, &p2[0], &p2[1], &p2[2]);
      trans[0] = p1[0] - p2[0];
      trans[1] = p1[1] - p2[1];
      trans[2] = p1[2] - p2[2];

      if (isd == GL_TRUE) translate(trans);
      else fprintf(stderr,"ERROR in GluUnProject\n");
/*      translate(trans);*/
/*      glMatrixMode(GL_PROJECTION);
      glTranslated(0.0,sh,0.0);
      glMatrixMode(GL_MODELVIEW);*/
    }
  }
}

void Do_key_lr(int iopt, int is_shift)  /* iopt=0 for left, 1 for right  */
{
  double sh=-5.0;
  GLint isd;
  XYZ trans;
  XYZ p1, p2;
  GLdouble modelMatrix[16], projMatrix[16];
  GLint viewport[4];

  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);

  if(iopt==1) sh=5.0;
  if(is_shift)
  {
    if (move_camera)
    {
      camera_x += sh * 0.2;
      set_scale();
    }
    else
    {
/*      move_molecule[0] += sh*5.0;*/
      isd = gluUnProject(0.0, 0.0, 0.0, modelMatrix, projMatrix, viewport, &p1[0], &p1[1], &p1[2]);
      isd = gluUnProject(sh*5.0, 0.0, 0.0, modelMatrix, projMatrix, viewport, &p2[0], &p2[1], &p2[2]);
      trans[0] = p2[0] - p1[0];
      trans[1] = p2[1] - p1[1];
      trans[2] = p2[2] - p1[2];

      if (isd == GL_TRUE) translate(trans);
      else fprintf(stderr,"ERROR in GluUnProject\n");
/*      glMatrixMode(GL_PROJECTION);
      glTranslated(sh*5.0,0.0,0.0);
      glMatrixMode(GL_MODELVIEW);*/
    }
  }
  else
  {
    if (move_camera)
    {
      camera_x += sh * 0.02;
      set_scale();
    }
    else
    {
/*      move_molecule[0] += sh;*/
    isd = gluUnProject(0.0, 0.0, 0.0, modelMatrix, projMatrix, viewport, &p1[0], &p1[1], &p1[2]);
    isd = gluUnProject(sh, 0.0, 0.0, modelMatrix, projMatrix, viewport, &p2[0], &p2[1], &p2[2]);
    trans[0] = p2[0] - p1[0];
    trans[1] = p2[1] - p1[1];
    trans[2] = p2[2] - p1[2];

    if (isd == GL_TRUE) translate(trans);
    else fprintf(stderr,"ERROR in GluUnProject\n");
/*      glMatrixMode(GL_PROJECTION);
      glTranslated(sh,0.0,0.0);
      glMatrixMode(GL_MODELVIEW);*/
    }
  }
}

void do_key_page(int iopt, int ishift)
{
  if (!m) return;
  if (m->ngrids)
  {
    if (iopt) luscus_gtk_select_orbital_down();
    else luscus_gtk_select_orbital_up();
    rerender_3d();
  }
  else if (ishift) glTranslated(0.0,0.0,iopt * 0.1);
  else if (m->n_selected == 1)
  {
    if (iopt) change_element_by_one(1);
    else change_element_by_one(-1);
    rerender_3d();
    redraw();
  }
  else if (m->n_selected == 2)
  {
    if (iopt)
      change_bond_type((get_bond_type(m->selected[0], m->selected[1])+1)%7);
    else
      change_bond_type((get_bond_type(m->selected[0], m->selected[1])-1)%7);
      rerender_3d();
  }
  /*else if (m->n_selected == 3)
  {
    change angle value();
  }*/
  else if (m->nvibr)
  {
    if (iopt)
    {
      ivib++;
      if (ivib >= m->nvibr) ivib = m->nvibr-1;
    }
    else
    {
      ivib--;
      if (ivib < 0) ivib = 0;
    }
    set_current_graph_data();
    rerender_3d();
    redraw();
  }
  else if (n_geometries)
  {
    if (iopt) callback_geo_forw(NULL, NULL);
    else callback_geo_back(NULL, NULL);
  } 
}

void do_key_insert(void)
{
  if (m->n_selected <= 1)
  {
    if (!last_fragment) add_atoms(0);
    else add_fragment(0);
  }
  draw_all_pixdata();
}

void luscus_gtk_move_light_state(void)
{
  if(Input_Data.move_light==0) glGetDoublev(GL_MODELVIEW_MATRIX,keepmvm);
  else glLoadMatrixd(keepmvm);

  Input_Data.move_light=!Input_Data.move_light;
}

/*----------------------*/

void redraw(void)
{
  GdkWindow *window;
  GtkAllocation allocation;
  window = gtk_widget_get_window(drawingarea);
  gtk_widget_get_allocation(drawingarea, &allocation);
  gdk_window_invalidate_rect(window, &allocation, FALSE);
}

double grayscale(double r, double g, double b)
{
  return r*0.299+g*0.587+b*0.114;
}

/*---- timers -----*/

void luscus_gtk_freq_timer(void)
{
  g_timeout_add(30, luscus_freq_timer, NULL);
}

gboolean luscus_freq_timer(gpointer data)
{
  /*select another list member*/
  rerender_3d();
  redraw();
  if (!m->nvibr) return FALSE;
  else return TRUE;
}

void start_check_for_changes(void)
{
  g_timeout_add_seconds(10, check_for_changes, NULL);
}

gboolean check_for_changes(gpointer data)
{
  check_if_file_changed();
  return TRUE;
}

/*render functions*/

void draw_molecule_for_multiview(void)
{
  int i, j;
  int ia1, ia2;
  double middle[3];
  double angle, norm0, norm1, lenbond;
  GLUquadricObj *quadObj;

  quadObj=gluNewQuadric();
  if (quadObj == 0) return;

  /*draw atoms - low resolution*/

  for(i = 0; i < m->natom; i++)
  {
    glPushMatrix();
    glTranslated(m->xyz[i][0], m->xyz[i][1], m->xyz[i][2]);
    SetAtomColor(i);
    luscus_gtk_draw_sphere(1, m->elem[i].vdw_rad/2., 7, 7);

    glPopMatrix();
  }

  /*draw bonds - low resolution*/

  for(i = 0; i < m->nbond; i++)
  {
    if (m->bond[i].bond_type == SINGLE_BOND || m->bond[i].bond_type == DOUBLE_BOND || m->bond[i].bond_type == TRIPLE_BOND || m->bond[i].bond_type == S_P_BOND)
    {
      ia1 = m->bond[i].iat1;
      ia2 = m->bond[i].iat2;
      for(j = 0; j < 3; j++) middle[j] = 0.5F * (m->xyz[ia1][j] + m->xyz[ia2][j]);

      CalcCyllinder(m->xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);

      glPushMatrix();
        glTranslated(m->xyz[ia1][0], m->xyz[ia1][1], m->xyz[ia1][2]);
        glRotated(angle, norm0, norm1, 0.0);
        SetAtomColor(ia1);
        gluCylinder(quadObj, 0.1, 0.1, lenbond, 7, 7);
      glPopMatrix();

      CalcCyllinder(middle, m->xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      glPushMatrix();
        glTranslated(middle[0], middle[1], middle[2]);
        glRotated(angle, norm0, norm1, 0.0);
        SetAtomColor(ia2);
        gluCylinder(quadObj, 0.1, 0.1, lenbond, 7, 7);
      glPopMatrix();
    }
  }
  gluDeleteQuadric(quadObj);
}


void SetAtomColor(int iatom)
{
  int i;
  double gray;
  if (iatom >= m->natom) return;

  for (i = 0; i < m->n_marked; i++)
  {
    if (m->marked[i] == iatom)
    {
      glColor3d(0.0f, 1.0f, 1.0f);
    }
  }

  if (Input_Data.bw)
  {
    gray = grayscale(m->elem[iatom].color[0], m->elem[iatom].color[1], m->elem[iatom].color[2]);
    glColor3d(gray, gray, gray);
  }
  else
    glColor3fv(m->elem[iatom].color);
}

void draw_vibrations(void)
{
  int i, j;
  XYZ *xyz;

  xyz = (XYZ*) malloc(sizeof(XYZ) * m->natom);

  if (is_is) rvibr += 0.1 * Input_Data.frequency_speed;
  else rvibr -= 0.1 * Input_Data.frequency_speed;
  if (rvibr > 1.0 || rvibr < -1.0) is_is = 1 - is_is;

  for(i = 0; i < m->natom; i++)
    for(j = 0; j < 3; j++)
      xyz[i][j] = m->xyz[i][j] + m->normal_mode[ivib][i][j] * rvibr * Input_Data.frequency_amplitude;

  draw_atoms_xyz(xyz);
  draw_bonds_xyz(xyz);

  vibr_count++;

  free(xyz);
}

void draw_atoms_xyz(XYZ *xyz)
{
  int i, j;
  double tmprad;
  for(i = 0; i < m->natom; i++)
  {
    if (Input_Data.style_atoms) tmprad = m->elem[i].vdw_rad/30.;
    else tmprad = m->elem[i].vdw_rad/3.;
    glPushMatrix();
    glTranslated(xyz[i][0], xyz[i][1], xyz[i][2]);
    SetAtomColor(i);
    for(j = 0; j < m->n_marked; j++)
      if (m->marked[j] == i) glColor3d(0.0f, 1.0f, 1.0f);

    switch(Input_Data.atomshape)
    {
      case 0:
        luscus_gtk_draw_sphere(1, tmprad * Input_Data.fatness_a, slices, stacks);
        for (j = 0; j < m->n_selected; j++)
        {
          if (m->selected[j] == i)
          {
            if (j == 0) glColor3d(0.0f, 0.0f, 1.0f);
            else glColor3d(1.0f, 0.0f, 1.0f);
            glDisable(GL_LIGHTING);
            luscus_gtk_draw_sphere(0, tmprad * Input_Data.fatness_a + 0.02, 20, 20);
            glEnable(GL_LIGHTING);
          }
        }
        break;
#ifdef GTK_GLEXT
      case 1:
        luscus_gtk_draw_torus(1, tmprad/2., tmprad, 30, 30);
        break;
      case 2:
        glFrontFace(GL_CW);
        luscus_gtk_draw_teapot(1, tmprad);
        glFrontFace(GL_CCW);
        break;
      case 3:
        luscus_gtk_draw_cube(1, tmprad);
        break;
#endif
      default:
        luscus_gtk_draw_sphere(1, tmprad * Input_Data.fatness_a, slices, stacks);
    }
    glPopMatrix();

  }
}

void draw_bonds_xyz(XYZ *xyz)
{
  int i, j;
  int i1;
  double csize=0.05*Input_Data.fatness_b;
  double middle[3];
  double angle, norm0, norm1, lenbond;
  double delta[3];
  double delta2[3];
  double plain[3];
  int ia1, ia2, itmpb;
  GLUquadricObj *quadObj;

  quadObj=gluNewQuadric();

  plain[0]=10.0;
  plain[1]=20.0;
  plain[2]=30.0;

  if (quadObj == 0) return;

  for(i = 0; i < m->nbond; i++)
  {
    ia1 = m->bond[i].iat1;
    ia2 = m->bond[i].iat2;
    itmpb = m->bond[i].bond_type;
    for(j = 0; j < 3; j++) middle[j] = 0.5F * (xyz[ia1][j] + xyz[ia2][j]);
    if (Input_Data.style_atoms && itmpb != NO_BOND) itmpb = LINE_BOND;

    if (itmpb == NO_BOND) continue;

    else if (itmpb == LINE_BOND)
    {
      glDisable(GL_LIGHTING);
      glBegin(GL_LINE_STRIP);
      SetAtomColor(ia1);
      glVertex3dv(xyz[ia1]);
      glVertex3dv(middle);
      SetAtomColor(ia2);
      glVertex3dv(middle);
      glVertex3dv(xyz[ia2]);
      glEnd();
      glEnable(GL_LIGHTING);
    }
    else if (itmpb == SINGLE_BOND)
    {
      CalcCyllinder(xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);
      
      glPushMatrix();
        glTranslated(xyz[ia1][0], xyz[ia1][1], xyz[ia1][2]);
        glRotated(angle, norm0, norm1, 0.0);
        SetAtomColor(ia1);
        gluCylinder(quadObj, csize, csize ,lenbond, slices, stacks);
      glPopMatrix();

      CalcCyllinder(middle, xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      glPushMatrix();
        glTranslated(middle[0], middle[1], middle[2]);
        glRotated(angle, norm0, norm1, 0.0);
        SetAtomColor(ia2);
        gluCylinder(quadObj, csize, csize, lenbond, slices, stacks);
      glPopMatrix();
    }
    else if (itmpb == PARTIAL_BOND)
    {
      CalcCyllinder(xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);

      for(j = 0; j < 3; j++) delta[j] = (middle[j] - xyz[ia1][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
	glPushMatrix();
          SetAtomColor(ia1);
          glTranslated(xyz[ia1][0]+delta[0]*j, xyz[ia1][1]+delta[1]*j, xyz[ia1][2]+delta[2]*j);
          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize, lenbond/20, slices, stacks);
        glPopMatrix();
      }

      CalcCyllinder(middle, xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      for(j = 0; j < 3; j++) delta[j] = (-middle[j] + xyz[ia2][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
        glPushMatrix();
	  glTranslated(middle[0]+delta[0]*j, middle[1]+delta[1]*j, middle[2]+delta[2]*j);
          glRotated(angle,norm0, norm1, 0.0);
          SetAtomColor(ia2);
          gluCylinder(quadObj, csize, csize, lenbond/20, slices, stacks);
        glPopMatrix();
      }
    }
    else if(itmpb == S_P_BOND)
    {
      /*first partial*/
      CalcCyllinder(xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);

      Calc2Cyllinder(xyz[ia1][0], xyz[ia1][1], xyz[ia1][2],
                     middle[0], middle[1], middle[2],
                     plain[0], plain[1], plain[2],
                     &delta[0],&delta[1],&delta[2]);

      for(j = 0; j < 3; j++) delta2[j] = (middle[j] - xyz[ia1][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
        glPushMatrix();
          SetAtomColor(ia1);
     
          glTranslated(xyz[ia1][0]+delta[0]+delta2[0]*j,
                       xyz[ia1][1]+delta[1]+delta2[1]*j,
                       xyz[ia1][2]+delta[2]+delta2[2]*j);

          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize ,lenbond/20,slices,stacks);
        glPopMatrix();
      }
      /*second partial*/
      CalcCyllinder(middle, xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      for(j = 0; j < 3; j++) delta2[j] = (-middle[j] + xyz[ia2][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
     	glPushMatrix();
          glTranslated(middle[0]+delta[0]+delta2[0]*j,
                       middle[1]+delta[1]+delta2[1]*j,
                       middle[2]+delta[2]+delta2[2]*j);
          glRotated(angle,norm0, norm1, 0.0);
          SetAtomColor(ia2);
          gluCylinder(quadObj, csize, csize,lenbond/20,slices,stacks);
        glPopMatrix();
      }

      /* First part of the full bond */

      CalcCyllinder(xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);
      glPushMatrix();
        SetAtomColor(ia1);
        glTranslated(xyz[ia1][0]-delta[0],
                     xyz[ia1][1]-delta[1],
                     xyz[ia1][2]-delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();

      /* Second part of the full bond */

      CalcCyllinder(xyz[ia2], middle, &angle, &norm0, &norm1, &lenbond);
      glPushMatrix();
        SetAtomColor(ia2);
        glTranslated(xyz[ia2][0]-delta[0],
                     xyz[ia2][1]-delta[1],
                     xyz[ia2][2]-delta[2]);
        glRotated(angle, norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();
    }
    else if (itmpb == DOUBLE_BOND || itmpb == TRIPLE_BOND)
    {
      plain[0]=10.0;
      plain[1]=20.0;
      plain[2]=30.0;
      /* find bounded atom */
      i1 = find_first_bounded(ia1, ia2);
      if (i1 >= 0) for(j = 0; j < 3; j++) plain[j] = xyz[i1][j];

      CalcCyllinder(xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);

      Calc2Cyllinder(xyz[ia1][0], xyz[ia1][1], xyz[ia1][2],
                     middle[0], middle[1], middle[2],
                     plain[0], plain[1], plain[2],
                     &delta[0], &delta[1], &delta[2]);
      if (delta[0] == 0 && delta[1] == 0 && delta[2] == 0) fprintf(stderr, "Bug, Can't draw\n");

      glPushMatrix();
        glTranslated(xyz[ia1][0]+delta[0], xyz[ia1][1]+delta[1], xyz[ia1][2]+delta[2]);
        glRotated(angle,norm0, norm1, 0.0);

        SetAtomColor(ia1);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();

      glPushMatrix();
        glTranslated(xyz[ia1][0]-delta[0], xyz[ia1][1]-delta[1], xyz[ia1][2]-delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();

      if (m->bond[i].bond_type == TRIPLE_BOND)
      {
        glPushMatrix();
          glTranslated(xyz[ia1][0], xyz[ia1][1], xyz[ia1][2]);
          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
        glPopMatrix();
      }

      CalcCyllinder(middle, xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      glPushMatrix();
        glTranslated(middle[0]+delta[0], middle[1]+delta[1], middle[2]+delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        SetAtomColor(ia2);
        gluCylinder(quadObj, csize, csize,lenbond,slices,stacks);
      glPopMatrix();

      glPushMatrix();
        glTranslated(middle[0]-delta[0], middle[1]-delta[1], middle[2]-delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize,lenbond,slices,stacks);
      glPopMatrix();

      if(m->bond[i].bond_type == TRIPLE_BOND)
      {
        glPushMatrix();
          glTranslated(middle[0], middle[1], middle[2]);
          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
         glPopMatrix();
      }
    }

    if (Input_Data.povray || Input_Data.ps)
      print_bond_eps(ia1, ia2, middle);

  }
  gluDeleteQuadric(quadObj);
}

void draw_atoms(void)
{
  int i, j;
/*  double diam;*/
  double tmprad;
/*  SymmetryElement X[8];*/

  GLfloat blackamb[4] = { 0.2F,0.1F,0.1F,0.5F };
  GLfloat blackdif[4] = { 0.2F,0.1F,0.0F,0.5F };
  GLfloat blackspec[4] = { 0.5F,0.5F,0.5F,0.5F };
  GLfloat whiteamb[4] = { 0.7F,0.7F,0.4F,0.5F };
  GLfloat whitedif[4] = { 0.8F,0.7F,0.4F,0.5F };
  GLfloat whitespec[4] = { 0.8F,0.7F,0.4F,0.5F };
  GLfloat copperamb[4] = { 0.24F,0.2F,0.07F,1.0F };
  GLfloat copperdif[4] = { 0.75F,0.61F,0.22F,1.0F };
  GLfloat copperspec[4] = { 0.32F,0.25F,0.17F,1.0F };
  GLfloat darkamb[4] = { 0.10F,0.10F,0.10F,1.0F };
  GLfloat darkdif[4] = { 0.6F,0.6F,0.6F,1.0F };
  GLfloat darkspec[4] = { 0.25F,0.25F,0.25F,1.0F };

  switch(Input_Data.ntexture)
  {
    case 1:
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,whitedif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,whiteamb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,whitespec);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,(GLfloat)40);
      break;
    case 2:
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,blackdif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,blackamb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,blackspec);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,(float)40.f);
      break;
    case 3:
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,darkdif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,darkamb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,darkspec);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,(float)100.f);
      break;
    case 4:
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,copperdif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,copperamb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,copperspec);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,(float)40.f);
      break;
    case 5:
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,whitedif);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,whiteamb);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,whitespec);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,(float)100.f);
  }

  for(i = 0; i < m->natom; i++)
  {
    if (Input_Data.style_atoms) tmprad = m->elem[i].vdw_rad/30.;
    else tmprad = m->elem[i].vdw_rad/3.;
    glPushMatrix();
    glTranslated(m->xyz[i][0], m->xyz[i][1], m->xyz[i][2]);
    SetAtomColor(i);
    for(j = 0; j < m->n_marked; j++)
      if (m->marked[j] == i) glColor3d(0.0f, 1.0f, 1.0f);

#ifdef GTK_GLEXT
    switch(Input_Data.atomshape)
    {
      case 0:
#endif
        luscus_gtk_draw_sphere(1, tmprad * Input_Data.fatness_a, slices, stacks);
        for (j = 0; j < m->n_selected; j++)
        {
          if (m->selected[j] == i)
          {
            if (j == 0) glColor3d(0.0f, 0.0f, 1.0f);
            else glColor3d(1.0f, 0.0f, 1.0f);
            glDisable(GL_LIGHTING);
            luscus_gtk_draw_sphere(0, tmprad * Input_Data.fatness_a + 0.02, 20, 20);
            glEnable(GL_LIGHTING);
          }
        }
        /* povray starts here */
        if (Input_Data.povray)
        {
          fprintf(fl_povray, "\nsphere {\n <%10.6f,%10.6f,%10.6f>, %f\n  pigment { color rgb <%7.4f %7.4f %7.4f>\n  }\n}\n",
                  m->xyz[i][0], m->xyz[i][1], m->xyz[i][2],
                  m->elem[i].vdw_rad/2.F, m->elem[i].color[0], m->elem[i].color[1], m->elem[i].color[2]);
        }
        /* povray ends here */
        /* gveps starts here */
        if(Input_Data.ps)
        {
          gveps_atom(fl_ps, SIZE_SMALL, m->elem[i].vdw_rad/2.F,
                     m->xyz[i][0], m->xyz[i][1], m->xyz[i][2],
                     m->elem[i].color[0], m->elem[i].color[1], m->elem[i].color[2], &ctl);
        }
        /* gveps ends here */
#ifdef GTK_GLEXT
        break;
      case 1:
        luscus_gtk_draw_torus(1, tmprad/2., tmprad, 30, 30);
        break;
      case 2:
        glFrontFace(GL_CW);
        luscus_gtk_draw_teapot(1, tmprad);
        glFrontFace(GL_CCW);
        break;
      case 3:
        luscus_gtk_draw_cube(1, tmprad);
        break;
      default:
        luscus_gtk_draw_sphere(1, tmprad * Input_Data.fatness_a, slices, stacks);
    }
#endif
    glPopMatrix();
  }
}

void draw_bonds(void)
{
  int i, j;
  int i1;
  double csize=0.05*Input_Data.fatness_b;
  double middle[3];
  double angle, norm0, norm1, lenbond, deltalen;
  double delta[3];
  double delta2[3];
  double plain[3];
  int ia1, ia2, itmpb;
  GLUquadricObj *quadObj;

  quadObj=gluNewQuadric();

/*  plain[0]=10.0;
  plain[1]=20.0;
  plain[2]=30.0;*/

  if (quadObj == 0) return;

  for(i = 0; i < m->nbond; i++)
  {
    ia1 = m->bond[i].iat1;
    ia2 = m->bond[i].iat2;
    itmpb = m->bond[i].bond_type;
    for(j = 0; j < 3; j++) middle[j] = 0.5F * (m->xyz[ia1][j] + m->xyz[ia2][j]);
    if (Input_Data.style_atoms && itmpb != NO_BOND) itmpb = LINE_BOND;

    if (itmpb == NO_BOND) continue;
    else if (itmpb == LINE_BOND)
    {
      glDisable(GL_LIGHTING);
      glBegin(GL_LINE_STRIP);
      SetAtomColor(ia1);
      glVertex3dv(m->xyz[ia1]);
      glVertex3dv(middle);
      SetAtomColor(ia2);
      glVertex3dv(middle);
      glVertex3dv(m->xyz[ia2]);
      glEnd();
      glEnable(GL_LIGHTING);
    }
    else if (itmpb == SINGLE_BOND)
    {
      CalcCyllinder(m->xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);
      
      glPushMatrix();
        glTranslated(m->xyz[ia1][0], m->xyz[ia1][1], m->xyz[ia1][2]);
        glRotated(angle, norm0, norm1, 0.0);
        SetAtomColor(ia1);
        gluCylinder(quadObj, csize, csize ,lenbond, slices, stacks);
      glPopMatrix();

      CalcCyllinder(middle, m->xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      glPushMatrix();
        glTranslated(middle[0], middle[1], middle[2]);
        glRotated(angle, norm0, norm1, 0.0);
        SetAtomColor(ia2);
        gluCylinder(quadObj, csize, csize, lenbond, slices, stacks);
      glPopMatrix();
    }
    else if (itmpb == PARTIAL_BOND)
    {
      CalcCyllinder(m->xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);

      for(j = 0; j < 3; j++) delta[j] = (middle[j] - m->xyz[ia1][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
	glPushMatrix();
          SetAtomColor(ia1);
          glTranslated(m->xyz[ia1][0]+delta[0]*j, m->xyz[ia1][1]+delta[1]*j, m->xyz[ia1][2]+delta[2]*j);
          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize, lenbond/20, slices, stacks);
        glPopMatrix();
      }

      CalcCyllinder(middle, m->xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      for(j = 0; j < 3; j++) delta[j] = (-middle[j] + m->xyz[ia2][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
        glPushMatrix();
	  glTranslated(middle[0]+delta[0]*j, middle[1]+delta[1]*j, middle[2]+delta[2]*j);
          glRotated(angle,norm0, norm1, 0.0);
          SetAtomColor(ia2);
          gluCylinder(quadObj, csize, csize, lenbond/20, slices, stacks);
        glPopMatrix();
      }
    }
    else if(itmpb == S_P_BOND)
    {
      /*first partial*/
      CalcCyllinder(m->xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);

      i1 = find_first_bounded(ia2, ia1);
      np3(m->xyz[ia1], m->xyz[ia2], m->xyz[i1], plain);
      np3(m->xyz[ia1], m->xyz[ia2], plain, delta);

      for(j = 0; j < 3; j++) delta[j] -= m->xyz[ia1][j];
      deltalen = 1.0 /(sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]) * 15.0);
      for(j = 0; j < 3; j++) delta[j] *= deltalen;
      if (delta[0] == 0 && delta[1] == 0 && delta[2] == 0) fprintf(stderr, "Bug, Can't draw\n");

/*      Calc2Cyllinder(m->xyz[ia1][0], m->xyz[ia1][1], m->xyz[ia1][2],
                     middle[0], middle[1], middle[2],
                     plain[0], plain[1], plain[2],
                     &delta[0],&delta[1],&delta[2]); */

      for(j = 0; j < 3; j++) delta2[j] = (middle[j] - m->xyz[ia1][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
        glPushMatrix();
          SetAtomColor(ia1);
     
          glTranslated(m->xyz[ia1][0]+delta[0]+delta2[0]*j,
                       m->xyz[ia1][1]+delta[1]+delta2[1]*j,
                       m->xyz[ia1][2]+delta[2]+delta2[2]*j);

          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize ,lenbond/20,slices,stacks);
        glPopMatrix();
      }
      /*second partial*/
      CalcCyllinder(middle, m->xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      for(j = 0; j < 3; j++) delta2[j] = (-middle[j] + m->xyz[ia2][j]) / 10.;

      for(j = 0; j < 10; j++)
      {
     	glPushMatrix();
          glTranslated(middle[0]+delta[0]+delta2[0]*j,
                       middle[1]+delta[1]+delta2[1]*j,
                       middle[2]+delta[2]+delta2[2]*j);
          glRotated(angle,norm0, norm1, 0.0);
          SetAtomColor(ia2);
          gluCylinder(quadObj, csize, csize,lenbond/20,slices,stacks);
        glPopMatrix();
      }

      /* First part of the full bond */

      CalcCyllinder(m->xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);
      glPushMatrix();
        SetAtomColor(ia1);
        glTranslated(m->xyz[ia1][0]-delta[0],
                     m->xyz[ia1][1]-delta[1],
                     m->xyz[ia1][2]-delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();

      /* Second part of the full bond */

      CalcCyllinder(m->xyz[ia2], middle, &angle, &norm0, &norm1, &lenbond);
      glPushMatrix();
        SetAtomColor(ia2);
        glTranslated(m->xyz[ia2][0]-delta[0],
                     m->xyz[ia2][1]-delta[1],
                     m->xyz[ia2][2]-delta[2]);
        glRotated(angle, norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();
    }
    else if (itmpb == DOUBLE_BOND || itmpb == TRIPLE_BOND)
    {
/*      plain[0]=10.0;
      plain[1]=20.0;
      plain[2]=30.0;*/
      /* find bounded atom */
      i1 = find_first_bounded(ia2, ia1);
/*      if (i1 >= 0) for(j = 0; j < 3; j++) plain[j] = m->xyz[i1][j];*/

      CalcCyllinder(m->xyz[ia1], middle, &angle, &norm0, &norm1, &lenbond);
      np3(m->xyz[ia1], m->xyz[ia2], m->xyz[i1], plain);
      np3(m->xyz[ia1], m->xyz[ia2], plain, delta);

/*      Calc2Cyllinder(m->xyz[ia1][0], m->xyz[ia1][1], m->xyz[ia1][2],
                     middle[0], middle[1], middle[2],
                     plain[0], plain[1], plain[2],
                     &delta[0], &delta[1], &delta[2]);*/
      for(j = 0; j < 3; j++) delta[j] -= m->xyz[ia1][j];
      deltalen = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
      if (deltalen < 1.0e-6)
      {
        /*orient bonds somehow*/
        plain[0] = plain[1] = plain[2] = 0.5774;
        np3(m->xyz[ia1], m->xyz[ia2], plain, delta);
        deltalen = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
      }

      for(j = 0; j < 3; j++) delta[j] /= (deltalen * 15.0);
      if (delta[0] == 0.0 && delta[1] == 0.0 && delta[2] == 0.0) fprintf(stderr, "Bug, Can't draw\n");

      glPushMatrix();
        glTranslated(m->xyz[ia1][0]+delta[0], m->xyz[ia1][1]+delta[1], m->xyz[ia1][2]+delta[2]);
        glRotated(angle,norm0, norm1, 0.0);

        SetAtomColor(ia1);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();

      glPushMatrix();
        glTranslated(m->xyz[ia1][0]-delta[0], m->xyz[ia1][1]-delta[1], m->xyz[ia1][2]-delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
      glPopMatrix();

      if (m->bond[i].bond_type == TRIPLE_BOND)
      {
        glPushMatrix();
          glTranslated(m->xyz[ia1][0], m->xyz[ia1][1], m->xyz[ia1][2]);
          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
        glPopMatrix();
      }

      CalcCyllinder(middle, m->xyz[ia2], &angle, &norm0, &norm1, &lenbond);

      glPushMatrix();
        glTranslated(middle[0]+delta[0], middle[1]+delta[1], middle[2]+delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        SetAtomColor(ia2);
        gluCylinder(quadObj, csize, csize,lenbond,slices,stacks);
      glPopMatrix();

      glPushMatrix();
        glTranslated(middle[0]-delta[0], middle[1]-delta[1], middle[2]-delta[2]);
        glRotated(angle,norm0, norm1, 0.0);
        gluCylinder(quadObj, csize, csize,lenbond,slices,stacks);
      glPopMatrix();

      if(m->bond[i].bond_type == TRIPLE_BOND)
      {
        glPushMatrix();
          glTranslated(middle[0], middle[1], middle[2]);
          glRotated(angle,norm0, norm1, 0.0);
          gluCylinder(quadObj, csize, csize ,lenbond,slices,stacks);
         glPopMatrix();
      }
    }

    if (Input_Data.povray || Input_Data.ps)
      print_bond_eps(ia1, ia2, middle);

  }
  gluDeleteQuadric(quadObj);
}

void draw_axes(void)
{
  GLUquadricObj *quadObj;
  GLdouble len, rad;

  len = 0.5 * Calc_Diameter();
  rad = 0.01 * Calc_Diameter();

  if (Input_Data.povray)
  {
    fprintf(fl_povray, "\ncone\n{\n <0.0, 0.0, 0.0>, %f\n  <%f, 0.0, 0.0>, %f\n pigment { color rgb <1.0  0.0  0.0>\n }\n}\n", rad, len, rad);
    fprintf(fl_povray, "\ncone\n{\n <%f, 0.0, 0.0>, %f\n  <%f, 0.0, 0.0>, 0.0\n pigment { color rgb <1.0  0.0  0.0>\n }\n}\n\n", len, 2.0*rad, len+rad);

    fprintf(fl_povray, "\ncone\n{\n <0.0, 0.0, 0.0>, %f\n  <0.0, %f, 0.0>, %f\n pigment { color rgb <0.0  1.0  0.0>\n }\n}\n", rad, len, rad);
    fprintf(fl_povray, "\ncone\n{\n <0.0, %f, 0.0>, %f\n  <0.0, %f, 0.0>, 0.0\n pigment { color rgb <0.0  1.0  0.0>\n }\n}\n\n", len, 2.0*rad, len+rad);

    fprintf(fl_povray, "\ncone\n{\n <0.0, 0.0, 0.0>, %f\n  <0.0, 0.0, %f>, %f\n pigment { color rgb <0.0  0.0  1.0>\n }\n}\n", rad, len, rad);
    fprintf(fl_povray, "\ncone\n{\n <0.0, 0.0, %f>, %f\n  <0.0, 0.0, %f>, 0.0\n pigment { color rgb <0.0  0.0  1.0>\n }\n}\n\n", len, 2.0*rad, len+rad);

    return;
  }

  glDisable(GL_LIGHTING);

  quadObj=gluNewQuadric();

  glPushMatrix();
    if (Input_Data.bw) glColor4d(0.9f,0.9f,0.9f,1.0f);
    else glColor4d(1.0f,0.0f,0.0f,1.0f);
    glRotated(90.0f, 0.0, 1.0, 0.0);
    gluCylinder(quadObj, rad, rad, len, slices, stacks);
    glTranslated(0.0, 0.0, len);
    gluCylinder(quadObj, 2.0*rad, 0.00, 10.0*rad, slices, stacks);
  glPopMatrix();

  glPushMatrix();
    if (Input_Data.bw) glColor4d(0.5f,0.5f,0.5f,1.0f);
    else glColor4d(0.0f,1.0f,0.0f,1.0f);
    glRotated(90.0f, -1.0, 0.0, 0.0);
    gluCylinder(quadObj, rad, rad, len, slices, stacks);
    glTranslated(0.0, 0.0, len);
    gluCylinder(quadObj, 2.0*rad, 0.00, 10.0*rad, slices, stacks);
  glPopMatrix();

  glPushMatrix();
    if (Input_Data.bw) glColor4d(0.0f,0.0f,0.0f,1.0f);
    else glColor4d(0.0f,0.0f,1.0f,1.0f);
    gluCylinder(quadObj, rad, rad, len, slices, stacks);
    glTranslated(0.0, 0.0, len);
    gluCylinder(quadObj, 2.0*rad, 0.00, 10.0*rad, slices, stacks);
  glPopMatrix();

  glEnable(GL_LIGHTING);
  gluDeleteQuadric(quadObj);
}

void draw_dipole(void)
{
  GLUquadricObj *quadObj;
  XYZ nul = {0.F, 0.F, 0.F};
  double angle, norm0, norm1, len;
  double csize = 0.025*Input_Data.fatness_b;

  if (Input_Data.povray)
  {
    fprintf(fl_povray, "\ncone\n{\n <%7.4f, %7.4f, %7.4f>, %7.4f\n  <%7.4f, %7.4f, %7.4f>, %7.4f\n pigment { color rgb <%7.4f %7.4f %7.4f>\n }\n}\n\n",
    0.8*m->dipole[0], 0.8*m->dipole[1], 0.8*m->dipole[2], csize,
    1.8*m->dipole[0], 1.8*m->dipole[1], 1.8*m->dipole[2], csize,
    Input_Data.extracolor[0], Input_Data.extracolor[1], Input_Data.extracolor[2]);

    fprintf(fl_povray, "\ncone\n{\n <%7.4f, %7.4f, %7.4f>, %7.4f\n  <%7.4f, %7.4f, %7.4f>, 0.0\n pigment { color rgb <%7.4f %7.4f %7.4f>\n }\n}\n\n",
    0.8*m->dipole[0], 0.8*m->dipole[1], 0.8*m->dipole[2], 1.65 * csize,
    1.8*m->dipole[0], 1.8*m->dipole[1], 1.8*m->dipole[2],
    Input_Data.extracolor[0], Input_Data.extracolor[1], Input_Data.extracolor[2]);
  }

  quadObj=gluNewQuadric();

  CalcCyllinder(nul, m->dipole, &angle, &norm0, &norm1, &len);

  if (Input_Data.bw)
  {
    double gray;
    gray = grayscale(Input_Data.extracolor[0], Input_Data.extracolor[1], Input_Data.extracolor[2]);
    glColor3d(gray, gray, gray);
  }
  else
    glColor3fv(Input_Data.extracolor);

  glPushMatrix();
    glRotated(angle, norm0, norm1, 0.0);
    gluCylinder(quadObj, csize, csize, 0.8*len, slices, stacks);
  glPopMatrix();

  glPushMatrix();
    glTranslated(0.8*m->dipole[0], 0.8*m->dipole[1], 0.8*m->dipole[2]);
    glRotated(angle, norm0, norm1, 0.0);
    gluCylinder(quadObj, csize*1.65, 0.F, 0.2*len, slices, stacks);

  glPopMatrix();

  gluDeleteQuadric(quadObj);
}

void draw_vectors(void)
{
  int i;
  double angle, norm0, norm1, len;
  GLUquadricObj *quadObj;

  if (!Input_Data.povray) quadObj = gluNewQuadric();

  for (i = 0; i < m->nvector; i++)
  {

    if (Input_Data.povray)
    {
      fprintf(fl_povray, "\ncone\n{\n <%7.4f, %7.4f, %7.4f>, %7.4f\n  <%7.4f, %7.4f, %7.4f>, %7.4f\n pigment { color rgbt <%7.4f %7.4f %7.4f %7.4f>\n }\n}\n\n",
      m->vector1[i][0], m->vector1[i][1], m->vector1[i][2], m->radius[i],
      m->vector1[i][0] + (1.F-m->sharpness[i])*(m->vector2[i][0] - m->vector1[i][0]),
      m->vector1[i][1] + (1.F-m->sharpness[i])*(m->vector2[i][1] - m->vector1[i][1]),
      m->vector1[i][2] + (1.F-m->sharpness[i])*(m->vector2[i][2] - m->vector1[i][2]), m->radius[i],
      m->vector_color[i][0], m->vector_color[i][1], m->vector_color[i][2], 1.0 - m->vector_color[i][3]);
  
      fprintf(fl_povray, "\ncone\n{\n <%7.4f, %7.4f, %7.4f>, %7.4f\n  <%7.4f, %7.4f, %7.4f>, 0.0\n pigment { color rgbt <%7.4f %7.4f %7.4f %7.4f>\n }\n}\n\n",
      m->vector1[i][0] + (1.F-m->sharpness[i])*(m->vector2[i][0] - m->vector1[i][0]),
      m->vector1[i][1] + (1.F-m->sharpness[i])*(m->vector2[i][1] - m->vector1[i][1]),
      m->vector1[i][2] + (1.F-m->sharpness[i])*(m->vector2[i][2] - m->vector1[i][2]), m->radius[i]*1.65,
      m->vector2[i][0], m->vector2[i][1], m->vector2[i][2],
      m->vector_color[i][0], m->vector_color[i][1], m->vector_color[i][2], 1.0 - m->vector_color[i][3]);
    }
    else
    {
      CalcCyllinder(m->vector1[i], m->vector2[i], &angle, &norm0, &norm1, &len);

      if (Input_Data.bw)
      {
        double gray;
        gray = grayscale((double) m->vector_color[i][0], (double) m->vector_color[i][1], (double) m->vector_color[i][2]);
        glColor4d(gray, gray, gray, (double) m->vector_color[i][3]);
      }
      else
        glColor4fv(m->vector_color[i]);
  
      glPushMatrix();
        glTranslated(m->vector1[i][0], m->vector1[i][1], m->vector1[i][2]);
        glRotated(angle, norm0, norm1, 0.0);
        gluCylinder(quadObj, m->radius[i], m->radius[i], (1.F-m->sharpness[i])*len, slices, stacks);
      glPopMatrix();
  
      glPushMatrix();
        glTranslated(m->vector1[i][0] + (1.F-m->sharpness[i])*(m->vector2[i][0] - m->vector1[i][0]),
                     m->vector1[i][1] + (1.F-m->sharpness[i])*(m->vector2[i][1] - m->vector1[i][1]),
                     m->vector1[i][2] + (1.F-m->sharpness[i])*(m->vector2[i][2] - m->vector1[i][2]));
        glRotated(angle, norm0, norm1, 0.0);
        gluCylinder(quadObj, 1.65 * m->radius[i], 0.F, m->sharpness[i]*len, slices, stacks);
      glPopMatrix();
    }
  }

  if (!Input_Data.povray) gluDeleteQuadric(quadObj);
}

void draw_triangles(void)
{
  int i;
  XYZ norm;
  for(i = 0; i < m->ntriangle; i++)
  {
    np3(m->triangle1[i], m->triangle2[i], m->triangle3[i], norm); /*DO test if np3 givs normal to two vectors!!!!!!!!*/
    if (Input_Data.povray)
    {
      fprintf(fl_povray, "\nsmooth_triangle\n{\n<%7.4f %7.4f %7.4f>, <%7.4f %7.4f %7.4f>\n<%7.4f %7.4f %7.4f>, <%7.4f %7.4f %7.4f>\n<%7.4f %7.4f %7.4f>, <%7.4f %7.4f %7.4f>\n  pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
              m->triangle1[i][0], m->triangle1[i][1], m->triangle1[i][2],  norm[0], norm[1], norm[2],
              m->triangle2[i][0], m->triangle2[i][1], m->triangle2[i][2],  norm[0], norm[1], norm[2],
              m->triangle3[i][0], m->triangle3[i][1], m->triangle3[i][2],  norm[0], norm[1], norm[2],
              m->triangle_color[i][0], m->triangle_color[i][1], m->triangle_color[i][2], 1.0 - m->triangle_color[i][3]);
    }
    else
    {
      if (Input_Data.bw)
      {
        double gray;
        gray = grayscale((double) m->triangle_color[i][0], (double) m->triangle_color[i][1], (double) m->triangle_color[i][2]);
        glColor4d(gray, gray, gray, (double) m->triangle_color[i][3]);
      }
      else
        glColor4fv(m->triangle_color[i]);
      glBegin(GL_POLYGON);
        glNormal3dv(norm);
        glVertex3dv(m->triangle1[i]);
        glVertex3dv(m->triangle2[i]);
        glVertex3dv(m->triangle3[i]);
      glEnd();
    }
  }
}

void draw_spheres(void)
{
  int i;

/*  glDepthMask(GL_FALSE);*/
  for(i = 0; i < m->nsphere; i++)
  {
    if(Input_Data.povray)
    {
      fprintf(fl_povray, "sphere {\n  <%7.4f %7.4f %7.4f> %7.4f\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
              m->sphere_center[i][0], m->sphere_center[i][1], m->sphere_center[i][2], m->sphere_radius[i], 
              m->sphere_color[i][0], m->sphere_color[i][1], m->sphere_color[i][2], 1.0 - m->sphere_color[i][3]);
    }
    else
    {
      glPushMatrix();
        if (Input_Data.bw)
        {
          double gray;
          gray = grayscale((double) m->sphere_color[i][0], (double) m->sphere_color[i][1], (double) m->sphere_color[i][2]);
          glColor4d(gray, gray, gray, (double) m->sphere_color[i][3]);
        }
        else
          glColor4fv(m->sphere_color[i]);
/*      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);*/
        glTranslated(m->sphere_center[i][0], m->sphere_center[i][1], m->sphere_center[i][2]);
        luscus_gtk_draw_sphere(1, m->sphere_radius[i], 2 * slices, 2* stacks);
      glPopMatrix();
    }
  }
/*  glDepthMask(GL_TRUE);*/
}

void draw_surfaces(void)
{
  int i;
  XYZ ut1, ut2, ut3, ut4;
  XYZ norm;

  for(i = 0; i < m->nsurf; i++)
  {
    pl3to4(m->surf1[i], m->surf2[i], m->surf3[i], ut1, ut2, ut3, ut4);
    np3(m->surf1[i], m->surf2[i], m->surf3[i], norm);
    if(Input_Data.povray)
    {
      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
              ut1[0], ut1[1], ut1[2], ut2[0], ut2[1], ut2[2], ut3[0], ut3[1], ut3[2], ut4[0], ut4[1], ut4[2],
              m->surf_color[i][0], m->surf_color[i][1], m->surf_color[i][2], 1.0 - m->surf_color[i][3]);
    }
    else
    {
      double gray;
      if(Input_Data.bw)
      {
        gray=grayscale((double) m->surf_color[i][0], (double) m->surf_color[i][1], (double) m->surf_color[i][2]);
        printf("gray = %f\n", gray);
        glColor4d(gray, gray, gray, (double) m->surf_color[i][3]);
      }
      else
        glColor4fv(m->surf_color[i]);

      glBegin(GL_QUADS);
        glNormal3dv(norm);
        glVertex3dv(ut1);
        glVertex3dv(ut2);
        glVertex3dv(ut3);
        glVertex3dv(ut4);
      glEnd();
    }
  }
}

void draw_cells(void)
{
  int i, j;
  XYZ in0, in1, in2, in3, in12, in13, in23, in123;
  XYZ no;


  for(i = 0; i < m->ncells; i++)
  {
    if (!Input_Data.povray)
    {
      if (Input_Data.bw)
      {
        double gray;
        gray = grayscale((double) m->cell_color[i][0], (double) m->cell_color[i][1], (double) m->cell_color[i][2]);
        glColor4d(gray, gray, gray, (double) m->cell_color[i][3]);
      }
      else
        glColor4fv(m->cell_color[i]);
    }
    for(j = 0; j < 3; j++)
    {
      in0[j] = m->cell1[i][j];
      in1[j] = m->cell2[i][j];
      in2[j] = m->cell3[i][j];
      in3[j] = m->cell4[i][j];
    }

    pl4to6(in0, in1, in2, in3, in12, in13, in23, in123);
    if (Input_Data.povray)
    {
      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
      in0[0], in0[1], in0[2], in1[0], in1[1], in1[2], in12[0], in12[1], in12[2], in2[0], in2[1], in2[2],
      m->cell_color[i][0], m->cell_color[i][1], m->cell_color[i][2], 1.0 - m->cell_color[i][3]);

      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
      in0[0], in0[1], in0[2], in2[0], in2[1], in2[2], in23[0], in23[1], in23[2], in3[0], in3[1], in3[2],
      m->cell_color[i][0], m->cell_color[i][1], m->cell_color[i][2], 1.0 - m->cell_color[i][3]);

      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
      in0[0], in0[1], in0[2], in3[0], in3[1], in3[2], in13[0], in13[1], in13[2], in1[0], in1[1], in1[2],
      m->cell_color[i][0], m->cell_color[i][1], m->cell_color[i][2], 1.0 - m->cell_color[i][3]);

      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
      in1[0], in1[1], in1[2], in12[0], in12[1], in12[2], in123[0], in123[1], in123[2], in13[0], in13[1], in13[2],
      m->cell_color[i][0], m->cell_color[i][1], m->cell_color[i][2], 1.0 - m->cell_color[i][3]);

      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
      in2[0], in2[1], in2[2], in23[0], in23[1], in23[2], in123[0], in123[1], in123[2], in12[0], in12[1], in12[2],
      m->cell_color[i][0], m->cell_color[i][1], m->cell_color[i][2], 1.0 - m->cell_color[i][3]);

      fprintf(fl_povray, "polygon { 4\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   <%7.4f %7.4f %7.4f>\n   pigment{ rgbt <%7.4f %7.4f %7.4f %7.4f> }\n}\n",
      in3[0], in3[1], in3[2], in23[0], in23[1], in23[2], in123[0], in123[1], in123[2], in13[0], in13[1], in13[2],
      m->cell_color[i][0], m->cell_color[i][1], m->cell_color[i][2], 1.0 - m->cell_color[i][3]);
    }
    else
    {
/*      glPushMatrix();
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(3.0);
        glBegin(GL_LINES);
          glVertex3dv(in0);
          glVertex3dv(in1);
        
          glVertex3dv(in1);
          glVertex3dv(in12);

          glVertex3dv(in12);
          glVertex3dv(in2);

          glVertex3dv(in2);
          glVertex3dv(in0);

          glVertex3dv(in2);
          glVertex3dv(in23);

          glVertex3dv(in23);
          glVertex3dv(in3);

          glVertex3dv(in0);
          glVertex3dv(in3);

          glVertex3dv(in3);
          glVertex3dv(in13);

          glVertex3dv(in13);
          glVertex3dv(in1);

          glVertex3dv(in12);
          glVertex3dv(in123);

          glVertex3dv(in123);
          glVertex3dv(in13);

          glVertex3dv(in23);
          glVertex3dv(in123);
        glEnd();
      glPopMatrix();*/ /*frame-cell*/

      glBegin(GL_QUADS);
        np3(in0, in1, in12, no);
        glNormal3dv(no);
        glVertex3dv(in0);
        glVertex3dv(in1);
        glVertex3dv(in12);
        glVertex3dv(in2);

        np3(in0, in2, in23, no);
        glNormal3dv(no);
        glVertex3dv(in0);
        glVertex3dv(in2);
        glVertex3dv(in23);
        glVertex3dv(in3);

        np3(in0, in3, in13, no);
        glNormal3dv(no);
        glVertex3dv(in0);
        glVertex3dv(in3);
        glVertex3dv(in13);
        glVertex3dv(in1);

        np3(in1, in12, in123, no);
        glNormal3dv(no);
        glVertex3dv(in1);
        glVertex3dv(in12);
        glVertex3dv(in123);
        glVertex3dv(in13);

        np3(in2, in23, in123, no);
        glNormal3dv(no);
        glVertex3dv(in2);
        glVertex3dv(in23);
        glVertex3dv(in123);
        glVertex3dv(in12);

        np3(in3, in23, in123, no);
        glNormal3dv(no);
        glVertex3dv(in3);
        glVertex3dv(in23);
        glVertex3dv(in123);
        glVertex3dv(in13);
      glEnd();
    }
  }
}

void draw_grid(void)
{
  glDepthMask(GL_FALSE);
  if (surf == NULL) return;

  /* Draw all surfaces */
  if (Input_Data.show_epot)
  {
    msrf_Draw1(All_Surfaces, Input_Data.electrostatic_poten_color, Input_Data.povray, fl_povray,
              Input_Data.ps, fl_ps, &ctl);
  }
  else
    msrf_Draw(All_Surfaces, Input_Data.neg_pos_color, Input_Data.povray, fl_povray,
              Input_Data.ps, fl_ps, &ctl);

  glDepthMask(GL_TRUE);
}

void draw_selected_go(void)
{
  int i;
  double angle, norm0, norm1, len;
  int slen, fontsize;
  XYZ in0, in1, in2, in3, in12, in13, in23, in123;
  XYZ cen;
  GLUquadricObj *quadObj = gluNewQuadric();
  gluQuadricDrawStyle(quadObj, GLU_LINE);
/*static int selected_go_num = 0;
static int selected_go_type = 0;*/
  glDisable(GL_LIGHTING);
  switch(selected_go_type)
  {
    case SPHERE:
      glPushMatrix();
        glColor3d(0.0, 1.0, 1.0);
        glTranslated(m->sphere_center[selected_go_num][0], m->sphere_center[selected_go_num][1], m->sphere_center[selected_go_num][2]);
        luscus_gtk_draw_sphere(0, m->sphere_radius[selected_go_num]+0.01, slices, stacks);
      glPopMatrix();
    break;
    case VECTOR:
      CalcCyllinder(m->vector1[selected_go_num], m->vector2[selected_go_num], &angle, &norm0, &norm1, &len);
      glColor3d(0.0, 1.0, 1.0);
      glPushMatrix();
        glTranslated(m->vector1[selected_go_num][0], m->vector1[selected_go_num][1], m->vector1[selected_go_num][2]);
        glRotated(angle, norm0, norm1, 0.0);
        gluCylinder(quadObj, m->radius[selected_go_num]+0.01, m->radius[selected_go_num]+0.01, (1.F-m->sharpness[selected_go_num])*len, slices, stacks);
      glPopMatrix();

      glPushMatrix();
        glTranslated(m->vector1[selected_go_num][0] + (1.F-m->sharpness[selected_go_num])*(m->vector2[selected_go_num][0] - m->vector1[selected_go_num][0]),
                     m->vector1[selected_go_num][1] + (1.F-m->sharpness[selected_go_num])*(m->vector2[selected_go_num][1] - m->vector1[selected_go_num][1]),
                     m->vector1[selected_go_num][2] + (1.F-m->sharpness[selected_go_num])*(m->vector2[selected_go_num][2] - m->vector1[selected_go_num][2]));
        glRotated(angle, norm0, norm1, 0.0);
        gluCylinder(quadObj, 1.65 * m->radius[selected_go_num], 0.F, m->sharpness[selected_go_num]*len, slices, stacks);
      glPopMatrix();
      break;
    case TRIANGLE:
      glColor3d(0.0, 1.0, 1.0);
      glBegin(GL_LINE_LOOP);
        glVertex3dv(m->triangle1[selected_go_num]);
        glVertex3dv(m->triangle2[selected_go_num]);
        glVertex3dv(m->triangle3[selected_go_num]);
      glEnd();
      break;
    case SURFACE:
      glColor3d(0.0, 1.0, 1.0);
      pl3to4(m->surf1[selected_go_num], m->surf2[selected_go_num], m->surf3[selected_go_num], in0, in1, in2, in3);
      glBegin(GL_LINES);
        glVertex3dv(in0);
        glVertex3dv(in2);

        glVertex3dv(in1);
        glVertex3dv(in3);
      glEnd();
      break;
    case CELL:
      glColor3d(0.0, 1.0, 1.0);

      for(i = 0; i < 3; i++)
      {
        in0[i] = m->cell1[selected_go_num][i];
        in1[i] = m->cell2[selected_go_num][i];
        in2[i] = m->cell3[selected_go_num][i];
        in3[i] = m->cell4[selected_go_num][i];
      }
      pl4to6(in0, in1, in2, in3, in12, in13, in23, in123); /*??????*/
      for(i = 0; i < 3; i++) cen[i] = 0.5 * (in0[i] + in123[i]);
      for(i = 0; i < 3; i++)
      {
        in0[i] = cen[i] + 1.01 * (in0[i] - cen[i]);
        in1[i] = cen[i] + 1.01 * (in1[i] - cen[i]);
        in2[i] = cen[i] + 1.01 * (in2[i] - cen[i]);
        in3[i] = cen[i] + 1.01 * (in3[i] - cen[i]);
        in12[i] = cen[i] + 1.01 * (in12[i] - cen[i]);
        in13[i] = cen[i] + 1.01 * (in13[i] - cen[i]);
        in23[i] = cen[i] + 1.01 * (in23[i] - cen[i]);
        in123[i] = cen[i] + 1.01 * (in123[i] - cen[i]);
      }

      glBegin(GL_LINE_LOOP);
        glVertex3dv(in0);
        glVertex3dv(in1);
        glVertex3dv(in12);
        glVertex3dv(in2);
      glEnd();
      glBegin(GL_LINE_LOOP);
        glVertex3dv(in3);
        glVertex3dv(in23);
        glVertex3dv(in123);
        glVertex3dv(in13);
      glEnd();
      glBegin(GL_LINES);
        glVertex3dv(in0);
        glVertex3dv(in3);

        glVertex3dv(in1);
        glVertex3dv(in13);

        glVertex3dv(in2);
        glVertex3dv(in23);

        glVertex3dv(in12);
        glVertex3dv(in123);
      glEnd();

      break;
/*    case TEXTBOX:
      slen = strlen(m->textboxes[selected_go_num].message);

      i = strlen(m->textboxes[selected_go_num].font);
      while(m->textboxes[selected_go_num].font[i] == 32) i--;
      while(m->textboxes[selected_go_num].font[i] != 32) i--;
      fontsize = atoi(m->textboxes[selected_go_num].font+i);

      wr = slen * (int) (fontsize * 12.0 / 16.0) + 5;
      hr = (int) (fontsize * 16.0 / 12.0) + 1;

      glPushMatrix();
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      glColor3d(0.0, 1.0, 1.0);

      glBegin(GL_LINES);
        glVertex3d((double) m->textboxes[selected_go_num].coord_x + 1, (double) m->textboxes[selected_go_num].coord_y + 1, 0.0);
        glVertex3d((double) m->textboxes[selected_go_num].coord_x + 1, (double) m->textboxes[selected_go_num].coord_y - hr - 1, 0.0);
        glVertex3d((double) m->textboxes[selected_go_num].coord_x + wr + 1, (double) m->textboxes[selected_go_num].coord_y - hr - 1, 0.0);
        glVertex3d((double) m->textboxes[selected_go_num].coord_x + wr + 1, (double) m->textboxes[selected_go_num].coord_y + 1, 0.0);
      glEnd();

      glPopMatrix();*/

    default: break;
  }
  glEnable(GL_LIGHTING);
  gluDeleteQuadric(quadObj);
}

void luscus_draw_3d(void)
{
  if (m)
  {
    glNewList(list_3d, GL_COMPILE_AND_EXECUTE);
    
    if (m->nvibr && !Input_Data.povray && !Input_Data.ps)
      draw_vibrations();
    else
    {
      if (m->natom && !Input_Data.hide_atoms) draw_atoms();

      if (m->nbond && !Input_Data.hide_bonds) draw_bonds();
    }

    if (Input_Data.show_axis) draw_axes();

    if (m->ishow & HAS_DIPOLE) draw_dipole();

    if (m->nvector) draw_vectors();

    if (m->ntriangle) draw_triangles();

    if (m->nsphere) draw_spheres();

    if (m->nsurf) draw_surfaces();

    if (m->ncells) draw_cells();

    if (m->ngrids) draw_grid();

    if (selected_go_type) draw_selected_go();

    glEndList();
  }
}

void luscus_draw_2d(void)
{
  gint w, h;

  gv_gtk_get_screen_size(&w, &h);

  if (shift_move == 1) draw_rectangle(w, h);
  if (!m) return;
  if (m->ntextboxes) 
    draw_textboxes(w, h);
  if (textbox_state == 2) draw_textbox_rectangle(w, h);

  if (!m->natom) return;
  if (m->pixdata)
    draw_labels(w, h);

}

void draw_rectangle(int w, int h)
{
  int t00, t01, t11, t10; /*coordinates*/
  GLdouble x00, x01, x11, x10;

  t00=MIN(shift_move_x_b, shift_move_x_e);
  t10=MAX(shift_move_x_b, shift_move_x_e);
  t01=MIN(shift_move_y_b, shift_move_y_e);
  t11=MAX(shift_move_y_b, shift_move_y_e);

  x00 = x_size * 2.F * (GLdouble) t00 / (GLdouble) w - x_size - camera_x;
  x10 = x_size * 2.F * (GLdouble) t10 / (GLdouble) w - x_size - camera_x;
  x01 = y_size * (1.F - 2.F * (GLdouble) t01 / (GLdouble) h) - camera_y;
  x11 = y_size * (1.F - 2.F * (GLdouble) t11 / (GLdouble) h) - camera_y;

  glPushMatrix();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glColor4d(0.0, 1.0, 1.0, 1.0);
  glBegin(GL_LINES);

  glVertex3d(x00, x11, 0.0); glVertex3d(x10, x11, 0.0);
  glVertex3d(x10, x11, 0.0); glVertex3d(x10, x01, 0.0);
  glVertex3d(x10, x01, 0.0); glVertex3d(x00, x01, 0.0);
  glVertex3d(x00, x01, 0.0); glVertex3d(x00, x11, 0.0);

  glEnd();
  glPopMatrix();
}

void draw_labels(int w, int h)
{
  int i, j;
  int show;
  float *plotwidth;

  GLdouble *screen_coord;

  GLdouble modelMatrix[16], projMatrix[16];
  GLint viewport[4];

  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);

  screen_coord = (GLdouble*) malloc(3 * m->natom * sizeof(GLdouble));

/*  if (w / daw != 1 && h / dah != 1)
  {

    -------render strings with large fonts!
  }*/
  plotwidth = (float*) malloc(m->natom * sizeof(float));


  for(i = 0; i < m->natom; i++)
  {
    plotwidth[i] = 0.0;
    gluProject(m->xyz[i][0],     m->xyz[i][1],       m->xyz[i][2],
               modelMatrix,      projMatrix,         viewport,
               screen_coord+i*3, screen_coord+i*3+1, screen_coord+i*3+2);

    screen_coord[i*3] =   x_size * 2.F * screen_coord[i*3]   / (GLdouble) w - x_size - camera_x;
    screen_coord[i*3+1] = y_size * 2.F * screen_coord[i*3+1] / (GLdouble) h - y_size - camera_y;
  }

  glPushMatrix();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
/*  glColor3fv(Input_Data.label_color);*/

  glPixelTransferf(GL_RED_SCALE, 0.0); /*extremly slow!!!*/
  glPixelTransferf(GL_GREEN_SCALE, 0.0);
  glPixelTransferf(GL_BLUE_SCALE, 0.0);
  glPixelTransferf(GL_RED_BIAS, Input_Data.label_color[0]);
  glPixelTransferf(GL_GREEN_BIAS, Input_Data.label_color[1]);
  glPixelTransferf(GL_BLUE_BIAS, Input_Data.label_color[2]);

  for(i = 0; i < m->natom; i++)
  {
    show = 1;

    if (screen_coord[i*3]+camera_x < -x_size) show = 0;
    if (screen_coord[i*3+1]+camera_y < -y_size) show = 0;
    if (screen_coord[i*3]+camera_x > x_size) show = 0;
    if (screen_coord[i*3+1]+camera_y > y_size) show = 0;

    if (show)
      for(j = 0; j < m->natom; j++)
        if (fabs(screen_coord[j*3]   - screen_coord[i*3]) < 0.5/scale &&
            fabs(screen_coord[j*3+1] - screen_coord[i*3+1]) < 0.5/scale &&
                 screen_coord[j*3+2] < screen_coord[i*3+2]) show = 0;
    if (show)
    {
      if (Input_Data.label_atoms & (1 << 1))
        if (m->pixdata[i].pix_index.pixels)
        {
          glRasterPos3d(screen_coord[i*3] + plotwidth[i], screen_coord[i*3+1], molecule_diameter);
          glDrawPixels(m->pixdata[i].pix_index.width, m->pixdata[i].pix_index.height, GL_ALPHA, GL_UNSIGNED_BYTE, m->pixdata[i].pix_index.pixels);
          plotwidth[i] += m->pixdata[i].pix_index.width * x_size * 2.F / (GLdouble) w;
        }

      if (Input_Data.label_atoms & (1 << 3))
        if (m->pixdata[i].pix_symm.pixels)
        {
          glRasterPos3d(screen_coord[i*3] + plotwidth[i], screen_coord[i*3+1], molecule_diameter);
          glDrawPixels(m->pixdata[i].pix_symm.width, m->pixdata[i].pix_symm.height, GL_ALPHA, GL_UNSIGNED_BYTE, m->pixdata[i].pix_symm.pixels);
          plotwidth[i] += m->pixdata[i].pix_symm.width * x_size * 2.F / (GLdouble) w;
        }

      if (Input_Data.label_atoms & (1 << 4) && m->ishow & HAS_ATOMNAMES)
        if (m->pixdata[i].pix_name.pixels)
        {
          glRasterPos3d(screen_coord[i*3] + plotwidth[i], screen_coord[i*3+1], molecule_diameter);
          glDrawPixels(m->pixdata[i].pix_name.width, m->pixdata[i].pix_name.height, GL_ALPHA, GL_UNSIGNED_BYTE, m->pixdata[i].pix_name.pixels);
          plotwidth[i] += m->pixdata[i].pix_name.width * x_size * 2.F / (GLdouble) w;
        }

      if (Input_Data.label_atoms & (1 << 2) && m->ishow & HAS_ATOMNUMS)
        if (m->pixdata[i].pix_addnum.pixels)
        {
          glRasterPos3d(screen_coord[i*3] + plotwidth[i], screen_coord[i*3+1], molecule_diameter);
          glDrawPixels(m->pixdata[i].pix_addnum.width, m->pixdata[i].pix_addnum.height, GL_ALPHA, GL_UNSIGNED_BYTE, m->pixdata[i].pix_addnum.pixels);
          plotwidth[i] += m->pixdata[i].pix_addnum.width * x_size * 2.F / (GLdouble) w;
        }

      if (Input_Data.label_atoms & (1 << 5) && m->ishow & HAS_MULLIKEN)
        if (m->pixdata[i].pix_charge_m.pixels)
        {
          glRasterPos3d(screen_coord[i*3] + plotwidth[i], screen_coord[i*3+1], molecule_diameter);
          glDrawPixels(m->pixdata[i].pix_charge_m.width, m->pixdata[i].pix_charge_m.height, GL_ALPHA, GL_UNSIGNED_BYTE, m->pixdata[i].pix_charge_m.pixels);
          plotwidth[i] += m->pixdata[i].pix_charge_m.width * x_size * 2.F / (GLdouble) w;
        }

      if (Input_Data.label_atoms & (1 << 6) && m->ishow & HAS_LOPROP)
        if (m->pixdata[i].pix_charge_l.pixels)
        {
          glRasterPos3d(screen_coord[i*3] + plotwidth[i], screen_coord[i*3+1], molecule_diameter);
          glDrawPixels(m->pixdata[i].pix_charge_l.width, m->pixdata[i].pix_charge_l.height, GL_ALPHA, GL_UNSIGNED_BYTE, m->pixdata[i].pix_charge_l.pixels);
          plotwidth[i] += m->pixdata[i].pix_charge_l.width * x_size * 2.F / (GLdouble) w;
        }

    }
  }

  glPopMatrix();


  free(plotwidth);
  free(screen_coord);
}

void draw_textboxes(int w, int h)
{
  int i/*, j, k*/;
  int test_x, test_y;
  int same_aspect_ratio = 1;
  double ratio_x, ratio_y;

  if (!m->textboxes) return;


  gv_gtk_get_screen_size(&test_x, &test_y);
  if (test_x - w > 5 && test_y - h < 5) same_aspect_ratio = 0;

  glPushMatrix();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);


  for(i = 0; i < m->ntextboxes; i++)
  {
    if (m->textboxes[i].pixtext.pixels)
    {
      glPixelTransferf(GL_RED_SCALE, 0.0); 
      glPixelTransferf(GL_GREEN_SCALE, 0.0);
      glPixelTransferf(GL_BLUE_SCALE, 0.0);
      glPixelTransferf(GL_RED_BIAS, m->textboxes[i].color[0]);
      glPixelTransferf(GL_GREEN_BIAS, m->textboxes[i].color[1]);
      glPixelTransferf(GL_BLUE_BIAS, m->textboxes[i].color[2]);

      glRasterPos3d(x_size * (2.F * m->textboxes[i].coord_x / (GLdouble) test_x - 1.0) - camera_x,
                    y_size * (1.F - 2.F * (m->textboxes[i].coord_y + m->textboxes[i].pixtext.height-2) / (GLdouble) test_y) - camera_y,
                    molecule_diameter);
      if (!same_aspect_ratio)
      {
        ratio_x = (GLdouble) w / (GLdouble) test_x;
        ratio_y = (GLdouble) h / (GLdouble) test_y;
        if (ratio_x < ratio_y) ratio_y = ratio_x;
        else ratio_x = ratio_y;

        glPixelZoom(ratio_x, ratio_y);
      }

      glDrawPixels(m->textboxes[i].pixtext.width, m->textboxes[i].pixtext.height,
                   GL_ALPHA, GL_UNSIGNED_BYTE, m->textboxes[i].pixtext.pixels);

      glPixelTransferf(GL_RED_SCALE, 1.0); 
      glPixelTransferf(GL_GREEN_SCALE, 1.0);
      glPixelTransferf(GL_BLUE_SCALE, 1.0);
      glPixelTransferf(GL_RED_BIAS, 0.0);
      glPixelTransferf(GL_GREEN_BIAS, 0.0);
      glPixelTransferf(GL_BLUE_BIAS, 0.0);
    }
  }
  glPopMatrix();
}

void draw_textbox_rectangle(int w, int h)
{
  int i;
  int wr, hr; /*coordinates*/
  int slen;
  int fontsize;

  int t00, t01, t11, t10; /*coordinates*/
  GLdouble x00, x01, x11, x10;

  if (!m) return;
  if (!m->ntextboxes) return;

  if (m->textboxes[m->ntextboxes-1].message == NULL) slen = 0;
  else slen = strlen(m->textboxes[m->ntextboxes-1].message);

  if (m->textboxes[m->ntextboxes-1].font == NULL) fontsize = 12;
  else
  {
    i = strlen(m->textboxes[m->ntextboxes-1].font);
    while(m->textboxes[m->ntextboxes-1].font[i] == 32) i--;
    while(m->textboxes[m->ntextboxes-1].font[i] != 32) i--;
    fontsize = atoi(m->textboxes[m->ntextboxes-1].font+i);
  }

/*  gc = gdk_gc_new(GDK_DRAWABLE(gld));
  gdk_gc_set_rgb_fg_color(gc, &rectangle_color);*/

  wr = slen * (int) (fontsize * 12.0 / 16.0) + 5;
  hr = (int) (fontsize * 16.0 / 12.0) + 1;

  t00 = m->textboxes[m->ntextboxes-1].coord_x - 1;
  t10 = t00 + m->textboxes[m->ntextboxes-1].pixtext.width+5;
  t01 = m->textboxes[m->ntextboxes-1].coord_y - 1;
  t11 = t01 + (m->textboxes[m->ntextboxes-1].pixtext.height ? m->textboxes[m->ntextboxes-1].pixtext.height : 12) +2;

  x00 = x_size * 2.F * (GLdouble) t00 / (GLdouble) w - x_size - camera_x;
  x10 = x_size * 2.F * (GLdouble) t10 / (GLdouble) w - x_size - camera_x;
  x01 = y_size * (1.F - 2.F * (GLdouble) t01 / (GLdouble) h) - camera_y;
  x11 = y_size * (1.F - 2.F * (GLdouble) t11 / (GLdouble) h) - camera_y;

  glPushMatrix();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glColor4d(0.0, 1.0, 1.0, 1.0);
  glBegin(GL_LINES);

  glVertex3d(x00, x11, 0.0); glVertex3d(x10, x11, 0.0);
  glVertex3d(x10, x11, 0.0); glVertex3d(x10, x01, 0.0);
  glVertex3d(x10, x01, 0.0); glVertex3d(x00, x01, 0.0);
  glVertex3d(x00, x01, 0.0); glVertex3d(x00, x11, 0.0);

  glEnd();
  glPopMatrix();

/*  gdk_gc_unref(gc);*/

}

/*----- OFFSCREEN RENDREING -----*/

/*void luscus_save_builtin_graphic_format(gchar* filename, gchar* type, gint width, gint height)
{
  printf("SHOULD BE RECODED!\n");
  return;
}*/

/*----- CALLBACK FUNCTIONS -----*/

void callback_maximize(GtkWidget *widget, gpointer data)
{
  fullscreen=!fullscreen;
  if(fullscreen) gtk_window_fullscreen(GTK_WINDOW(window));
  else gtk_window_unfullscreen(GTK_WINDOW(window));
}

static gboolean luscus_gtk_mouse_button_press(GtkWidget *da, GdkEventButton *event, gpointer data)
{
  GLint viewport[4], viewport_r[4];
  GtkAllocation allocation;
  gdouble w, h;

  int i;
  int col_sel;
  int ir2, rrr2;
  int float_select=-1;
  double rrr_select;
  int ifound;
  int pressed_center[2];
  double ff;

BEGIN_OGL_STUFF;
#ifdef GTK_GLEXT
  GdkWindow *window;
  window = gtk_widget_get_window(drawingarea);
#endif

  gtk_widget_get_allocation(drawingarea, &allocation);
  w = (double) allocation.width;
  h = (double) allocation.height;

  if (!m) return FALSE;

  if (textbox_state == 3)
  {
    m->textboxes[itextbox].coord_x = (int) event->x;
    m->textboxes[itextbox].coord_y = (int) event->y;
    textbox_state = 0;
    redraw();
    append_backup();
    return TRUE;
  }
  if (textbox_state == 2)
  {
    textbox_state = 0;
    draw_pixdata_textbox(m->ntextboxes-1);
    redraw();
    append_backup();
    return TRUE;
  }
  else if (textbox_state == 1)
  {

    gdk_window_set_cursor(window, NULL);
    if (m->ntextboxes)
    {
      m->textboxes[m->ntextboxes-1].coord_x = (int) event->x;
      m->textboxes[m->ntextboxes-1].coord_y = (int) event->y;
    }
    textbox_state = 2;
    luscus_gtk_update_3Dobject_info();
  }

  if (event->button == 2)
  {
    if (m->n_selected) unselect_all();
    else if (m->n_marked) unmark_all();
/*    status_sel();*/
    redraw();
  }

  if (event->button == 1 || event->button == 3)
  {
    x_prev = 2. * event->x / w - 1.;
    y_prev = 2. * (h - event->y) / h - 1.;

    pressed_center[0]= (int) event->x ;
    pressed_center[1] = h - event->y;

    rrr_select=3000.F;

    for(i = 0; i < m->natom; i++)
    {
      glPushMatrix();
      glRasterPos3d(m->xyz[i][0], m->xyz[i][1], m->xyz[i][2]);
      glGetIntegerv(GL_CURRENT_RASTER_POSITION,&viewport[0]);
      glPopMatrix();

      ff=m->elem[i].vdw_rad/2.;

      glPushMatrix();
      glRasterPos3d(m->xyz[i][0]+ff, m->xyz[i][1]+ff, m->xyz[i][2]+ff);
      glGetIntegerv(GL_CURRENT_RASTER_POSITION,&viewport_r[0]);
      glPopMatrix();

      rrr2=(viewport[0]-pressed_center[0])*(viewport[0]-pressed_center[0])+
           (viewport[1]-pressed_center[1])*(viewport[1]-pressed_center[1]);

      ir2=(viewport_r[0]-viewport[0])*(viewport_r[0]-viewport[0])+
          (viewport_r[1]-viewport[1])*(viewport_r[1]-viewport[1]);

      if (rrr2 < ir2 && rrr2 < rrr_select)
      {
        rrr_select=rrr2;
        float_select=i;
      }
    }
  }

  if (event->button == 1)
  {
    if (m->n_selected < 4 && float_select > -1)
    {
      col_sel=1;
      for(i = 0; i < 4; i++) if(m->selected[i] == float_select) col_sel=0;

      if(col_sel==1)
      {
        m->selected[m->n_selected] = float_select;
        m->n_selected++;
        rerender_3d();
      }
      if (m->n_selected == 1) luscus_xyz_editor_select_row(m->selected[0]);
      luscus_gtk_update_upon_select_or_mark();
    }
  }
  else if (event->button == 3)
  {
    if(float_select >= 0)
    {
      ifound=-1;
      for(i = 0; i < m->n_marked; i++)
        if (m->marked[i]==float_select) ifound=i;

      if(ifound==-1)
      {
        mark_atom(float_select);
        rerender_3d();
/*        greyStek[ngreyStek]=float_select;
        ngreyStek++;*/
/*        on_unselect();*/
      }
      else
      {
        unmark_atom(ifound);
        rerender_3d();
/*        for(i=ifound; i<ngreyStek-1; i++)
           greyStek[i]=greyStek[i+1];
        ngreyStek--;*/
/*        on_unselect();*/
      }

    }
  }

  redraw();

END_OGL_STUFF
  return TRUE;
}

static gboolean luscus_gtk_mouse_button_release(GtkWidget *da, GdkEventButton *event, gpointer data)
{
  shift_move = 0;
/*  luscus_gtk_interac_with_gtk_part(myCountSelect, ngreyStek, selected_internal_coord_value);*/
  redraw();
  return TRUE;
}

static gboolean luscus_gtk_mouse_move(GtkWidget *da, GdkEventMotion *event, gpointer data)
{
  gdouble w, h;
  gdouble x;
  gdouble y;
  double x_ax, y_ax, z_ax;
  double angle;
  double dx, dy;
  GtkAllocation allocation;

BEGIN_OGL_STUFF;
  gtk_widget_get_allocation(drawingarea, &allocation);
  w = (double) allocation.width;
  h = (double) allocation.height;
  x = 2. * event->x / w - 1.;
  y = 2. * (h - event->y) / h - 1.;

  if (event->state & GDK_BUTTON1_MASK)
  {
    if (control_pressed)
    {
      dx=x-x_prev;
      dy=y-y_prev;

      if (2.0*fabs(dy) < fabs(dx)) 
      {
        expand_or_retract = dx / 10.0;
        do_expand_or_retract();
      }
    }
    else
    {
      x_ax = (y_prev - y);
      y_ax = (x - x_prev);
      z_ax = (x_prev * y - x * y_prev);

      angle = sqrt((x_ax*x_ax+y_ax*y_ax+z_ax*z_ax) / (x*x+y*y+1) * (x_prev*x_prev+y_prev*y_prev+1));
      angle*=(180.0e0/M_PI);

      /*cumulative angle*/
  
      tot_ang += fabs(angle);
      if (tot_ang > 45.0)
      {
        tot_ang = 0.;
        if (m->ngrids)
          luscus_gtk_resort_surfaces();
      }

      glMatrixMode(GL_MODELVIEW);
      glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
      glLoadIdentity();
      glRotated(angle,x_ax,y_ax,z_ax);
      glMultMatrixd(mvm);
    }
  }
  else if (event->state & GDK_BUTTON2_MASK)
  {
    dx=x-x_prev;
    dy=y-y_prev;

    if (2.0*fabs(dy) < fabs(dx)) 
    {
      if (dx > 0.0) scale*=1.02;
      else scale/=1.02;
      set_scale();
    }
    else
    {
/*      if (dy > 0.0) Do_key_plus();
      else Do_key_minus();*/
    }
  }
  else if (event->state & GDK_BUTTON3_MASK)
  {
    if (shift_move == 0)
    {
      shift_move_x_b=(int) event->x;
      shift_move_y_b=(int) event->y;
      shift_move=1;
      return TRUE;
    }
    if (shift_move==1)
    {
      int t00,t01,t11,t10;
      int ifound,i,j;
      GLint viewport[4];
      shift_move_x_e=(int) event->x;
      shift_move_y_e=(int) event->y;
      t00=MIN(shift_move_x_b, shift_move_x_e);
      t10=MAX(shift_move_x_b, shift_move_x_e);
      t01=MIN(shift_move_y_b, shift_move_y_e);
      t11=MAX(shift_move_y_b, shift_move_y_e);

      for(i = 0; i < m->natom; i++)
      {
        glPushMatrix();
        glRasterPos3d(m->xyz[i][0], m->xyz[i][1], m->xyz[i][2]);
        glGetIntegerv(GL_CURRENT_RASTER_POSITION,viewport);
        if(viewport[0]>t00 && viewport[0] < t10
           && h - viewport[1] > t01
           && h - viewport[1] < t11)
        {
          ifound=-1;
          for (j = 0;j < m->n_marked; j++)
            if (m->marked[j]==i) ifound=j;
  
          if (ifound==-1) mark_atom(i);
        }
      }
    }
  }

  redraw();
END_OGL_STUFF

  x_prev = x;
  y_prev = y;
  return TRUE;
}

void luscus_gtk_mouse_scroll(GtkWidget *widget, GdkEventScroll  *event, gpointer data)
{
  if (event->direction == GDK_SCROLL_UP) do_key_page(1,0);
  else if (event->direction == GDK_SCROLL_DOWN) do_key_page(0,0);
}

static gboolean luscus_draw_init(GtkWidget *da, gpointer data)
{
  GLfloat lpos1[4]={-.5,-.5,1,0};
  GLfloat defaultam[4] = { 0.1F,0.1F,0.1F,1.0F };
  GLfloat defaultdif[4] = { 0.7F,0.7F,0.7F,1.0F };
  GLfloat defaultspec[4] = { 0.2F,0.2F,0.2F,1.0F };
  float tmp[]={1.0, 1.0, 1.0};
  gint first = GPOINTER_TO_INT(data);
BEGIN_OGL_STUFF;
  scale = 1.0;
  expand_or_retract = 1.0;
  camera_x = camera_y = 0.0;
  list_3d = glGenLists(1);
  redraw_3d = 1;

  glLightfv(GL_LIGHT0,GL_POSITION,lpos1);
  glLightfv(GL_LIGHT0,GL_AMBIENT,defaultam);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,defaultdif);
  glLightfv(GL_LIGHT0,GL_SPECULAR,defaultspec);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);

  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, tmp);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, (float) 64); /*for atoms it is better 64 (more shiny)*/
  
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glDisable(GL_CULL_FACE);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
END_OGL_STUFF
  return TRUE;
}

#ifdef GTK3
static gboolean luscus_draw_reshape(GtkWidget *da, GdkEvent *event, gpointer data)
#endif
#ifdef GTK2
static gboolean luscus_draw_reshape(GtkWidget *da, GdkEventConfigure *event, gpointer data)
#endif
{
  gint w, h;
  GtkAllocation allocation;
BEGIN_OGL_STUFF;
  gtk_widget_get_allocation(da, &allocation);
  w = allocation.width;
  h = allocation.height;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  screen_width = (double) w;
  screen_height = (double) h;

  if (w < h)
  {
    x_size = molecule_diameter / (scale * 2.F);
    y_size = molecule_diameter * (double) h / ((double) w * scale * 2.F);
  }
  else
  {
    x_size = molecule_diameter * (double) w / ((double) h * scale * 2.F);
    y_size = molecule_diameter / (scale * 2.F);
  }

  glOrtho(-x_size - camera_x, x_size - camera_x,
          -y_size - camera_y, y_size - camera_y,
          -2.0*molecule_diameter, 2.0*molecule_diameter);

  glMatrixMode(GL_MODELVIEW);

  glViewport(0, 0, w, h);

END_OGL_STUFF
return TRUE;
}

#ifdef GTK3
static gboolean luscus_draw_display(GtkWidget *da, cairo_t * cr, gpointer data)
#endif
#ifdef GTK2
static gboolean luscus_draw_display(GtkWidget *da, GdkEventExpose *event, gpointer data)
#endif
{
  double bw;
  GLfloat lpos[4]={0.,0.,1.,0};
BEGIN_OGL_STUFF;
  glMatrixMode(GL_MODELVIEW);

  if (Input_Data.bw)
  {
    bw=grayscale(Input_Data.background_color[0], Input_Data.background_color[1], Input_Data.background_color[2]);
    glClearColor(bw, bw, bw, 0.F);
  }
  else
    glClearColor(Input_Data.background_color[0], Input_Data.background_color[1],
                 Input_Data.background_color[2], Input_Data.background_color[3]); 

  glClearDepth(1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (Input_Data.move_light)
  {
    glPushMatrix();
    glLightfv(GL_LIGHT0, GL_POSITION, lpos);
    glLoadMatrixd(keepmvm);
  }

  if (redraw_3d)
  {
    luscus_draw_3d();
    redraw_3d = 0;
  }
  else if(glIsList(list_3d)) glCallList(list_3d);

  luscus_draw_2d();

  if (Input_Data.move_light) glPopMatrix();

  glFlush();
SWAP_OGL_BUFFERS
END_OGL_STUFF
  return TRUE;
}

void Init_Gui(int argc, char **argv)
{
  GtkWidget *vbox;
  GtkWidget *hbox;
  GtkWidget *menubar;
  GtkWidget *notebook;
  GtkWidget *eventbox;
  GdkColor white;
  GtkWidget *window;
#ifdef GTK_GLEXT
  GdkGLConfig *glconfig;
#else
  XVisualInfo *xvisual;
  GdkScreen *screen;

  Window root;
  int xscreen;

#ifdef GTK2
  GdkColormap *colormap;
#endif

  int attributes[] ={GLX_RGBA, 
                     GLX_ALPHA_SIZE, 1,
                     GLX_RED_SIZE, 1, 
                     GLX_GREEN_SIZE, 1, 
                     GLX_BLUE_SIZE, 1, 
                     GLX_DOUBLEBUFFER, True, 
                     GLX_DEPTH_SIZE, 12, 
                     None};
#endif

  molecule_diameter = 5.F;

  gtk_disable_setlocale();
  gtk_init(&argc, &argv);
#ifdef GTK_GLEXT
  gtk_gl_init(&argc, &argv);
  glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH | GDK_GL_MODE_DOUBLE);
#else
  tmpglxvis = &glxvis;
#endif

  gdk_color_parse("white", &white);

  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  g_signal_connect(G_OBJECT(window), "delete_event", G_CALLBACK(Kill_Gui), NULL);
  pix = gdk_pixbuf_new_from_xpm_data(molcas015002_xpm);
  gtk_window_set_icon(GTK_WINDOW(window), pix);
  g_signal_connect(G_OBJECT(window), "key_press_event", G_CALLBACK(luscus_gtk_key_press), NULL);
  g_signal_connect(G_OBJECT(window), "key_release_event", G_CALLBACK(luscus_gtk_key_release), NULL);

#ifndef GTK_GLEXT
  glxvis.display = gdk_x11_get_default_xdisplay();
  xscreen = DefaultScreen(glxvis.display);
  screen = gdk_screen_get_default();
  xvisual = glXChooseVisual(glxvis.display, xscreen, attributes); /*xscreen => glxvis.id */
  glxvis.visual = gdk_x11_screen_lookup_visual(screen, xvisual->visualid); /*xscreen => glxvis.id */
  root = RootWindow(glxvis.display, xscreen);
  glxvis.xcolormap = XCreateColormap(glxvis.display, root, xvisual->visual, AllocNone);
  gtk_widget_set_visual(window, glxvis.visual);
  glxvis.context = glXCreateContext(glxvis.display, xvisual, NULL, TRUE);
  XFree(xvisual);
#endif

#ifdef GTK2
  vbox = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox), FALSE);
#endif
  gtk_container_add(GTK_CONTAINER(window), vbox);

  menubar = make_menubar(window);
  gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, FALSE, 0);
  gtk_widget_show(menubar);

  eventbox = gtk_event_box_new();
  gtk_box_pack_start(GTK_BOX(vbox), eventbox, FALSE, FALSE, 0);
  gtk_widget_modify_base(eventbox, GTK_STATE_NORMAL, &white);
  gtk_widget_modify_fg(eventbox, GTK_STATE_NORMAL, &white);
  gtk_widget_modify_bg(eventbox, GTK_STATE_NORMAL, &white);

  statusbar1 = gtk_statusbar_new();
  gtk_container_add(GTK_CONTAINER(eventbox), statusbar1);
  gtk_widget_show(statusbar1);
  gtk_widget_show(eventbox);

  eventbox = gtk_event_box_new();
  gtk_box_pack_start(GTK_BOX(vbox), eventbox, FALSE, FALSE, 0);
  gtk_widget_modify_base(eventbox, GTK_STATE_NORMAL, &white);
  gtk_widget_modify_fg(eventbox, GTK_STATE_NORMAL, &white);
  gtk_widget_modify_bg(eventbox, GTK_STATE_NORMAL, &white);

  statusbar2 = gtk_statusbar_new();
  gtk_container_add(GTK_CONTAINER(eventbox), statusbar2);
  gtk_widget_show(statusbar2);
  gtk_widget_show(eventbox);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

  drawingarea = gtk_drawing_area_new();
  gtk_widget_set_double_buffered(drawingarea, FALSE);
#ifdef GTK_GLEXT
  gtk_widget_set_gl_capability(drawingarea, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE);
#endif

//#ifdef GTK2
//  gtk_widget_set_colormap(drawingarea, colormap);
//#endif

  gtk_box_pack_start(GTK_BOX(hbox), drawingarea, TRUE, TRUE, 0);
  gtk_widget_set_size_request(GTK_WIDGET(drawingarea), Input_Data.init_screen_size, Input_Data.init_screen_size);
  gtk_widget_add_events(GTK_WIDGET(drawingarea), GDK_BUTTON1_MOTION_MASK |
                                     GDK_BUTTON2_MOTION_MASK |
                                     GDK_BUTTON3_MOTION_MASK |
                                     GDK_SCROLL_MASK |
                                     GDK_BUTTON_PRESS_MASK |
                                     GDK_BUTTON_RELEASE_MASK |
                                     GDK_EXPOSURE_MASK);
  g_signal_connect_after(G_OBJECT(drawingarea), "realize", G_CALLBACK(luscus_draw_init), NULL); 
  g_signal_connect(G_OBJECT(drawingarea), "configure_event", G_CALLBACK(luscus_draw_reshape), NULL);
#ifdef GTK3
  g_signal_connect(G_OBJECT(drawingarea), "draw", G_CALLBACK(luscus_draw_display), NULL);
#elif GTK2
  g_signal_connect(G_OBJECT(drawingarea), "expose_event", G_CALLBACK(luscus_draw_display), NULL); 
#endif
  g_signal_connect(G_OBJECT(drawingarea), "motion_notify_event", G_CALLBACK(luscus_gtk_mouse_move), NULL);
  g_signal_connect(G_OBJECT(drawingarea), "button_press_event", G_CALLBACK(luscus_gtk_mouse_button_press), NULL);
  g_signal_connect(G_OBJECT(drawingarea), "button_release_event", G_CALLBACK(luscus_gtk_mouse_button_release), NULL);
  g_signal_connect(G_OBJECT(drawingarea), "scroll-event", G_CALLBACK(luscus_gtk_mouse_scroll), NULL);

  gtk_widget_show(drawingarea);

  notebook = make_notebook();
  gtk_box_pack_start(GTK_BOX(hbox), notebook, FALSE, FALSE, 0);
  gtk_widget_show(notebook);

  gtk_widget_show(hbox);
  gtk_widget_show(vbox);
  gtk_widget_show(window);

  textbox_state = 0;

}

void Start_Gui(void)
{
  print_plugin_warnings();

  gtk_main();

  if (glIsList(list_3d)) glDeleteLists(list_3d, 1);
#ifndef GTK_GLEXT
  glXDestroyContext(tmpglxvis->display, tmpglxvis->context);
  XFreeColormap(tmpglxvis->display, tmpglxvis->xcolormap);
  g_object_unref(G_OBJECT(tmpglxvis->visual));
#endif
}

void Kill_Gui(void)
{
  remove_backup();
  gtk_main_quit();
}

