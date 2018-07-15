/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdlib.h>
#include<string.h>
#include<gtk/gtk.h>

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
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_gtk.h"

int tgaSave(char *filename, int w, int h)
{
  unsigned char cGarbage = 0, type, mode, aux;
  short int iGarbage = 0;
  int i;
  int rc;
  unsigned char pixelDepth=32;
  FILE *file;
  int x, y;
  unsigned char *imageData;
#ifdef TRACE
    puts("Trace: tgaSave");
#endif
    printf("TGA SAVE!\n");

/*  x = GuiState.win_state.x;
  y = GuiState.win_state.y;*/
  imageData = (unsigned char *) malloc(sizeof(unsigned char) * w * h * 4);

  glReadPixels(0, 0, w, h, GL_RGBA,GL_UNSIGNED_BYTE, (GLvoid *) imageData);

  file=fopen(filename, "wb");
  if(file == NULL)
  {
    make_warning("Error: Can't create file");
  /*  sprintf(LOG,"Can't create file %s\n",filename);
    printlog(LOG);*/
    return 1;
  }
  mode = pixelDepth / 8;
  if ((pixelDepth == 24) || (pixelDepth == 32)) type = 2;
  else type = 3;

/* write the header  */

  rc=fwrite(&cGarbage, sizeof(unsigned char), 1, file);
  rc=fwrite(&cGarbage, sizeof(unsigned char), 1, file);

  rc=fwrite(&type, sizeof(unsigned char), 1, file);

  rc=fwrite(&iGarbage, sizeof(short int), 1, file);
  rc=fwrite(&iGarbage, sizeof(short int), 1, file);
  rc=fwrite(&cGarbage, sizeof(unsigned char), 1, file);
  rc=fwrite(&iGarbage, sizeof(short int), 1, file);
  rc=fwrite(&iGarbage, sizeof(short int), 1, file);

  rc=fwrite(&w, sizeof(short int), 1, file);
  rc=fwrite(&h, sizeof(short int), 1, file);
  rc=fwrite(&pixelDepth, sizeof(unsigned char), 1, file);

  rc=fwrite(&cGarbage, sizeof(unsigned char), 1, file);

  /* convert the image data from RGB(a) to BGR(A) */
  if (mode >= 3)
    for (i=0; i < w * h * mode ; i+= mode)
    {
      aux = imageData[i];
      imageData[i] = imageData[i+2];
      imageData[i+2] = aux;
    }

  rc=fwrite(imageData, sizeof(unsigned char),
           (unsigned int)(w * h * mode), file);
  fclose(file);
  free(imageData);
  luscus_gtk_pop_message_from_statusbar2();
  luscus_gtk_push_message_to_statusbar2("Bitmap saved");
  printf("TGA SAVED!\n");

  return 0;
}
#ifdef GTK_GLEXT
void pixbuf_save(char *filename, char *type, int width, int height)
{
  unsigned char *imageData;
  int rowstride;
  GdkPixbuf *pixbuf;
  GError *error = NULL;
  unsigned char *imageRow;
  int i;

  imageData = (unsigned char *) malloc(sizeof(unsigned char) * width * height * 4);
  imageRow = (unsigned char *) malloc(sizeof(unsigned char) * width * 4);

  glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid *) imageData);

  for (i = 0; i < height/2; i++)
  {
    memcpy(imageRow, (imageData+i*width*4), width*4);
    memcpy((imageData+i*width*4), (imageData+(height-i-1)*width*4), width*4);
    memcpy((imageData+(height-i-1)*width*4), imageRow, width*4);
  }
  pixbuf = gdk_pixbuf_new_from_data(imageData, GDK_COLORSPACE_RGB, TRUE, 8, width, height, width*4, NULL, NULL);

  if (pixbuf) gdk_pixbuf_save(pixbuf, filename, type, &error, NULL);
  g_object_unref(G_OBJECT(pixbuf));

  free(imageData);
}

static gboolean example_screen_init(GtkWidget *widget, gpointer data)
{
  GLfloat lpos1[4] = {-0.5, -0.5, 1,0};
  GLfloat defaultam[4] = {0.1F, 0.1F, 0.1F, 1.0F};
  GLfloat defaultdif[4] = {0.7F, 0.7F, 0.7F, 1.0F};
  GLfloat defaultspec[4] = {0.2F, 0.2F, 0.2F, 1.0F};
  float tmp[] = {1.0, 1.0, 1.0};
  gint w, h;
  GtkAllocation allocation;
  double x_size, y_size;
  double scale;
  double camera_x, camera_y;
  double molecule_diameter;

  GLdouble *mvm_main;

  GdkWindow *window = gtk_widget_get_window(widget); \
  GdkGLContext *glcontext = gtk_widget_get_gl_context(widget); \
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget); \
  if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) return FALSE;

  gtk_widget_get_allocation(widget, &allocation);
  w = allocation.width;
  h = allocation.height;

  glLightfv(GL_LIGHT0,GL_POSITION,lpos1);
  glLightfv(GL_LIGHT0,GL_AMBIENT,defaultam);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,defaultdif);
  glLightfv(GL_LIGHT0,GL_SPECULAR,defaultspec);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);

  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, tmp);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, (float) 64); /*for atoms itis better 64 (more shiny)*/
  
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glDisable(GL_CULL_FACE);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  get_sizes_data(&x_size, &y_size, &scale, &camera_x, &camera_y, &molecule_diameter);
  glOrtho(-x_size - camera_x, x_size - camera_x,
          -y_size - camera_y, y_size - camera_y,
          -2.0*molecule_diameter, 2.0*molecule_diameter);

  glMatrixMode(GL_MODELVIEW);
  glViewport(0, 0, w, h);

  mvm_main = get_main_mvm();
  glLoadMatrixd(mvm_main);

  glClearColor(Input_Data.background_color[0], Input_Data.background_color[1], Input_Data.background_color[2], Input_Data.background_color[3]);
  glClearDepth(1.0);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gdk_gl_drawable_gl_end(gldrawable);

  return TRUE;
}

#ifdef GTK2
static gboolean example_screen_draw(GtkWidget* widget, GdkEventExpose* event, gpointer data)
#endif
#ifdef GTK3
static gboolean example_screen_draw(GtkWidget* widget, cairo_t* cr, gpointer data)
#endif
{
  double bw;
  GLuint list_3d;
  gdouble wid, hei;
#ifndef GTK_OLD
  GtkAllocation allocation;
#endif

  GdkWindow *window = gtk_widget_get_window(widget);
  GdkGLContext *glcontext = gtk_widget_get_gl_context(widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);
#ifdef GTK_OLD
  wid = widget->allocation.width;
  hei = widget->allocation.height;
#else
  gtk_widget_get_allocation(widget, &allocation);

  wid = allocation.width;
  hei = allocation.height;
#endif

  if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) return FALSE;

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

  list_3d = (GLuint) get_list_3d();

  if(m)
  {
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

    if (m->ntextboxes) draw_textboxes(wid, hei);
  }

  glFlush();
  gdk_gl_drawable_gl_end(gldrawable);

  return TRUE;
}

void luscus_save_builtin_graphic_format(gchar* filename, gchar* type, gint width, gint height)
{
  GtkWidget *dialog;
  GtkWidget *sw, *ca;
  GtkWidget *daw;
  GdkPixbuf *pixbuf;
  GdkWindow *wind;
  GtkAllocation allocation;
  gint result;

/*  GError *error = NULL;*/

  GdkGLConfig *glconfig;

  dialog = gtk_dialog_new_with_buttons("Screenshot", NULL, GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                                           GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                           GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                           NULL);
  sw = gtk_scrolled_window_new(NULL, NULL);
  gtk_widget_set_size_request(GTK_WIDGET(sw), 500, 500);
  ca = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
  gtk_box_pack_start(GTK_BOX(ca), sw, TRUE, TRUE, 0);

  glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH | GDK_GL_MODE_SINGLE);

  daw = gtk_drawing_area_new();
  gtk_widget_set_size_request(GTK_WIDGET(daw), width, height);
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(sw), daw);
  gtk_widget_set_double_buffered(daw, FALSE);

  gtk_widget_set_gl_capability(daw, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE);

  g_signal_connect_after(G_OBJECT(daw), "realize", G_CALLBACK(example_screen_init), NULL);
#ifdef GTK2
  g_signal_connect(G_OBJECT(daw), "expose_event", G_CALLBACK(example_screen_draw), NULL);
#endif
#ifdef GTK3
  g_signal_connect(G_OBJECT(daw), "draw", G_CALLBACK(example_screen_draw), NULL);
#endif

  gtk_widget_show(daw);
  gtk_widget_show(sw);

  result = gtk_dialog_run(GTK_DIALOG(dialog));

  if (result == GTK_RESPONSE_ACCEPT)
  {
/*save image to file*/
#ifdef EBUG
    printf("SAVING SCREENSHOT\n");
#endif
    if (strcmp(type, "tga") == 0) tgaSave(filename, width, height);
    else pixbuf_save(filename, type, width, height);

/*    {
      wind = gtk_widget_get_window(daw);
      gtk_widget_get_allocation(daw, &allocation);
      printf("allocation = %d %d\n", allocation.width, allocation.height);
      gdk_window_invalidate_rect(wind, &allocation, FALSE);
      gdk_window_process_updates(wind, TRUE);
#ifdef GTK2
      pixbuf = gdk_pixbuf_get_from_drawable(NULL, wind, gdk_colormap_get_system(), 0, 0, 0, 0, width, height);
#endif
#ifdef GTK3
      pixbuf = gdk_pixbuf_get_from_window(wind, 0, 0, width, height);
#endif
      if (pixbuf) gdk_pixbuf_save(pixbuf, filename, type, &error, NULL);
      g_object_unref(G_OBJECT(pixbuf));
    }*/
  }

  gtk_widget_destroy(dialog);
}
#else /*non-gtk-glext*/
void luscus_save_builtin_graphic_format(gchar* filename, gchar* type, GdkPixbuf *pixbuf, gint width, gint height)
{
  GError *error = NULL;

  if (strcmp(type, "tga") == 0) tgaSave(filename, width, height);
  else if (pixbuf) gdk_pixbuf_save(pixbuf, filename, type, &error, NULL);
}
#endif

