/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<gtk/gtk.h>
#ifdef GTK_GLEXT
#include<gtk/gtkgl.h>
#endif
#include<math.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_gtk.h"

#define GEOM_MODE 1
#define VIB_MODE 2

/*MAKE SURE THAT LINE DATA IS DEALLOCATED!!!*/

typedef struct line
{
  double x_b; /*begin*/
  double y_b;
  double x_e; /*end*/
  double y_e;
} LINE;

typedef struct struct_tics
{
  double tics_min;
  double tics_max;
  double tics_diff;
} S_TICS;

gboolean graph_exist = FALSE;
LINE g_scale;
LINE *g_lines;
int nlines;
static int *freq_n =0;
static int graph_mode;

void luscus_gtk_make_graph(GtkWidget*, char*);
#ifdef GTK_OLD
void graph_expose_event(GtkWidget*, gpointer);
#elif GTK2
void graph_expose_event(GtkWidget*, gpointer);
#elif GTK3
gboolean graph_expose_event(GtkWidget*, cairo_t*, gpointer);
#endif
/*int freq_count(void);*/
int freq_index(int);
int current_data = 0;
GtkWidget *drawarea;

void set_current_graph_data(void)
{
#ifndef GTK_OLD
  GdkWindow *window;
  GtkAllocation allocation;
#endif
  if (graph_mode == GEOM_MODE) current_data = igeo;
  else if (graph_mode == VIB_MODE) current_data = ivib;
  else return;
  if (graph_exist)
  {
#ifdef GTK_OLD
    gdk_window_invalidate_rect(drawarea->window, &drawarea->allocation, FALSE);
#else
    window = gtk_widget_get_window(drawarea);
    gtk_widget_get_allocation(drawarea, &allocation);
    gdk_window_invalidate_rect(window, &allocation, FALSE);
#endif
  }
}

S_TICS find_tics(double *n_min, double *n_max)
{
  int i;
  double num_min, num_max;
  S_TICS n_tics;

  num_min = *n_min;
  num_max = *n_max;

  /*find order of magnitude*/
  for(i = 6; i > -10 && (num_max - num_min) < pow(10.F,(double) i); i--);

  n_tics.tics_diff = pow(10.F,(double) i);
  if (num_min > 0.)
  {
    for(n_tics.tics_min = pow(10.F,(double) i) * floor(num_min/pow(10.F,(double) i)) - 10.F * n_tics.tics_diff; n_tics.tics_min < num_min; n_tics.tics_min += n_tics.tics_diff);
    n_tics.tics_min-=n_tics.tics_diff;
  }
  else
    for(n_tics.tics_min = pow(10.F,(double) i) * floor(num_min/pow(10.F,(double) i)) + 10.F * n_tics.tics_diff; n_tics.tics_min > num_min; n_tics.tics_min -= n_tics.tics_diff);

  if (num_max > 0.)
  {
    for(n_tics.tics_max = pow(10.F,(double) i) * floor(num_max/pow(10.F,(double) i)) + n_tics.tics_diff * 10.F; n_tics.tics_max >= num_max; n_tics.tics_max -= n_tics.tics_diff);
    n_tics.tics_max += n_tics.tics_diff;
  }
  else
    for(n_tics.tics_max = pow(10.F,(double) i) * floor(num_max/pow(10.F,(double) i)) -n_tics.tics_diff * 10.F; n_tics.tics_max < num_max; n_tics.tics_max += n_tics.tics_diff);

  return n_tics;
}

/**
@return Count of current normal modes.
*/
/*int freq_count()
{
  int i, count = 0, result = 0;

  for (i=0;  i< m->nvibr;  ++i)
    if (m->ir_intensity[i] > 0.001)
      ++result;

  g_free(freq_n);
#ifdef TRACE_MALLOC
  puts("allocate freq_n");
#endif
  freq_n = (int*) g_malloc (result * sizeof(int));

  for (i=0;  i < m->nvibr;  ++i)
    if (m->ir_intensity[i] > 0.001)
    {
      freq_n[count] = i;
      ++count;
    }

  return result;
}*/

int freq_index (int i)
{
  return freq_n [i];
}

void kill_graph_dialog(GtkWidget *dialog, gpointer data)
{
  graph_exist = FALSE;
  gtk_widget_destroy(dialog);
  g_free(g_lines);
  graph_mode = 0;
}

void luscus_gtk_show_graphical_energies(GtkWidget *button, gpointer data)
{
  int i;
  /*prepare data for writing a graph*/

  graph_mode = GEOM_MODE;
  nlines = n_geometries;
  current_data = igeo;
  g_lines = (LINE*) g_malloc(sizeof(LINE) * nlines);
  g_scale.x_b = 0.F;
  g_scale.y_b = 1.0E+6;
  g_scale.x_e = (double) n_geometries;
  g_scale.y_e = -1.0E+6;

  for(i = 0; i < n_geometries; i++)
  {
    if (fabs(mol[i+1].geo_energy) > 1.0e-7)
    {
      if (mol[i].geo_energy < g_scale.y_b) g_scale.y_b = mol[i].geo_energy;
      if (mol[i].geo_energy > g_scale.y_e) g_scale.y_e = mol[i].geo_energy;
    }

    if (i < n_geometries - 1)
    {
      g_lines[i].x_b = (double) i;
      g_lines[i].y_b = mol[i].geo_energy;
      if (fabs(mol[i+1].geo_energy) < 1.0e-7)
      {
        g_lines[i].y_e = mol[i].geo_energy;
        g_lines[i].x_e = (double) i;
      }
      else
      {
        g_lines[i].y_e = mol[i+1].geo_energy;
        g_lines[i].x_e = (double) i+1;
      }
    }
   /* else if (fabs(mol[i].geo_energy) < 1.0e-7)*/ /*this will ensure that if the datapoint is missing, there will be no 0.0 E 2014.01.27*/
/*    {
      g_lines[i].x_b = (double) i;
      g_lines[i].x_e = (double) i;
      g_lines[i].y_b = mol[i].geo_energy;
      g_lines[i].y_e = mol[i].geo_energy;
    }*/
    else
    {
      g_lines[i].x_b = (double) i;
      g_lines[i].x_e = (double) i;
      g_lines[i].y_b = mol[i].geo_energy;
      g_lines[i].y_e = mol[i].geo_energy;
    }
  }

  luscus_gtk_make_graph(button, "Energies");
}

void luscus_gtk_show_vib_spectrum(GtkWidget *button, gpointer data)
{
  int i;
  /*prepare data for writing a graph*/
  
  graph_mode = VIB_MODE;
  current_data = ivib;
  g_lines = (LINE*) g_malloc(sizeof(LINE) * m->nvibr);
  g_scale.x_b = 1.0E+6;
  g_scale.y_b = 1.0E+6;
  g_scale.x_e = -1.0E+6;
  g_scale.y_e = -1.0E+6;

  nlines = m->nvibr;
  for(i = 0; i < m->nvibr; i++)
  {
    if (m->ir_intensity[i] < g_scale.y_b) g_scale.y_b = m->ir_intensity[i];
    if (m->ir_intensity[i] > g_scale.y_e) g_scale.y_e = m->ir_intensity[i];
    if (m->freq[i] < g_scale.x_b) g_scale.x_b = m->freq[i];
    if (m->freq[i] > g_scale.x_e) g_scale.x_e = m->freq[i];
    g_lines[i].x_b = g_lines[i].x_e = m->freq[i];
    g_lines[i].y_b = 0.F;
    g_lines[i].y_e = m->ir_intensity[i];
  }

  /*make some space before min. value and max. value*/
  g_scale.x_e += 50.;

  luscus_gtk_make_graph(button, "Vibrational spectrum");

}

#ifdef GTK2
void graph_expose_event(GtkWidget *drawingarea, gpointer data)
#endif
#ifdef GTK3
gboolean graph_expose_event(GtkWidget *drawingarea, cairo_t *cr, gpointer data)
#endif
{
  int i;
/*  cairo_surface_t *coord, *lines;*/
  gchar *tmp;

  double tmpd;
#ifdef GTK_OLD
  double w = (double) drawingarea->allocation.width;
  double h = (double) drawingarea->allocation.height;
#else
  double w, h;
  GdkWindow *window;
  GtkAllocation allocation;
#endif
#ifndef GTK3
  cairo_t *cr;
#endif

  double pad_x = 35., pad_y = 15.; /*x and y padding used for coordinate axis*/
  double x, y; /*local x and y coordinates*/

  S_TICS x_tics;
  S_TICS y_tics;

  x_tics = find_tics(&g_scale.x_b, &g_scale.x_e);
  y_tics = find_tics(&g_scale.y_b, &g_scale.y_e);

#ifdef GTK_OLD
  cr = gdk_cairo_create(drawingarea->window);
#elif GTK2
  window = gtk_widget_get_window(drawingarea);
  gtk_widget_get_allocation(drawarea, &allocation);
  w = allocation.width;
  h = allocation.height;
  cr = gdk_cairo_create(window);
#elif GTK3
  window = gtk_widget_get_window(drawingarea);
  gtk_widget_get_allocation(drawarea, &allocation);
  w = allocation.width;
  h = allocation.height;
#endif

  cairo_set_source_rgb(cr, 1.F, 1.F, 1.F);
  cairo_paint(cr);

  cairo_set_line_width (cr, 1.0);
  cairo_set_line_cap (cr, CAIRO_LINE_CAP_SQUARE);
  cairo_set_source_rgb(cr, 0, 0, 0);

  cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 12.0);

  /*draw axes*/

  cairo_move_to(cr, pad_x, (double) h - pad_y);
  cairo_line_to(cr, (double) w, (double) h - pad_y); /*x axis*/
  cairo_move_to(cr, pad_x, (double) h - pad_y);
  cairo_line_to(cr, pad_x, 0.); /*y axis*/

  /*draw x tics*/
  for (tmpd = x_tics.tics_min; tmpd < x_tics.tics_max + 0.00001F; tmpd += x_tics.tics_diff)
  {
    x = pad_x + (tmpd - x_tics.tics_min) * (w - 2*pad_x) / (x_tics.tics_max - x_tics.tics_min);
/*     	g_scale.x_b) * (w - pad_x) / (g_scale.x_e - g_scale.x_b); */
    y = h - pad_y;
    cairo_move_to(cr, x, y);
    y = h - pad_y + 5;
    cairo_line_to(cr, x, y);
  }

  tmp = g_strdup_printf("%.1f", x_tics.tics_min);
  cairo_move_to(cr, pad_x - 4, h - 2.0);
  cairo_show_text(cr, tmp);
  g_free(tmp);

  tmp = g_strdup_printf("%.1f", x_tics.tics_max);
  cairo_move_to(cr, x - 15.0, h - 2.0);
  cairo_show_text(cr, tmp);
  g_free(tmp);

  /*draw y tics*/

  for(tmpd = y_tics.tics_min; tmpd < y_tics.tics_max + 0.00001F; tmpd += y_tics.tics_diff)
  {
    y = pad_y + (tmpd - y_tics.tics_min) * (h - 2*pad_y) / (y_tics.tics_max - y_tics.tics_min);
/*    y = pad_y + (tmpd - g_scale.x_b) * (h - pad_y) / (g_scale.y_e - g_scale.y_b);*/
    x = pad_x;
    cairo_move_to(cr, x, h - y);
    x = pad_x - 5;
    cairo_line_to(cr, x, h - y);
  }

  tmp = g_strdup_printf("%.1f", y_tics.tics_min);
  cairo_move_to(cr, 2.0, h - pad_y + 4.);
  cairo_show_text(cr, tmp);
  g_free(tmp);

  tmp = g_strdup_printf("%.1f", y_tics.tics_max);
  cairo_move_to(cr, 2.0, h - y + 4.);
  cairo_show_text(cr, tmp);
  g_free(tmp);

  cairo_stroke(cr);
/*  cairo_stroke_preserve(cr);
  cairo_new_sub_path(cr);*/

  
  /*draw lines*/
  cairo_set_source_rgb(cr, 0.F, 0.F, 0.95F);
  for (i = 0; i < nlines; i++)
  {
    x = pad_x + (g_lines[i].x_b - x_tics.tics_min) * (w - 2*pad_x) / (x_tics.tics_max - x_tics.tics_min);
    y = pad_y + (g_lines[i].y_b - y_tics.tics_min) * (h - 2*pad_y) / (y_tics.tics_max - y_tics.tics_min);
/*    cairo_move_to (cr, g_lines[i].x_b, g_lines[i].y_b);
    cairo_line_to (cr, g_lines[i].x_e, g_lines[i].y_e);*/
    cairo_move_to(cr, x, h - y);

    x = pad_x + (g_lines[i].x_e - x_tics.tics_min) * (w - 2*pad_x) / (x_tics.tics_max - x_tics.tics_min);
    y = pad_y + (g_lines[i].y_e - y_tics.tics_min) * (h - 2*pad_y) / (y_tics.tics_max - y_tics.tics_min);

    cairo_line_to(cr, x, h - y);
    cairo_stroke(cr);
  }

  /**/
  cairo_new_sub_path(cr);
  cairo_set_source_rgb(cr, 1.F, 0.F, 0.F);
  x = pad_x + (g_lines[current_data].x_b - x_tics.tics_min) * (w - 2*pad_x) / (x_tics.tics_max - x_tics.tics_min);
  y = pad_y + (g_lines[current_data].y_b - y_tics.tics_min) * (h - 2*pad_y) / (y_tics.tics_max - y_tics.tics_min);
  cairo_arc(cr, x, h-y, 3.0, 0.0, 2.0*M_PI);
  /**/
  cairo_stroke(cr);

#ifdef OLD_GTK
  cairo_destroy(cr);
#endif
#ifdef GTK2
  cairo_destroy(cr);
#endif
#ifdef GTK3
  return FALSE;
#endif
}

void luscus_gtk_make_graph(GtkWidget *widget, char *title)
{
  GtkWidget *graph_dialog;

  gint run;

  if (graph_exist) return; /*only one graph should be active; otherwise errors could occur*/
  graph_exist = TRUE;
  graph_dialog = gtk_dialog_new_with_buttons(title, GTK_WINDOW(gtk_widget_get_toplevel(widget)), GTK_DIALOG_DESTROY_WITH_PARENT, GTK_STOCK_OK, GTK_RESPONSE_ACCEPT, NULL);

  g_signal_connect(G_OBJECT(graph_dialog), "response", G_CALLBACK(kill_graph_dialog), NULL);

/*  g_signal_connect_swapped(graph_dialog, "response", G_CALLBACK(gtk_widget_destroy), graph_dialog);*/

  drawarea = gtk_drawing_area_new();
  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(graph_dialog))), drawarea, TRUE, TRUE, 0);
  gtk_widget_set_size_request(GTK_WIDGET(drawarea), 420, 250);

#ifdef GTK2
  g_signal_connect(G_OBJECT(drawarea), "expose_event", G_CALLBACK(graph_expose_event), NULL);
#endif
#ifdef GTK3
  g_signal_connect(G_OBJECT(drawarea), "draw", G_CALLBACK(graph_expose_event), NULL);
#endif

  gtk_widget_show(drawarea);

  gtk_widget_show(graph_dialog);

}

