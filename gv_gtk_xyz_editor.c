/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdlib.h>
#include<string.h>
#include<gtk/gtk.h>
#include<gdk/gdkkeysyms.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_gtk.h"

GtkWidget *atom_treeview = NULL;

void luscus_destroy_xyz_edit_dialog(GtkWidget *dialog, gpointer data)
{
  gtk_widget_destroy(atom_treeview);
  gtk_widget_destroy(dialog);
  atom_treeview = NULL;
}

void callback_xyz_edit(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *content_area;
  GtkWidget *sw;

  GtkCellRenderer *renderer;
  GtkTreeViewColumn *column;

  dialog = gtk_dialog_new_with_buttons("xyz editor",
                                        GTK_WINDOW(window),
                                        GTK_DIALOG_DESTROY_WITH_PARENT,
                                        GTK_STOCK_OK,
                                        GTK_RESPONSE_NONE,
                                        NULL);

  gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 400);

  g_signal_connect(dialog, "response", G_CALLBACK(luscus_destroy_xyz_edit_dialog), NULL);
/*  g_signal_connect_swapped(dialog, "response",
                           G_CALLBACK(gtk_widget_destroy), dialog);*/

  content_area = gtk_dialog_get_content_area(GTK_DIALOG (dialog));

  /*expand!!!!!!*/

  sw = gtk_scrolled_window_new(NULL, NULL);
#ifdef GTK3
  gtk_widget_set_vexpand(sw, TRUE);
#endif
  gtk_container_add(GTK_CONTAINER(content_area), sw);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);


  atom_treeview = gtk_tree_view_new_with_model(GTK_TREE_MODEL(atom_list));
  g_signal_connect(G_OBJECT(atom_treeview), "row-activated", G_CALLBACK(luscus_select_atom), NULL);
  gtk_container_add(GTK_CONTAINER (sw), atom_treeview);

/*column number*/
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes("#", renderer, "text", COLUMN_NUMBER, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

/*column atom symbol*/
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes("Symbol", renderer, "text", COLUMN_ATNAME, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

/*column X coordinate*/
  renderer = gtk_cell_renderer_text_new();
  g_object_set(renderer, "editable", TRUE, NULL);
  g_signal_connect(renderer, "edited", G_CALLBACK(luscus_change_x_coordinate), NULL); 
  column = gtk_tree_view_column_new_with_attributes("X", renderer, "text", COLUMN_X, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

/*column Y coordinate*/
  renderer = gtk_cell_renderer_text_new();
  g_object_set(renderer, "editable", TRUE, NULL);
  g_signal_connect(renderer, "edited", G_CALLBACK(luscus_change_y_coordinate), NULL); 
  column = gtk_tree_view_column_new_with_attributes("Y", renderer, "text", COLUMN_Y, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

/*column Z coordinate*/
  renderer = gtk_cell_renderer_text_new();
  g_object_set(renderer, "editable", TRUE, NULL);
  g_signal_connect(renderer, "edited", G_CALLBACK(luscus_change_z_coordinate), NULL);
  column = gtk_tree_view_column_new_with_attributes("Z", renderer, "text", COLUMN_Z, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

/*column numeration*/
  renderer = gtk_cell_renderer_text_new();
  g_object_set(renderer, "editable", TRUE, NULL);
  g_signal_connect(renderer, "edited", G_CALLBACK(luscus_change_atom_numbering), NULL);
  column = gtk_tree_view_column_new_with_attributes("numeration", renderer, "text", COLUMN_NUMER, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

/*column atom symbol*/
  renderer = gtk_cell_renderer_text_new();
  g_object_set(renderer, "editable", TRUE, NULL);
  g_signal_connect(renderer, "edited", G_CALLBACK(luscus_change_atom_name), NULL);
  column = gtk_tree_view_column_new_with_attributes("Name", renderer, "text", COLUMN_NAME, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(atom_treeview), column);

  gtk_widget_show(atom_treeview);
  
  gtk_widget_show(sw);

  gtk_widget_show(dialog);
}

void luscus_select_atom(GtkTreeView *tree_view, GtkTreePath *path, GtkTreeViewColumn *column, gpointer data)
{
  gint iatom;
  iatom = gtk_tree_path_get_indices(path)[0];

  if (iatom > m->natom) return;
  unselect_all();
  m->n_selected++;
  m->selected[0] = iatom;
  redraw();
}

void luscus_change_x_coordinate(GtkCellRendererText* cell, const gchar *path_string, const gchar *new_text, gpointer data)
{
  GtkTreePath *path = gtk_tree_path_new_from_string(path_string);
  GtkTreeIter iter;
  gint iatom;

  iatom = gtk_tree_path_get_indices(path)[0];
  gtk_tree_model_get_iter(GTK_TREE_MODEL(atom_list), &iter, path);

  m->xyz[iatom][0] = g_strtod(new_text, NULL);
  gtk_list_store_set(atom_list, &iter, COLUMN_X, m->xyz[iatom][0], -1);
  rerender_3d();
  redraw();

  gtk_tree_path_free (path);
}

void luscus_change_y_coordinate(GtkCellRendererText* cell, const gchar *path_string, const gchar *new_text, gpointer data)
{
  GtkTreePath *path = gtk_tree_path_new_from_string (path_string);
  GtkTreeIter iter;
  gint iatom;

  iatom = gtk_tree_path_get_indices(path)[0];
  gtk_tree_model_get_iter (GTK_TREE_MODEL(atom_list), &iter, path);

  m->xyz[iatom][1] = g_strtod(new_text, NULL);
  gtk_list_store_set(atom_list, &iter, COLUMN_Y, m->xyz[iatom][1], -1);
  rerender_3d();
  redraw();

  gtk_tree_path_free (path);
}

void luscus_change_z_coordinate(GtkCellRendererText* cell, const gchar *path_string, const gchar *new_text, gpointer data)
{
  GtkTreePath *path = gtk_tree_path_new_from_string (path_string);
  GtkTreeIter iter;
  gint iatom;

  iatom = gtk_tree_path_get_indices(path)[0];
  gtk_tree_model_get_iter (GTK_TREE_MODEL(atom_list), &iter, path);

  m->xyz[iatom][2] = g_strtod(new_text, NULL);
  gtk_list_store_set(atom_list, &iter, COLUMN_Z, m->xyz[iatom][2], -1);
  rerender_3d();
  redraw();

  gtk_tree_path_free (path);
}

void luscus_change_atom_numbering(GtkCellRendererText* cell, const gchar *path_string, const gchar *new_text, gpointer data)
{
  GtkTreePath *path = gtk_tree_path_new_from_string (path_string);
  GtkTreeIter iter;
  gint iatom;

  iatom = gtk_tree_path_get_indices(path)[0];
  gtk_tree_model_get_iter (GTK_TREE_MODEL(atom_list), &iter, path);
  printf("Changing numbering of atom #%d = |%s|\n", iatom, new_text);

  if (m->additional_numeration)
  {
    m->additional_numeration[iatom] = atoi(new_text);
    gtk_list_store_set(atom_list, &iter, COLUMN_NUMER, m->additional_numeration[iatom], -1);
/*    redraw();*/
  }
  if (!(m->ishow & HAS_ATOMNUMS))
  {
    m->ishow ^= HAS_ATOMNUMS;
    luscus_gtk_menubar_show_or_hide_widgets();
  }

  gtk_tree_path_free (path);
}

void luscus_change_atom_name(GtkCellRendererText* cell, const gchar *path_string, const gchar *new_text, gpointer data)
{
  GtkTreePath *path = gtk_tree_path_new_from_string (path_string);
  GtkTreeIter iter;
  gint iatom;

  iatom = gtk_tree_path_get_indices(path)[0];
  gtk_tree_model_get_iter (GTK_TREE_MODEL(atom_list), &iter, path);

  if (m->name)
  {
    if (m->name[iatom]) free(m->name[iatom]);
    m->name[iatom] = strdup(new_text);
    gtk_list_store_set(atom_list, &iter, COLUMN_NAME, new_text, -1);
    /*redraw();*/
  }
  if (!(m->ishow & HAS_ATOMNAMES))
  {
    m->ishow ^= HAS_ATOMNAMES;
    luscus_gtk_menubar_show_or_hide_widgets();
  }

  gtk_tree_path_free (path);
}

void luscus_xyz_editor_select_row(int iat)
{
  GtkTreePath *path;
  GtkTreeSelection *selection;

  if (atom_treeview == NULL) return;

#ifdef EBUG
  gchar *ctmp;
  gdouble xtmp, ytmp, ztmp;
  GtkTreeIter iter;
  gboolean test;
#endif

  selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(atom_treeview));
  path = gtk_tree_path_new_from_indices(iat, -1);

#ifdef EBUG
  test = gtk_tree_model_get_iter(GTK_TREE_MODEL(atom_list), &iter, path);
  gtk_tree_model_get(GTK_TREE_MODEL(atom_list), &iter, COLUMN_ATNAME, &ctmp, COLUMN_X, &xtmp, COLUMN_Y, &ytmp, COLUMN_Z, &ztmp, -1);
  printf("selected: %s %f %f %f\n", ctmp, xtmp, ytmp, ztmp);
#endif

  gtk_tree_selection_select_path(selection, path);

}

