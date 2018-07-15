/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*handle atom lists*/
#include<gtk/gtk.h>
#include"luscus.h"
#include"gv.h"
#include"gv_gtk.h"

void init_atom_list(void)
{
  atom_list = gtk_list_store_new(N_COLUMNS, G_TYPE_UINT, G_TYPE_STRING, G_TYPE_DOUBLE, G_TYPE_DOUBLE, G_TYPE_DOUBLE, G_TYPE_UINT, G_TYPE_STRING);
}

void deallocate_atom_list(void)
{
  gboolean test;
  GtkTreeIter iter;
  GtkTreePath *path;

  if (atom_list)
  {
    path = gtk_tree_path_new_first();
    test = gtk_tree_model_get_iter(GTK_TREE_MODEL(atom_list), &iter, path);
    gtk_tree_path_free(path);

    if (test)
      while(gtk_list_store_remove(atom_list, &iter));
  }

#ifdef EBUG
  printf("ALL atoms removed from list!\n"); fflush(stdout);
#endif
}

void insert_atom_into_list(int iat)
{
  GtkTreeIter iter;
  gtk_list_store_append(atom_list, &iter);
  gtk_list_store_set(atom_list, &iter,
                     COLUMN_NUMBER, iat+1,
                     COLUMN_ATNAME, m->elem[iat].name,
                     COLUMN_X, m->xyz[iat][0],
                     COLUMN_Y, m->xyz[iat][1],
                     COLUMN_Z, m->xyz[iat][2],
                     COLUMN_NUMER, m->additional_numeration[iat],
                     COLUMN_NAME, m->name[iat],
                     -1);
#ifdef EBUG
  printf("Adding atom #%d to list!\n", iat); fflush(stdout);
#endif
}

void insert_all_atoms_into_list(void)
{
  int i;
  GtkTreeIter iter;

  for (i = 0; i < m->natom; i++)
  {
    gtk_list_store_append(atom_list, &iter);
    gtk_list_store_set(atom_list, &iter,
                       COLUMN_NUMBER, i+1,
                       COLUMN_ATNAME, m->elem[i].name,
                       COLUMN_X, m->xyz[i][0],
                       COLUMN_Y, m->xyz[i][1], 
                       COLUMN_Z, m->xyz[i][2], 
                       COLUMN_NUMER, m->additional_numeration[i],
                       COLUMN_NAME, m->name[i],
                       -1);
#ifdef EBUG
  printf("Adding atom #%d to list!\n", i); fflush(stdout);
#endif
  }
}

void remove_atom_from_list(int iat)
{
  int i;
/*  int iiat =iat;*/
  GtkTreeIter iter;
  GtkTreePath *path;
  gboolean test;
#ifdef EBUG
  gchar *tmp;
#endif

 /*find path from atom number*/
  path = gtk_tree_path_new_from_indices(iat, -1);
  test = gtk_tree_model_get_iter(GTK_TREE_MODEL(atom_list), &iter, path);

#ifdef SMECE
  /*move iter to the begining of the list*/
  gtk_tree_model_get_iter_first(GTK_TREE_MODEL(atom_list), &iter);
#ifdef EBUG
  if (! gtk_list_store_iter_is_valid(atom_list, &iter))
  {
    printf("ERROR: starting iter is invalid!\n");
  }
#endif
  /*move iter down the list iat times*/
  for(i = 0; i < iat; i++)
  {
#ifdef EBUG
    gtk_tree_model_get(GTK_TREE_MODEL(atom_list), &iter, COLUMN_ATNAME, &tmp, -1);
    printf("listing atom: |%s|\n", tmp);
#endif
    test = gtk_tree_model_iter_next(GTK_TREE_MODEL(atom_list), &iter);
    if (!test) break;
  }
#endif

  /*remove list element pointed by iter*/
  if (test)
    test = gtk_list_store_remove(atom_list, &iter);

#ifdef EBUG
  printf("Removing atom #%d from the list!\n", iat); fflush(stdout);
#endif

  if (test)
  {
    gtk_list_store_set(atom_list, &iter, COLUMN_NUMBER, (iat++)+1, -1);
    /*change indices of all atoms after iat*/
    while(gtk_tree_model_iter_next(GTK_TREE_MODEL(atom_list), &iter))
      gtk_list_store_set(atom_list, &iter, COLUMN_NUMBER, (iat++)+1, -1);

      /*make something with NUMER and NAME*/
  }

}

void change_atom_parameters_in_list(int iat)
{
  int i;
  GtkTreeIter iter;
  gboolean test = FALSE;

  /*move iter to the begining of the list*/
  gtk_tree_model_get_iter_first(GTK_TREE_MODEL(atom_list), &iter);

  for(i = 0; i < iat; i++)
  {
    test = gtk_tree_model_iter_next(GTK_TREE_MODEL(atom_list), &iter);
    if (!test) break;
  }

  if (test)
    gtk_list_store_set(atom_list, &iter,
                       COLUMN_NUMBER, iat+1,
                       COLUMN_ATNAME, m->elem[i].name,
                       COLUMN_X, m->xyz[i][0],
                       COLUMN_Y, m->xyz[i][1], 
                       COLUMN_Z, m->xyz[i][2], 
                       COLUMN_NUMER, m->additional_numeration[i],
                       COLUMN_NAME, m->name[i],
                       -1);
}

