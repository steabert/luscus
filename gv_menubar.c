/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<string.h>
#include<gtk/gtk.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_gtk.h"

GtkWidget *menu_item_save;
GtkWidget *menu_undo;
GtkWidget *menu_item_xyz_edit;
GtkWidget *menu_editable;
GtkWidget *menu_label_numeration;
GtkWidget *menu_label_names;
GtkWidget *menu_label_muliken;
GtkWidget *menu_label_loprop;
GtkWidget *menu_automatic_bonding;
GtkWidget *menu_optimize_mechanics;
GtkWidget *menu_optimize_semiempirical;
GtkWidget *menu_optimize_ab_initio;
GtkWidget *menu_orbital_calculation;
int n_calc_defs;
DEF_CALC_T *calc_defs;

GtkWidget *make_menubar(GtkWidget *window)
{
  int i;
  int n_calc_menus;
  GtkWidget *menubar;
  GtkWidget *menu;
  GtkWidget *submenu;
  GtkWidget *menu_item;
  GtkAccelGroup *ag;

  menubar = gtk_menu_bar_new();
  menu = gtk_menu_new();

  ag = gtk_accel_group_new();
  gtk_window_add_accel_group(GTK_WINDOW(window), ag);
  
   /*File*/
  /*-Open-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_OPEN, ag);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_open_file), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Save-*/
  menu_item_save = gtk_image_menu_item_new_from_stock(GTK_STOCK_SAVE, ag);
  g_signal_connect(G_OBJECT(menu_item_save), "activate", G_CALLBACK(callback_save_file), (gpointer) 0);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item_save);
  gtk_widget_show(menu_item_save);

  /*-Save as-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SAVE_AS, ag);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_save_file), (gpointer) 1);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Close-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_CLOSE, ag);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_close_file), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Quit-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_QUIT, ag);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(Kill_Gui), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_menu_item_new_with_mnemonic("_File");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

  /*Edit*/

  menu = gtk_menu_new();
  /*Undo*/

  menu_undo = gtk_image_menu_item_new_from_stock(GTK_STOCK_UNDO, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_undo), "Undo");
  g_signal_connect(G_OBJECT(menu_undo), "activate", G_CALLBACK(callback_do_undo), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_undo);
  gtk_widget_show(menu_undo);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item_xyz_edit = gtk_menu_item_new_with_label("xyz editor");
  g_signal_connect(G_OBJECT(menu_item_xyz_edit), "activate", G_CALLBACK(callback_xyz_edit), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item_xyz_edit);
  gtk_widget_show(menu_item_xyz_edit);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_editable = gtk_check_menu_item_new_with_label("Editable");
  if (m)
  {
#ifdef EBUG
    printf("M defined!\n");
#endif
    if (m->editable) gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_editable), TRUE);
    else gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_editable), FALSE);
  }
  else
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_editable), TRUE);
  g_signal_connect(G_OBJECT(menu_editable), "activate", G_CALLBACK(callback_make_editable), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_editable);
  gtk_widget_show(menu_editable);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Background color-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Background Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_background), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-label color-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Label Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_label), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-orbital color-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Orbital negative Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_neg_orbital), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Orbital positive Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_pos_orbital), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-electrostatic surface color-*/
  
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Electrostatic negative Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_neg_epot), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Electrostatic positive Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_pos_epot), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);
  
  /*-plain color-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_COLOR, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Draw Color");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_plain), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SELECT_FONT, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Change font");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_change_font), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-atom properties-*/
  menu_item = gtk_menu_item_new_with_label("Atom properties");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_adjust_atom_properties), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-automatic bonding-*/
  menu_automatic_bonding = gtk_check_menu_item_new_with_label("Determine bonding");
  g_signal_connect(G_OBJECT(menu_automatic_bonding), "activate", G_CALLBACK(callback_change_automatic_bonding), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_automatic_bonding);
  gtk_widget_show(menu_automatic_bonding);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Save Settings-*/
/*  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_SAVE, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Save settings");*/
  menu_item = gtk_menu_item_new_with_label("Save settings");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_save_settings), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_menu_item_new_with_mnemonic("_Edit");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

    /*View*/

  menu = gtk_menu_new();

  /*-move light source-*/
  menu_item = gtk_menu_item_new_with_label("Move light");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_move_light), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Create viewpoint-*/
  menu_item = gtk_menu_item_new_with_label("Create viewpoint");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_create_viewpoint), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-Restore viewpoint-*/
  menu_item = gtk_menu_item_new_with_label("Restore viewpoint");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_restore_viewpoint), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

     /*-Color/grayscale-*/
  menu_item = gtk_check_menu_item_new_with_label("Grayscale");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), FALSE);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_grayscale), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-animate-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_MEDIA_PLAY, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Animate");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_animate), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_check_menu_item_new_with_label("Show atoms");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_show_atoms), NULL);
  gtk_widget_show(menu_item);

  menu_item = gtk_check_menu_item_new_with_label("Show bonds");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_show_bonds), NULL);
  gtk_widget_show(menu_item);

  menu_item = gtk_check_menu_item_new_with_label("Show axes");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), FALSE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_show_axes), NULL);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_check_menu_item_new_with_label("Dot/line");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_callback_show_dot_line), NULL);
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), FALSE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_FULLSCREEN, ag);
  gtk_menu_item_set_label(GTK_MENU_ITEM(menu_item), "Maximize");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_maximize), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_menu_item_new_with_mnemonic("_View");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

  menu = gtk_menu_new();

  /*calculation*/

  /*-MM-*/
  menu_optimize_mechanics = gtk_menu_item_new_with_label("molecular mechanics");
/*  g_signal_connect(G_OBJECT(menu_optimize_mechanics), "activate", G_CALLBACK(callback_move_light), NULL);*/
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_optimize_mechanics);
  gtk_widget_show(menu_optimize_mechanics);

  submenu = gtk_menu_new();
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_optimize_mechanics), submenu);

  for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
  {
    if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"mm"))
    {
      menu_item = gtk_menu_item_new_with_label(calc_defs[i].plugin_name);
      gtk_menu_shell_append(GTK_MENU_SHELL(submenu), menu_item);
      g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_start_calculation), GINT_TO_POINTER(i));
      n_calc_menus++;
      gtk_widget_show(menu_item);
    }
  }
  if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_mechanics), FALSE);

  /*-semiempirical-*/
  menu_optimize_semiempirical = gtk_menu_item_new_with_label("semiempirecal");
/*  g_signal_connect(G_OBJECT(menu_optimize_semiempirical), "activate", G_CALLBACK(callback_move_light), NULL);*/
  submenu = gtk_menu_new();
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_optimize_semiempirical), submenu);

  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_optimize_semiempirical);
  gtk_widget_show(menu_optimize_semiempirical);

  for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
  {
    if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"semiem"))
    {
      menu_item = gtk_menu_item_new_with_label(calc_defs[i].plugin_name);
      gtk_menu_shell_append(GTK_MENU_SHELL(submenu), menu_item);
      g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_start_calculation), GINT_TO_POINTER(i));
      n_calc_menus++;
      gtk_widget_show(menu_item);
    }
  }
  if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_semiempirical), FALSE);

  /*-ab-initio-*/
  menu_optimize_ab_initio = gtk_menu_item_new_with_label("Ab initio");
/*  g_signal_connect(G_OBJECT(menu_optimize_ab_initio), "activate", G_CALLBACK(menu_optimize_ab_initio), NULL);*/
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_optimize_ab_initio);
  gtk_widget_show(menu_optimize_ab_initio);

  submenu = gtk_menu_new();
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_optimize_ab_initio), submenu);

  for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
  {
    if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"abinitio"))
    {
      menu_item = gtk_menu_item_new_with_label(calc_defs[i].plugin_name);
      gtk_menu_shell_append(GTK_MENU_SHELL(submenu), menu_item);
      g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_start_calculation), GINT_TO_POINTER(i));
      n_calc_menus++;
      gtk_widget_show(menu_item);
    }
  }
  if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_ab_initio), FALSE);

  /*- -*/

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*- orbital calculation -*/
  menu_orbital_calculation = gtk_menu_item_new_with_label("Orbital calculations");
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_orbital_calculation);
  gtk_widget_show(menu_orbital_calculation);

  submenu = gtk_menu_new();
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_orbital_calculation), submenu);

  for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
  {
#ifdef EBUG
    printf("i = %d type = |%s| name = |%s|\n", i, calc_defs[i].type, calc_defs[i].plugin_name);
#endif
    if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type, "orbital"))
    {
#ifdef EBUG
      printf("plugin name = |%s|\n", calc_defs[i].plugin_name);
#endif
      menu_item = gtk_menu_item_new_with_label(calc_defs[i].description);
      gtk_menu_shell_append(GTK_MENU_SHELL(submenu), menu_item);
      g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_start_calculation), GINT_TO_POINTER(i));
      n_calc_menus++;
      gtk_widget_show(menu_item);
    }
  }
  if (m)
  {
    if (m->ngrids)
    {
      if (n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), TRUE);
      else gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), FALSE);
    }
    else
      gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), FALSE);
  }
  else
    gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), FALSE);

  menu_item = gtk_menu_item_new_with_mnemonic("_Calculation");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

  menu = gtk_menu_new();

  /*labels*/

  menu_item = gtk_check_menu_item_new_with_label("Show indices");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), FALSE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_show_indices), NULL);
  gtk_widget_show(menu_item);

  menu_label_numeration = gtk_check_menu_item_new_with_label("Show numeration");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_label_numeration), FALSE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_label_numeration);
  g_signal_connect(G_OBJECT(menu_label_numeration), "activate", G_CALLBACK(luscus_gtk_show_numeration), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(menu_label_numeration), FALSE);
  gtk_widget_show(menu_label_numeration);

  menu_item = gtk_check_menu_item_new_with_label("Show atom symbols");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), FALSE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_show_symbols), NULL);
  gtk_widget_show(menu_item);

  menu_label_names = gtk_check_menu_item_new_with_label("Show atom names");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_label_names), FALSE);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_label_names);
  g_signal_connect(G_OBJECT(menu_label_names), "activate", G_CALLBACK(luscus_gtk_show_names), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(menu_label_names), FALSE);
  gtk_widget_show(menu_label_names);

  menu_label_muliken = gtk_check_menu_item_new_with_label("Show Mulliken charges");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_label_muliken), FALSE);
  g_signal_connect(G_OBJECT(menu_label_muliken), "activate", G_CALLBACK(luscus_gtk_show_mulliken), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_label_muliken);
  gtk_widget_set_sensitive(GTK_WIDGET(menu_label_muliken), FALSE);
  gtk_widget_show(menu_label_muliken);

  menu_label_loprop = gtk_check_menu_item_new_with_label("Show Loprop charges");
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_label_loprop), FALSE);
  g_signal_connect(G_OBJECT(menu_label_loprop), "activate", G_CALLBACK(luscus_gtk_show_loprop), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_label_loprop);
  gtk_widget_set_sensitive(GTK_WIDGET(menu_label_loprop), FALSE);
  gtk_widget_show(menu_label_loprop);

  menu_item = gtk_menu_item_new_with_mnemonic("_Labels");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

  menu = gtk_menu_new();

  /*Screenshot*/
  menu_item = gtk_menu_item_new_with_label("Screenshot");
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(callback_screenshot), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_menu_item_new_with_mnemonic("_Screenshot");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

    /*Help*/
  menu = gtk_menu_new();

  /*-Help-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_HELP, ag);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_callback_help), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  menu_item = gtk_separator_menu_item_new();
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);

  /*-About-*/
  menu_item = gtk_image_menu_item_new_from_stock(GTK_STOCK_ABOUT, ag);
  g_signal_connect(G_OBJECT(menu_item), "activate", G_CALLBACK(luscus_gtk_callback_about), NULL);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu), menu_item);
  gtk_widget_show(menu_item);
  
  menu_item = gtk_menu_item_new_with_mnemonic("_Help");
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_item), menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(menubar), menu_item);
  gtk_widget_show(menu_item);

  gtk_widget_show(menubar);
  return menubar;
}

void luscus_gtk_menubar_show_or_hide_widgets(void)
{
  int i;
  int filetype;
  int n_calc_menus;

#ifdef EBUG
  printf("CHECKING WIDGETS IN MENUBAR TO SHOW/HIDE\n");
#endif

  filetype = get_filetype(get_input_filename());
  if (input_filetypes[filetype].backward) gtk_widget_set_sensitive(GTK_WIDGET(menu_item_save), TRUE);
  else gtk_widget_set_sensitive(GTK_WIDGET(menu_item_save), FALSE);

  if (m->ishow & HAS_ATOMNUMS) gtk_widget_set_sensitive(GTK_WIDGET(menu_label_numeration), TRUE);
  else gtk_widget_set_sensitive(GTK_WIDGET(menu_label_numeration), FALSE);

  if (m->ishow & HAS_ATOMNAMES) gtk_widget_set_sensitive(GTK_WIDGET(menu_label_names), TRUE);
  else gtk_widget_set_sensitive(GTK_WIDGET(menu_label_names), FALSE);

  if (m->ishow & HAS_MULLIKEN) gtk_widget_set_sensitive(GTK_WIDGET(menu_label_muliken), TRUE);
  else gtk_widget_set_sensitive(GTK_WIDGET(menu_label_muliken), FALSE);

  if (m->ishow & HAS_LOPROP) gtk_widget_set_sensitive(GTK_WIDGET(menu_label_loprop), TRUE);
  else gtk_widget_set_sensitive(GTK_WIDGET(menu_label_loprop), FALSE);

  if (Input_Data.automatic_rebonding)
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_automatic_bonding), TRUE);
  else
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_automatic_bonding), FALSE);

  if (m->editable)
  {
    gtk_widget_set_sensitive(GTK_WIDGET(menu_item_xyz_edit), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(menu_undo), TRUE);
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_editable), TRUE);

    /*molecular mechanics*/
    for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
    {
      if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"mm"))
        n_calc_menus++;
    }
    if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_mechanics), FALSE);
  
    /*semiempirical*/
    for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
    {
      if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"semiem"))
        n_calc_menus++;
    }
    if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_semiempirical), FALSE);
  
    /*ab-initio*/
    for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
    {
      if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"abinitio"))
        n_calc_menus++;
    }
    if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_ab_initio), FALSE);
  }
  else
  {
    gtk_widget_set_sensitive(GTK_WIDGET(menu_item_xyz_edit), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(menu_undo), FALSE);
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_editable), FALSE);

    gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_mechanics), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_semiempirical), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(menu_optimize_ab_initio), FALSE);
  }
  if (m->ngrids)
  {
    /*orbital calculations*/
    for(i = 0, n_calc_menus=0; i < n_calc_defs; i++)
    {
      if (calc_defs[i].plugin_name && strcasestr(calc_defs[i].type,"orbital"))
        n_calc_menus++;
    }
    if (! n_calc_menus) gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), FALSE);
    else gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), TRUE);
  }
  else
  {
    gtk_widget_set_sensitive(GTK_WIDGET(menu_orbital_calculation), FALSE);
  }
  
}

