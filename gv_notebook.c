/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<gtk/gtk.h>
#include<string.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_gtk.h"
/*#include"gv_gtk_fragments.h"*/
#include "center16.xpm"
#include "icon_revert20.xpm"

/*#define NONE 0
#define SPHERE 1
#define VECTOR 2
#define TRIANGLE 3
#define SURFACE 4
#define CELL 5

typedef struct gtk_3d_description
{
  GtkWidget *button;
  char i3dtype;
  int inum;
} GTK_3D_DESC;*/

int n_fragments;
FRAG_DATA *frag;
GtkWidget **freq_buttons = NULL;
GtkWidget **freq_labels = NULL;
int n_vibs = 0;
int n_3D_desc = 0;
GTK_3D_DESC *gtk_3d_desc = NULL;

static gboolean luscus_gtk_find_selected_orbital(gint, gint);


char *bond_description[NUM_OF_BONDS] =
{
  "No bond",
  "Single bond",
  "double bond",
  "triple bond",
  "partial bond",
  "One and half bond",
  "Line"
};

char *orbital_description_short[] = 
{
  "u", "f", "i", "1", "2", "3", "s", "d"
};

char *angle_description[NUM_OF_ANGLES] =
{
  "60°",
  "90°",
  "109.47°",
  "120°",
  "180°"
};

char *orbital_description[NUM_ORB_TYPES] =
{
  "undefined",
  "frozen",
  "inactive",
  "RAS 1",
  "RAS 2",
  "RAS 3",
  "secondary",
  "deleted"
};

MOL *m;
INPUT_DATA Input_Data;

GtkWidget *notebook = NULL;

GtkWidget *vbox_watched;
int n_label_watched;
GtkWidget **labels_watched = NULL;
GtkWidget *vbox_vibration;
GtkWidget *vbox_edit;
GtkWidget *table_fragments;
GtkWidget *vbox_symmetry;
GtkWidget *vbox_draw;
GtkWidget *vbox_orbitals;

GtkTreeIter *selected_orbital;
GtkListStore *store;
GtkTreeSelection *selection;
GtkWidget *treeview;
GtkWidget *geo_separator;
GtkWidget *geo_button;
GtkWidget *geo_play_button;
GtkWidget *geo_hbox;
GtkWidget *geo_slider;
GtkWidget *label_igeo;
GtkWidget *spin_bond;
GtkWidget *spin_angle;
GtkWidget *spin_tors;
/*GtkWidget *spin_play_speed;*/
GtkWidget *hbox_xyz, *spin_x, *spin_y, *spin_z;
GtkWidget *grid_expander;
GtkWidget *geo_expander;
GtkWidget *grid_frame_transp;
GtkWidget *grid_frame_isosurf;
GtkWidget *all_orbitals_button;
GtkWidget *epot_button;
/*GtkWidget *dens_diff_button;*/
GtkWidget *vib_spectrum_button;
GtkWidget *vbox_freqs;
GtkWidget *vbox_3Dobject;

GtkWidget *combo_button;

GtkWidget *button_add_atom;
GtkWidget *button_remove_atom;
GtkWidget *button_change_atom;
GtkWidget *button_undo;

/*GtkWidget *label_bond_angle_torsion;
GtkWidget *label_unit;*/
GtkWidget *button_change_bond;
GtkWidget *button_change_angle;
GtkWidget *button_add_dummy;
GtkWidget *button_remove_selection;
GtkWidget *button_unmark;
GtkWidget *button_renumber;
GtkWidget *button_mark_H, *button_mark_reverse, *button_mark_element, *button_neighbor;
GtkWidget *button_watch;
GtkWidget *button_unwatch;
GtkWidget *button_symmetry;
GtkWidget *button_average_symm;
GtkWidget *button_trans_symm;
GtkWidget *button_inversion_symm;
GtkWidget *button_mirror_symm;
GtkWidget *spin_transl_vect;
GtkWidget *hbox_transl;
GtkWidget *button_c2;
GtkWidget *spin_axis;
GtkWidget *hbox_rot_axis;

GtkWidget *button_arrow;
GtkWidget *button_sphere;
GtkWidget *button_plain;
GtkWidget *button_triangle;
GtkWidget *button_cell;
GtkWidget *button_clear_drawed;

GtkWidget **button_custom_grid = NULL;
gint n_custom_grid = 0;
/*GtkWidget *element_dialog;*/
#ifdef GTK2
  GtkObject *tranparency_level, *isosurface_level, *vibration_speed, *vibration_amplitude;
  GtkObject *adj_bond;
  GtkObject *adj_angle;
  GtkObject *adj_tors;
  GtkObject *translation_magnitude, *rotation_axis;
  GtkObject *adj_geoms, *adj_play_speed;
#elif GTK3
  GtkAdjustment *tranparency_level, *isosurface_level, *vibration_speed, *vibration_amplitude;
  GtkAdjustment *adj_bond;
  GtkAdjustment *adj_angle;
  GtkAdjustment *adj_tors;
  GtkAdjustment *translation_magnitude, *rotation_axis;
  GtkAdjustment *adj_geoms, *adj_play_speed;
#endif

GtkWidget *make_notebook(void)
{
  GtkWidget *vbox_outer, *hbox, *vbox_inner;
  GtkWidget *label;
  GtkWidget *scale;
  GtkWidget *frame;
  GtkWidget *table;
  GtkWidget *button;
  GtkWidget *icon;
  GtkWidget *separator;
  GtkWidget *scrolled_window;
  GtkWidget *entry;
  GtkWidget *expander;

/*  GList *combo_list = NULL;*/

  GdkPixbuf *pixb;
  gint i;
  gint ix, iy;
  gchar *textbuf;

#ifdef GTK2
  vbox_outer = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_outer = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_outer), FALSE);
#endif

  notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(vbox_outer), notebook, TRUE, TRUE, 0);

  /*--------frequency--------*/

#ifdef GTK2
  vbox_vibration = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_vibration = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_vibration), FALSE);
#endif
  frame = gtk_frame_new("vibration speed");
  gtk_box_pack_start(GTK_BOX(vbox_vibration), frame, FALSE, FALSE, 0);

  vibration_speed = gtk_adjustment_new(0.2, 0.F, 1.F, 0.01F, 0.1F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(vibration_speed), Input_Data.frequency_speed);
  g_signal_connect(G_OBJECT(vibration_speed), "value-changed", G_CALLBACK(luscus_gtk_change_vibration_speed), NULL);
#ifdef GTK2
  scale = gtk_hscale_new(GTK_ADJUSTMENT(vibration_speed));
#elif GTK3
  scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, vibration_speed);
#endif

  gtk_scale_set_digits(GTK_SCALE(scale), 2);
  gtk_scale_set_draw_value(GTK_SCALE(scale), FALSE);
  gtk_container_add(GTK_CONTAINER(frame), scale);
  gtk_widget_show(scale);

  gtk_widget_show(frame);

 /* amplitude */

  frame = gtk_frame_new("amplitude");
  gtk_box_pack_start(GTK_BOX(vbox_vibration), frame, FALSE, FALSE, 0);

  vibration_amplitude = gtk_adjustment_new(0.F, 0.F, 2.F, 1.F, 1.F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(vibration_amplitude), Input_Data.frequency_amplitude);
  g_signal_connect(G_OBJECT(vibration_amplitude), "value-changed", G_CALLBACK(luscus_gtk_change_vibration_amplitude), NULL);
#ifdef GTK2
  scale = gtk_hscale_new(GTK_ADJUSTMENT(vibration_amplitude));
#elif GTK3
  scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, vibration_amplitude);
#endif
  gtk_scale_set_digits(GTK_SCALE(scale), 0);
  gtk_scale_set_draw_value(GTK_SCALE(scale), FALSE);
  gtk_container_add(GTK_CONTAINER(frame), scale);
  gtk_widget_show(scale);

  gtk_widget_show(frame);

  /*frequency buttons*/
  scrolled_window = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
  gtk_box_pack_start(GTK_BOX(vbox_vibration), scrolled_window, TRUE, TRUE, 0);
  
#ifdef GTK2
  vbox_freqs = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_freqs = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_freqs), FALSE);
#endif
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window), vbox_freqs);

  gtk_widget_hide(vbox_freqs);
  gtk_widget_show(scrolled_window);

  vib_spectrum_button = gtk_button_new_with_label("vibrational spectrum");
  gtk_box_pack_start(GTK_BOX(vbox_vibration), vib_spectrum_button, FALSE, FALSE, 5);
  g_signal_connect(G_OBJECT(vib_spectrum_button), "clicked", G_CALLBACK(luscus_gtk_show_vib_spectrum), NULL);
  gtk_widget_hide(vib_spectrum_button);  

  gtk_widget_hide(vbox_vibration);

  label = gtk_label_new("vibrations");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox_vibration, label);

  /*--------edit--------*/

#ifdef GTK2
  vbox_edit = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_edit = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_edit), FALSE);
#endif

  label = gtk_label_new("edit molecule");
  gtk_box_pack_start(GTK_BOX(vbox_edit), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

#ifdef GTK2
  hbox = gtk_hbox_new(TRUE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), TRUE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  button_add_atom = gtk_button_new_from_stock(GTK_STOCK_ADD);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_add_atom), "Add an atom"); /*add bond*/
  gtk_box_pack_start(GTK_BOX(hbox), button_add_atom, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_add_atom), "clicked", G_CALLBACK(luscus_gtk_add_fragment), NULL);
  gtk_widget_show(button_add_atom);

  button_remove_atom = gtk_button_new_from_stock(GTK_STOCK_REMOVE);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_remove_atom), "delete atom"); /*delete bond*/
  gtk_box_pack_start(GTK_BOX(hbox), button_remove_atom, FALSE, FALSE, 0);
  gtk_widget_show(button_remove_atom);
  g_signal_connect(G_OBJECT(button_remove_atom), "clicked", G_CALLBACK(luscus_gtk_remove_fragment), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_remove_atom), FALSE);
  gtk_widget_show(hbox);

  button_change_atom = gtk_button_new_with_label("Change");
  gtk_box_pack_start(GTK_BOX(hbox), button_change_atom, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_change_atom), "clicked", G_CALLBACK(callback_change_atom), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_change_atom), FALSE);
  gtk_widget_show(button_change_atom);

  button_undo = gtk_button_new_from_stock(GTK_STOCK_UNDO);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_undo), "Undo last action");
  gtk_box_pack_start(GTK_BOX(hbox), button_undo, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_undo), "clicked", G_CALLBACK(callback_do_undo), NULL);
  gtk_widget_show(button_undo);

  adj_bond = gtk_adjustment_new(0.000F, 0.F, 1000.F, 0.01F, 0.1F, 0.F);
  adj_angle = gtk_adjustment_new(0.000F, 0.F, 180.F, 0.1F, 1.0F, 0.F);
  adj_tors = gtk_adjustment_new(0.000F, -360.0F, 360.0F, 0.2F, 1.0F, 0.F);

#ifdef GTK2
  spin_bond = gtk_spin_button_new(GTK_ADJUSTMENT(adj_bond), 0.01, 2);
  spin_angle = gtk_spin_button_new(GTK_ADJUSTMENT(adj_angle), 0.01, 2);
  spin_tors = gtk_spin_button_new(GTK_ADJUSTMENT(adj_tors), 0.01, 2);

  hbox = gtk_hbox_new(FALSE, 0);
#elif GTK3
  spin_bond = gtk_spin_button_new(adj_bond, 0.01, 2);
  spin_angle = gtk_spin_button_new(adj_angle, 0.01, 2);
  spin_tors = gtk_spin_button_new(adj_tors, 0.01, 2);

  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  label = gtk_label_new("bond: ");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  g_signal_connect(G_OBJECT(spin_bond), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_bond), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_bond), "value-changed", G_CALLBACK(luscus_gtk_change_bond_value), NULL);

  gtk_box_pack_start(GTK_BOX(hbox), spin_bond, FALSE, FALSE, 0);
  gtk_widget_set_sensitive(GTK_WIDGET(spin_bond), FALSE);
  gtk_widget_show(spin_bond);

  label = gtk_label_new("\303\205");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  gtk_widget_show(hbox);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  label = gtk_label_new("angle: ");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  g_signal_connect(G_OBJECT(spin_angle), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_angle), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_angle), "value-changed", G_CALLBACK(luscus_gtk_change_angle_value), NULL);

/*  g_signal_connect(G_OBJECT(spin_bond_ang_tor), "value-changed", G_CALLBACK(luscus_gtk_change_bond_angle_tors_value), NULL);*/

  gtk_box_pack_start(GTK_BOX(hbox), spin_angle, FALSE, FALSE, 0);
  gtk_widget_set_sensitive(GTK_WIDGET(spin_angle), FALSE);
  gtk_widget_show(spin_angle);

  label = gtk_label_new("\313\232");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  gtk_widget_show(hbox);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  label = gtk_label_new("dihedral: ");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  g_signal_connect(G_OBJECT(spin_tors), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_tors), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_tors), "value-changed", G_CALLBACK(luscus_gtk_change_tors_value), NULL);

/*  g_signal_connect(G_OBJECT(spin_bond_ang_tor), "value-changed", G_CALLBACK(luscus_gtk_change_bond_angle_tors_value), NULL);*/

  gtk_box_pack_start(GTK_BOX(hbox), spin_tors, FALSE, FALSE, 0);
  gtk_widget_set_sensitive(GTK_WIDGET(spin_tors), FALSE);
  gtk_widget_show(spin_tors);

  label = gtk_label_new("\313\232");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  gtk_widget_show(hbox);

#ifdef GTK2
  hbox_xyz = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox_xyz = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox_xyz), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox_xyz, FALSE, FALSE, 0);
 
  label = gtk_label_new("x:");
  gtk_box_pack_start(GTK_BOX(hbox_xyz), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  spin_x = gtk_spin_button_new_with_range(-999.0, 999.0, 0.01);
  gtk_spin_button_set_increments(GTK_SPIN_BUTTON(spin_x), 0.01, 0.1);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin_x), 3);
  gtk_box_pack_start(GTK_BOX(hbox_xyz), spin_x, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(spin_x), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_x), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_x), "value-changed", G_CALLBACK(luscus_gtk_change_xyz_coordinate), GINT_TO_POINTER(0));
  gtk_widget_show(spin_x);

  label = gtk_label_new("y:");
  gtk_box_pack_start(GTK_BOX(hbox_xyz), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  spin_y = gtk_spin_button_new_with_range(-999.0, 999.0, 0.01);
  gtk_spin_button_set_increments(GTK_SPIN_BUTTON(spin_x), 0.01, 0.1);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin_y), 3);
  gtk_box_pack_start(GTK_BOX(hbox_xyz), spin_y, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(spin_y), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_y), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_y), "value-changed", G_CALLBACK(luscus_gtk_change_xyz_coordinate), GINT_TO_POINTER(1));
  gtk_widget_show(spin_y);

  label = gtk_label_new("z:");
  gtk_box_pack_start(GTK_BOX(hbox_xyz), label, FALSE, FALSE, 0);
  gtk_widget_show(label);
  
  spin_z = gtk_spin_button_new_with_range(-999.0, 999.0, 0.01);
  gtk_spin_button_set_increments(GTK_SPIN_BUTTON(spin_x), 0.01, 0.1);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin_z), 3);
  gtk_box_pack_start(GTK_BOX(hbox_xyz), spin_z, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(spin_z), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_z), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_z), "value-changed", G_CALLBACK(luscus_gtk_change_xyz_coordinate), GINT_TO_POINTER(2));
  gtk_widget_show(spin_z);

  gtk_widget_hide(hbox_xyz);

#ifdef GTK_OLD
  button_change_bond = gtk_combo_box_new_text();
  for(i = 0; i < NUM_OF_BONDS; i++)
    gtk_combo_box_append_text(GTK_COMBO_BOX(button_change_bond), bond_description[i]);
#else
  button_change_bond = gtk_combo_box_text_new();
  for(i = 0; i < NUM_OF_BONDS; i++)
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(button_change_bond), bond_description[i]);
#endif
  g_signal_connect(G_OBJECT(button_change_bond), "changed", G_CALLBACK(luscus_gtk_change_bond_callback), NULL);
  gtk_box_pack_start(GTK_BOX(vbox_edit), button_change_bond, FALSE, FALSE, 0);
  gtk_widget_hide(button_change_bond);

#ifdef GTK_OLD
  button_change_angle = gtk_combo_box_new_text();
  for(i = 0; i < NUM_OF_ANGLES; i++)
    gtk_combo_box_append_text(GTK_COMBO_BOX(button_change_angle), angle_description[i]);
  gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_angle), -1);
#else
  button_change_angle = gtk_combo_box_text_new();
  for(i = 0; i < NUM_OF_ANGLES; i++)
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(button_change_angle), angle_description[i]);
  gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_angle), -1);
#endif
  g_signal_connect(G_OBJECT(button_change_angle), "changed", G_CALLBACK(luscus_gtk_change_angle_callback), NULL);
  gtk_box_pack_start(GTK_BOX(vbox_edit), button_change_angle, FALSE, FALSE, 0);
  gtk_widget_hide(button_change_angle);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif

  gtk_box_pack_start(GTK_BOX(vbox_edit), separator, FALSE, FALSE, 5);
  gtk_widget_show(separator);

  button_remove_selection = gtk_button_new_with_label("unselect");
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_remove_selection), "add dummy atoms (reference points) on the direction of axis");
  gtk_box_pack_start(GTK_BOX(vbox_edit), button_remove_selection, FALSE, FALSE, 0);
  gtk_widget_set_sensitive(GTK_WIDGET(button_remove_atom), FALSE);
  g_signal_connect(G_OBJECT(button_remove_selection), "clicked", G_CALLBACK(luscus_gtk_remove_selection), NULL);
  gtk_widget_show(button_remove_selection);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), separator, FALSE, FALSE, 5);
  gtk_widget_show(separator);

  label = gtk_label_new("mark atoms");
  gtk_box_pack_start(GTK_BOX(vbox_edit), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

#ifdef GTK2
  hbox = gtk_hbox_new(TRUE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  icon = gtk_image_new_from_stock(GTK_STOCK_CLEAR, GTK_ICON_SIZE_BUTTON);
  button_unmark = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button_unmark), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_unmark), "unmark all marked atoms");
  gtk_box_pack_start(GTK_BOX(hbox), button_unmark, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_unmark), "clicked", G_CALLBACK(luscus_gtk_unmark), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_unmark), FALSE);
  gtk_widget_show(button_unmark);

  icon = gtk_image_new_from_stock(GTK_STOCK_SORT_ASCENDING, GTK_ICON_SIZE_BUTTON);
  button_renumber = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button_renumber), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_renumber), "renumber atoms and put marked atoms first");
  gtk_box_pack_start(GTK_BOX(hbox), button_renumber, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_renumber), "clicked", G_CALLBACK(luscus_gtk_sort_mark), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_renumber), FALSE);
  gtk_widget_show(button_renumber);

  textbuf = g_markup_escape_text("<span size=xx-large>H</span>", -1); /* g_markup_printf_escaped("<>H<>");*/
  label = gtk_label_new(NULL);
  gtk_label_set_markup(GTK_LABEL(label), textbuf);
  g_free(textbuf);  

  button_mark_H = gtk_button_new_with_label("H");
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_mark_H), "Mark hydrogen atoms in the molecule ");
  gtk_box_pack_start(GTK_BOX(hbox), button_mark_H, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_mark_H), "clicked", G_CALLBACK(luscus_gtk_mark_H_atoms), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_mark_H), TRUE);
  gtk_widget_show(button_mark_H);
  gtk_widget_show(label);

  pixb = gdk_pixbuf_new_from_xpm_data(icon_revert_xpm);
  icon = gtk_image_new_from_pixbuf(pixb);
  button_mark_reverse = gtk_button_new();
  gtk_container_add(GTK_CONTAINER(button_mark_reverse), icon);
/*  button_mark_reverse = gtk_button_new_with_label("reverse");*/
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_mark_reverse), "Reverse marked atoms");
  gtk_box_pack_start(GTK_BOX(hbox), button_mark_reverse, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_mark_reverse), "clicked", G_CALLBACK(luscus_gtk_mark_reverse), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_mark_reverse), TRUE);
  gtk_widget_show(button_mark_reverse);
  gtk_widget_show(icon);

  icon = gtk_image_new_from_stock(GTK_STOCK_SELECT_ALL, GTK_ICON_SIZE_BUTTON);
  button_mark_element = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button_mark_element), icon);
/*  button_mark_element = gtk_button_new_with_label("mark element");*/
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_mark_element), "Mark all atoms which are the same element as selected atom");
  gtk_box_pack_start(GTK_BOX(hbox), button_mark_element, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_mark_element), "clicked", G_CALLBACK(luscus_gtk_mark_element), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_mark_element), FALSE);
  gtk_widget_show(button_mark_element);

  button_neighbor = gtk_button_new_with_label("\342\212\266");
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_neighbor), "Mark atoms bonded to the selected atom");
  gtk_box_pack_start(GTK_BOX(hbox), button_neighbor, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_neighbor), "clicked", G_CALLBACK(luscus_gtk_mark_neighbor), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_neighbor), FALSE);
  gtk_widget_show(button_neighbor);

  gtk_widget_show(hbox);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), separator, FALSE, FALSE, 5);
  gtk_widget_show(separator);

  label = gtk_label_new("watch values");
  gtk_box_pack_start(GTK_BOX(vbox_edit), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

#ifdef GTK2
  hbox = gtk_hbox_new(TRUE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), TRUE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  button_watch = gtk_button_new_from_stock(GTK_STOCK_ADD);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_watch), "Watch the value of the bond/angle/torsion");
  gtk_box_pack_start(GTK_BOX(hbox), button_watch, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_watch), "clicked", G_CALLBACK(luscus_gtk_watch_value), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_watch), FALSE);
  gtk_widget_show(button_watch);

  button_unwatch = gtk_button_new_from_stock(GTK_STOCK_REMOVE);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_unwatch), "Remove all watched values");
  gtk_box_pack_start(GTK_BOX(hbox), button_unwatch, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_unwatch), "clicked", G_CALLBACK(luscus_gtk_unwatch_values), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_unwatch), FALSE);
  gtk_widget_show(button_unwatch);

  gtk_widget_show(hbox);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), separator, FALSE, FALSE, 5);
  gtk_widget_show(separator);

  label = gtk_label_new("dummy atoms");
  gtk_box_pack_start(GTK_BOX(vbox_edit), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_edit), hbox, FALSE, FALSE, 0);

  button_add_dummy = gtk_button_new_from_stock(GTK_STOCK_ADD);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_add_dummy), "add dummy atoms");
  gtk_box_pack_start(GTK_BOX(hbox), button_add_dummy, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_add_dummy), "clicked", G_CALLBACK(luscus_gtk_add_dummy_atoms), NULL);
  gtk_widget_show(button_add_dummy);

  button = gtk_button_new_from_stock(GTK_STOCK_REMOVE);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "remove dummy atoms");
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_remove_dummy_atoms), NULL);
  gtk_widget_show(button);

  gtk_widget_show(hbox);

  gtk_widget_show(vbox_edit);
  label = gtk_label_new("edit");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox_edit, label);

  /*-------Fragments-------*/

  table_fragments = gtk_table_new(5, 4, TRUE);

  for (i=0; i < n_fragments; i++)
  {
    ix = i%4;
    iy = (i - ix)/4;

    button = gtk_button_new();
/*    pix = gdk_pixbuf_new_from_xpm_data(a_xpms[i]);*/
/*    icon = gtk_image_new_from_pixbuf(pix);*/
#ifdef WINDOWS
    textbuf = g_strconcat(rcdir, "\\", frag[i].image_file_name, NULL);
#else
    textbuf = g_strconcat(rcdir, "/", frag[i].image_file_name, NULL);
#endif
    icon = gtk_image_new_from_file(textbuf);
    g_free(textbuf);
    gtk_container_add(GTK_CONTAINER(button), icon);
    gtk_table_attach(GTK_TABLE(table_fragments), button, ix, ix+1, iy, iy+1, GTK_FILL, GTK_FILL, 1, 1);

    g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_select_fragment), GINT_TO_POINTER(i));
    gtk_widget_show(button);
    gtk_widget_show(icon);
  }

  gtk_widget_show(table_fragments);
  label = gtk_label_new("fragments");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), table_fragments, label);
  gtk_widget_show(table_fragments);

   /*symmetry*/

#ifdef GTK2
  vbox_symmetry = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_symmetry = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_symmetry), FALSE);
#endif

#ifdef HAS_MSYM
  /*here goes MSYM stuff*/
#else
  button_symmetry = gtk_button_new_with_label("Symmetry"); /*inversion*/
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_symmetry), "analyze the symmetry of the molecule and display symmetry elements");
  gtk_box_pack_start(GTK_BOX(vbox_symmetry), button_symmetry, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_symmetry), "clicked", G_CALLBACK(luscus_gtk_alpply_symmetry), NULL);
/*  gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);*/
  gtk_widget_show(button_symmetry);

  button_average_symm = gtk_button_new_with_label("Average symmetry");
  gtk_box_pack_start(GTK_BOX(vbox_symmetry), button_average_symm, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_average_symm), "clicked", G_CALLBACK(luscus_gtk_force_symmetry), NULL);
/*  gtk_widget_set_sensitive(GTK_WIDGET(button_average_symm), FALSE);*/
  gtk_widget_show(button_average_symm);
#endif

  button_inversion_symm = gtk_button_new_with_label("apply inversion symmetry");
  gtk_box_pack_start(GTK_BOX(vbox_symmetry), button_inversion_symm, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_inversion_symm), "clicked", G_CALLBACK(luscus_gtk_alpply_symmetry), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_inversion_symm), FALSE);
  gtk_widget_show(button_inversion_symm);

#ifdef GTK2
  hbox_transl = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox_transl = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox_transl), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_symmetry), hbox_transl, FALSE, FALSE, 0);

  button_trans_symm = gtk_button_new_with_label("apply translation");
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_trans_symm), "Translate molecule in the direction of selected atoms");
  gtk_box_pack_start(GTK_BOX(hbox_transl), button_trans_symm, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_trans_symm), "clicked", G_CALLBACK(luscus_gtk_apply_translation), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_trans_symm), FALSE);
  gtk_widget_show(button_trans_symm);

  label = gtk_label_new(" r=");
  gtk_box_pack_start(GTK_BOX(hbox_transl), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  translation_magnitude = gtk_adjustment_new(0.F, 0.F, 100.F, 1.F, 1.F, 0.F);
#ifdef GTK2
  spin_transl_vect = gtk_spin_button_new(GTK_ADJUSTMENT(translation_magnitude), 0.01, 2);
#elif GTK3
  spin_transl_vect = gtk_spin_button_new(translation_magnitude, 0.01, 2);
#endif
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin_transl_vect), 2);
  g_signal_connect(G_OBJECT(spin_transl_vect), "value-changed", G_CALLBACK(change_translation_magnitude), NULL);
  g_signal_connect(G_OBJECT(spin_transl_vect), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_transl_vect), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  gtk_box_pack_start(GTK_BOX(hbox_transl), spin_transl_vect, FALSE, FALSE, 0);
  gtk_widget_set_sensitive(GTK_WIDGET(spin_transl_vect), FALSE);
  gtk_widget_show(spin_transl_vect);

  gtk_widget_show(hbox_transl);

#ifdef GTK2
  hbox_rot_axis = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox_rot_axis = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox_transl), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_symmetry), hbox_rot_axis, FALSE, FALSE, 0);

  button_c2 = gtk_button_new_with_label("apply Cn symmetry");
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_c2), "apply Cn symmetry operation around axis defined by selected atoms");
  gtk_box_pack_start(GTK_BOX(hbox_rot_axis), button_c2, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_c2), "clicked", G_CALLBACK(luscus_gtk_alpply_rot_symmetry), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_c2), FALSE);
  gtk_widget_show(button_c2);

  label = gtk_label_new(" n=");
  gtk_box_pack_start(GTK_BOX(hbox_rot_axis), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  rotation_axis = gtk_adjustment_new(2.F, 2.F, 100.F, 1.F, 1.F, 0.F);
#ifdef GTK2
  spin_axis = gtk_spin_button_new(GTK_ADJUSTMENT(rotation_axis), 0.01, 2);
#elif GTK3
  spin_axis = gtk_spin_button_new(rotation_axis, 0.01, 2);
#endif
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(spin_axis), 0);
  g_signal_connect(G_OBJECT(spin_axis), "value-changed", G_CALLBACK(change_rotation_axis_order), NULL);
  g_signal_connect(G_OBJECT(spin_axis), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(spin_axis), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  gtk_box_pack_start(GTK_BOX(hbox_rot_axis), spin_axis, FALSE, FALSE, 0);
  gtk_widget_set_sensitive(GTK_WIDGET(spin_axis), FALSE);
  gtk_widget_show(spin_axis);

  gtk_widget_show(hbox_rot_axis);

  button_mirror_symm = gtk_button_new_with_label("apply mirror symmetry");
  gtk_box_pack_start(GTK_BOX(vbox_symmetry), button_mirror_symm, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_mirror_symm), "clicked", G_CALLBACK(luscus_gtk_alpply_symmetry), NULL);
  gtk_widget_set_sensitive(GTK_WIDGET(button_mirror_symm), FALSE);
  gtk_widget_show(button_mirror_symm);

  label = gtk_label_new("Symmetry");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox_symmetry, label);
  gtk_widget_show(vbox_symmetry);

  /*objects*/

#ifdef GTK2
  vbox_draw = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_draw = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_draw), FALSE);
#endif

  button_arrow = gtk_button_new_with_label("draw arrow");
  gtk_box_pack_start(GTK_BOX(vbox_draw), button_arrow, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_arrow), "clicked", G_CALLBACK(luscus_add_arrow), NULL);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_arrow), "select two atoms that define vector");
  gtk_widget_set_sensitive(GTK_WIDGET(button_arrow), FALSE);
  gtk_widget_show(button_arrow);

  button_sphere = gtk_button_new_with_label("draw sphere");
  gtk_box_pack_start(GTK_BOX(vbox_draw), button_sphere, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_sphere), "clicked", G_CALLBACK(luscus_add_sphere), NULL);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_sphere), "select atom that defines center of sphere");
  gtk_widget_set_sensitive(GTK_WIDGET(button_sphere), FALSE);
  gtk_widget_show(button_sphere);

  button_plain = gtk_button_new_with_label("draw plain");
  gtk_box_pack_start(GTK_BOX(vbox_draw), button_plain, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_plain), "clicked", G_CALLBACK(luscus_add_plain), NULL);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_plain), "select three atoms that define plain");
  gtk_widget_set_sensitive(GTK_WIDGET(button_plain), FALSE);
  gtk_widget_show(button_plain);

  button_triangle = gtk_button_new_with_label("draw triangle");
  gtk_box_pack_start(GTK_BOX(vbox_draw), button_triangle, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_triangle), "clicked", G_CALLBACK(luscus_add_triangle), NULL);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_triangle), "select three atoms that define triangle");
  gtk_widget_set_sensitive(GTK_WIDGET(button_triangle), FALSE);
  gtk_widget_show(button_triangle);

  button_cell = gtk_button_new_with_label("draw cell");
  gtk_box_pack_start(GTK_BOX(vbox_draw), button_cell, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button_cell), "clicked", G_CALLBACK(luscus_add_cell), NULL);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_cell), "select four atoms that define a cell");
  gtk_widget_set_sensitive(GTK_WIDGET(button_cell), FALSE);
  gtk_widget_show(button_cell);

  button = gtk_button_new_with_label("insert text");
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_insert_textbox), NULL);
  gtk_box_pack_start(GTK_BOX(vbox_draw), button, FALSE, FALSE, 0);
  gtk_widget_show(button);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_draw), separator, FALSE, FALSE, 5);
  gtk_widget_show(separator);

  button_clear_drawed = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
  gtk_box_pack_start(GTK_BOX(vbox_draw), button_clear_drawed, FALSE, FALSE, 0);
  g_signal_connect(GTK_WIDGET(button_clear_drawed), "clicked", G_CALLBACK(luscus_clear_drawings), NULL);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button_clear_drawed), "Clear all vectors, spheres, triangles, plains and cells");
  gtk_widget_set_sensitive(GTK_WIDGET(button_clear_drawed), FALSE);
  gtk_widget_show(button_clear_drawed);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_draw), separator, FALSE, FALSE, 5);
  gtk_widget_show(separator);

  scrolled_window = gtk_scrolled_window_new(NULL, NULL);
  gtk_box_pack_start(GTK_BOX(vbox_draw), scrolled_window, TRUE, TRUE, 0);
#ifdef GTK2
  vbox_3Dobject = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_3Dobject = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_3Dobject), FALSE);
#endif
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window), vbox_3Dobject);
  gtk_widget_show(vbox_3Dobject);
  gtk_widget_show(scrolled_window);

  /*color_button, transparency, radius, */

  label = gtk_label_new("draw");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox_draw, label);
  gtk_widget_show(vbox_draw);

  /*--------grid--------*/

#ifdef GTK2
  vbox_orbitals = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_orbitals = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_orbitals), FALSE);
#endif

  /*orbitals*/

  label = gtk_label_new("orbitals");
  gtk_widget_show(vbox_orbitals);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox_orbitals, label);

  /* treeview */

  scrolled_window = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
  gtk_box_pack_start(GTK_BOX(vbox_orbitals), scrolled_window, TRUE, TRUE, 0);

  treeview = luscus_gtk_make_treeview();
  gtk_container_add(GTK_CONTAINER(scrolled_window), treeview);
  gtk_widget_show(treeview);

  gtk_widget_show(scrolled_window);

#ifdef GTK2
  hbox = gtk_hbox_new(TRUE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_orbitals), hbox, FALSE, FALSE, 0);

  button = gtk_button_new_from_stock(GTK_STOCK_DELETE);
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(delete_orbital_from_list), NULL);
  gtk_widget_show(button);

  button = gtk_button_new_from_stock(GTK_STOCK_UNDELETE);
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(undelete_all_orbitals_from_list), NULL);
  gtk_widget_show(button);

/*  button = gtk_button_new_with_label("show density");
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(show_density_callback), NULL);
  gtk_widget_show(button);*/

  gtk_widget_show(hbox);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_orbitals), hbox, FALSE, FALSE, 0);

  label = gtk_label_new("filter");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(entry), "focus-in-event", G_CALLBACK(luscus_gtk_disable_keys), NULL);
  g_signal_connect(G_OBJECT(entry), "focus-out-event", G_CALLBACK(luscus_gtk_enable_keys), NULL);
  g_signal_connect(G_OBJECT(entry), "activate", G_CALLBACK(orbital_filter_callback), NULL);
  gtk_widget_show(entry);

  gtk_widget_show(hbox);

  gtk_widget_hide(vbox_orbitals);

  /*general*/

/*  vbox = gtk_vbox_new(FALSE, 0);*/

  expander = gtk_expander_new("zoom/translation");
  gtk_expander_set_expanded(GTK_EXPANDER(expander), FALSE);
  gtk_box_pack_start(GTK_BOX(vbox_outer), expander, FALSE, FALSE, 0);
  gtk_widget_show(expander);

#ifdef GTK2
  vbox_inner = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_inner = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_inner), FALSE);
#endif
  gtk_container_add(GTK_CONTAINER(expander), vbox_inner);

#ifdef GTK_OLD
  combo_button = gtk_combo_box_new_text();
  gtk_combo_box_append_text(GTK_COMBO_BOX(combo_button), "Move molecule");
  gtk_combo_box_append_text(GTK_COMBO_BOX(combo_button), "Move camera");
#else
  combo_button = gtk_combo_box_text_new();
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_button), "Move molecule");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_button), "Move camera");
#endif
  gtk_box_pack_start(GTK_BOX(vbox_inner), combo_button, FALSE, FALSE, 0);
  
/*  combo_list = g_list_append(combo_list, "Move molecule");
  combo_list = g_list_append(combo_list, "Move camera");
  gtk_combo_set_popdown_strings(GTK_COMBO(combo_button), combo_list);
  g_list_free(combo_list);*/

  if (move_camera)
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo_button), 1);
  else
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo_button), 0);
  g_signal_connect(G_OBJECT(combo_button), "changed", G_CALLBACK(luscus_toggle_movement_mode), NULL);
  gtk_widget_show(combo_button);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_inner), hbox, FALSE, FALSE, 5);
  gtk_widget_show(hbox);

  table = gtk_table_new(3, 3, TRUE);
  gtk_box_pack_start(GTK_BOX(hbox), table, TRUE, FALSE, 5);

  icon = gtk_image_new_from_stock(GTK_STOCK_GO_UP, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "move up");
  gtk_table_attach(GTK_TABLE(table), button, 1,2, 0,1, GTK_SHRINK, GTK_SHRINK, 3, 3);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_move_up), NULL);
  gtk_widget_show(button);

  icon = gtk_image_new_from_stock(GTK_STOCK_GO_FORWARD, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "move right");
  gtk_table_attach(GTK_TABLE(table), button, 2,3, 1,2, GTK_SHRINK, GTK_SHRINK, 3, 3);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_move_forw), NULL);
  gtk_widget_show(button);

  icon = gtk_image_new_from_stock(GTK_STOCK_GO_BACK, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "move left");
  gtk_table_attach(GTK_TABLE(table), button, 0,1, 1,2, GTK_SHRINK, GTK_SHRINK, 3, 3);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_move_back), NULL);
  gtk_widget_show(button);

  icon = gtk_image_new_from_stock(GTK_STOCK_GO_DOWN, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "move down");
  gtk_table_attach(GTK_TABLE(table), button, 1,2, 2,3, GTK_SHRINK, GTK_SHRINK, 3, 3);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_move_down), NULL);
  gtk_widget_show(button);

  pixb = gdk_pixbuf_new_from_xpm_data(center_xpm);
  icon = gtk_image_new_from_pixbuf(pixb);
  button = gtk_button_new();
  gtk_container_add(GTK_CONTAINER(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "center");
  gtk_table_attach(GTK_TABLE(table), button, 1,2, 1,2, GTK_SHRINK, GTK_SHRINK, 3, 3);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_center), NULL);
  gtk_widget_show(button);
  gtk_widget_show(icon);

  gtk_widget_show(table);
  gtk_widget_show(hbox);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_inner), hbox, FALSE, FALSE, 5);

  icon = gtk_image_new_from_stock(GTK_STOCK_ZOOM_OUT, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "zoom out");
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_zoom), (gpointer) -1);
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, FALSE, 0);
  gtk_widget_show(button);

  icon = gtk_image_new_from_stock(GTK_STOCK_ZOOM_100, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "Normal size");
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_zoom), (gpointer) 0);
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, FALSE, 0);
  gtk_widget_show(button);

  icon = gtk_image_new_from_stock(GTK_STOCK_ZOOM_IN, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "zoom in");
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_zoom), (gpointer) 1);
  gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, FALSE, 0);
  gtk_widget_show(button);

  gtk_widget_show(hbox);
  gtk_widget_show(vbox_inner);

/*multiple geometries options*/

  geo_expander = gtk_expander_new("geometries");
  gtk_expander_set_expanded(GTK_EXPANDER(geo_expander), FALSE);
  gtk_box_pack_start(GTK_BOX(vbox_outer), geo_expander, FALSE, FALSE, 0);
  gtk_widget_hide(geo_expander);

#ifdef GTK2
  vbox_inner = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_inner = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_inner), FALSE);
#endif

  gtk_container_add(GTK_CONTAINER(geo_expander), vbox_inner);

#ifdef GTK2
  geo_separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  geo_separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_inner), geo_separator, FALSE, FALSE, 5);
  gtk_widget_show(geo_separator);

  geo_button = gtk_button_new_with_label("graphical info");
  gtk_box_pack_start(GTK_BOX(vbox_inner), geo_button, FALSE, FALSE, 5);
  g_signal_connect(G_OBJECT(geo_button), "clicked", G_CALLBACK(luscus_gtk_show_graphical_energies), NULL);
  gtk_widget_show(geo_button);

#ifdef GTK2
  geo_hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  geo_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(geo_hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox_inner), geo_hbox, FALSE, FALSE, 5);

  adj_play_speed = gtk_adjustment_new(0.00F, 0.F, 1.F, 0.01F, 0.1F, 0.F);
 /* g_signal_connect();*/

  /*
#ifdef GTK2
  spin_play_speed = gtk_spin_button_new(GTK_ADJUSTMENT(adj_play_speed), 0.01, 2);
#endif
#ifdef GTK3
  spin_play_speed = gtk_spin_button_new(adj_play_speed, 0.01, 2);
#endif

  gtk_box_pack_start(GTK_BOX(geo_hbox), spin_play_speed, TRUE, FALSE, 0);
  gtk_widget_show(spin_play_speed);
  */


  icon = gtk_image_new_from_stock(GTK_STOCK_GOTO_FIRST, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "first geometry");
  gtk_box_pack_start(GTK_BOX(geo_hbox), button, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(callback_geo_first), NULL);
  gtk_widget_show(button);

  icon = gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY, GTK_ICON_SIZE_BUTTON);
  geo_play_button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(geo_play_button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(geo_play_button), "Play trajectory");
  gtk_box_pack_start(GTK_BOX(geo_hbox), geo_play_button, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(geo_play_button), "clicked", G_CALLBACK(callback_geo_play), NULL);
  gtk_widget_show(geo_play_button);

/*  icon = gtk_image_new_from_stock(GTK_STOCK_MEDIA_PAUSE, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "pause");
  gtk_box_pack_start(GTK_BOX(geo_hbox), button, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(callback_geo_forw), NULL);
  gtk_widget_show(button);*/

  icon = gtk_image_new_from_stock(GTK_STOCK_GOTO_LAST, GTK_ICON_SIZE_BUTTON);
  button = gtk_button_new();
  gtk_button_set_image(GTK_BUTTON(button), icon);
  gtk_widget_set_tooltip_text(GTK_WIDGET(button), "last geometry");
  gtk_box_pack_start(GTK_BOX(geo_hbox), button, TRUE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(callback_geo_last), NULL);
  gtk_widget_show(button);

  gtk_widget_show(geo_hbox);

  /*slider*/
  adj_geoms = gtk_adjustment_new(0.F, 0.F, (gdouble) n_geometries, 1.F, 1.F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_geoms), 0);
  g_signal_connect(G_OBJECT(adj_geoms), "value-changed", G_CALLBACK(callback_adj_geometry), NULL);
#ifdef GTK2
  geo_slider = gtk_hscale_new(GTK_ADJUSTMENT(adj_geoms));
#endif
#ifdef GTK3
  geo_slider = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adj_geoms);
#endif
  gtk_scale_set_digits(GTK_SCALE(geo_slider), 0);
  gtk_scale_set_draw_value(GTK_SCALE(geo_slider), TRUE);
  gtk_box_pack_start(GTK_BOX(vbox_inner), geo_slider, FALSE, FALSE, 0);
  gtk_widget_show(geo_slider);

  label_igeo = gtk_label_new(NULL);
  gtk_box_pack_start(GTK_BOX(vbox_inner), label_igeo, FALSE, FALSE, 0);
  gtk_widget_show(label_igeo);

  gtk_widget_show(vbox_inner);

  /*grid info*/

  grid_expander = gtk_expander_new("surface");
  gtk_expander_set_expanded(GTK_EXPANDER(grid_expander), FALSE);
  gtk_box_pack_start(GTK_BOX(vbox_outer), grid_expander, FALSE, FALSE, 0);
  gtk_widget_hide(grid_expander);

#ifdef GTK2
  vbox_inner = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_inner = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_inner), FALSE);
#endif
  gtk_container_add(GTK_CONTAINER(grid_expander), vbox_inner);

  grid_frame_transp = gtk_frame_new("transparency");
  gtk_box_pack_start(GTK_BOX(vbox_inner), grid_frame_transp, FALSE, FALSE, 0);

  tranparency_level = gtk_adjustment_new(0.0, 0.F, 1.F, 0.01F, 0.1F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(tranparency_level), Input_Data.neg_pos_color[0][3]);
#ifdef GTK2
  scale = gtk_hscale_new(GTK_ADJUSTMENT(tranparency_level));
#elif GTK3
  scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, tranparency_level);
#endif
  gtk_scale_set_digits(GTK_SCALE(scale), 2);
  g_signal_connect(G_OBJECT(tranparency_level), "value-changed", G_CALLBACK(luscus_gtk_change_transparency_level), NULL);
  gtk_container_add(GTK_CONTAINER(grid_frame_transp), scale);
  gtk_widget_show(scale);

  gtk_widget_hide(grid_frame_transp);

  grid_frame_isosurf = gtk_frame_new("isosurface");

  gtk_box_pack_start(GTK_BOX(vbox_inner), grid_frame_isosurf, FALSE, FALSE, 0);
  isosurface_level = gtk_adjustment_new(0.0, 0.F, 0.2F, 0.001F, 0.01F, 0.F); 
  gtk_adjustment_set_value(GTK_ADJUSTMENT(isosurface_level), Input_Data.lev);

#ifdef GTK2
  scale = gtk_hscale_new(GTK_ADJUSTMENT(isosurface_level));
#elif GTK3
  scale = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, isosurface_level);
#endif
  gtk_scale_set_digits(GTK_SCALE(scale), 3);
  g_signal_connect(G_OBJECT(isosurface_level), "value-changed", G_CALLBACK(luscus_gtk_change_isosurface_level), NULL);
  gtk_container_add(GTK_CONTAINER(grid_frame_isosurf), scale);
  gtk_widget_show(scale);

  gtk_widget_hide(grid_frame_isosurf);

  all_orbitals_button = gtk_button_new_with_label("View all orbitals");
  gtk_box_pack_start(GTK_BOX(vbox_inner), all_orbitals_button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(all_orbitals_button), "clicked", G_CALLBACK(luscus_gtk_show_all_orbitals), NULL);
  gtk_widget_hide(all_orbitals_button);

  epot_button = gtk_button_new_with_label("Show electrostatic potential");
  gtk_box_pack_start(GTK_BOX(vbox_inner), epot_button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(epot_button), "clicked", G_CALLBACK(luscus_gtk_show_electrostatic_potential), NULL);
  gtk_widget_hide(epot_button);

  gtk_widget_show(vbox_inner);

  /*density difference*/
/*  dens_diff_button = gtk_button_new_with_label("Calculate density difference");
  gtk_box_pack_start(GTK_BOX(vbox_outer), dens_diff_button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(dens_diff_button), "clicked", G_CALLBACK(luscus_calculate_dens_difference), NULL);
  gtk_widget_hide(dens_diff_button);*/

  expander = gtk_expander_new("whatched values");
  gtk_expander_set_expanded(GTK_EXPANDER(expander), FALSE);
  gtk_box_pack_start(GTK_BOX(vbox_outer), expander, TRUE, TRUE, 0);
  gtk_widget_show(expander);

  scrolled_window = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);

#ifdef GTK2
  vbox_watched = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox_watched = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox_watched), FALSE);
#endif
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window), vbox_watched);
  gtk_widget_hide(vbox_watched);

  gtk_container_add(GTK_CONTAINER(expander), scrolled_window);
/*  gtk_box_pack_start(GTK_BOX(vbox_outer), scrolled_window, TRUE, TRUE, 0);*/
  gtk_widget_show(scrolled_window);

 /*list of all geometrical elements*/

/*  label = gtk_label_new("general");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label);
  gtk_widget_show(vbox);*/
  gtk_widget_show(notebook);

  gtk_widget_show(vbox_outer);
  return vbox_outer;
}

void luscus_gtk_diseable_editing(void)
{
/*GtkWidget *vbox_vibration;
GtkWidget *vbox_edit;
GtkWidget *vbox_fragments;
GtkWidget *vbox_symmetry;
GtkWidget *vbox_draw;
GtkWidget *vbox_orbitals;*/
  gtk_widget_set_sensitive(GTK_WIDGET(vbox_edit), FALSE);
  gtk_widget_set_sensitive(GTK_WIDGET(table_fragments), FALSE);
  gtk_widget_set_sensitive(GTK_WIDGET(vbox_symmetry), FALSE);
  gtk_widget_set_sensitive(GTK_WIDGET(vbox_draw), FALSE);

  luscus_gtk_menubar_show_or_hide_widgets();
}

void luscus_gtk_enable_editing(void)
{
  gtk_widget_set_sensitive(GTK_WIDGET(vbox_edit), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(table_fragments), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(vbox_symmetry), TRUE);
  gtk_widget_set_sensitive(GTK_WIDGET(vbox_draw), TRUE);
/*  gtk_widget_set_sensitive(GTK_WIDGET(notebook), TRUE);*/
  luscus_gtk_menubar_show_or_hide_widgets();
}

GtkWidget *luscus_gtk_make_treeview(void)
{
  GtkTreeIter iter_combo;
  GtkCellRenderer *renderer;
  GtkTreeViewColumn *column;
  GtkListStore *list_store_combo;
  GtkTreeSortable *sortable;
  gint i;
/*  gint data_sym, data_orbnum;
  gdouble data_energy, data_occu;
  gint orb_num_type;*/

  /*-------creating a view component-------*/

  treeview = gtk_tree_view_new();

  list_store_combo = gtk_list_store_new(1, G_TYPE_STRING);

  for(i = 0; i < NUM_ORB_TYPES; i++)
  {
    gtk_list_store_append(list_store_combo, &iter_combo);
    gtk_list_store_set(list_store_combo, &iter_combo, 0, orbital_description[i], -1);
  }

  renderer = gtk_cell_renderer_combo_new();
  g_object_set(renderer, "text-column", 0, "model", list_store_combo, "editable", TRUE, "has-entry", FALSE, NULL);
  g_signal_connect(renderer, "edited", G_CALLBACK(orbital_combo_edited), (gpointer) list_store_combo);
  gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "orbital type", renderer, "text", ORBITAL_TYPE, NULL);

  renderer = gtk_cell_renderer_text_new();
  gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "symmetry", renderer, "text", ORBITAL_SYMM, NULL);

  renderer = gtk_cell_renderer_text_new();
  gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "index", renderer, "text", ORBITAL_NUM, NULL);

  renderer = gtk_cell_renderer_text_new();
  gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "energy", renderer, "text", ORBITAL_ENERG, NULL);

  renderer = gtk_cell_renderer_text_new();
  gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, "occupancy", renderer, "text", ORBITAL_OCC, NULL);

  renderer = gtk_cell_renderer_text_new();
  gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview), -1, " ", renderer, "text", ORBITAL_EDITED, NULL);

  /*-------creating store(model)-------*/

  store = gtk_list_store_new(NUM_COLUMNS, G_TYPE_STRING, G_TYPE_INT, G_TYPE_INT, G_TYPE_DOUBLE, G_TYPE_DOUBLE, G_TYPE_STRING);
  sortable = GTK_TREE_SORTABLE(store);

  gtk_tree_sortable_set_sort_func(sortable, ORBITAL_SYMM, luscus_gtk_orbital_compare_func, GINT_TO_POINTER(ORBITAL_SYMM), NULL);
  gtk_tree_sortable_set_sort_func(sortable, ORBITAL_NUM,  luscus_gtk_orbital_compare_func, GINT_TO_POINTER(ORBITAL_NUM), NULL);
  gtk_tree_sortable_set_sort_func(sortable, ORBITAL_ENERG, luscus_gtk_orbital_compare_func, GINT_TO_POINTER(ORBITAL_ENERG), NULL);
  gtk_tree_sortable_set_sort_func(sortable, ORBITAL_OCC,  luscus_gtk_orbital_compare_func, GINT_TO_POINTER(ORBITAL_OCC), NULL);
  column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), ORBITAL_SYMM);
  gtk_tree_view_column_set_sort_column_id(column, ORBITAL_SYMM);
  column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), ORBITAL_NUM);
  gtk_tree_view_column_set_sort_column_id(column, ORBITAL_NUM);
  column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), ORBITAL_ENERG);
  gtk_tree_view_column_set_sort_column_id(column, ORBITAL_ENERG);
  column = gtk_tree_view_get_column(GTK_TREE_VIEW(treeview), ORBITAL_OCC);
  gtk_tree_view_column_set_sort_column_id(column, ORBITAL_OCC);

  selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));
  gtk_tree_selection_set_mode(selection, GTK_SELECTION_SINGLE | GTK_SELECTION_NONE);
  if (selected_orbital) gtk_tree_selection_select_iter(selection, selected_orbital);
  gtk_tree_selection_set_select_function(selection, luscus_gtk_select_orbital_from_list, NULL, NULL);

  return treeview;
}

gint luscus_gtk_orbital_compare_func(GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b, gpointer data)
{
  gint sort_col = GPOINTER_TO_INT(data);
  gint int_data1, int_data2;
  gdouble double_data1, double_data2;
  gint res = 0;

  switch(sort_col)
  {
    case ORBITAL_SYMM:
      gtk_tree_model_get(model, a, ORBITAL_SYMM, &int_data1, -1);
      gtk_tree_model_get(model, b, ORBITAL_SYMM, &int_data2, -1);
      res = (int_data1 > int_data2) ? -1 : 1;
      return res;
    case ORBITAL_NUM:
      gtk_tree_model_get(model, a, ORBITAL_NUM, &int_data1, -1);
      gtk_tree_model_get(model, b, ORBITAL_NUM, &int_data2, -1);
      res = (int_data1 > int_data2) ? -1 : 1;
      return res;
    case ORBITAL_ENERG:
      gtk_tree_model_get(model, a, ORBITAL_ENERG, &double_data1, -1);
      gtk_tree_model_get(model, b, ORBITAL_ENERG, &double_data2, -1);
      res = (double_data1 > double_data2) ? 1 : -1;
      return res;
    case ORBITAL_OCC:
      gtk_tree_model_get(model, a, ORBITAL_OCC, &double_data1, -1);
      gtk_tree_model_get(model, b, ORBITAL_OCC, &double_data2, -1);
      res = (double_data1 > double_data2) ? -1 : 1;
      return res;
  }
  return res;
}

gboolean luscus_gtk_select_orbital_from_list(GtkTreeSelection *selection, GtkTreeModel *model, GtkTreePath *path,
                                         gboolean path_selected, gpointer data)
{
/*  gint test_sym, test_orbnum;*/
  GtkTreeIter iter;
  char lin[82];
  int orb_symm, orb_index;
  double orb_occ, orb_energ;
  gchar *orb_type, *orb_edit;

  if (gtk_tree_model_get_iter(model, &iter, path))
  {
    if (!gtk_tree_selection_iter_is_selected(selection, &iter))
    {
/*      gtk_tree_model_get(GTK_TREE_MODEL(store) , &iter, ORBITAL_SYMM, &test_sym, ORBITAL_NUM, &test_orbnum, -1);
      luscus_gtk_find_selected_orbital(test_sym, test_orbnum);*/

      gtk_tree_model_get(GTK_TREE_MODEL(store), &iter,
                                                      ORBITAL_TYPE, &orb_type,
                                                      ORBITAL_SYMM, &orb_symm, 
                                                      ORBITAL_NUM,  &orb_index,
                                                      ORBITAL_ENERG,&orb_energ,
                                                      ORBITAL_OCC,  &orb_occ,
                                                      ORBITAL_EDITED,&orb_edit,  -1);
      snprintf(lin, 82, "%1s orbital sym: %2d index: %2d energy: %12.4f occ: %5.2f type: %s", orb_edit, orb_symm, orb_index, orb_energ, orb_occ, orb_type);
      luscus_gtk_pop_message_from_statusbar1();
      luscus_gtk_push_message_to_statusbar1(lin);

      luscus_gtk_find_selected_orbital(orb_symm, orb_index);
    }
  }

  return TRUE;
}

void luscus_gtk_show_or_hide_widgets(void)
{
  int i;

  for(i = 0; i < n_custom_grid; i++)
    if (GTK_IS_WIDGET(button_custom_grid[i]))
      gtk_widget_destroy(button_custom_grid[i]);
  n_custom_grid = 0;
  if (button_custom_grid) g_free(button_custom_grid);
  button_custom_grid = NULL;

  if (!m)
  {
    gtk_widget_hide(vbox_vibration);
    gtk_widget_hide(vbox_orbitals);
    gtk_widget_hide(geo_expander);
/*    gtk_widget_hide(geo_button);
    gtk_widget_hide(geo_hbox);
    gtk_widget_hide(label_igeo);*/
    gtk_widget_hide(grid_frame_transp);
    gtk_widget_hide(grid_frame_isosurf);
    gtk_widget_hide(all_orbitals_button);
    gtk_widget_hide(grid_expander);
    gtk_widget_hide(epot_button);
/*    gtk_widget_hide(dens_diff_button);*/
    gtk_widget_hide(vib_spectrum_button);
    gtk_widget_hide(vbox_freqs);
/*    gtk_widget_hide(vbox_3Dobject);*/
    gtk_widget_set_sensitive(button_remove_atom, FALSE);
    gtk_widget_set_sensitive(button_change_atom, FALSE);
    gtk_widget_hide(button_change_bond);
    gtk_widget_set_sensitive(button_add_dummy, FALSE);
    gtk_widget_set_sensitive(button_remove_selection, FALSE);
    gtk_widget_set_sensitive(button_mark_H, FALSE);
    gtk_widget_set_sensitive(button_mark_reverse, FALSE);
    gtk_widget_set_sensitive(button_mark_element, FALSE);
    gtk_widget_set_sensitive(button_neighbor, FALSE);
    return;
  }
  if (m->nvibr)
  {
    gtk_widget_show(vbox_vibration);
    gtk_widget_show(vib_spectrum_button);
    gtk_widget_show(vbox_freqs);
  }
  else
  {
    gtk_widget_hide(vbox_vibration);
    gtk_widget_hide(vib_spectrum_button);
    gtk_widget_hide(vbox_freqs);
  }
  if (n_geometries > 1)
  {
    gtk_widget_show(geo_expander);
    gtk_adjustment_set_upper(GTK_ADJUSTMENT(adj_geoms), (gdouble) n_geometries-1);
/*    gtk_widget_show(geo_button);
    gtk_widget_show(geo_hbox);
    gtk_widget_show(label_igeo);*/
  }
  else
  {
    gtk_widget_hide(geo_expander);
/*    gtk_widget_hide(geo_button);
    gtk_widget_hide(geo_hbox);
    gtk_widget_hide(label_igeo);*/
  }
  if (m->nvector + m->ntriangle + m->nsphere + m->nsurf + m->ncells)
  {
    gtk_widget_show(vbox_3Dobject);
  }
  if (m->ngrids)
  {
    gtk_widget_show(grid_frame_transp);
    gtk_widget_show(grid_frame_isosurf);
    gtk_widget_show(all_orbitals_button);
    gtk_widget_show(grid_expander);
/*    gtk_widget_show(dens_diff_button);*/
    gtk_widget_show(vbox_orbitals);

    for(i = 0; i < m->ngrids; i++)
    {
      if (m->grid_type[i] == CUSTOM)
      {
        n_custom_grid++;
        button_custom_grid = (GtkWidget**) g_realloc(button_custom_grid, sizeof(GtkWidget*) * n_custom_grid);
        button_custom_grid[n_custom_grid - 1] = gtk_button_new_with_label(m->grid.titlesArr[i]);
        g_signal_connect(G_OBJECT(button_custom_grid[n_custom_grid - 1]), "clicked", G_CALLBACK(show_custom_grid), GINT_TO_POINTER(i));
        gtk_box_pack_start(GTK_BOX(vbox_orbitals), button_custom_grid[n_custom_grid - 1], FALSE, FALSE, 0);
        gtk_widget_show(button_custom_grid[n_custom_grid - 1]);
      }
    }
  }
  else
  {
    gtk_widget_hide(grid_expander);
    gtk_widget_hide(grid_frame_transp);
    gtk_widget_hide(grid_frame_isosurf);
    gtk_widget_hide(all_orbitals_button);
    gtk_widget_hide(vbox_orbitals);
  }
/*  if (m->ngrids) gtk_widget_show(vbox_orbitals);
  else gtk_widget_hide(vbox_orbitals);*/
  if (m->ishow & HAS_E_POT) gtk_widget_show(epot_button);
  else gtk_widget_hide(epot_button);
  load_grid_data();
}

void luscus_gtk_update_upon_select_or_mark(void)
{
  int i;
  gchar *tmp;

  if (m->n_selected == 0)
  {
    gtk_widget_set_sensitive(button_change_atom, FALSE);
    gtk_widget_hide(button_change_bond);
    gtk_widget_hide(button_change_angle);
    gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_angle), -1);
    gtk_widget_hide(hbox_xyz);

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), TRUE);
    luscus_gtk_pop_message_from_statusbar2();
    gtk_widget_set_sensitive(GTK_WIDGET(button_average_symm), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_trans_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_trans_symm), "select two atoms that define the translation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_transl_vect), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_c2), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_c2), "select two atoms that define the rotation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_axis), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_inversion_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_inversion_symm), "select an atom that defines the center of inversion");
    gtk_widget_set_sensitive(GTK_WIDGET(button_mirror_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_mirror_symm), "select three atoms that define mirror plane");

    gtk_widget_set_sensitive(GTK_WIDGET(button_add_atom), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_remove_selection), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_add_dummy), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_H), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_reverse), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_element), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_neighbor), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_watch), FALSE);
    gtk_button_set_label(GTK_BUTTON(button_symmetry), "apply symmetry");

    gtk_widget_set_sensitive(GTK_WIDGET(spin_bond), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_angle), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_tors), FALSE);

/*    gtk_label_set_text(GTK_LABEL(label_bond_angle_torsion), NULL);
    gtk_label_set_text(GTK_LABEL(label_unit), NULL);
    gtk_widget_hide(spin_bond_ang_tor);*/

    gtk_widget_set_sensitive(GTK_WIDGET(button_arrow), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_arrow), "select two atoms that define vector");
    gtk_widget_set_sensitive(GTK_WIDGET(button_sphere), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_sphere), "select one or two that defines center (and radius) of sphere");
    gtk_widget_set_sensitive(GTK_WIDGET(button_plain), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_plain), "select three atoms that define plain");
    gtk_widget_set_sensitive(GTK_WIDGET(button_triangle), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_triangle), "select three atoms that define triangle");
    gtk_widget_set_sensitive(GTK_WIDGET(button_cell), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_cell), "select four atoms that define a cell");
  }
  else if (m->n_selected == 1)
  {
    gtk_widget_set_sensitive(button_change_atom, TRUE);
    gtk_widget_hide(button_change_bond);
    gtk_widget_hide(button_change_angle);
    gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_angle), -1);
    gtk_widget_show(hbox_xyz);

    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_x), m->xyz[m->selected[0]][0]);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_y), m->xyz[m->selected[0]][1]);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_z), m->xyz[m->selected[0]][2]);

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_average_symm), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_trans_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_trans_symm), "select two atoms that define the translation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_transl_vect), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_c2), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_c2), "select two atoms that define the rotation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_axis), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_inversion_symm), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_inversion_symm), "apply inversion symmetry");
    gtk_widget_set_sensitive(GTK_WIDGET(button_mirror_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_mirror_symm), "select three atoms that define mirror plane");

    gtk_widget_set_sensitive(GTK_WIDGET(button_add_atom), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_remove_selection), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_add_dummy), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_H), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_reverse), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_element), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_neighbor), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_watch), FALSE);

    gtk_widget_set_sensitive(GTK_WIDGET(spin_bond), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_angle), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_tors), FALSE);

/*    gtk_label_set_text(GTK_LABEL(label_bond_angle_torsion), "atom: ");
    gtk_label_set_text(GTK_LABEL(label_unit), m->elem[m->selected[0]].name);
    gtk_widget_hide(spin_bond_ang_tor);*/

    gtk_widget_set_sensitive(GTK_WIDGET(button_arrow), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_arrow), "select two atoms that define vector");
    gtk_widget_set_sensitive(GTK_WIDGET(button_sphere), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_sphere), "draw sphere around selected atom");
    gtk_widget_set_sensitive(GTK_WIDGET(button_plain), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_plain), "select three atoms that define plain");
    gtk_widget_set_sensitive(GTK_WIDGET(button_triangle), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_triangle), "select three atoms that define triangle");
    gtk_widget_set_sensitive(GTK_WIDGET(button_cell), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_cell), "select four atoms that define a cell");    

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);
    tmp = g_strdup_printf("%s %d (%f %f %f)", m->elem[m->selected[0]].name, m->selected[0]+1, m->xyz[m->selected[0]][0],  m->xyz[m->selected[0]][1], m->xyz[m->selected[0]][2]);
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2(tmp);
    g_free(tmp);
  }
  else if (m->n_selected == 2)
  {
    gtk_widget_set_sensitive(button_change_atom, FALSE);
    gtk_widget_show(button_change_bond);
    gtk_widget_hide(button_change_angle);
    gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_angle), -1);
    gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_bond), get_bond_type(m->selected[0], m->selected[1]));
    gtk_adjustment_set_value(GTK_ADJUSTMENT(translation_magnitude), get_selected_bond_length());  /*OPTION: use default value from Input_Data*/
    gtk_widget_hide(hbox_xyz);

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_average_symm), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_trans_symm), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_trans_symm), "Translate molecule in the direction of selected atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_transl_vect), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_c2), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_c2), "apply Cn symmetry operation around axis defined by selected atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_axis), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_inversion_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_inversion_symm), "select an atom that defines the center of inversion");
    gtk_widget_set_sensitive(GTK_WIDGET(button_mirror_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_mirror_symm), "select three atoms that define mirror plane");

    gtk_widget_set_sensitive(GTK_WIDGET(button_add_atom), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_remove_selection), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_H), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_reverse), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_element), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_neighbor), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_watch), TRUE);

    gtk_widget_set_sensitive(GTK_WIDGET(spin_bond), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_angle), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_tors), FALSE);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_bond), get_selected_bond_length());

/*    gtk_label_set_text(GTK_LABEL(label_bond_angle_torsion), "bond");
    gtk_label_set_text(GTK_LABEL(label_unit), "\303\205");
    gtk_adjustment_set_upper(GTK_ADJUSTMENT(adj_bond_angle_tors), 100.F);
    gtk_adjustment_set_lower(GTK_ADJUSTMENT(adj_bond_angle_tors), 0.01F);
    gtk_adjustment_set_step_increment(GTK_ADJUSTMENT(adj_bond_angle_tors),0.01F);
    gtk_adjustment_set_page_increment(GTK_ADJUSTMENT(adj_bond_angle_tors),0.1F);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_bond_angle_tors), get_selected_bond_length());
    gtk_widget_show(spin_bond_ang_tor);*/
    gtk_widget_set_sensitive(GTK_WIDGET(button_add_dummy), TRUE);

    gtk_widget_set_sensitive(GTK_WIDGET(button_arrow), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_arrow), "draw vector between selecter atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_sphere), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_sphere), "draw sphere around selected atom");
    gtk_widget_set_sensitive(GTK_WIDGET(button_plain), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_plain), "select three atoms that define plain");
    gtk_widget_set_sensitive(GTK_WIDGET(button_triangle), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_triangle), "select three atoms that define triangle");
    gtk_widget_set_sensitive(GTK_WIDGET(button_cell), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_cell), "select four atoms that define a cell");

/*    if (m->n_selected != m->natom)
    {
      gtk_button_set_label(GTK_BUTTON(button_symmetry), "apply translation symmetry");
    }
    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);*/

    tmp = g_strdup_printf("bond: %s %d-%s %d  %f", m->elem[m->selected[0]].name, m->selected[0]+1,
                                                   m->elem[m->selected[1]].name, m->selected[1]+1, get_selected_bond_length());
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2(tmp);
    g_free(tmp);
  }
  else if (m->n_selected == 3)
  {
    gtk_widget_set_sensitive(button_change_atom, FALSE);
    gtk_widget_hide(button_change_bond);
    gtk_widget_show(button_change_angle);
    gtk_widget_hide(hbox_xyz);

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_average_symm), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_trans_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_trans_symm), "select two atoms that define the translation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_transl_vect), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_c2), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_c2), "select two atoms that define the rotation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_axis), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_inversion_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_inversion_symm), "select an atom that defines the center of inversion");
    gtk_widget_set_sensitive(GTK_WIDGET(button_mirror_symm), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_mirror_symm), "apply mirror symmetry");

    gtk_widget_set_sensitive(GTK_WIDGET(button_add_atom), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_remove_selection), TRUE);

    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_H), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_reverse), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_element), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_neighbor), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_watch), TRUE);


    gtk_widget_set_sensitive(GTK_WIDGET(spin_bond), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_angle), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_tors), FALSE);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_angle), get_selected_angle_value());

/*    gtk_label_set_text(GTK_LABEL(label_bond_angle_torsion), "angle");
    gtk_label_set_text(GTK_LABEL(label_unit), "\313\232");
    gtk_adjustment_set_upper(GTK_ADJUSTMENT(adj_bond_angle_tors), 180.F);
    gtk_adjustment_set_lower(GTK_ADJUSTMENT(adj_bond_angle_tors), 0.01F);
    gtk_adjustment_set_step_increment(GTK_ADJUSTMENT(adj_bond_angle_tors),0.2F);
    gtk_adjustment_set_page_increment(GTK_ADJUSTMENT(adj_bond_angle_tors),1.0F);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_bond_angle_tors), get_selected_angle_value());*/

    gtk_widget_set_sensitive(GTK_WIDGET(button_arrow), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_arrow), "draw vector between selecter atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_sphere), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_sphere), "select atom that defines center of sphere");
    gtk_widget_set_sensitive(GTK_WIDGET(button_plain), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_plain), "draw plain through selected atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_triangle), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_triangle), "draw triangle through selected atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_cell), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_cell), "select four atoms that define a cell");    

/*    if (m->n_selected != m->natom)
    {
      gtk_widget_show(spin_bond_ang_tor);
      gtk_button_set_label(GTK_BUTTON(button_symmetry), "apply mirror symmetry");
    }
    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), TRUE);*/

    tmp = g_strdup_printf("angle: %s %d-%s %d-%s %d  %f", m->elem[m->selected[0]].name, m->selected[0]+1,
                                                          m->elem[m->selected[1]].name, m->selected[1]+1,
                                                          m->elem[m->selected[2]].name, m->selected[2]+1,
                                                          get_selected_angle_value());
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2(tmp);
    g_free(tmp);
  }
  else if (m->n_selected == 4)
  {
    gtk_widget_set_sensitive(button_change_atom, FALSE);
    gtk_widget_hide(button_change_bond);
    gtk_widget_hide(button_change_angle);
    gtk_combo_box_set_active(GTK_COMBO_BOX(button_change_angle), -1);    
    gtk_widget_hide(hbox_xyz);

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_average_symm), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_trans_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_trans_symm), "select two atoms that define the translation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_transl_vect), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_c2), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_c2), "select two atoms that define the rotation vector");
    gtk_widget_set_sensitive(GTK_WIDGET(spin_axis), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_inversion_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_inversion_symm), "select an atom that defines the center of inversion");
    gtk_widget_set_sensitive(GTK_WIDGET(button_mirror_symm), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_mirror_symm), "select three atoms that define mirror plane");

    gtk_widget_set_sensitive(GTK_WIDGET(button_add_atom), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_remove_selection), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_H), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_reverse), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_mark_element), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_neighbor), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_watch), TRUE);

    gtk_widget_set_sensitive(GTK_WIDGET(spin_bond), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_angle), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(spin_tors), TRUE);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_tors), get_selected_dihedral_value());

/*    gtk_label_set_text(GTK_LABEL(label_bond_angle_torsion), "dihedral");
    gtk_label_set_text(GTK_LABEL(label_unit), "\313\232");
    gtk_adjustment_set_upper(GTK_ADJUSTMENT(adj_bond_angle_tors), 360.F);
    gtk_adjustment_set_lower(GTK_ADJUSTMENT(adj_bond_angle_tors),-360.F);
    gtk_adjustment_set_step_increment(GTK_ADJUSTMENT(adj_bond_angle_tors),0.2F);
    gtk_adjustment_set_page_increment(GTK_ADJUSTMENT(adj_bond_angle_tors),1.0F);
    gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_bond_angle_tors), get_selected_dihedral_value());
    gtk_widget_show(spin_bond_ang_tor);*/
    gtk_widget_set_sensitive(GTK_WIDGET(button_add_dummy), FALSE);

    gtk_widget_set_sensitive(GTK_WIDGET(button_arrow), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_arrow), "draw vector between selecter atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_sphere), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_sphere), "select atom that defines center of sphere");
    gtk_widget_set_sensitive(GTK_WIDGET(button_plain), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_plain), "draw plain through selected atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_triangle), FALSE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_triangle), "draw triangle through selected atoms");
    gtk_widget_set_sensitive(GTK_WIDGET(button_cell), TRUE);
    gtk_widget_set_tooltip_text(GTK_WIDGET(button_cell), "draw cell through selected atoms");

    gtk_widget_set_sensitive(GTK_WIDGET(button_symmetry), FALSE);

    tmp = g_strdup_printf("dihedral:%s %d-%s %d-%s %d-%s %d  %f",
                          m->elem[m->selected[0]].name, m->selected[0]+1,
                          m->elem[m->selected[1]].name, m->selected[1]+1,
                          m->elem[m->selected[2]].name, m->selected[2]+1, 
                          m->elem[m->selected[3]].name, m->selected[3]+1, get_selected_dihedral_value());
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2(tmp);
    g_free(tmp);
  }
  if (m->n_marked)
  {
    gtk_widget_set_sensitive(GTK_WIDGET(button_unmark), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_renumber), TRUE);
  }
  else
  {
    gtk_widget_set_sensitive(GTK_WIDGET(button_unmark), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(button_renumber), FALSE);
  }

  if (m->n_selected == 1 || m->n_selected == 2 || m->n_marked) gtk_widget_set_sensitive(GTK_WIDGET(button_remove_atom), TRUE);
  else gtk_widget_set_sensitive(GTK_WIDGET(button_remove_atom), FALSE);

  if (m->nvector || m->ntriangle || m->nsphere || m->nsurf || m->ncells)
    gtk_widget_set_sensitive(GTK_WIDGET(button_clear_drawed), TRUE);
  else
    gtk_widget_set_sensitive(GTK_WIDGET(button_clear_drawed), FALSE);

  /**/
  if (n_watched)
  {
    gtk_widget_set_sensitive(GTK_WIDGET(button_unwatch), TRUE);
    gtk_widget_show(vbox_watched);

    if (n_watched != n_label_watched)
    {
      for(i = 0; i < n_label_watched; i++)
       	gtk_widget_destroy(labels_watched[i]);
      g_free(labels_watched);
      labels_watched = (GtkWidget**) g_malloc(sizeof(GtkWidget*) * n_watched);
      for(i = 0; i < n_watched; i++)
      {
        switch(watch[i].coord_type)
        {
          case C_BOND:
            tmp = g_strdup_printf("bond: %d-%d: %f", watch[i].atom[0], watch[i].atom[1],
                  get_bond_length(watch[i].atom[0] - 1, watch[i].atom[1] - 1));
          break;
          case C_ANGLE:
            tmp = g_strdup_printf("angle: %d-%d-%d: %f", watch[i].atom[0],
                  watch[i].atom[1], watch[i].atom[2],
                  get_angle_value(watch[i].atom[0] - 1, watch[i].atom[1]-1, watch[i].atom[2]-1));
          break;
          case C_DIHEDRAL:
            tmp = g_strdup_printf("dihedral: %d-%d-%d-%d: %f", watch[i].atom[0],
            watch[i].atom[1], watch[i].atom[2], watch[i].atom[3],
            get_dihedral_value(watch[i].atom[0]-1, watch[i].atom[1]-1, watch[i].atom[2]-1, watch[i].atom[3]-1));
          break;
          default:
          tmp = NULL;
        }
        labels_watched[i] = gtk_label_new(tmp);
        gtk_box_pack_start(GTK_BOX(vbox_watched), labels_watched[i], FALSE, FALSE, 0);
        gtk_widget_show(labels_watched[i]);
        g_free(tmp);
      }
      n_label_watched = n_watched;
    }
  }
  else
  {
    gtk_widget_set_sensitive(GTK_WIDGET(button_unwatch), FALSE);
    if (n_label_watched)
    {
      for(i = 0; i < n_label_watched; i++)
       	gtk_widget_destroy(labels_watched[i]);
      g_free(labels_watched);
      labels_watched = NULL;
      n_label_watched = 0;
    }
  }

/*GtkWidget *vbox_watched;
int n_label_watched;
GtkWidget **labels_watched;*/

}

void change_watched_data(void)
{
  int i;
  char *tmp;
  if (n_label_watched != n_watched) return;
  for(i = 0; i < n_label_watched; i++)
  {
    switch(watch[i].coord_type)
    {
      case C_BOND:
        tmp = g_strdup_printf("bond: %d-%d: %f", watch[i].atom[0], watch[i].atom[1],
              get_bond_length(watch[i].atom[0]-1, watch[i].atom[1]-1));
      break;
      case C_ANGLE:
        tmp = g_strdup_printf("angle: %d-%d-%d: %f", watch[i].atom[0],
              watch[i].atom[1], watch[i].atom[2],
              get_angle_value(watch[i].atom[0]-1, watch[i].atom[1]-1, watch[i].atom[2]-1));
      break;
      case C_DIHEDRAL:
        tmp = g_strdup_printf("dihedral: %d-%d-%d-%d: %f", watch[i].atom[0],
        watch[i].atom[1], watch[i].atom[2], watch[i].atom[3],
        get_dihedral_value(watch[i].atom[0]-1, watch[i].atom[1]-1, watch[i].atom[2]-1, watch[i].atom[3]-1));
      break;
      default:
      tmp = NULL;
    }
    gtk_label_set_text(GTK_LABEL(labels_watched[i]), tmp);
    g_free(tmp);
  }
}

void luscus_gtk_update_geo_info(void)
{
  gint i;
  gchar *textbuf;
  gchar *tmpbuf;
  gchar *txt_freq = NULL;
  gchar *txt_irint = NULL;
  gchar *txt_ramint = NULL;

  if (m->ishow & HAS_ENERGY && m->ishow & HAS_RMS_GRAD && m->ishow & HAS_MAX_GRAD)
    textbuf = g_strdup_printf("Geometry: %d\nEnergy: %f\nRMS grad: %f\nmax. grad: %f", igeo, m->geo_energy, m->rms_grad, m->max_grad);
  else if (m->ishow & HAS_ENERGY && m->ishow & HAS_RMS_GRAD) 
    textbuf = g_strdup_printf("Geometry: %d\nEnergy: %f\nRMS grad: %f", igeo, m->geo_energy, m->rms_grad);
  else if (m->ishow & HAS_ENERGY && m->ishow & HAS_MAX_GRAD)
    textbuf = g_strdup_printf("Geometry: %d\nEnergy: %f\nmax. grad: %f", igeo, m->geo_energy, m->max_grad);
  else if (m->ishow & HAS_RMS_GRAD && m->ishow & HAS_MAX_GRAD)
    textbuf = g_strdup_printf("Geometry: %d\nRMS grad: %f\nmax. grad: %f", igeo, m->rms_grad, m->max_grad);
  else if (m->ishow & HAS_ENERGY) textbuf = g_strdup_printf("Geometry: %d\nEnergy: %f", igeo, m->geo_energy);
  else if (m->ishow & HAS_RMS_GRAD) textbuf = g_strdup_printf("Geometry: %d\nRMS grad: %f", igeo, m->rms_grad);
  else if (m->ishow & HAS_MAX_GRAD) textbuf = g_strdup_printf("Geometry: %d\nmax. grad: %f", igeo, m->max_grad);
  else textbuf = g_strdup_printf("Geometry: %d", igeo);
  gtk_label_set_text(GTK_LABEL(label_igeo), textbuf);
  g_free(textbuf);

  /*delete buttons from the last geom*/
  deallocate_vib_buttons();

  n_vibs = m->nvibr;
  freq_buttons = (GtkWidget**) g_malloc(sizeof(GtkWidget*) * m->nvibr);
  freq_labels = (GtkWidget**) g_malloc(sizeof(GtkWidget*) * m->nvibr);

  for(i = 0; i < m->nvibr; i++)
  {
    /*this allocation/deallocation of strings is bizarre!*/
    txt_freq = g_markup_printf_escaped("%.2f cm<sup>-1</sup>", m->freq[i]);
    textbuf = g_strdup(txt_freq);
    
    if (m->ishow & HAS_IRINT)
    {
      txt_irint = g_markup_printf_escaped(" IR int.=%.2f", m->ir_intensity[i]);
      g_free(textbuf);
      textbuf = g_strconcat(txt_freq, txt_irint, NULL);
      g_free(txt_irint);
    }
    if (m->ishow & HAS_RAMAN)
    {
      txt_ramint = g_markup_printf_escaped(" Raman int.=%.2f", m->raman_intensity[i]);
      tmpbuf = g_strconcat(textbuf, txt_ramint, NULL);
      g_free(textbuf);
      textbuf = tmpbuf;
/*      textbuf = g_strconcat(txt_freq, txt_irint, NULL);*/
    }

    freq_labels[i] = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(freq_labels[i]), textbuf);
    freq_buttons[i] = gtk_button_new();
    gtk_container_add(GTK_CONTAINER(freq_buttons[i]), freq_labels[i]);
    gtk_box_pack_start(GTK_BOX(vbox_freqs), freq_buttons[i], FALSE, FALSE, 0);
    g_signal_connect(G_OBJECT(freq_buttons[i]), "clicked", G_CALLBACK(luscus_gtk_goto_freq), GINT_TO_POINTER(i));

    gtk_widget_show(freq_buttons[i]);
    gtk_widget_show(freq_labels[i]);
    g_free(txt_freq);
    g_free(textbuf);
  }
  /*if vibrations exists, start frequency timer*/
  if (m->nvibr) luscus_gtk_freq_timer();

  if (m->editable) luscus_gtk_enable_editing();
  else luscus_gtk_diseable_editing();
}

void deallocate_vib_buttons(void)
{
  int i;
  /*delete buttons from the last geom*/

  if (freq_labels)
    for(i = 0; i < n_vibs; i++)
      gtk_widget_destroy(freq_labels[i]);

  if (freq_buttons)
    for(i = 0; i < n_vibs; i++)
      gtk_widget_destroy(freq_buttons[i]);

  g_free(freq_buttons);
  g_free(freq_labels);
  freq_buttons = NULL;
  freq_labels = NULL;
}

/*ALLOCATE ADDITIONAL MEMORY FOR NORMAL MODES*/

void load_grid_data(void)
{
  int i;
  int nnc_grid = 0;
  GtkTreeIter iter;

  if (m->ngrids == 0) return;
  for(i = 0; i < m->ngrids; i++)
  {
    if (m->grid_type[i] == ORBITAL)
    {
      nnc_grid++;
      gtk_list_store_append(store, &iter);

      if (i == iorb) selected_orbital = &iter;

      gtk_list_store_set(store, &iter,
                         ORBITAL_TYPE, orbital_description[m->orbital_type[i]],
                         ORBITAL_SYMM, m->grid_symmetry[i],
                         ORBITAL_NUM, m->grid_index[i],
                         ORBITAL_ENERG, m->grid_energy[i],
                         ORBITAL_OCC, m->grid_occ[i],
                         ORBITAL_EDITED, " ",
                        -1);
    }
  }

  if (nnc_grid)
  {
    gtk_tree_view_set_model(GTK_TREE_VIEW(treeview), GTK_TREE_MODEL(store));
    luscus_gtk_search_orbital_in_list(m->grid_symmetry[iorb], m->grid_index[iorb]);
  }
/*  gtk_tree_selection_select_iter(selection, selected_orbital);*/
}

void delete_all_orbitals_from_list(void)
{
  GtkTreeIter iter;
  GtkTreeModel *model;
  gboolean next = TRUE;

  if (!m) return;
  if (m->ngrids <= 0) return;

  model = gtk_tree_view_get_model(GTK_TREE_VIEW(treeview));
  gtk_tree_view_set_model(GTK_TREE_VIEW(treeview), NULL);

/*  while(next)
  {
    printf("REMOVING ");
    next = gtk_tree_model_get_iter_first(model, &iter);
    if (next)
    {
      gboolean irem;
      irem = gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
      printf("REMOVING ITER irem = %d!\n", irem);
    }
  }*/
  gtk_list_store_clear(GTK_LIST_STORE(model));
}

void delete_orbital_from_list(GtkWidget *button, gpointer data)
{
  GtkTreeIter iter;
  GtkTreeModel *model;

  if (gtk_tree_selection_get_selected(selection, &model, &iter))
    gtk_list_store_remove(store, &iter);
}

void undelete_all_orbitals_from_list(GtkWidget *widget, gpointer data)
{
  GtkTreeIter iter;
  GtkTreeModel *model;
  gboolean next = TRUE;
  gint i;
/*  gint data_sym, data_orbnum;
  gdouble data_energy, data_occu;
  gint orb_num_type;*/

    /*detach model from the treeview*/

  model = gtk_tree_view_get_model(GTK_TREE_VIEW(treeview));
  g_object_ref(model);
  gtk_tree_view_set_model(GTK_TREE_VIEW(treeview), NULL);

    /*remove all iters (data)*/

  while(next)
  {
    next = gtk_tree_model_get_iter_first(model, &iter);
    if (next) gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
  }

  /*add new data*/

  selected_orbital = NULL;
  for (i = 0; i < m->ngrids; i++)
  {
    if (m->grid_type[i] == ORBITAL)
    {
    /*acquireing tree iter*/

      gtk_list_store_append(store, &iter);
      if (i == iorb)
      {
       	selected_orbital = &iter;
      }

      if (m->edited[i])
      {
        gtk_list_store_set(store, &iter,
                           ORBITAL_TYPE, orbital_description[m->orbital_type[i]],
                           ORBITAL_SYMM, m->grid_symmetry[i],
                           ORBITAL_NUM, m->grid_index[i],
                           ORBITAL_ENERG, m->grid_energy[i],
                           ORBITAL_OCC, m->grid_occ[i],
                           ORBITAL_EDITED, "*",
                          -1);
      }
      else
      {
        gtk_list_store_set(store, &iter,
                           ORBITAL_TYPE, orbital_description[m->orbital_type[i]],
                           ORBITAL_SYMM, m->grid_symmetry[i],
                           ORBITAL_NUM, m->grid_index[i],
                           ORBITAL_ENERG, m->grid_energy[i],
                           ORBITAL_OCC, m->grid_occ[i],
                           ORBITAL_EDITED, " ",
                          -1);
      }
    }
  }

  /*reattach the model to the treeview*/

  gtk_tree_view_set_model(GTK_TREE_VIEW(treeview), model);
  g_object_unref(model);
}

gint luscus_gtk_change_orbital_type(gint sym, gint num)
{
  gint i;
  gint data_sym = 0, data_num = 0, data_type;
  /*changes orbital type return new orbital type*/
  GtkTreeModel *model;
  GtkTreeIter iter;

  gchar *tmptext;
  gboolean next = TRUE;

  model = GTK_TREE_MODEL(store);
  gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

/*  gint iiorb = 0;*/

  while(next)
  {
    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, ORBITAL_SYMM, &data_sym, -1);
    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, ORBITAL_NUM, &data_num, -1);

    if (next && sym == data_sym && num == data_num) next = FALSE;
    else
    {
      next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
/*      iiorb++;*/
    }
  }

  gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, ORBITAL_TYPE, &tmptext, -1);

  for(data_type = 0; g_strcmp0(tmptext, orbital_description[data_type]) != 0; data_type++);

  data_type++;
  if (data_type > 7) data_type = 1;


  gtk_list_store_set(store, &iter, ORBITAL_TYPE, orbital_description[data_type], -1); 
  gtk_list_store_set(store, &iter, ORBITAL_TYPE, orbital_description[data_type], ORBITAL_EDITED, "*", -1); 

  for(i = 0; i < m->ngrids; i++)
    if (m->grid_type[i] == ORBITAL)
      if (m->grid_symmetry[i] == sym && m->grid_index[i] == num)
      {
        m->orbital_type[i] = data_type;
        m->edited[i] = 1;
        luscus_gtk_multiview_redraw_list(i, data_type);
      }

  luscus_show_substaces();
  return data_type;
}

void orbital_combo_edited(GtkCellRendererText *cell, gchar *path_string, gchar *new_text, GtkListStore *list_store_combo)
{
  GtkTreeModel *model;
  GtkTreeIter iter;
  gboolean obtain_iter;
  gint i;
  gint osym, oind;
  gint oldtype, newtype;

  model = GTK_TREE_MODEL(store);
  obtain_iter = gtk_tree_model_get_iter_from_string(model, &iter, path_string);

  gtk_tree_model_get(model, &iter, ORBITAL_SYMM, &osym, ORBITAL_NUM, &oind, -1);
  for(i = 0; i < m->ngrids && m->grid_symmetry[i] != osym || m->grid_index[i] != oind; i++);
  oldtype = m->grid_type[i];
  m->edited[i] = 1;

  newtype = 0;
  while(g_strcmp0(new_text, orbital_description[newtype])) newtype++;
  m->orbital_type[i] = newtype;
 
  gtk_list_store_set(store, &iter, ORBITAL_TYPE, new_text, ORBITAL_EDITED, "*", -1);

  m->grid.iniIndex[oldtype]--;
  m->grid.iniIndex[newtype]++;
  luscus_show_substaces();
}

void change_orbital_type(int itype)
{
  GtkTreeModel *model;
  GtkTreeIter iter;
  int isym, iindex;
  int iorb;

  if (gtk_tree_selection_get_selected(selection, &model, &iter))
  {
    gtk_list_store_set(store, &iter, ORBITAL_TYPE, orbital_description[itype], ORBITAL_EDITED, "*", -1); 
    gtk_tree_model_get(model, &iter, ORBITAL_SYMM, &isym, ORBITAL_NUM, &iindex, -1);
    iorb = get_orbital(isym, iindex);
    if (iorb >= 0) m->orbital_type[iorb] = itype;
    m->edited[iorb] = 1;
    luscus_gtk_multiview_redraw_list(iorb, itype);
  }
  else printf("NO orbital selected!\n");
  luscus_show_substaces();
}

void luscus_show_substaces(void)
{
  int i;
  int test[NUM_ORB_TYPES];
  gchar tmp[512];
  gchar *text;
  tmp[0] = 0;
  g_strlcat(tmp, "Subspaces:", 512);
  for(i = 0; i < NUM_ORB_TYPES; i++) test[i] = 0;
  for(i = 0; i < m->ngrids; i++)
    if (strstr(m->grid.titlesArr[i], "Orbital"))
      test[m->orbital_type[i]]++;

  for(i = 0; i < NUM_ORB_TYPES; i++)
    if (test[i])
    {
      text = g_strdup_printf(" %s:%d;", orbital_description_short[i], test[i]);
      g_strlcat(tmp, text, 512);
      g_free(text);
    }

  luscus_gtk_pop_message_from_statusbar2();
  luscus_gtk_push_message_to_statusbar2(tmp);
}

void orbital_filter_callback(GtkWidget *entry, gpointer data)
{
  gchar *string;
  int isym, inum;
  string = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));


  if (string[0] != '#') gtk_entry_set_text(GTK_ENTRY(entry), "ERROR: wrong filter format!");
  else
  {
    if (string[1] == ':')
    {
      luscus_gtk_get_filter_str(string);
/*      get_filter_str(string);
      validate_filter();*/
    }
    else
    {
      sscanf(string+1, "%d %d", &isym, &inum);
/*      printf("searching orbital %d %d\n", isym, inum);*/
      if (luscus_gtk_find_selected_orbital(isym, inum)) gtk_entry_set_text(GTK_ENTRY(entry), "ERROR: No such orbital!");
      luscus_gtk_search_orbital_in_list(isym, inum);
    }
  }

  g_free(string);
}

void luscus_gtk_get_filter_str(gchar *string) /*this function is replacement for the "get_filter_str"*/
{                                         /*this function operates on the GTK storage of orbitals*/
  gint i = 1, j, k;
  gint slen = strlen((char*) string);
  gchar *tmp = NULL;
  gboolean filter_sym = FALSE;
  gboolean filter_occu = FALSE;
  gboolean filter_index = FALSE;
  gboolean filter_energ = FALSE;
  gchar testchar;
  gchar *use_sym = NULL;
  gdouble min_energ, max_energ, min_occu, max_occu;
  gchar *use_index = NULL;
/*  GtkTreeModel *model;*/
  GtkTreeIter iter, olditer;
  gboolean next;
  gboolean do_erase;

  gint int_data;
  gdouble double_data;
  gchar *char_data;

  use_sym = g_malloc(sizeof(gchar));
  use_sym[0] = 0;
  use_index = g_malloc(sizeof(gchar));
  use_index[0] = 0;

  while(i < slen)
  {
    while(string[i] != 's' && string[i] != 'o' && string[i] != 'e' && string[i] != 'i' && i < slen) i++;
    if (i == slen) return;

    testchar = string[i];
    for(j = 1; string[i+j] != ';' && i+j <= slen; j++)
    {
      tmp = g_realloc(tmp, (j+1) * sizeof(gchar));
      tmp[j-1] = string[i+j];
      tmp[j] = 0;
    }
    i += j;

    switch(testchar)
    {
      case 's':
        filter_sym = TRUE;
        k = 1;
        for(j = 0; tmp[j] != 0; j++)
          if (g_ascii_isdigit(tmp[j]))
          {
            use_sym = g_realloc(use_sym, ++k);
            use_sym[k-2] = tmp[j];
            use_sym[k-1] = 0;
          }
        break;
      case 'o':
        filter_occu = TRUE;
        min_occu = g_ascii_strtod(tmp, NULL);
/*      max_occu = g_ascii_strtod(index(tmp,':')+1, NULL);*/
        max_occu = g_ascii_strtod(strstr(tmp,":")+1, NULL);
        break;
      case 'e':
        filter_energ = TRUE;
        min_energ = g_ascii_strtod(tmp, NULL);
/*      max_energ = g_ascii_strtod(index(tmp,':')+1, NULL);*/
        max_energ = g_ascii_strtod(strstr(tmp,":")+1, NULL);
        break;
      case 'i':
        filter_index = TRUE;
        k = 1;
        for(j = 0; tmp[j] != 0; j++)
          if (tmp[j] == 'F' || tmp[j] == 'I' || tmp[j] == '1' || tmp[j] == '2' || tmp[j] == '3' || tmp[j] == 'S' || tmp[j] == 'D' || tmp[j] == 'U')
          {
            use_index = g_realloc(use_index, ++k);
            use_index[k-2] = tmp[j];
            use_index[k-1] = 0;
          }
        break;
      default: break;
    }
    g_free(tmp);
    tmp = NULL;
  }

  /*Now delete orbitals that doesn't match criteria*/

  if (filter_sym) 
  {
    if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) return;
    next = TRUE;

    while(next)
    {
      olditer = iter;
      next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter); /*first move to the next list element in order not to delete current element (current iter might be lost)*/
      
      gtk_tree_model_get(GTK_TREE_MODEL(store), &olditer, ORBITAL_SYMM, &int_data, -1);
      do_erase = TRUE;

      for(i = 0; use_sym[i] != 0; i++)
      {
        if (use_sym[i] - 48 == int_data) /* this are ASCII representations! */
          do_erase = FALSE;
      }

      if (do_erase) gtk_list_store_remove(store, &olditer); 
    }
  }
  if (filter_occu)
  {
    if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) return;
    next = TRUE;

    while(next)
    {
      olditer = iter;
      next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);

      gtk_tree_model_get(GTK_TREE_MODEL(store), &olditer, ORBITAL_OCC, &double_data, -1);

      if (double_data < max_occu && double_data > min_occu) do_erase = FALSE;
      else do_erase = TRUE;

      if (do_erase) gtk_list_store_remove(store, &olditer);
    }
  }
  if (filter_index)
  {
    if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) return;
    next = TRUE;

    while(next)
    {
      olditer = iter;
      next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);

      gtk_tree_model_get(GTK_TREE_MODEL(store), &olditer, ORBITAL_TYPE, &char_data, -1);
      do_erase = TRUE;

      for(i = 0; use_index[i] != 0; i++)
      {
        switch(use_index[i])
	{
	  case 'F': if (g_strcmp0(char_data, "frozen") == 0) do_erase = FALSE; break;
	  case 'I': if (g_strcmp0(char_data, "inactive") == 0) do_erase = FALSE; break;
	  case '1': if (g_strcmp0(char_data, "RAS 1") == 0) do_erase = FALSE; break;
	  case '2': if (g_strcmp0(char_data, "RAS 2") == 0) do_erase = FALSE; break;
	  case '3': if (g_strcmp0(char_data, "RAS 3") == 0) do_erase = FALSE; break;
	  case 'S': if (g_strcmp0(char_data, "secondary") == 0) do_erase = FALSE; break;
	  case 'D': if (g_strcmp0(char_data, "deleted") == 0) do_erase = FALSE; break;
	  case 'U': if (g_strcmp0(char_data, "unknown") == 0) do_erase = FALSE; break;
	}
      }

      if (do_erase) gtk_list_store_remove(store, &olditer);
    }
  }
  if (filter_energ)
  {
    if (!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter)) return;
    next = TRUE;

    while(next)
    {
      olditer = iter;
      next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);

      gtk_tree_model_get(GTK_TREE_MODEL(store), &olditer, ORBITAL_ENERG, &double_data, -1);

      if (double_data < max_energ && double_data > min_energ) do_erase = FALSE;
      else do_erase = TRUE;

      if (do_erase) gtk_list_store_remove(store, &olditer);
    }
  }

  g_free(use_sym);
  g_free(use_index);
  return;
}

void luscus_gtk_search_orbital_in_list(gint isym, gint inum)
{
  gboolean next = TRUE;
  gint tsym, tnum;
  GtkTreeIter iter;

  gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

  while(next)
  {
    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, ORBITAL_SYMM, &tsym, -1);
    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, ORBITAL_NUM, &tnum, -1);

    if (next && isym == tsym && inum == tnum)
    {
      gtk_tree_selection_select_iter(selection, &iter);
      next = FALSE;
    }	

    next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
  }

}

void luscus_gtk_select_orbital_up(void)
{
  GtkTreeModel *mod;
  GtkTreeIter iter;
  gboolean next;
  char lin[82];
  int orb_symm, orb_index;
  double orb_occ, orb_energ;
  gchar *orb_type, *orb_edit;
  if (!gtk_tree_selection_get_selected(selection, &mod, NULL))
  { /*none selected!*/
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    gtk_tree_selection_select_iter(selection, &iter);
    return;
  }

  if (gtk_tree_selection_get_selected(selection, &mod, &iter))
    next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
  if (next)
    gtk_tree_selection_select_iter(selection, &iter);
  else /*last orbital is selected; goto first*/
  {
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    gtk_tree_selection_select_iter(selection, &iter);
  }
  gtk_tree_model_get(GTK_TREE_MODEL(store), &iter,
                                                  ORBITAL_TYPE, &orb_type,
                                                  ORBITAL_SYMM, &orb_symm, 
                                                  ORBITAL_NUM,  &orb_index,
                                                  ORBITAL_ENERG,&orb_energ,
                                                  ORBITAL_OCC,  &orb_occ,
                                                  ORBITAL_EDITED,&orb_edit,  -1);
  snprintf(lin, 82, "%1s orbital sym: %2d index: %2d energy: %12.4f occ: %5.2f type: %s", orb_edit, orb_symm, orb_index, orb_energ, orb_occ, orb_type);

  /*print orbital data*/
  luscus_gtk_pop_message_from_statusbar1();
  luscus_gtk_push_message_to_statusbar1(lin);
}

void luscus_gtk_select_orbital_down(void)
{
  GtkTreeIter iter, last_iter;
  GtkTreeModel *mod;
  GtkTreePath *path;
  gboolean next;
  char lin[82];
  int orb_symm, orb_index;
  double orb_occ, orb_energ;
  gchar *orb_type, *orb_edit;
 /* gint nrow;*/

  if (gtk_tree_selection_get_selected(selection, &mod, NULL))
  {/*there is selected row*/
    gtk_tree_selection_get_selected(selection, &mod, &iter);
    path = gtk_tree_model_get_path(GTK_TREE_MODEL(store), &iter);
    if (path) next = gtk_tree_path_prev(path);

    if (next)
    {
      gtk_tree_model_get_iter(GTK_TREE_MODEL(store), &iter, path);
      gtk_tree_path_free(path);
    }
    else /*first orbital is selected -> goto last*/
    {
      gtk_tree_path_free(path);

      next = TRUE;
      while(next) /*rewind to last it might be CPU intensive if the list is very large!*/
      {
        last_iter = iter;
        next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
      }
      iter = last_iter;
    }
  }
  else
  { /*select last row*/
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter); /*goto first*/

    next = TRUE;
    while(next) /*rewind to last it might be CPU intensive if the list is very large!*/
    {
      last_iter = iter;
      next = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    }
    iter = last_iter;
  }
  gtk_tree_selection_select_iter(selection, &iter);

  /*print orbital data*/
  gtk_tree_model_get(GTK_TREE_MODEL(store), &iter,
                                                  ORBITAL_TYPE, &orb_type,
                                                  ORBITAL_SYMM, &orb_symm, 
                                                  ORBITAL_NUM,  &orb_index,
                                                  ORBITAL_ENERG,&orb_energ,
                                                  ORBITAL_OCC,  &orb_occ,
                                                  ORBITAL_EDITED,&orb_edit,  -1);
  snprintf(lin, 82, "%1s orbital sym: %2d index: %2d energy: %12.4f occ: %5.2f type: %s", orb_edit, orb_symm, orb_index, orb_energ, orb_occ, orb_type);

  /*print orbital data*/
  luscus_gtk_pop_message_from_statusbar1();
  luscus_gtk_push_message_to_statusbar1(lin);
}

static gboolean luscus_gtk_find_selected_orbital(gint test_sym, gint test_orbnum)
{
  gint iiorb = 0;
/*  gint data_sym = 0, data_orbnum = 0, data_type;
  gdouble data_energy, data_occu;*/
  gboolean next = TRUE;

  while(next)
  {
    if (m->grid_symmetry[iiorb] == test_sym && m->grid_index[iiorb] == test_orbnum) next = FALSE;
    else iiorb++;
    if (iiorb >= m->ngrids) return TRUE;
  }

  iorb = iiorb;
  do_load_grid();

  rerender_3d();

  redraw();
  return FALSE;
}

void show_custom_grid(GtkWidget *button, gpointer data)
{
  iorb = GPOINTER_TO_INT(data);
  gtk_tree_selection_unselect_all(selection);
  do_load_grid();
  rerender_3d();
  redraw();
}

int get_orbital(int isym, int iindex)
{
  int i;

  for(i = 0; i < m->ngrids; i++)
    if (isym == m->grid_symmetry[i] && iindex == m->grid_index[i]) return i;

  /*orbital is not found!*/
  return -1;
}

void luscus_gtk_update_3Dobject_info(void)
{
  int i, j;
  gchar buff1[20];
  GdkColor color;

  /*destroy all gtk widgets that describe 3D objects*/
  for (i = 0; i < n_3D_desc; i++)
    if (GTK_IS_WIDGET(gtk_3d_desc[i].button))
      gtk_widget_destroy(gtk_3d_desc[i].button);
  g_free(gtk_3d_desc);

  n_3D_desc = m->nsphere + m->nvector + m->ntriangle + m->nsurf + m->ncells + m->ntextboxes;
  gtk_3d_desc = g_malloc(n_3D_desc * sizeof(GTK_3D_DESC));

  /*sphere*/
  j = 0;
  for(i = 0; i < m->nsphere; i++, j++)
  {
    g_snprintf(buff1, 20, "sphere #%d", i+1);
    gtk_3d_desc[j].i3dtype=SPHERE;
    gtk_3d_desc[j].inum = i;
    gtk_3d_desc[j].button = gtk_button_new_with_label(buff1);
    gtk_box_pack_start(GTK_BOX(vbox_3Dobject), gtk_3d_desc[j].button, FALSE, FALSE, 0);
    color.red = (guint16) (m->sphere_color[i][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->sphere_color[i][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->sphere_color[i][2] * (gfloat) G_MAXUINT16);
    gtk_widget_modify_bg(GTK_WIDGET(gtk_3d_desc[j].button), GTK_STATE_NORMAL, &color);
    g_signal_connect(G_OBJECT(gtk_3d_desc[j].button), "clicked", G_CALLBACK(luscus_gtk_modify_3Dobject), (gpointer) (&gtk_3d_desc[j]));

    gtk_widget_show(gtk_3d_desc[j].button);
  }
  /*vector*/
  for(i = 0; i < m->nvector; i++, j++)
  {
    g_snprintf(buff1, 20, "vector #%d", i+1);
    gtk_3d_desc[j].i3dtype=VECTOR;
    gtk_3d_desc[j].inum = i;
    gtk_3d_desc[j].button = gtk_button_new_with_label(buff1);
    gtk_box_pack_start(GTK_BOX(vbox_3Dobject), gtk_3d_desc[j].button, FALSE, FALSE, 0);
    color.red = (guint16) (m->vector_color[i][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->vector_color[i][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->vector_color[i][2] * (gfloat) G_MAXUINT16);
    gtk_widget_modify_bg(GTK_WIDGET(gtk_3d_desc[j].button), GTK_STATE_NORMAL, &color);
    g_signal_connect(G_OBJECT(gtk_3d_desc[j].button), "clicked", G_CALLBACK(luscus_gtk_modify_3Dobject), (gpointer) (&gtk_3d_desc[j]));
    gtk_widget_show(gtk_3d_desc[j].button);
  }
  /*triangle*/
  for(i = 0; i < m->ntriangle; i++, j++)
  {
    g_snprintf(buff1, 20, "triangle #%d", i+1);
    gtk_3d_desc[j].i3dtype=TRIANGLE;
    gtk_3d_desc[j].inum = i;
    gtk_3d_desc[j].button = gtk_button_new_with_label(buff1);
    gtk_box_pack_start(GTK_BOX(vbox_3Dobject), gtk_3d_desc[j].button, FALSE, FALSE, 0);
    color.red = (guint16) (m->triangle_color[i][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->triangle_color[i][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->triangle_color[i][2] * (gfloat) G_MAXUINT16);
    gtk_widget_modify_bg(GTK_WIDGET(gtk_3d_desc[j].button), GTK_STATE_NORMAL, &color);
    g_signal_connect(G_OBJECT(gtk_3d_desc[j].button), "clicked", G_CALLBACK(luscus_gtk_modify_3Dobject), (gpointer) (&gtk_3d_desc[j]));
    gtk_widget_show(gtk_3d_desc[j].button);
  }
  /*surface*/
  for(i = 0; i < m->nsurf; i++, j++)
  {
    g_snprintf(buff1, 20, "surface #%d", i+1);
    gtk_3d_desc[j].i3dtype=SURFACE;
    gtk_3d_desc[j].inum = i;
    gtk_3d_desc[j].button = gtk_button_new_with_label(buff1);
    gtk_box_pack_start(GTK_BOX(vbox_3Dobject), gtk_3d_desc[j].button, FALSE, FALSE, 0);
    color.red = (guint16) (m->surf_color[i][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->surf_color[i][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->surf_color[i][2] * (gfloat) G_MAXUINT16);
    gtk_widget_modify_bg(GTK_WIDGET(gtk_3d_desc[j].button), GTK_STATE_NORMAL, &color);
    g_signal_connect(G_OBJECT(gtk_3d_desc[j].button), "clicked", G_CALLBACK(luscus_gtk_modify_3Dobject), (gpointer) (&gtk_3d_desc[j]));
    gtk_widget_show(gtk_3d_desc[j].button);
  }
  /*cell*/
  for(i = 0; i < m->ncells; i++, j++)
  {
    printf("PUTTING CELL BUTTON #%d\n", i); fflush(stdout);
    gtk_3d_desc[j].i3dtype=CELL;
    gtk_3d_desc[j].inum = i;
    g_snprintf(buff1, 20, "cell #%d", i+1);
    gtk_3d_desc[j].button = gtk_button_new_with_label(buff1);
    gtk_box_pack_start(GTK_BOX(vbox_3Dobject), gtk_3d_desc[j].button, FALSE, FALSE, 0);
    color.red = (guint16) (m->cell_color[i][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->cell_color[i][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->cell_color[i][2] * (gfloat) G_MAXUINT16);
    gtk_widget_modify_bg(GTK_WIDGET(gtk_3d_desc[j].button), GTK_STATE_NORMAL, &color);
    g_signal_connect(G_OBJECT(gtk_3d_desc[j].button), "clicked", G_CALLBACK(luscus_gtk_modify_3Dobject), (gpointer) (&gtk_3d_desc[j]));
    gtk_widget_show(gtk_3d_desc[j].button);
  }

#ifdef EBUG
  printf("checking textboxes\n");
#endif
  /*textboxes*/
  for(i = 0; i < m->ntextboxes; i++, j++)
  {
#ifdef EBUG
    printf("textbox #%d\n", i);
#endif
    gtk_3d_desc[j].i3dtype=TEXTBOX;
    gtk_3d_desc[j].inum = i;
    g_snprintf(buff1, 20, "textbox #%d", i+1);
    gtk_3d_desc[j].button = gtk_button_new_with_label(buff1);
    gtk_box_pack_start(GTK_BOX(vbox_3Dobject), gtk_3d_desc[j].button, FALSE, FALSE, 0);

    color.red = (guint16) (m->textboxes[i].color[0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->textboxes[i].color[1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->textboxes[i].color[2] * (gfloat) G_MAXUINT16);

    gtk_widget_modify_text(GTK_WIDGET(gtk_3d_desc[j].button), GTK_STATE_NORMAL, &color);
    g_signal_connect(G_OBJECT(gtk_3d_desc[j].button), "clicked", G_CALLBACK(luscus_gtk_modify_3Dobject), (gpointer) (&gtk_3d_desc[j]));
    gtk_widget_show(gtk_3d_desc[j].button);
    if (!m->textboxes[i].pixtext.pixels)
      draw_pixdata_textbox(i);
  }

/*int n_3D_desc = 0;
GTK_3D_DESC *gtk_3d_desc;*/
}

void geo_play_button_play(int mode)
{
  GtkWidget *image;
  if (mode == 0) image = gtk_image_new_from_stock(GTK_STOCK_MEDIA_PAUSE, GTK_ICON_SIZE_BUTTON);
  else image = gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY, GTK_ICON_SIZE_BUTTON);
  gtk_button_set_image(GTK_BUTTON(geo_play_button), image);
}

/*typedef struct gtk_3d_description
{
  GtkWidget *button;
  char i3dtype;
  int inum;
} GTK_3D_DESC;*/

