/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdlib.h>
#include<string.h>
#include<gtk/gtk.h>
#ifdef GTK_GLEXT
#include<gtk/gtkgl.h>
#endif
//#include<gdk-pixbuf/gdk-pixbuf.h>
#include<gdk/gdkkeysyms.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_gtk.h"

/*int showing_symbols = 0;
int showing_indices = 0;
int showing_numeration = 0;
int showing_names = 0;
int showing_mulliken = 0;
int showing_loprop = 0;*/

int rot_axis_order = 2;
int n_input_types;
int elem_selected = 1;
int control_pressed;
INPUT_FORMAT *input_filetypes;
/*INPUT_DATA Input_Data;*/
GtkWidget *element_dialog;
GtkWidget *atom_prop_label;
GtkWidget *atom_prop_color_button;
GtkWidget *atom_prop_atom_size_spin;
GtkWidget *atom_prop_bond_size_spin;
GtkWidget *atom_prop_atom_valency_spin;

ELEM_DATA *e;
ELEM_DATA *tmp_elem;

gint calc_timer;
char calc_alive = 0;

double angle_values[] = 
{
  60., 90., 109.45, 120., 180.
};

void luscus_check_writable_pic_files(GdkPixbufFormat *data, GSList **list)
{
  if (gdk_pixbuf_format_is_writable(data))
    *list = g_slist_prepend (*list, data);
}


void luscus_toggle_movement_mode(GtkWidget *combo, gpointer data)
{
  luscus_gtk_pop_message_from_statusbar2();
  if (gtk_combo_box_get_active(GTK_COMBO_BOX(combo)) == 1)
  {
    move_camera=1;
    luscus_gtk_push_message_to_statusbar2("moving camera");
  }
  else
  {
    move_camera=0;
    luscus_gtk_push_message_to_statusbar2("moving molecule");
  }
}

void luscus_gtk_move_up(GtkWidget *button, gpointer data)
{
  Do_key_ud(1, 0);
  redraw(); 
}

void luscus_gtk_move_down(GtkWidget *button, gpointer data)
{
  Do_key_ud(0, 0);
  redraw(); 
}

void luscus_gtk_move_forw(GtkWidget *button, gpointer data)
{
  Do_key_lr(1, 0);
  redraw();
}

void luscus_gtk_move_back(GtkWidget *button, gpointer data)
{
  Do_key_lr(0, 0);
  redraw();
}

void luscus_gtk_center(GtkWidget *button, gpointer data)
{
  Do_center();
/*  set_origin_molecule();*/
  redraw();
}

void luscus_gtk_zoom(GtkWidget *button, gint in_out)
{
  change_scale(in_out);
  redraw();
}

void callback_geo_forw(GtkWidget *button, gpointer data)
{
  if (!m) return;
  if (igeo == n_geometries-1) return;
  igeo++;
  if (igeo >= n_geometries) igeo = n_geometries - 1;
  get_new_section();
  luscus_gtk_update_geo_info();
  set_current_graph_data();
  set_scale();
  rerender_3d();
  redraw();
}

void callback_geo_back(GtkWidget *button, gpointer data)
{
  if (!m) return;
  if (igeo == 0) return;
  igeo--;
  if (igeo <= 0) igeo = 0;
  get_new_section();
  luscus_gtk_update_geo_info();
  set_current_graph_data();
  set_scale();
  rerender_3d();
  redraw();
}

void callback_geo_first(GtkWidget *button, gpointer data)
{
  if (!m) return;
  if (igeo == 0) return;
  igeo = 0;
  get_new_section();
  luscus_gtk_update_geo_info();
  set_current_graph_data();
  rerender_3d();
  redraw();
}

void callback_geo_last(GtkWidget *button, gpointer data)
{
  if (!m) return;
  if (igeo == n_geometries) return;
  igeo = n_geometries - 1;
  get_new_section();
  luscus_gtk_update_geo_info();
  set_current_graph_data();
  rerender_3d();
  redraw();
}

#ifdef GTK3
void callback_adj_geometry(GtkAdjustment *adj, gpointer data)
#endif
#ifdef GTK2
void callback_adj_geometry(GtkObject *adj, gpointer data)
#endif
{
  if (!m) return;
  igeo = (int) gtk_adjustment_get_value(GTK_ADJUSTMENT(adj));
  if (igeo >= n_geometries-1) igeo = n_geometries-1;
  get_new_section();
  luscus_gtk_update_geo_info();
  set_current_graph_data();
  rerender_3d();
  redraw();
}

void callback_geo_play(GtkWidget *button, gpointer data)
{
  if (Input_Data.animate)
  {
    geo_play_button_play(1);
    Input_Data.animate = FALSE;
  }
  else
  {
    geo_play_button_play(0);
    Input_Data.animate = TRUE;
    luscus_gtk_start_animation(0);
  }
  redraw();
}

/*grid*/

void luscus_gtk_change_transparency_level(GtkWidget* adjustment, gpointer data)
{
  Input_Data.electrostatic_poten_color[0][3] =
  Input_Data.electrostatic_poten_color[1][3] =
  Input_Data.neg_pos_color[0][3] = Input_Data.neg_pos_color[1][3] =
  gtk_adjustment_get_value(GTK_ADJUSTMENT(adjustment));
  rerender_3d();
  redraw();
}

void luscus_gtk_change_isosurface_level(GtkWidget* adjustment, gpointer data)
{
  Input_Data.lev = gtk_adjustment_get_value(GTK_ADJUSTMENT(adjustment));
  make_surfaces(); 
  rerender_3d();
  redraw();
}

void luscus_gtk_show_all_orbitals(GtkWidget* toggle, gpointer data)
{
  luscus_gtk_show_multiorb_window();
}

void luscus_gtk_show_electrostatic_potential(GtkWidget* button, gpointer data)
{
  if (Input_Data.show_epot) Input_Data.show_epot = 0;
  else Input_Data.show_epot = 1;
  rerender_3d();
  redraw();
}

/*vibrations*/

void luscus_gtk_change_vibration_speed(GtkWidget *adjustment, gpointer data)
{
  Input_Data.frequency_speed = gtk_adjustment_get_value(GTK_ADJUSTMENT(adjustment));
  redraw();
}

void luscus_gtk_change_vibration_amplitude(GtkWidget *adjustment, gpointer data)
{
  Input_Data.frequency_amplitude = gtk_adjustment_get_value(GTK_ADJUSTMENT(adjustment));
  redraw();
}

/*fragments*/

void luscus_gtk_add_fragment(GtkWidget *button, gpointer data)
{
  if (m->n_selected == 2) change_bond_type(SINGLE_BOND);
  else do_key_insert();
  printf("INSERTING FRAGMENT\n");
  rerender_3d();
  redraw();
}

void luscus_gtk_remove_fragment(GtkWidget* widget, gpointer data)
{
  if (m->n_selected == 2) change_bond_type(NO_BOND);
  else delete_coord(0);
  rerender_3d();
  redraw();
}

void luscus_gtk_change_bond_callback(GtkWidget *combo_bond, gpointer data)
{
  change_bond_type(gtk_combo_box_get_active(GTK_COMBO_BOX(combo_bond)));
  rerender_3d();
  redraw();
}

void luscus_gtk_change_angle_callback(GtkWidget *combo_angle, gpointer data)
{
  int iang;
  iang = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_angle));
  luscus_gtk_move_coord(angle_values[iang]);
  rerender_3d();
  redraw();
}

void callback_change_atom(GtkWidget* widget, gpointer data)
{
  /*PERIODIC SYSTEM*/
  GtkWidget *vbox, *hbox;
  GtkWidget *persys; /*GtkTable with buttons arreanged in form of periodic system*/
  GtkWidget *expander;

  GtkWidget *label;
#ifdef GTK2
  GtkObject *adj_atom_size, *adj_bond_size, *adj_valency;
#else
  GtkAdjustment *adj_atom_size, *adj_bond_size, *adj_valency;
#endif
  GdkColor color;
  gint response;

  GtkWidget *custom_atom_color_button;
  GtkWidget *custom_atom_valency_spin;
  GtkWidget *custom_atom_size_spin;
  GtkWidget *custom_bond_size_spin;

  printf("CHANGE ATOM\n");

  element_dialog = gtk_dialog_new_with_buttons("Choose element", GTK_WINDOW(gtk_widget_get_toplevel(widget)),
                                               GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                               GTK_STOCK_CLOSE, GTK_RESPONSE_REJECT,
                                               NULL);

    /*periodic system*/
  persys = luscus_gtk_make_periodic_system();
  vbox = gtk_dialog_get_content_area(GTK_DIALOG(element_dialog));
  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(element_dialog))), persys, FALSE, FALSE, 0);
   /*custom data change*/
  expander = gtk_expander_new("custom data");

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_container_add(GTK_CONTAINER(expander), hbox);

  color.red = (guint16) (m->elem[m->selected[0]].color[0] * (gfloat) G_MAXUINT16);
  color.green = (guint16) (m->elem[m->selected[0]].color[1] * (gfloat) G_MAXUINT16);
  color.blue = (guint16) (m->elem[m->selected[0]].color[2] * (gfloat) G_MAXUINT16);

  custom_atom_color_button = gtk_color_button_new_with_color(&color);
  gtk_box_pack_start(GTK_BOX(hbox), custom_atom_color_button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(custom_atom_color_button), "color-set", G_CALLBACK(custom_atom_color_changed_callback), NULL);
  gtk_widget_show(custom_atom_color_button);
  /* -atom size- */
  label = gtk_label_new("Van der Waals radius:");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  adj_atom_size = gtk_adjustment_new(0.0, 0.F, 3.F, 0.01F, 0.1F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_atom_size), m->elem[m->selected[0]].vdw_rad);

#ifdef GTK2
  custom_atom_size_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj_atom_size), 0.01, 2);
#else
  custom_atom_size_spin = gtk_spin_button_new(adj_atom_size, 0.01, 2);
#endif
  g_signal_connect(G_OBJECT(custom_atom_size_spin), "value-changed", G_CALLBACK(custom_atom_size_changed_callback), NULL);
  gtk_box_pack_start(GTK_BOX(hbox), custom_atom_size_spin, FALSE, FALSE, 3);
  gtk_widget_show(custom_atom_size_spin);

   /* -bond size- */
  label = gtk_label_new("bond size:");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  adj_bond_size = gtk_adjustment_new(0.0, 0.F, 3.F, 0.01F, 0.1F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_bond_size), m->elem[m->selected[0]].bond_rad);

#ifdef GTK2
  custom_bond_size_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj_bond_size), 0.01, 2);
#else
  custom_bond_size_spin = gtk_spin_button_new(adj_bond_size, 0.01, 2);
#endif
  g_signal_connect(G_OBJECT(custom_bond_size_spin), "value-changed", G_CALLBACK(custom_bond_size_changed_callback), NULL);

  gtk_box_pack_start(GTK_BOX(hbox), custom_bond_size_spin, FALSE, FALSE, 0);
  gtk_widget_show(custom_bond_size_spin);

  label = gtk_label_new("valency:");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

    /* -valency- */

  adj_valency = gtk_adjustment_new(0.0, 0.F, 30.F, 1.0F, 1.0F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_valency), m->elem[m->selected[0]].valency);

#ifdef GTK2
  custom_atom_valency_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj_valency), 0.01, 2);
#else
  custom_atom_valency_spin = gtk_spin_button_new(adj_valency, 0.01, 2);
#endif
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(custom_atom_valency_spin), 0);
  g_signal_connect(G_OBJECT(custom_atom_valency_spin), "value-changed", G_CALLBACK(custom_valency_changed_callback), NULL);

  gtk_box_pack_start(GTK_BOX(hbox), custom_atom_valency_spin, FALSE, FALSE, 0);
  gtk_widget_show(custom_atom_valency_spin);

  gtk_widget_show(hbox);
  gtk_widget_show(expander);
  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(element_dialog))), expander, FALSE, FALSE, 0);

  response = gtk_dialog_run(GTK_DIALOG(element_dialog));
  if (GTK_IS_WIDGET(element_dialog)) gtk_widget_destroy(element_dialog);
  rerender_3d();
  redraw();
}

GtkWidget *luscus_gtk_make_periodic_system(void)
{
  GtkWidget *table;
  GtkWidget *button;
  gint i;

  gint maxx = 0, maxy = 0;

  for (i = 1; i < number_of_elements; i++)
  {
    if (e[i].periodic_pos_x > maxx) maxx = e[i].periodic_pos_x;
    if (e[i].periodic_pos_y > maxy) maxy = e[i].periodic_pos_y;
  }

  table = gtk_table_new(maxy+1, maxx+1, TRUE);

  for(i = 1; i < number_of_elements; i++)
  {
    button = gtk_button_new_with_label(e[i].name);
    g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_select_element), GINT_TO_POINTER(i));
    gtk_table_attach(GTK_TABLE(table), button, e[i].periodic_pos_x, e[i].periodic_pos_x+1, e[i].periodic_pos_y, e[i].periodic_pos_y+1, GTK_FILL, GTK_FILL, 0, 0);
    gtk_widget_show(button);
  }
  gtk_widget_show(table);

  return table;
}

gboolean luscus_gtk_select_element(GtkWidget *button, gint iatom)
{
  int i;
  for(i = 0; i < 3; i++)
    m->elem[m->selected[0]].color[i] = e[iatom].color[i];
  m->elem[m->selected[0]].vdw_rad = e[iatom].vdw_rad;
  m->elem[m->selected[0]].bond_rad = e[iatom].bond_rad;
  m->elem[m->selected[0]].valency = e[iatom].valency;
  free(m->elem[m->selected[0]].name);
  m->elem[m->selected[0]].name = strdup(e[iatom].name);
  change_atom_parameters_in_list(m->selected[0]);
  append_backup();
  gtk_widget_destroy(element_dialog);
  return FALSE;
}

void custom_atom_color_changed_callback(GtkWidget *col_button, gpointer data)
{
  GdkColor color;

  gtk_color_button_get_color(GTK_COLOR_BUTTON(col_button), &color);

  if (m->selected[0] < 0) return;

  m->elem[m->selected[0]].color[0] = (double) color.red / (double) G_MAXUINT16;
  m->elem[m->selected[0]].color[1] = (double) color.green / (double) G_MAXUINT16;
  m->elem[m->selected[0]].color[2] = (double) color.blue / (double) G_MAXUINT16;
  m->ishow ^= HAS_COLOUR;
}

void custom_atom_size_changed_callback(GtkWidget *spin, gpointer data)
{
  if (m->selected[0] < 0) return;
  m->elem[m->selected[0]].vdw_rad = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));
}

void custom_bond_size_changed_callback(GtkWidget *spin, gpointer data)
{
  if (m->selected[0] < 0) return;
  m->elem[m->selected[0]].bond_rad = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));
}

void custom_valency_changed_callback(GtkWidget *spin, gpointer data)
{
  if (m->selected[0] < 0) return;
  m->elem[m->selected[0]].valency = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));
}

gboolean luscus_gtk_disable_keys(GtkWidget *widget, GdkEvent *event, gpointer data)
{
  accept_keys = 0;
  return FALSE;
}

gboolean luscus_gtk_enable_keys(GtkWidget *widget, GdkEvent *event, gpointer user_data)
{
  accept_keys = 1;
  return FALSE;
}

/*this function is replaced with three functions (below) ->0.3.6*/
/*void luscus_gtk_change_bond_angle_tors_value(GtkWidget* adj, gpointer data)
{
  gdouble value = gtk_spin_button_get_value(GTK_SPIN_BUTTON(adj));
  luscus_gtk_move_coord(value);
  redraw();
}*/

void luscus_gtk_change_bond_value(GtkWidget* adj, gpointer data)
{
  gdouble value = gtk_spin_button_get_value(GTK_SPIN_BUTTON(adj));
  luscus_gtk_move_bond(value);
  rerender_3d();
  redraw();
}

void luscus_gtk_change_angle_value(GtkWidget* adj, gpointer data)
{
  gdouble value = gtk_spin_button_get_value(GTK_SPIN_BUTTON(adj));
  luscus_gtk_move_angle(value);
  rerender_3d();
  redraw();
}

void luscus_gtk_change_tors_value(GtkWidget* adj, gpointer data)
{
  gdouble value = gtk_spin_button_get_value(GTK_SPIN_BUTTON(adj));
  luscus_gtk_move_torsion(value);
  rerender_3d();
  redraw();
}

void luscus_gtk_change_xyz_coordinate(GtkWidget* spin, gpointer data)
{
  int icoord=GPOINTER_TO_INT(data);
  double value;

  if (m->n_selected < 1 || m->n_selected > 4) return;
  value = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));
  accumulate_motion(m->xyz[m->selected[0]][icoord] - value);
  m->xyz[m->selected[0]][icoord] = value;
  rerender_3d();

  redraw();
  change_watched_data();
  if (Input_Data.automatic_rebonding) rebond();
  change_atom_parameters_in_list(m->selected[0]);
}

void luscus_gtk_add_dummy_atoms(GtkWidget* widget, gpointer data)
{
  add_atoms(1);
  rerender_3d();
  redraw();
}

void luscus_gtk_remove_dummy_atoms(GtkWidget* widget, gpointer data)
{
  delete_dummy_atoms();
  rerender_3d();
  redraw();
}

void luscus_gtk_remove_selection(GtkWidget* widget, gpointer data)
{
  unselect_all();
  redraw();
}

void luscus_gtk_unmark(GtkWidget* widget, gpointer data)
{
  unmark_all();
  rerender_3d();
  redraw();
}

void luscus_gtk_sort_mark(GtkWidget* widget, gpointer data)
{
  renumber_marked();
  rerender_3d();
  redraw();
}

void luscus_gtk_mark_H_atoms(GtkWidget* widget, gpointer data)
{
  mark_H();
  rerender_3d();
  redraw();
}

void luscus_gtk_mark_reverse(GtkWidget* widget, gpointer data)
{
  reverse_marked();
  rerender_3d();
  redraw();
}

void luscus_gtk_mark_element(GtkWidget* widget, gpointer data)
{
  mark_element_as_selected();
  rerender_3d();
  redraw();
}

void luscus_gtk_mark_neighbor(GtkWidget* widget, gpointer data)
{
  mark_neighbor();
  rerender_3d();
  redraw();
}

void luscus_gtk_watch_value(GtkWidget* widget, gpointer data)
{
  add_watched_coord();
  luscus_gtk_update_upon_select_or_mark();
}

void luscus_gtk_unwatch_values(GtkWidget* widget, gpointer data)
{
  remove_watched_coords();
  luscus_gtk_update_upon_select_or_mark(); 
}

void luscus_gtk_select_fragment(GtkWidget *adjustment, gpointer data)
{
  add_fragment(GPOINTER_TO_INT(data));
  rerender_3d();
  redraw();
}

#ifndef HAS_MSYM
void luscus_gtk_alpply_symmetry(GtkWidget* button, gpointer data)
{
  luscus_gtk_set_vec_rot(0);
  do_symmetry();
  rerender_3d();
  redraw();
}
#endif

void luscus_gtk_apply_translation(GtkWidget* button, gpointer data)
{
  luscus_gtk_set_vec_rot(0);
  do_symmetry();
  rerender_3d();
  redraw();
}

void change_translation_magnitude(GtkWidget *sb, gpointer data)
{
  Input_Data.symmetry_translation = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sb));
}

void luscus_gtk_alpply_rot_symmetry(GtkWidget* button, gpointer data)
{
  luscus_gtk_set_vec_rot(rot_axis_order);
  do_symmetry();
  rerender_3d();
  redraw();
}

void change_rotation_axis_order(GtkWidget *sb, gpointer data)
{
  rot_axis_order = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(sb));
}

#ifndef HAS_MSYM

void luscus_gtk_force_symmetry(GtkWidget* button, gpointer data)
{
  Input_Data.force_symmetry = 1;
  do_symmetry();
  rerender_3d();
  redraw();
}
#endif

void luscus_add_arrow(GtkWidget* button, gpointer data)
{
  if (m->n_selected != 2) return;
  add_vector(m, m->xyz[m->selected[0]], m->xyz[m->selected[1]]);
  rerender_3d();
  redraw();
  append_backup();
}

void luscus_add_sphere(GtkWidget* button, gpointer data)
{
  if (m->n_selected == 0) return;
  else if (m->n_selected == 1) add_sphere(m, m->xyz[m->selected[0]], 2.0);
  else if (m->n_selected == 2) add_sphere(m, m->xyz[m->selected[0]], rr_dist(m->xyz[m->selected[0]], m->xyz[m->selected[1]]));
  else return;
  rerender_3d();
  redraw();
  append_backup();
}

void luscus_add_plain(GtkWidget *button, gpointer data)
{
  if (m->n_selected != 3) return;
  add_surface(m, m->xyz[m->selected[0]], m->xyz[m->selected[1]], m->xyz[m->selected[2]]);
  rerender_3d();
  redraw();
  append_backup();
}

void luscus_add_triangle(GtkWidget *button, gpointer data)
{
  if (m->n_selected != 3) return;
  add_triangle(m, m->xyz[m->selected[0]], m->xyz[m->selected[1]], m->xyz[m->selected[2]]);
  rerender_3d();
  redraw();
  append_backup();
}

void luscus_add_cell(GtkWidget *button, gpointer data)
{
  if (m->n_selected != 4) return;
  add_cell(m, m->xyz[m->selected[0]], m->xyz[m->selected[1]], m->xyz[m->selected[2]], m->xyz[m->selected[3]]);
  rerender_3d();
  redraw();
  append_backup();
}

void luscus_insert_textbox(GtkWidget *button, gpointer data)
{
  if (!textbox_state) /*this way shit won't happen if user clicks "insert text" twice*/
  {
    allocate_textbox(m);
    m->textboxes[m->ntextboxes-1].coord_x = 0;
    m->textboxes[m->ntextboxes-1].coord_y = 0;

    m->textboxes[m->ntextboxes-1].color[0] = 0.f;
    m->textboxes[m->ntextboxes-1].color[1] = 0.f;
    m->textboxes[m->ntextboxes-1].color[2] = 0.f;
    m->textboxes[m->ntextboxes-1].color[3] = 1.0;

    m->textboxes[m->ntextboxes-1].font = g_strdup(Input_Data.font);
    m->textboxes[m->ntextboxes-1].message = NULL;
    m->textboxes[m->ntextboxes-1].pixtext.width = 0;
    m->textboxes[m->ntextboxes-1].pixtext.height = 0;
    m->textboxes[m->ntextboxes-1].pixtext.pixels = NULL;

    textbox_state = 1;
    luscus_set_cursor_text();
  }
}

void luscus_clear_drawings(GtkWidget *button, gpointer data)
{
  deallocate_vectors(m);
  deallocate_triangles(m);
  deallocate_spheres(m);
  deallocate_surfaces(m);
  deallocate_cells(m);
  deallocate_textboxes(m);
  luscus_gtk_update_upon_select_or_mark();
  luscus_gtk_update_3Dobject_info();
  rerender_3d();
  redraw();
  append_backup();
}

void luscus_gtk_modify_3Dobject(GtkWidget *button, GTK_3D_DESC* desc)
{
  int i;
  GtkWidget *dialog;
  GtkWidget *vbox, *hbox;
  GtkWidget *label;
  GtkWidget *color_button;
/*  GtkWidget *hscale;*/
  GtkWidget *spin_button;
  GtkWidget *separator;
  GtkWidget *entry_00, *entry_01, *entry_02;
  GtkWidget *entry_10, *entry_11, *entry_12;
  GtkWidget *entry_20, *entry_21, *entry_22;
  GtkWidget *entry_30, *entry_31, *entry_32;
  gint response;
  gchar title[20];
  gchar buffer[30];
  GdkColor color;
#ifdef GTK2
  GtkObject *adj_transparency;
#endif
#ifdef GTK3
  GtkAdjustment *adj_transparency;
#endif
  double new_value;
  GtkWidget *font_button;
  gint width, height;
  GtkWidget *button_1;

  set_go_selected(desc->i3dtype, desc->inum);
  rerender_3d();
  redraw();
    
  if (desc->i3dtype == SPHERE)
  {
    g_snprintf(title, 20, "sphere #%d", desc->inum+1);
    color.red = (guint16) (m->sphere_color[desc->inum][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->sphere_color[desc->inum][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->sphere_color[desc->inum][2] * (gfloat) G_MAXUINT16);    
    adj_transparency = gtk_adjustment_new(m->sphere_color[desc->inum][3], 0.0, 1.0, 0.01, 0.1, 0.0);
  }
  else if (desc->i3dtype == VECTOR)
  {
    g_snprintf(title, 20, "vector #%d", desc->inum+1);
    color.red = (guint16) (m->vector_color[desc->inum][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->vector_color[desc->inum][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->vector_color[desc->inum][2] * (gfloat) G_MAXUINT16);
    adj_transparency = gtk_adjustment_new(m->vector_color[desc->inum][3], 0.0, 1.0, 0.01, 0.1, 0.0);
  }
  else if (desc->i3dtype == TRIANGLE)
  {
    g_snprintf(title, 20, "triangle #%d", desc->inum+1);
    color.red = (guint16) (m->triangle_color[desc->inum][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->triangle_color[desc->inum][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->triangle_color[desc->inum][2] * (gfloat) G_MAXUINT16);
    adj_transparency = gtk_adjustment_new(m->triangle_color[desc->inum][3], 0.0, 1.0, 0.01, 0.1, 0.0);
  }
  else if (desc->i3dtype == SURFACE)
  {
    g_snprintf(title, 20, "surface #%d", desc->inum+1);
    color.red = (guint16) (m->surf_color[desc->inum][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->surf_color[desc->inum][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->surf_color[desc->inum][2] * (gfloat) G_MAXUINT16);
    adj_transparency = gtk_adjustment_new(m->surf_color[desc->inum][3], 0.0, 1.0, 0.01, 0.1, 0.0);
  }
  else if (desc->i3dtype == CELL)
  {
    g_snprintf(title, 20, "cell #%d", desc->inum+1);
    color.red = (guint16) (m->cell_color[desc->inum][0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->cell_color[desc->inum][1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->cell_color[desc->inum][2] * (gfloat) G_MAXUINT16);
    adj_transparency = gtk_adjustment_new(m->cell_color[desc->inum][3], 0.0, 1.0, 0.01, 0.1, 0.0);
  }
  else if (desc->i3dtype == TEXTBOX)
  {
    g_snprintf(title, 20, "textbox #%d", desc->inum+1);
    color.red = (guint16) (m->textboxes[desc->inum].color[0] * (gfloat) G_MAXUINT16);
    color.green = (guint16) (m->textboxes[desc->inum].color[1] * (gfloat) G_MAXUINT16);
    color.blue = (guint16) (m->textboxes[desc->inum].color[2] * (gfloat) G_MAXUINT16);

  }

  dialog = gtk_dialog_new_with_buttons (title, NONE, GTK_DIALOG_MODAL,
                                        GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                        GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                        GTK_STOCK_DELETE, GTK_RESPONSE_YES, NULL);

  vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

  if (desc->i3dtype != TEXTBOX)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("color");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    color_button = gtk_color_button_new_with_color(&color);
    gtk_box_pack_start(GTK_BOX(hbox), color_button, TRUE, TRUE, 0);
    /*g_signal_connect*/
    gtk_widget_show(color_button);

#ifdef GTK2
    separator = gtk_vseparator_new();
#endif
#ifdef GTK3
    separator = gtk_separator_new(GTK_ORIENTATION_VERTICAL);
#endif
    gtk_box_pack_start(GTK_BOX(hbox), separator, TRUE, TRUE, 0);
    gtk_widget_show(separator);
  
    label = gtk_label_new("transparency");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    spin_button = gtk_spin_button_new(GTK_ADJUSTMENT(adj_transparency), 0.01, 2);
    gtk_box_pack_start(GTK_BOX(hbox), spin_button, TRUE, TRUE, 0);
    gtk_widget_show(spin_button);

    gtk_widget_show(hbox);
  }

  if (desc->i3dtype == SPHERE)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("center");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->sphere_center[desc->inum][0]);
    entry_00 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_00, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_00), buffer);
    gtk_widget_show(entry_00);

    g_snprintf(buffer, 30, "%15.8f", m->sphere_center[desc->inum][1]);
    entry_01 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_01, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_01), buffer);
    gtk_widget_show(entry_01);

    g_snprintf(buffer, 30, "%15.8f", m->sphere_center[desc->inum][2]);
    entry_02 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_02, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_02), buffer);
    gtk_widget_show(entry_02);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("radius");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->sphere_radius[desc->inum]);
    entry_10 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_10, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_10), buffer);
    gtk_widget_show(entry_10);

    gtk_widget_show(hbox);
  }
  else if (desc->i3dtype == VECTOR)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("base");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->vector1[desc->inum][0]);
    entry_00 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_00, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_00), buffer);
    gtk_widget_show(entry_00);

    g_snprintf(buffer, 30, "%15.8f", m->vector1[desc->inum][1]);
    entry_01 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_01, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_01), buffer);
    gtk_widget_show(entry_01);

    g_snprintf(buffer, 30, "%15.8f", m->vector1[desc->inum][2]);
    entry_02 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_02, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_02), buffer);
    gtk_widget_show(entry_02);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("tip");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->vector2[desc->inum][0]);
    entry_10 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_10, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_10), buffer);
    gtk_widget_show(entry_10);

    g_snprintf(buffer, 30, "%15.8f", m->vector2[desc->inum][1]);
    entry_11 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_11, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_11), buffer);
    gtk_widget_show(entry_11);

    g_snprintf(buffer, 30, "%15.8f", m->vector2[desc->inum][2]);
    entry_12 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_12, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_12), buffer);
    gtk_widget_show(entry_12);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("radius");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->radius[desc->inum]);
    entry_20 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_20, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_20), buffer);
    gtk_widget_show(entry_20);

    label = gtk_label_new("sharpness");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->sharpness[desc->inum]);
    entry_21 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_21, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_21), buffer);
    gtk_widget_show(entry_21);

    gtk_widget_show(hbox);
  }
  else if (desc->i3dtype == TRIANGLE)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 1");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->triangle1[desc->inum][0]);
    entry_00 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_00, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_00), buffer);
    gtk_widget_show(entry_00);

    g_snprintf(buffer, 30, "%15.8f", m->triangle1[desc->inum][1]);
    entry_01 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_01, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_01), buffer);
    gtk_widget_show(entry_01);

    g_snprintf(buffer, 30, "%15.8f", m->triangle1[desc->inum][2]);
    entry_02 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_02, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_02), buffer);
    gtk_widget_show(entry_02);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 2");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->triangle2[desc->inum][0]);
    entry_10 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_10, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_10), buffer);
    gtk_widget_show(entry_10);

    g_snprintf(buffer, 30, "%15.8f", m->triangle2[desc->inum][1]);
    entry_11 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_11, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_11), buffer);
    gtk_widget_show(entry_11);

    g_snprintf(buffer, 30, "%15.8f", m->triangle2[desc->inum][2]);
    entry_12 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_12, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_12), buffer);
    gtk_widget_show(entry_12);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 3");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->triangle3[desc->inum][0]);
    entry_20 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_20, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_20), buffer);
    gtk_widget_show(entry_20);

    g_snprintf(buffer, 30, "%15.8f", m->triangle3[desc->inum][1]);
    entry_21 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_21, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_21), buffer);
    gtk_widget_show(entry_21);

    g_snprintf(buffer, 30, "%15.8f", m->triangle3[desc->inum][2]);
    entry_22 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_22, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_22), buffer);
    gtk_widget_show(entry_22);

    gtk_widget_show(hbox);
  }
  else if (desc->i3dtype == SURFACE)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("point 1");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->surf1[desc->inum][0]);
    entry_00 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_00, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_00), buffer);
    gtk_widget_show(entry_00);

    g_snprintf(buffer, 30, "%15.8f", m->surf1[desc->inum][1]);
    entry_01 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_01, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_01), buffer);
    gtk_widget_show(entry_01);

    g_snprintf(buffer, 30, "%15.8f", m->surf1[desc->inum][2]);
    entry_02 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_02, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_02), buffer);
    gtk_widget_show(entry_02);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("point 2");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->surf2[desc->inum][0]);
    entry_10 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_10, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_10), buffer);
    gtk_widget_show(entry_10);

    g_snprintf(buffer, 30, "%15.8f", m->surf2[desc->inum][1]);
    entry_11 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_11, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_11), buffer);
    gtk_widget_show(entry_11);

    g_snprintf(buffer, 30, "%15.8f", m->surf2[desc->inum][2]);
    entry_12 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_12, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_12), buffer);
    gtk_widget_show(entry_12);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("point 3");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->surf3[desc->inum][0]);
    entry_20 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_20, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_20), buffer);
    gtk_widget_show(entry_20);

    g_snprintf(buffer, 30, "%15.8f", m->surf3[desc->inum][1]);
    entry_21 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_21, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_21), buffer);
    gtk_widget_show(entry_21);

    g_snprintf(buffer, 30, "%15.8f", m->surf3[desc->inum][2]);
    entry_22 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_22, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_22), buffer);
    gtk_widget_show(entry_22);

    gtk_widget_show(hbox);
  }
  else if (desc->i3dtype == CELL)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 1");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->cell1[desc->inum][0]);
    entry_00 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_00, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_00), buffer);
    gtk_widget_show(entry_00);

    g_snprintf(buffer, 30, "%15.8f", m->cell1[desc->inum][1]);
    entry_01 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_01, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_01), buffer);
    gtk_widget_show(entry_01);

    g_snprintf(buffer, 30, "%15.8f", m->cell1[desc->inum][2]);
    entry_02 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_02, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_02), buffer);
    gtk_widget_show(entry_02);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 2");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->cell2[desc->inum][0]);
    entry_10 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_10, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_10), buffer);
    gtk_widget_show(entry_10);

    g_snprintf(buffer, 30, "%15.8f", m->cell2[desc->inum][1]);
    entry_11 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_11, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_11), buffer);
    gtk_widget_show(entry_11);

    g_snprintf(buffer, 30, "%15.8f", m->cell2[desc->inum][2]);
    entry_12 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_12, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_12), buffer);
    gtk_widget_show(entry_12);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 3");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->cell3[desc->inum][0]);
    entry_20 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_20, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_20), buffer);
    gtk_widget_show(entry_20);

    g_snprintf(buffer, 30, "%15.8f", m->cell3[desc->inum][1]);
    entry_21 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_21, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_21), buffer);
    gtk_widget_show(entry_21);

    g_snprintf(buffer, 30, "%15.8f", m->cell3[desc->inum][2]);
    entry_22 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_22, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_22), buffer);
    gtk_widget_show(entry_22);

    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("vertex 4");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);

    g_snprintf(buffer, 30, "%15.8f", m->cell4[desc->inum][0]);
    entry_30 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_30, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_30), buffer);
    gtk_widget_show(entry_30);

    g_snprintf(buffer, 30, "%15.8f", m->cell4[desc->inum][1]);
    entry_31 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_31, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_31), buffer);
    gtk_widget_show(entry_31);

    g_snprintf(buffer, 30, "%15.8f", m->cell4[desc->inum][2]);
    entry_32 = gtk_entry_new();
    gtk_box_pack_start(GTK_BOX(hbox), entry_32, TRUE, TRUE, 0);
    gtk_entry_set_text(GTK_ENTRY(entry_32), buffer);
    gtk_widget_show(entry_32);

    gtk_widget_show(hbox);
  }
  else if (desc->i3dtype == TEXTBOX)
  {
#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("text:");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    entry_00 = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(entry_00), m->textboxes[desc->inum].message);
    gtk_box_pack_start(GTK_BOX(hbox), entry_00, FALSE, FALSE, 0);
    gtk_widget_show(entry_00);
    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("font:");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    font_button = gtk_font_button_new_with_font(m->textboxes[desc->inum].font);
    gtk_box_pack_start(GTK_BOX(hbox), font_button, FALSE, FALSE, 0);
    gtk_widget_show(font_button);
    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("color:");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    color_button = gtk_color_button_new();
    gtk_color_button_set_color(GTK_COLOR_BUTTON(color_button), &color);
    gtk_box_pack_start(GTK_BOX(hbox), color_button, FALSE, FALSE, 0);
    gtk_widget_show(color_button);
    gtk_widget_show(hbox);

#ifdef GTK2
    hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
    hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);    
#endif
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

    label = gtk_label_new("coordinates:");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    gtk_widget_show(label);

    entry_10 = gtk_entry_new();
    g_snprintf(buffer, 20, "%d", m->textboxes[desc->inum].coord_x);
    gtk_entry_set_text(GTK_ENTRY(entry_10), buffer);
    gtk_box_pack_start(GTK_BOX(hbox), entry_10, FALSE, FALSE, 0);
    gtk_widget_show(entry_10);

    entry_11 = gtk_entry_new();
    g_snprintf(buffer, 20, "%d", m->textboxes[desc->inum].coord_y);
    gtk_entry_set_text(GTK_ENTRY(entry_11), buffer);
    gtk_box_pack_start(GTK_BOX(hbox), entry_11, FALSE, FALSE, 0);
    gtk_widget_show(entry_11);
    gtk_widget_show(hbox);

    button_1 = gtk_button_new_with_label("click on the display for new position");
    gtk_box_pack_start(GTK_BOX(vbox), button_1, FALSE, FALSE, 0);
    g_signal_connect(G_OBJECT(button_1), "clicked", G_CALLBACK(luscus_move_textbox_by_click), (gpointer) desc);
    g_signal_connect_swapped(G_OBJECT(button_1), "clicked", G_CALLBACK(gtk_widget_destroy), (gpointer) dialog);

    gtk_widget_show(button_1);
  }

  response = gtk_dialog_run(GTK_DIALOG(dialog));
  if (response == GTK_RESPONSE_ACCEPT)
  {
    gtk_color_button_get_color(GTK_COLOR_BUTTON(color_button), &color);
    if (desc->i3dtype != TEXTBOX)
      new_value = gtk_adjustment_get_value(GTK_ADJUSTMENT(adj_transparency));
    if (desc->i3dtype == SPHERE)
    {
      m->sphere_color[desc->inum][0] = (float) color.red / (float) G_MAXUINT16;
      m->sphere_color[desc->inum][1] = (float) color.green / (float) G_MAXUINT16;
      m->sphere_color[desc->inum][2] = (float) color.blue / (float) G_MAXUINT16;
      m->sphere_color[desc->inum][3] = new_value;
      m->sphere_center[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_00)));
      m->sphere_center[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_01)));
      m->sphere_center[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_02)));
      m->sphere_radius[desc->inum] = atof(gtk_entry_get_text(GTK_ENTRY(entry_10)));
    }
    else if (desc->i3dtype == VECTOR)
    {
      m->vector_color[desc->inum][0] = (float) color.red / (float) G_MAXUINT16;
      m->vector_color[desc->inum][1] = (float) color.green / (float) G_MAXUINT16;
      m->vector_color[desc->inum][2] = (float) color.blue / (float) G_MAXUINT16;
      m->vector_color[desc->inum][3] = new_value;

      m->vector1[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_00)));
      m->vector1[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_01)));
      m->vector1[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_02)));
      m->vector2[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_10)));
      m->vector2[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_11)));
      m->vector2[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_12)));
      m->radius[desc->inum] = atof(gtk_entry_get_text(GTK_ENTRY(entry_20)));
      m->sharpness[desc->inum] = atof(gtk_entry_get_text(GTK_ENTRY(entry_21)));
    }
    else if (desc->i3dtype == TRIANGLE)
    {
      m->triangle_color[desc->inum][0] = (float) color.red / (float) G_MAXUINT16;
      m->triangle_color[desc->inum][1] = (float) color.green / (float) G_MAXUINT16;
      m->triangle_color[desc->inum][2] = (float) color.blue / (float) G_MAXUINT16;
      m->triangle_color[desc->inum][3] = new_value;

      m->triangle1[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_00)));
      m->triangle1[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_01)));
      m->triangle1[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_02)));
      m->triangle2[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_10)));
      m->triangle2[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_11)));
      m->triangle2[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_12)));
      m->triangle3[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_20)));
      m->triangle3[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_21)));
      m->triangle3[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_22)));
    }
    else if (desc->i3dtype == SURFACE)
    {
      m->surf_color[desc->inum][0] = (float) color.red / (float) G_MAXUINT16;
      m->surf_color[desc->inum][1] = (float) color.green / (float) G_MAXUINT16;
      m->surf_color[desc->inum][2] = (float) color.blue / (float) G_MAXUINT16;
      m->surf_color[desc->inum][3] = new_value;

      m->surf1[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_00)));
      m->surf1[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_01)));
      m->surf1[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_02)));
      m->surf2[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_10)));
      m->surf2[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_11)));
      m->surf2[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_12)));
      m->surf3[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_20)));
      m->surf3[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_21)));
      m->surf3[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_22)));
    }
    else if (desc->i3dtype == CELL)
    {
      m->cell_color[desc->inum][0] = (float) color.red / (float) G_MAXUINT16;
      m->cell_color[desc->inum][1] = (float) color.green / (float) G_MAXUINT16;
      m->cell_color[desc->inum][2] = (float) color.blue / (float) G_MAXUINT16;
      m->cell_color[desc->inum][3] = new_value;

      m->cell1[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_00)));
      m->cell1[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_01)));
      m->cell1[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_02)));
      m->cell2[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_10)));
      m->cell2[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_11)));
      m->cell2[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_12)));
      m->cell3[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_20)));
      m->cell3[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_21)));
      m->cell3[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_22)));
      m->cell4[desc->inum][0] = atof(gtk_entry_get_text(GTK_ENTRY(entry_30)));
      m->cell4[desc->inum][1] = atof(gtk_entry_get_text(GTK_ENTRY(entry_31)));
      m->cell4[desc->inum][2] = atof(gtk_entry_get_text(GTK_ENTRY(entry_32)));
    }
    else if (desc->i3dtype == TEXTBOX)
    {
      gv_gtk_get_screen_size(&width, &height);

      m->textboxes[desc->inum].coord_x = atoi(gtk_entry_get_text(GTK_ENTRY(entry_10)));
      if (m->textboxes[desc->inum].coord_x > width) m->textboxes[desc->inum].coord_x = width;
      m->textboxes[desc->inum].coord_y = atoi(gtk_entry_get_text(GTK_ENTRY(entry_11)));
      if (m->textboxes[desc->inum].coord_y > height) m->textboxes[desc->inum].coord_y = height;

      m->textboxes[desc->inum].color[0] = (float) color.red/G_MAXUINT16;
      m->textboxes[desc->inum].color[1] = (float) color.green/G_MAXUINT16;
      m->textboxes[desc->inum].color[2] = (float) color.blue/G_MAXUINT16;
      m->textboxes[desc->inum].color[3] = 1.0;
     
      if (m->textboxes[desc->inum].font) free(m->textboxes[desc->inum].font);
      m->textboxes[desc->inum].font = g_strdup(gtk_font_button_get_font_name(GTK_FONT_BUTTON(font_button)));
      if (m->textboxes[desc->inum].message) free(m->textboxes[desc->inum].message);
      m->textboxes[desc->inum].message = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry_00)));
      draw_pixdata_textbox(desc->inum);
    }
  }
  else if (response == GTK_RESPONSE_YES)
  {
    if (desc->i3dtype == SPHERE)
    {
      for(i = desc->inum; i < m->nsphere-1; i++)
      {
	m->sphere_center[i][0] = m->sphere_center[i+1][0];
        m->sphere_center[i][1] = m->sphere_center[i+1][1];
        m->sphere_center[i][2] = m->sphere_center[i+1][2];
        m->sphere_radius[i] = m->sphere_radius[i+1];
        m->sphere_color[i][0] = m->sphere_color[i+1][0];
        m->sphere_color[i][1] = m->sphere_color[i+1][1];
        m->sphere_color[i][2] = m->sphere_color[i+1][2];
        m->sphere_color[i][3] = m->sphere_color[i+1][3];
      }
      m->nsphere--;
      m->sphere_center = (XYZ*) realloc(m->sphere_center, sizeof(XYZ) * m->nsphere);
      m->sphere_radius = (double*) realloc(m->sphere_radius, sizeof(double) * m->nsphere);
      m->sphere_color = (color_t*) realloc(m->sphere_color, sizeof(color_t) * m->nsphere);
    }
    else if (desc->i3dtype == VECTOR)
    {
      for(i = desc->inum; i < m->nvector-1; i++)
      {
	m->vector1[i][0] = m->vector1[i+1][0];
        m->vector1[i][1] = m->vector1[i+1][1];
        m->vector1[i][2] = m->vector1[i+1][2];
	m->vector2[i][0] = m->vector2[i+1][0];
        m->vector2[i][1] = m->vector2[i+1][1];
        m->vector2[i][2] = m->vector2[i+1][2];
        m->radius[i] = m->radius[i+1];
        m->sharpness[i] = m->sharpness[i+1];
        m->vector_color[i][0] = m->vector_color[i+1][0];
        m->vector_color[i][1] = m->vector_color[i+1][1];
        m->vector_color[i][2] = m->vector_color[i+1][2];
        m->vector_color[i][3] = m->vector_color[i+1][3];
      }
      m->nvector--;
      m->vector1 = (XYZ*) realloc(m->vector1, sizeof(XYZ) * m->nvector);
      m->vector2 = (XYZ*) realloc(m->vector2, sizeof(XYZ) * m->nvector);
      m->radius = (double*) realloc(m->radius, sizeof(double) * m->nvector);
      m->sharpness = (double*) realloc(m->sharpness, sizeof(double) * m->nvector);
      m->vector_color = (color_t*) realloc(m->vector_color, sizeof(color_t) * m->nvector);
    }
    else if (desc->i3dtype == TRIANGLE)
    {
      for(i = desc->inum; i < m->ntriangle-1; i++)
      {
	m->triangle1[i][0] = m->triangle1[i+1][0];
        m->triangle1[i][1] = m->triangle1[i+1][1];
        m->triangle1[i][2] = m->triangle1[i+1][2];
	m->triangle2[i][0] = m->triangle2[i+1][0];
        m->triangle2[i][1] = m->triangle2[i+1][1];
        m->triangle2[i][2] = m->triangle2[i+1][2];
	m->triangle3[i][0] = m->triangle3[i+1][0];
        m->triangle3[i][1] = m->triangle3[i+1][1];
        m->triangle3[i][2] = m->triangle3[i+1][2];
        m->triangle_color[i][0] = m->triangle_color[i+1][0];
        m->triangle_color[i][1] = m->triangle_color[i+1][1];
        m->triangle_color[i][2] = m->triangle_color[i+1][2];
        m->triangle_color[i][3] = m->triangle_color[i+1][3];
      }
      m->ntriangle--;
      m->triangle1 = (XYZ*) realloc(m->triangle1, sizeof(XYZ) * m->ntriangle);
      m->triangle2 = (XYZ*) realloc(m->triangle2, sizeof(XYZ) * m->ntriangle);
      m->triangle3 = (XYZ*) realloc(m->triangle3, sizeof(XYZ) * m->ntriangle);
      m->triangle_color = (color_t*) realloc(m->triangle_color, sizeof(color_t) * m->ntriangle);
    }
    else if (desc->i3dtype == SURFACE)
    {
      for(i = desc->inum; i < m->nsurf-1; i++)
      {
	m->surf1[i][0] = m->surf1[i+1][0];
        m->surf1[i][1] = m->surf1[i+1][1];
        m->surf1[i][2] = m->surf1[i+1][2];
	m->surf2[i][0] = m->surf2[i+1][0];
        m->surf2[i][1] = m->surf2[i+1][1];
        m->surf2[i][2] = m->surf2[i+1][2];
	m->surf3[i][0] = m->surf3[i+1][0];
        m->surf3[i][1] = m->surf3[i+1][1];
        m->surf3[i][2] = m->surf3[i+1][2];
        m->surf_color[i][0] = m->surf_color[i+1][0];
        m->surf_color[i][1] = m->surf_color[i+1][1];
        m->surf_color[i][2] = m->surf_color[i+1][2];
        m->surf_color[i][3] = m->surf_color[i+1][3];
      }
      m->nsurf--;
      m->surf1 = (XYZ*) realloc(m->surf1, sizeof(XYZ) * m->nsurf);
      m->surf2 = (XYZ*) realloc(m->surf2, sizeof(XYZ) * m->nsurf);
      m->surf3 = (XYZ*) realloc(m->surf3, sizeof(XYZ) * m->nsurf);
      m->surf_color = (color_t*) realloc(m->surf_color, sizeof(color_t) * m->nsurf);
    }
    else if (desc->i3dtype == CELL)
    {
      for(i = desc->inum; i < m->ncells-1; i++)
      {
    	m->cell1[i][0] = m->cell1[i+1][0];
        m->cell1[i][1] = m->cell1[i+1][1];
        m->cell1[i][2] = m->cell1[i+1][2];
        m->cell2[i][0] = m->cell2[i+1][0];
        m->cell2[i][1] = m->cell2[i+1][1];
        m->cell2[i][2] = m->cell2[i+1][2];
        m->cell3[i][0] = m->cell3[i+1][0];
        m->cell3[i][1] = m->cell3[i+1][1];
        m->cell3[i][2] = m->cell3[i+1][2];
        m->cell4[i][0] = m->cell4[i+1][0];
        m->cell4[i][1] = m->cell4[i+1][1];
        m->cell4[i][2] = m->cell4[i+1][2];
        m->cell_color[i][0] = m->cell_color[i+1][0];
        m->cell_color[i][1] = m->cell_color[i+1][1];
        m->cell_color[i][2] = m->cell_color[i+1][2];
        m->cell_color[i][3] = m->cell_color[i+1][3];    
      }
      m->ncells--;
      m->cell1 = (XYZ*) realloc(m->cell1, sizeof(XYZ) * m->ncells);
      m->cell2 = (XYZ*) realloc(m->cell2, sizeof(XYZ) * m->ncells);
      m->cell3 = (XYZ*) realloc(m->cell3, sizeof(XYZ) * m->ncells);
      m->cell4 = (XYZ*) realloc(m->cell4, sizeof(XYZ) * m->ncells);
      m->cell_color = (color_t*) realloc(m->cell_color, sizeof(color_t) * m->ncells);
    }
    else if (desc->i3dtype == TEXTBOX)
    {
      delete_ith_textbox(m, desc->inum);
    }
  }
  if (GTK_IS_WIDGET(dialog)) gtk_widget_destroy(dialog);

  unset_go_selected();
  rerender_3d();

  if (response != GTK_RESPONSE_REJECT)
  {
    luscus_gtk_update_3Dobject_info();
  }

  redraw();
  append_backup();
}

void luscus_move_textbox_by_click(GtkWidget* button, GTK_3D_DESC* desc)
{
  luscus_set_textbox_number(desc->inum);
  textbox_state = 3;
}

/*menubar functions*/

void callback_open_file(GtkWidget *widget, gpointer data)
{
  int i;
  gchar *tmp;
  GtkWidget *dialog;
  GtkFileFilter *ffilter;
  char *filename;
  char *file_description;

  dialog = gtk_file_chooser_dialog_new("open file", NULL, GTK_FILE_CHOOSER_ACTION_OPEN,
                                       GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                       GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                       NULL);
  if (get_current_directory())
    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), get_current_directory());
  else
    gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), "./");

/*  ffilter = gtk_file_filter_new();
  gtk_file_filter_set_name(ffilter, "all files");
  gtk_file_filter_add_pattern(ffilter, "*");
  gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), ffilter);*/
  
   /*filter for native luscus format*/
  ffilter = gtk_file_filter_new();
  tmp = g_strconcat("*.", input_filetypes[0].extension, NULL);
  gtk_file_filter_set_name(ffilter, input_filetypes[0].description);
  gtk_file_filter_add_pattern(ffilter, tmp);
  g_free(tmp);
  gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), ffilter);

  for(i = 0; i < n_input_types; i++)
    if (input_filetypes[i].forward)
    {
      ffilter = gtk_file_filter_new();
      tmp = g_strconcat("*.", input_filetypes[i].extension, NULL);
      gtk_file_filter_set_name(ffilter, input_filetypes[i].description);
      gtk_file_filter_add_pattern(ffilter, tmp);
      g_free(tmp);
      gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), ffilter);
    }

  if(gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
  {
    filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
    ffilter = gtk_file_chooser_get_filter(GTK_FILE_CHOOSER(dialog));
    file_description = (char*) gtk_file_filter_get_name(GTK_FILE_FILTER(ffilter));

    deallocate_vib_buttons();

   /**/

    if (strcmp(file_description, "all files") == 0) open_file(filename, NULL);
    else open_file(filename, file_description);
    g_free(filename);
  }

  gtk_widget_destroy(dialog);
  set_scale();
  rerender_3d();
  redraw();
}

void callback_save_file(GtkWidget *widget, gpointer data)
{
  gint i;
  gint mode;
  GtkWidget *dialog;
  gchar *filename;
  GtkFileFilter *ffilter;

  gchar *filetype;
  gchar *msg;
  mode = GPOINTER_TO_INT(data);

  filename = g_strdup(get_input_filename());
  if (filename == NULL) mode = 1;

  if (mode == 0)
  {
    luscus_gtk_do_save_file(filename, get_filetype(filename));
  /*  filename = ...*/
    /*close input*/
    /*luscus_gtk_do_save_file*/
  }
  else
  {
    dialog = gtk_file_chooser_dialog_new("Save file", NULL, GTK_FILE_CHOOSER_ACTION_SAVE,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                         GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
                                         NULL);

    gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);
    set_save_file_filename(dialog, mode);

    if(gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      ffilter = gtk_file_chooser_get_filter(GTK_FILE_CHOOSER(dialog));
      filetype = (gchar*) gtk_file_filter_get_name(ffilter);

      for (i = 0; i < n_input_types; i++)
        if (g_strcmp0(filetype, input_filetypes[i].description) == 0) break;

      g_free(filename);
      filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      luscus_gtk_do_save_file(filename, i);
    }

    gtk_widget_destroy(dialog);
  }

  msg = g_strdup_printf("File saved to: %s", filename);
  luscus_gtk_pop_message_from_statusbar2();
  luscus_gtk_push_message_to_statusbar2(msg);
  g_free(msg);
  g_free(filename);
}

void set_save_file_filename(GtkWidget *dialog, gint mode)
{
  gint i;
  gchar *fname;
  GtkFileFilter *ffilter;

  /*filter for luscus file*/
  ffilter = gtk_file_filter_new();
  gtk_file_filter_set_name(ffilter, input_filetypes[0].description);
  gtk_file_filter_add_pattern(ffilter, input_filetypes[0].extension);
  gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), ffilter);

  /*filters for other files*/
  for(i = 0; i < n_input_types; i++)
    if (input_filetypes[i].backward)
    {
      ffilter = gtk_file_filter_new();
      gtk_file_filter_set_name(ffilter, input_filetypes[i].description);
      gtk_file_filter_add_pattern(ffilter, input_filetypes[i].extension);
      gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), ffilter);
    }

  fname = nextfname();

  gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER(dialog), "./");
  if (fname == NULL) gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER(dialog), "newfile");
  else gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER(dialog), fname);

  free(fname);
}

void callback_close_file(GtkWidget* widget, gpointer data)
{
  close_file();
  deallocate_vib_buttons();
  rerender_3d();
  redraw();
  luscus_gtk_update_3Dobject_info();
}

void callback_do_undo(GtkWidget *widget, gpointer data)
{
  luscus_gtk_pop_message_from_statusbar2();
  get_backup();
/*  set_scale();*/
  rerender_3d();
  redraw();
  luscus_gtk_update_3Dobject_info();
}

void callback_make_editable(GtkWidget *widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget)))
  {
    m->editable = 1;
    luscus_gtk_enable_editing();
  }
  else
  {
    m->editable = 0;
    luscus_gtk_diseable_editing();
  }


}

void callback_change_background(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *colorsel;
  GdkColor *background_color;

  background_color = (GdkColor*) g_malloc(sizeof(GdkColor));
/*convert color_t to the GdkColor*/
  color_t_2_GdkColor(background_color, Input_Data.background_color);

  dialog = gtk_color_selection_dialog_new("Background color");
/*  colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(dialog)->colorsel);*/
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), background_color);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), background_color);
    GdkColor_2_color_t(background_color, Input_Data.background_color);
  }
  gtk_widget_destroy (dialog);
  redraw();

  g_free(background_color);
}

/* conversion of color_t in OpenGL format to the GTK format (GdkColor) */
void color_t_2_GdkColor(GdkColor *out_color, color_t in_color)
{
  out_color->red = (guint16) G_MAXUINT16 * ((double) in_color[0]);
  out_color->green = (guint16) G_MAXUINT16 * ((double) in_color[1]);
  out_color->blue = (guint16) G_MAXUINT16 * ((double) in_color[2]);
  out_color->pixel = (guint32) G_MAXUINT32 * ((double) in_color[3]);
}
/* conversion of the GTK format (GdkColor) to the color_t OpenGL[4] */

void GdkColor_2_color_t(GdkColor *in_color, color_t out_color)
{
  out_color[0] = (double) in_color->red / (double) G_MAXUINT16;
  out_color[1] = (double) in_color->green / (double) G_MAXUINT16;
  out_color[2] = (double) in_color->blue / (double) G_MAXUINT16;
  out_color[3] = (double) in_color->pixel / (double) G_MAXUINT32;
}

void callback_change_label(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *colorsel;
  GdkColor *label_color;

  label_color = (GdkColor*) g_malloc(sizeof(GdkColor));
  color_t_2_GdkColor(label_color, Input_Data.label_color);

  dialog = gtk_color_selection_dialog_new("Label color");
/*  colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(dialog)->colorsel);*/
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), label_color);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), label_color);
    GdkColor_2_color_t(label_color, Input_Data.label_color);
  }
  gtk_widget_destroy (dialog);
  redraw();

  g_free(label_color);
}

void callback_change_pos_orbital(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  /*GtkColorSelection*/ GtkWidget *colorsel;
  GdkColor *orbital_color;

  orbital_color = (GdkColor*) g_malloc(sizeof(GdkColor));
  color_t_2_GdkColor(orbital_color, Input_Data.neg_pos_color[0]);

  dialog = gtk_color_selection_dialog_new("Orbital color");
 /* colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(dialog)->colorsel);*/
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), orbital_color);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), orbital_color);
    GdkColor_2_color_t(orbital_color, Input_Data.neg_pos_color[0]);
  }
  gtk_widget_destroy(dialog);
  rerender_3d();
  redraw();

  g_free(orbital_color);
}

void callback_change_neg_orbital(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  /*GtkColorSelection*/ GtkWidget *colorsel;
  GdkColor *orbital_color;

  orbital_color = (GdkColor*) g_malloc(sizeof(GdkColor));
  color_t_2_GdkColor(orbital_color, Input_Data.neg_pos_color[1]);

  dialog = gtk_color_selection_dialog_new("Orbital color");
 /* colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(dialog)->colorsel);*/
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), orbital_color);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), orbital_color);
    GdkColor_2_color_t(orbital_color, Input_Data.neg_pos_color[1]);
  }
  gtk_widget_destroy(dialog);
  rerender_3d();
  redraw();

  g_free(orbital_color);
}

void callback_change_pos_epot(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *colorsel;
  GdkColor *epot_color;

  epot_color = (GdkColor*) g_malloc(sizeof(GdkColor));
  color_t_2_GdkColor(epot_color, Input_Data.electrostatic_poten_color[0]);

  dialog = gtk_color_selection_dialog_new("Electrostatic potential color");
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), epot_color);

  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), epot_color);
    GdkColor_2_color_t(epot_color, Input_Data.electrostatic_poten_color[0]);
  }
  gtk_widget_destroy(dialog);
  rerender_3d();
  redraw();

  g_free(epot_color);
}

void callback_change_neg_epot(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *colorsel;
  GdkColor *epot_color;

  epot_color = (GdkColor*) g_malloc(sizeof(GdkColor));
  color_t_2_GdkColor(epot_color, Input_Data.electrostatic_poten_color[1]);

  dialog = gtk_color_selection_dialog_new("Electrostatic potential color");
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), epot_color);

  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), epot_color);
    GdkColor_2_color_t(epot_color, Input_Data.electrostatic_poten_color[1]);
  }
  gtk_widget_destroy(dialog);
  rerender_3d();
  redraw();

  g_free(epot_color);
}

void callback_change_plain(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  /*GtkColorSelection*/ GtkWidget *colorsel;
  GdkColor *plain_color;

  plain_color = (GdkColor*) g_malloc(sizeof(GdkColor));
  color_t_2_GdkColor(plain_color, Input_Data.extracolor);

  dialog = gtk_color_selection_dialog_new("Plain color");
 /* colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(dialog)->colorsel);*/
  colorsel = gtk_color_selection_dialog_get_color_selection(GTK_COLOR_SELECTION_DIALOG(dialog));
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(colorsel), plain_color);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(colorsel), plain_color);
    GdkColor_2_color_t(plain_color, Input_Data.extracolor);
  }
  gtk_widget_destroy (dialog);
  rerender_3d();
  redraw();

  g_free(plain_color);
}

void callback_change_font(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;

#ifdef GTK2
  dialog = gtk_font_selection_dialog_new("Select font");
  gtk_font_selection_dialog_set_font_name(GTK_FONT_SELECTION_DIALOG(dialog), Input_Data.font);
#endif
#ifdef GTK3
  dialog = gtk_font_chooser_dialog_new("Select font", NULL);
  gtk_font_chooser_set_font(GTK_FONT_CHOOSER(dialog), Input_Data.font);
#endif

  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
  {
#ifdef GTK2
    g_strlcpy(Input_Data.font, gtk_font_selection_dialog_get_font_name(GTK_FONT_SELECTION_DIALOG(dialog)), 100);
#endif
#ifdef GTK3
    g_strlcpy(Input_Data.font, gtk_font_chooser_get_preview_text(GTK_FONT_CHOOSER(dialog)), 100);
#endif
    draw_all_pixdata();
  }

  gtk_widget_destroy(dialog);
}

void callback_adjust_atom_properties(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *hbox;
  GtkWidget *label;
  GtkWidget *persys;
  GtkWidget *separator;
  GdkColor color;
  gint response;
  gint i, j;


#ifdef GTK2
  GtkObject *adj_atom_size, *adj_bond_size, *adj_valency;
#endif
#ifdef GTK3
  GtkAdjustment *adj_atom_size, *adj_bond_size, *adj_valency;
#endif

  dialog = gtk_dialog_new_with_buttons("Atom properties", GTK_WINDOW(gtk_widget_get_toplevel(widget)),
                                       GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
				       GTK_STOCK_APPLY, GTK_RESPONSE_APPLY,
                                       GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
                                       GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                       NULL);

  tmp_elem = (ELEM_DATA*) g_malloc(sizeof(ELEM_DATA) * number_of_elements);

  /*copy data to the newly allocated array; Direct data changing will be on temporary array; this way, user can change properties on more than one atom*/
  for(i = 0; i < number_of_elements; i++)
  {
    tmp_elem[i].name = g_strdup(e[i].name);
    for(j = 0; j < 4; j++)
      tmp_elem[i].color[j] = e[i].color[j];
    tmp_elem[i].vdw_rad = e[i].vdw_rad;
    tmp_elem[i].bond_rad = e[i].bond_rad;
    tmp_elem[i].valency = e[i].valency;
    tmp_elem[i].periodic_pos_x = e[i].periodic_pos_x;
    tmp_elem[i].periodic_pos_y = e[i].periodic_pos_y;
  }

  persys = luscus_gtk_make_periodic_system1(tmp_elem);

  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(dialog))), persys, TRUE, TRUE, 0);

#ifdef GTK2
  separator = gtk_vseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(dialog))), separator, TRUE, TRUE, 0);
  gtk_widget_show(separator);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);  
#endif

  atom_prop_label = gtk_label_new(tmp_elem[elem_selected].name);
  gtk_box_pack_start(GTK_BOX(hbox), atom_prop_label, FALSE, FALSE, 0);
  gtk_widget_show(atom_prop_label);

  color_t_2_GdkColor(&color, tmp_elem[elem_selected].color);
  atom_prop_color_button = gtk_color_button_new();
  gtk_color_button_set_color(GTK_COLOR_BUTTON(atom_prop_color_button), &color);
  gtk_box_pack_start(GTK_BOX(hbox), atom_prop_color_button, FALSE, FALSE, 5);
  gtk_widget_show(atom_prop_color_button);

  label = gtk_label_new("Van der Waals radius:");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
  gtk_widget_show(label);

  adj_atom_size = gtk_adjustment_new(0.0, 0.F, 3.F, 0.01F, 0.1F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_atom_size), tmp_elem[elem_selected].vdw_rad);

#ifdef GTK2
  atom_prop_atom_size_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj_atom_size), 0.01, 2);
#endif
#ifdef GTK3
  atom_prop_atom_size_spin = gtk_spin_button_new(adj_atom_size, 0.01, 2);
#endif
  g_signal_connect(G_OBJECT(atom_prop_atom_size_spin), "value-changed", G_CALLBACK(atom_prop_atom_size_spin_callback), NULL);
  gtk_box_pack_start(GTK_BOX(hbox), atom_prop_atom_size_spin, FALSE, FALSE, 0);
  gtk_widget_show(atom_prop_atom_size_spin);

  label = gtk_label_new("bonding radius:");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
  gtk_widget_show(label);

  adj_bond_size = gtk_adjustment_new(0.0, 0.F, 3.F, 0.01F, 0.1F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_bond_size), tmp_elem[elem_selected].bond_rad);

#ifdef GTK2
  atom_prop_bond_size_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj_bond_size), 0.01, 2);
#endif
#ifdef GTK3
  atom_prop_bond_size_spin = gtk_spin_button_new(adj_bond_size, 0.01, 2);
#endif
  g_signal_connect(G_OBJECT(atom_prop_bond_size_spin), "value-changed", G_CALLBACK(atom_prop_bond_size_spin_callback), NULL);
  gtk_box_pack_start(GTK_BOX(hbox), atom_prop_bond_size_spin, FALSE, FALSE, 0);
  gtk_widget_show(atom_prop_bond_size_spin);

  label = gtk_label_new("valency:");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
  gtk_widget_show(label);

  adj_valency = gtk_adjustment_new(0.0, 0.F, 30.F, 1.0F, 1.0F, 0.F);
  gtk_adjustment_set_value(GTK_ADJUSTMENT(adj_valency), tmp_elem[elem_selected].valency);

#ifdef GTK2
  atom_prop_atom_valency_spin = gtk_spin_button_new(GTK_ADJUSTMENT(adj_valency), 0.01, 2);
#endif
#ifdef GTK3
  atom_prop_atom_valency_spin = gtk_spin_button_new(adj_valency, 0.01, 2);
#endif
  g_signal_connect(G_OBJECT(atom_prop_atom_valency_spin), "value-changed", G_CALLBACK(atom_prop_atom_valency_spin_callback), NULL);
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(atom_prop_atom_valency_spin), 0);
  gtk_box_pack_start(GTK_BOX(hbox), atom_prop_atom_valency_spin, FALSE, FALSE, 0);
  gtk_widget_show(atom_prop_atom_valency_spin);

  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(dialog))), hbox, TRUE, TRUE, 0);
  gtk_widget_show(hbox);

#ifdef GTK2
  separator = gtk_hseparator_new();
#endif
#ifdef GTK3
  separator = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
#endif
  gtk_box_pack_start(GTK_BOX(gtk_dialog_get_content_area(GTK_DIALOG(dialog))), separator, TRUE, TRUE, 0);
  gtk_widget_show(separator);

  response = gtk_dialog_run(GTK_DIALOG(dialog));

  if (response == GTK_RESPONSE_ACCEPT || response == GTK_RESPONSE_APPLY)
  {
    gtk_color_button_get_color(GTK_COLOR_BUTTON(atom_prop_color_button), &color);
    GdkColor_2_color_t(&color, e[elem_selected].color);
    e[elem_selected].vdw_rad = tmp_elem[elem_selected].vdw_rad;
    e[elem_selected].bond_rad = tmp_elem[elem_selected].bond_rad;
    e[elem_selected].valency = tmp_elem[elem_selected].valency;

    /*copy data to the newly allocated array; Direct data changing will be on temporary array; this way, user can change properties on more than one atom*/
/*    for(i = 0; i < number_of_elements; i++)
    {
      for(j = 0; j < 4; j++)
        e[i].color[j] = tmp_elem[i].color[j];
      e[i].vdw_rad = tmp_elem[i].vdw_rad;
      e[i].bond_rad = tmp_elem[i].bond_rad;
      e[i].valency = tmp_elem[i].valency;
    }*/
    for (i = 0; i < m->natom; i++)
    {
      for (j = 0; j < number_of_elements && g_strcmp0(e[j].name, m->elem[i].name) != 0; j++); /*searching for correct element*/
      m->elem[i].color[0] = e[j].color[0];
      m->elem[i].color[1] = e[j].color[1];
      m->elem[i].color[2] = e[j].color[2];
      m->elem[i].color[3] = e[j].color[3];
      m->elem[i].vdw_rad = e[j].vdw_rad;
      m->elem[i].bond_rad = e[j].bond_rad;
      m->elem[i].valency = e[j].valency;
    }

    if (response == GTK_RESPONSE_ACCEPT) save_atom_data();
  }
  if (GTK_IS_WIDGET(dialog)) gtk_widget_destroy(dialog);

  /*free temporary storage*/

  for(i = 0; i < number_of_elements; i++) g_free(tmp_elem[i].name);
  g_free(tmp_elem);
  rerender_3d();
  redraw();
}

GtkWidget *luscus_gtk_make_periodic_system1(ELEM_DATA* tmp_e)
{
  GtkWidget *table;
  GtkWidget *button;
  gint i;
  gint maxr = 0, maxc = 0;
  
  for(i = 1; i < number_of_elements; i++)
  {
    if (e[i].periodic_pos_x > maxc) maxc = e[i].periodic_pos_x;
    if (e[i].periodic_pos_y > maxr) maxr = e[i].periodic_pos_y;
  }

  table = gtk_table_new(maxr + 1, maxc + 1, TRUE);

  for(i = 1; i < number_of_elements; i++)
  {
    button = gtk_button_new_with_label(e[i].name);
    g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(luscus_gtk_atom_prop_select_element), (gpointer) GINT_TO_POINTER(i));
    gtk_table_attach(GTK_TABLE(table), button, e[i].periodic_pos_x, e[i].periodic_pos_x+1, e[i].periodic_pos_y, e[i].periodic_pos_y+1, GTK_FILL, GTK_FILL, 0, 0);
    gtk_widget_show(button);
  }
  gtk_widget_show(table);

  return table;
}

void callback_change_automatic_bonding(GtkWidget *widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget)))
  {
    Input_Data.automatic_rebonding = 1;
    rebond();
    rerender_3d();
    redraw();
  }
  else
    Input_Data.automatic_rebonding = 0;
}

void luscus_gtk_atom_prop_select_element(GtkWidget *widget, gpointer data)
{
  gint iatom = GPOINTER_TO_INT(data);
  GdkColor color;
  elem_selected = iatom; /*atoms are numbered according to the numeration in the atoms.c -> H is #1*/

  color_t_2_GdkColor(&color, tmp_elem[iatom].color);

  gtk_label_set_label(GTK_LABEL(atom_prop_label), tmp_elem[iatom].name);
  gtk_color_button_set_color(GTK_COLOR_BUTTON(atom_prop_color_button), &color);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(atom_prop_atom_size_spin), tmp_elem[iatom].vdw_rad);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(atom_prop_bond_size_spin), tmp_elem[iatom].bond_rad);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(atom_prop_atom_valency_spin), tmp_elem[iatom].valency);

}

void atom_prop_atom_size_spin_callback(GtkWidget *spin, gpointer data)
{
  tmp_elem[elem_selected].vdw_rad = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));
}

void atom_prop_bond_size_spin_callback(GtkWidget *spin, gpointer data)
{
  tmp_elem[elem_selected].bond_rad = gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));
}

void atom_prop_atom_valency_spin_callback(GtkWidget *spin, gpointer data)
{
  tmp_elem[elem_selected].valency = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin));
}

void callback_save_settings(GtkWidget *widget, gpointer data)
{
  save_settings_data();
}

void callback_move_light(GtkWidget *checkmenu, gpointer data)
{
  luscus_gtk_move_light_state();
}

void callback_create_viewpoint(GtkWidget *widget, gpointer data)
{
  ViewPoint(0);
}

void callback_restore_viewpoint(GtkWidget *widget, gpointer data)
{
  ViewPoint(1);
}

void callback_grayscale(GtkWidget *checkmenu, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(checkmenu))) Input_Data.bw = 1;
  else Input_Data.bw = 0;
  rerender_3d();
  redraw();
}

void callback_animate(GtkWidget *menuitem, gpointer data)
{
  if (Input_Data.animate)
  {
    Input_Data.animate = FALSE;
    gtk_image_menu_item_set_image(GTK_IMAGE_MENU_ITEM(menuitem), gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY, GTK_ICON_SIZE_SMALL_TOOLBAR));
    gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Animate");
  }
  else
  {
    gtk_image_menu_item_set_image(GTK_IMAGE_MENU_ITEM(menuitem), gtk_image_new_from_stock(GTK_STOCK_MEDIA_STOP, GTK_ICON_SIZE_SMALL_TOOLBAR));
    gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Stop animation");
    Input_Data.animate = TRUE;
    luscus_gtk_start_animation(0);
  }
  redraw();
}

void luscus_gtk_start_animation(int i)
{
  g_timeout_add((guint) (1000.0 * Input_Data.snapshot_delay), luscus_start_animation1, NULL);
}

gboolean luscus_start_animation1(gpointer data)
{
  animate(0);
  if (Input_Data.animate) return TRUE;
  else return FALSE;
}

void luscus_gtk_show_atoms(GtkWidget* menu_item, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menu_item))) Input_Data.hide_atoms = 0;
  else Input_Data.hide_atoms = 1;
  rerender_3d();
  redraw();
}

void luscus_gtk_show_bonds(GtkWidget* menu_item, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menu_item))) Input_Data.hide_bonds = 0;
  else Input_Data.hide_bonds = 1;
  rerender_3d();
  redraw();
}

void luscus_gtk_show_axes(GtkWidget* menu_item, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menu_item))) Input_Data.show_axis = 1;
  else Input_Data.show_axis = 0;
  rerender_3d();
  redraw();
}

void luscus_gtk_callback_show_dot_line(GtkWidget* menu_item, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menu_item))) Input_Data.style_atoms = 1;
  else Input_Data.style_atoms = 0;
  rerender_3d();
  redraw();
}

void luscus_gtk_show_labels(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 0;
  else Input_Data.label_atoms &= ~(1 << 0);
  redraw();
}

void luscus_gtk_show_indices(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 1;
  else Input_Data.label_atoms &= ~(1 << 1);
  redraw();
}

void luscus_gtk_show_numeration(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 2;
  else Input_Data.label_atoms &= ~(1 << 2);
  redraw();
}

void luscus_gtk_show_symbols(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 3;
  else Input_Data.label_atoms &= ~(1 << 3);
  redraw();
}

void luscus_gtk_show_names(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 4;
  else Input_Data.label_atoms &= ~(1 << 4);
  redraw();
}

void luscus_gtk_show_mulliken(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 5;
  else Input_Data.label_atoms &= ~(1 << 5);
  redraw();
}

void luscus_gtk_show_loprop(GtkWidget* widget, gpointer data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(widget))) Input_Data.label_atoms |= 1 << 6;
  else Input_Data.label_atoms &= ~(1 << 6);
  redraw();
}

void do_screenshot(int out)
{
  gchar *newfilename = NULL;
  int ifile = 0;

/*  newfilename = nextfname();*/
  do
  {
    printf("G_FREE(newfilename)0\n"); fflush(stdout);
    g_free(newfilename);
    if (out == 0) newfilename = g_strdup_printf("new_file_%.3d.tga", ifile);
    else if (out == 1) newfilename = g_strdup_printf("new_file_%.3d.eps", ifile);
    else if (out == 2) newfilename = g_strdup_printf("new_file_%.3d.eps", ifile);
    else if (out == 3) newfilename = g_strdup_printf("new_file_%.3d.pov", ifile);
    ifile++;
  } while(get_file_exist(newfilename) && ifile < 1000);

  luscus_gtk_push_message_to_statusbar2("screenshot done");

  if (out == 0) /*luscus_save_builtin_graphic_format(newfilename, "tga", Input_Data.init_screen_size, Input_Data.init_screen_size);*/ tgaSave(newfilename, Input_Data.init_screen_size, Input_Data.init_screen_size);
  else if (out == 1) initPS(0, newfilename);
  else if (out == 2) initPS(1, newfilename);
  else if (out == 3) initPovray(newfilename);

  printf("G_FREE(newfilename)1\n"); fflush(stdout);
  g_free(newfilename);
}

void callback_screenshot(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *combo;
  GtkWidget *hbox, *vbox, *ihbox;
  GtkWidget *label;
  GtkWidget *spin_button;
#ifdef GTK_GLEXT
#ifdef GTK2
  GtkObject *adj_w, *adj_h;
#elif GTK3
  GtkAdjustment *adj_w, *adj_h;
#endif
#else
  GdkWindow *wind;
  GdkPixbuf *pixbuf;
  GtkAllocation allocation;
#endif
  gint result;
  gchar *filename;
  gint num_file_saved;
  gint width, height;

/* ----- check writable formats ----- */

  GSList *formats = gdk_pixbuf_get_formats();
  GSList *writable_formats = NULL;
  GSList *lwritable_formats = NULL;
  GdkPixbufFormat *gpixformtmp;

  g_slist_foreach(formats, (GFunc) luscus_check_writable_pic_files, &writable_formats);
  g_slist_free(formats);

  lwritable_formats = writable_formats;

/*-----------------------------------*/

/* ----- check screen size ----- */

  gv_gtk_get_screen_size(&width, &height);

/* ----------------------------- */

      /*adj_w, *adj_h;*/
#ifdef GTK_GLEXT
  adj_w = gtk_adjustment_new((gdouble) width,  600.0, 4000.0, 1.0, 10.0, 0.0);
  adj_h = gtk_adjustment_new((gdouble) height, 600.0, 4000.0, 1.0, 10.0, 0.0);
#else
  /*GET PIXBUF BEFORE THE DIALOG IS INITIATED*/
  wind = gtk_widget_get_window(drawingarea);
  gtk_widget_get_allocation(drawingarea, &allocation);
  width = allocation.width;
  height = allocation.height;
#ifdef GTK2
  pixbuf = gdk_pixbuf_get_from_drawable(NULL, GDK_DRAWABLE(wind), NULL, 0, 0, 0, 0, width, height);
#endif
#ifdef GTK3
  pixbuf = gdk_pixbuf_get_from_window(wind, 0, 0, width, height);
#endif
#endif

  dialog = gtk_file_chooser_dialog_new("Save screenshot image", NULL, GTK_FILE_CHOOSER_ACTION_SAVE,
                                       GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL, GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT, NULL);

  /*  luscus_set_filename(NULL, dialog);*/
  gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), FALSE);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);  
#endif

#ifdef GTK_GLEXT
#ifdef GTK2
  vbox = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 0);

#ifdef GTK2
  ihbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  ihbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(ihbox), FALSE);  
#endif
  gtk_box_pack_start(GTK_BOX(vbox), ihbox, FALSE, FALSE, 0);

  label = gtk_label_new("Width:");
  gtk_box_pack_start(GTK_BOX(ihbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

#ifdef GTK2
  spin_button = gtk_spin_button_new(GTK_ADJUSTMENT(adj_w), 1.0, 0);
#elif GTK3
  spin_button = gtk_spin_button_new(adj_w, 1.0, 0);
#endif
  gtk_box_pack_start(GTK_BOX(ihbox), spin_button, FALSE, FALSE, 0);
  gtk_widget_show(spin_button);

  gtk_widget_show(ihbox);


#ifdef GTK2
  ihbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  ihbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(ihbox), FALSE);  
#endif

  gtk_box_pack_start(GTK_BOX(vbox), ihbox, FALSE, FALSE, 0);

  label = gtk_label_new("Height:");
  gtk_box_pack_start(GTK_BOX(ihbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

#ifdef GTK2
  spin_button = gtk_spin_button_new(GTK_ADJUSTMENT(adj_h), 1.0, 0);
#elif GTK3
  spin_button = gtk_spin_button_new(adj_h, 1.0, 0);
#endif
  gtk_box_pack_start(GTK_BOX(ihbox), spin_button, FALSE, FALSE, 0);
  gtk_widget_show(spin_button);


  label = gtk_label_new("pixels");
  gtk_box_pack_start(GTK_BOX(ihbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);

  gtk_widget_show(ihbox);

  gtk_widget_show(vbox);
#endif

#ifdef GTK_OLD
  combo = gtk_combo_box_new_text();
#else
  combo = gtk_combo_box_text_new_with_entry();
#endif
  gtk_box_pack_start(GTK_BOX(hbox), combo, TRUE, TRUE, 0);

  gtk_file_chooser_set_preview_widget(GTK_FILE_CHOOSER(dialog), FALSE);
  gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(dialog), hbox);
  gtk_widget_show(hbox);
  gtk_widget_show(combo);
#ifdef GTK_OLD
  gtk_combo_box_append_text(GTK_COMBO_BOX(combo), "TarGA (*.tga)");
  gtk_combo_box_append_text(GTK_COMBO_BOX(combo), "postscript (*.eps)");
  gtk_combo_box_append_text(GTK_COMBO_BOX(combo), "postscript level 2 (*.eps)");
  gtk_combo_box_append_text(GTK_COMBO_BOX(combo), "povray (*.pov)");
  while(lwritable_formats)
  {
    gpixformtmp = lwritable_formats->data;
    gtk_combo_box_append_text(GTK_COMBO_BOX(combo), gdk_pixbuf_format_get_description(gpixformtmp));
    lwritable_formats = lwritable_formats->next;
  }
#else
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), "TarGA (*.tga)");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), "postscript (*.eps)");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), "postscript level 2 (*.eps)");
  gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), "povray (*.pov)");
  while(lwritable_formats)
  {
    gpixformtmp = lwritable_formats->data;
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo), gdk_pixbuf_format_get_description(gpixformtmp));
    lwritable_formats = lwritable_formats->next;
  }
#endif

  gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);

/*  gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), fname);*/

  gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), "./");
/*  luscus_set_filename(combo, dialog);*/
/*  g_signal_connect(G_OBJECT(combo), "changed", G_CALLBACK(luscus_set_filename), (gpointer) dialog);*/
/*-----------------ERROR CALLING LUSCUS_SET_FILENAME-----------SEGMENTATION FAULT!!!-------------*/

  result = gtk_dialog_run(GTK_DIALOG(dialog));

  if(result == GTK_RESPONSE_ACCEPT)
  {

#ifdef GTK_GLEXT
    width = (gint) gtk_adjustment_get_value(GTK_ADJUSTMENT(adj_w));
    height = (gint) gtk_adjustment_get_value(GTK_ADJUSTMENT(adj_h));
#endif

    filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
#ifdef EBUG
    printf("filename = |%s|\n", filename); fflush(stdout);
#endif
    num_file_saved = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));

    if (num_file_saved > 3)
    {
      lwritable_formats = writable_formats;
      while(lwritable_formats)
      {
        gpixformtmp = lwritable_formats->data;
        /*gtk_combo_box_get_active_text has been deprecated since version 2.24 ... see manual...*/
#ifdef GTK2
#ifdef GTK_OLD
        if (g_strcmp0(gtk_combo_box_get_active_text(GTK_COMBO_BOX(combo)), gdk_pixbuf_format_get_description(gpixformtmp)) == 0)
#else
        if (g_strcmp0(gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo)), gdk_pixbuf_format_get_description(gpixformtmp)) == 0)
#endif
#endif
#ifdef GTK3
        if (g_strcmp0(gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo)), gdk_pixbuf_format_get_description(gpixformtmp)) == 0)       
#endif
        {
          gtk_widget_destroy(dialog);
          dialog = NULL;
#ifdef GTK_GLEXT
          luscus_save_builtin_graphic_format(filename, gdk_pixbuf_format_get_name(gpixformtmp), width, height);
#else
          luscus_save_builtin_graphic_format(filename, gdk_pixbuf_format_get_name(gpixformtmp), pixbuf, width, height);
#endif
          break;
        }
        lwritable_formats = lwritable_formats->next;
      }

    /*pronadji dodatne formate file-ova!*/

    }
    else
    {
      switch(num_file_saved)
      {
        case SCSHOT_TGA:
          gtk_widget_destroy(dialog);
          dialog = NULL;
#ifdef GTK_GLEXT
          luscus_save_builtin_graphic_format(filename, "tga", width, height);/*tgaSave(filename);*/ break;
#else
          luscus_save_builtin_graphic_format(filename, "tga", pixbuf, width, height); break;
#endif
        case SCSHOT_PS: initPS(0, filename); break;
        case SCSHOT_PS2: initPS(1, filename); break;
        case SCSHOT_POV: initPovray(filename); break;
#ifdef GTK_GLEXT
          luscus_save_builtin_graphic_format(filename, "tga", width, height);/*tgaSave(filename);*/ break;
#else
          luscus_save_builtin_graphic_format(filename, "tga", pixbuf, width, height); break;
#endif
      }
    }
    g_free(filename);
  }
  g_slist_free(writable_formats);
#ifndef GTK_GLEXT
  g_object_unref(G_OBJECT(pixbuf));
#endif
  if (dialog && GTK_IS_WIDGET(dialog)) gtk_widget_destroy(dialog);
}

void luscus_set_filename(GtkWidget *combo, GtkWidget *dialog)
{
  gchar *fname = NULL;
  int ifile = 0;
  gint num_file_saved; /*i-th screenshot file format*/

  GSList *formats = gdk_pixbuf_get_formats();
  GSList *writable_formats = NULL;
  GSList *lwritable_formats = NULL;
  GdkPixbufFormat *gpixformtmp;

  num_file_saved = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
   
/*  newfilename = nextfname();*/
  do
  {
      printf("_FREE(FNAME)\n"); fflush(stdout);
      g_free(fname);
      fname = g_strdup_printf("new_file_%.3d", ifile);
      printf("num_file_saved = %d\n", num_file_saved); fflush(stdout);

      if(num_file_saved > 3)
      {

/* find out extension here! */

        g_slist_foreach(formats, (GFunc) luscus_check_writable_pic_files, &writable_formats);
/*        g_slist_free(formats);*/

        lwritable_formats = writable_formats;

        while(lwritable_formats)
        {
          gpixformtmp = lwritable_formats->data;
          /*gtk_combo_box_get_active_text has been deprecated since version 2.24 ... see manual...*/
#ifdef GTK_OLD          
          if (g_strcmp0(gtk_combo_box_get_active_text(GTK_COMBO_BOX(combo)), gdk_pixbuf_format_get_description(gpixformtmp)) == 0)
#else
          if (g_strcmp0(gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo)), gdk_pixbuf_format_get_description(gpixformtmp)) == 0)
#endif
          {
            fname = g_strconcat(fname, ".", gdk_pixbuf_format_get_name(gpixformtmp), NULL);
            break;
          }
          lwritable_formats = lwritable_formats->next;
        }

        g_slist_free(writable_formats); 
      }
      else
      {
        switch(num_file_saved)
        {
          case SCSHOT_TGA: 
            fname = g_strconcat(fname, ".tga", NULL);
            break;
          case SCSHOT_PS:
            fname = g_strconcat(fname, ".eps", NULL);
            break;
          case SCSHOT_PS2:
            fname = g_strconcat(fname, ".eps", NULL);
            break;
          case SCSHOT_POV:
            fname = g_strconcat(fname, ".pov", NULL);
            break;
          default:
            fname = g_strconcat(fname, ".tga", NULL);
        }
      }
      printf("checking file |%s|\n", fname); fflush(stdout);
  
      ifile++;
  } while(get_file_exist(fname) && ifile < 1000);

  g_slist_free(formats);

  gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), fname);
  printf("last fname = |%s|\n", fname); fflush(stdout);
}

void luscus_gtk_goto_freq(GtkWidget *widget, gpointer data)
{
  gint tvib = GPOINTER_TO_INT(data);
  ivib = tvib;
  set_current_graph_data();
  redraw();
}

void luscus_start_calculation(GtkWidget* menu, gint icalc)
{
  int iarg = 0;
  int ic = 0;
  char num[4];
  char *dirname;

  gboolean ret;
  gchar *calc_filename = NULL;
  gchar *extra_filename = NULL;
  gchar *gargv[7];
  GPid gpid;
  GError *gerror = NULL;
  gchar *tmp = NULL;
  gint run;
  GtkWidget *dialog;
  GtkWidget *actionarea;
  GtkWidget *hbox;
  GtkWidget *label;
  GtkWidget *file_chooser;
  GtkWidget *entry;
  GtkFileFilter *ffilter;

  /*1. determine file name*/

  dirname = get_current_directory();
  do
  {
    if (calc_filename) g_free(calc_filename);
    ic++;
    sprintf(num, "%03d", ic);
#ifdef WINDOWS
    calc_filename = g_strconcat("\"", g_get_current_dir(), "\\","newfile_", num,".lus", "\"", NULL);
#else
    calc_filename = g_strconcat(g_get_current_dir(), "/newfile_", num,".lus", NULL);
#endif
  }
  while(get_file_exist(calc_filename));

#ifdef EBUG
  printf("the final name is |%s|\n", calc_filename);
#endif

  /*1a. setup plugin name*/
  gargv[0] = (gchar*) calc_defs[icalc].plugin_name;
  iarg++;
  /*2. determine name of the second file and/or argument if necessary*/


  if (calc_defs[icalc].extrafile || calc_defs[icalc].extraargsym && calc_defs[icalc].extraargval)
  {
    dialog = gtk_dialog_new_with_buttons("Calculation parameters", NULL,
                                         GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                         GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                         GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT, NULL);

    actionarea = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    if (calc_defs[icalc].extrafile)
    {
      file_chooser = gtk_file_chooser_button_new("compare density with this file", GTK_FILE_CHOOSER_ACTION_OPEN);
      gtk_box_pack_start(GTK_BOX(actionarea), file_chooser, FALSE, FALSE, 0);
      if (calc_defs[icalc].ext_out)
      {
        ffilter = gtk_file_filter_new();
        tmp = g_strconcat("*.", input_filetypes[calc_defs[icalc].iout].extension, NULL);
        gtk_file_filter_set_name(ffilter, input_filetypes[calc_defs[icalc].iout].description);
        gtk_file_filter_add_pattern(ffilter, tmp);
        g_free(tmp);
        tmp = NULL;
        gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(file_chooser), ffilter);
      }
      gtk_widget_show(file_chooser);
    }

    if (calc_defs[icalc].extraargsym && calc_defs[icalc].extraargval)
    {
#ifdef GTK2
      hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
      hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
      gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);      
#endif
      gtk_box_pack_start(GTK_BOX(actionarea), hbox, FALSE, FALSE, 0);

      if (calc_defs[icalc].extraargdes)
      {
        label = gtk_label_new(calc_defs[icalc].extraargdes);
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
        gtk_widget_show(label);
      }

      entry = gtk_entry_new();
      gtk_entry_set_text(GTK_ENTRY(entry), calc_defs[icalc].extraargval);
      gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
      gtk_widget_show(entry);

      gtk_widget_show(hbox);
    }

    run = gtk_dialog_run(GTK_DIALOG(dialog));

    if (run==GTK_RESPONSE_ACCEPT)
    {
      if (calc_defs[icalc].extraargsym && calc_defs[icalc].extraargval)
      {
        gargv[iarg++] = calc_defs[icalc].extraargsym;

        tmp = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
#ifdef EBUG
        printf("getting text from entry = |%s|\n", tmp);
#endif
        gargv[iarg++] = tmp;
      }
      if (calc_defs[icalc].extrafile)
      {
        extra_filename=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
        if (extra_filename == NULL)
        {
          make_warning("No file is selected, calculation can not proceed!");
          gtk_widget_destroy(dialog);
          if (tmp) g_free(tmp);
          return;
        }
        gargv[iarg++]=extra_filename;
      }
    }
    else if (run==GTK_RESPONSE_REJECT)
    {
      gtk_widget_destroy(dialog);
      if (tmp) g_free(tmp);
      return;
    }

    gtk_widget_destroy(dialog);
  }

  /*1. save file*/
#ifdef EBUG
  printf("saving file: |%s|\n", calc_filename);
#endif
  save_gv_file(calc_filename);
#ifdef EBUG
  printf("file |%s| saved\n", calc_filename); fflush(stdout);
#endif

  /*2. execute script*/

  gargv[iarg++] = (gchar*) calc_filename;
  gargv[iarg++] = (gchar*) NULL;

#ifdef EBUG
  printf("libpath:\n");
  printf("|%s|\n", calc_defs[icalc].libpath);
  printf("arguments:\n");
  printf("|%s|\n", gargv[0]);
  printf("|%s|\n", gargv[1]);
  printf("|%s|\n", gargv[2]);
  printf("|%s|\n", gargv[3]);
  printf("|%s|\n", gargv[4]);
  printf("|%s|\n", gargv[5]); fflush(stdout);
#endif

  ret = g_spawn_async(calc_defs[icalc].libpath, gargv, NULL, G_SPAWN_DO_NOT_REAP_CHILD, NULL, NULL, &gpid, &gerror);
  if (!ret)
  {
    make_warning(gerror->message);
    g_free(calc_filename);
    return;
  }
  else calc_alive = 1;
#ifdef EBUG
  printf("childs pid = %d\n", gpid);
#endif
  if (tmp) g_free(tmp);

  g_free(extra_filename);
  /*4. check if calculation is running*/
  close_file();

  printf("OPENING FILE: \n");

  open_file(calc_filename, NULL);
  
  calc_timer = g_timeout_add(1000, (GSourceFunc) luscus_check_if_calc_runs, (gpointer) calc_filename); /*watch the file every second*/
  g_child_watch_add(gpid, (GChildWatchFunc) luscus_calculation_over, NULL);

  luscus_gtk_push_message_to_statusbar2("Running calculation...");
}

gboolean luscus_check_if_calc_runs(char *calc_filename)
{
#ifdef EBUG
  printf("DEBUG: checking calculation! opening file: |%s|\n", calc_filename); fflush(stdout);
#endif
  check_if_file_changed();
  if (calc_alive)
  {
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2("Calculation running...");
    return TRUE; /*continue the timer*/
  }
  else
  {
    g_free(calc_filename); /*POSIBLLE BUG -check if calc_filename is needed after this function!!!*/
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2("Calculation completed");
    return FALSE; /*break the timer*/
  }
}

void luscus_calculation_over(GPid gpid, gint status, gpointer data)
{
#ifdef EBUG
  printf("DEBUG: calculation DONE!\n");fflush(stdout);
#endif
  calc_alive = 0;
  g_spawn_close_pid(gpid);
}

gboolean luscus_gtk_key_press(GtkWidget *da, GdkEventKey *event, gpointer data)
{
  gint i;
  gint x, y;
  GdkModifierType mask;


  if (textbox_state == 1) return TRUE;
  if (textbox_state == 2)
  {
    switch(event->keyval)
    {
#ifdef GTK_OLD
      case GDK_BackSpace:
#else
      case GDK_KEY_BackSpace:
#endif
      textbox_delete_last_char(m);
      draw_pixdata_textbox(m->ntextboxes-1);
      break;

#ifdef GTK_OLD
      case GDK_Escape:
#else
      case GDK_KEY_Escape:
#endif
      delete_last_textbox(m);

#ifdef GTK_OLD
      case 10:
      case GDK_KP_Enter:
      case GDK_ISO_Enter:
      case GDK_3270_Enter:
      case GDK_RockerEnter:
      case GDK_Return:
#else
      case 10:
      case GDK_KEY_KP_Enter:
      case GDK_KEY_ISO_Enter:
      case GDK_KEY_3270_Enter:
      case GDK_KEY_RockerEnter:
      case GDK_KEY_Return:
#endif
      textbox_state = 0;
      luscus_set_cursor_default();
      luscus_gtk_update_3Dobject_info();
      append_backup();
      break;
#ifdef GTK_OLD    /*fix characters from numerical keyboard 2014.01.30*/
      case GDK_KP_0:
#else
      case GDK_KEY_KP_0:
#endif
      textbox_delete_insert_char(m, '0');
      break;
#ifdef GTK_OLD
      case GDK_KP_1:
#else
      case GDK_KEY_KP_1:
#endif
      textbox_delete_insert_char(m, '1');
      break;
#ifdef GTK_OLD
      case GDK_KP_2:
#else
      case GDK_KEY_KP_2:
#endif
      textbox_delete_insert_char(m, '2');
      break;
#ifdef GTK_OLD
      case GDK_KP_3:
#else
      case GDK_KEY_KP_3:
#endif
      textbox_delete_insert_char(m, '3');
      break;
#ifdef GTK_OLD
      case GDK_KP_4:
#else
      case GDK_KEY_KP_4:
#endif
      textbox_delete_insert_char(m, '4');
      break;
#ifdef GTK_OLD
      case GDK_KP_5:
#else
      case GDK_KEY_KP_5:
#endif
      textbox_delete_insert_char(m, '5');
      break;
#ifdef GTK_OLD
      case GDK_KP_6:
#else
      case GDK_KEY_KP_6:
#endif
      textbox_delete_insert_char(m, '6');
      break;
#ifdef GTK_OLD
      case GDK_KP_7:
#else
      case GDK_KEY_KP_7:
#endif
      textbox_delete_insert_char(m, '7');
      break;
#ifdef GTK_OLD
      case GDK_KP_8:
#else
      case GDK_KEY_KP_8:
#endif
      textbox_delete_insert_char(m, '8');
      break;
#ifdef GTK_OLD
      case GDK_KP_9:
#else
      case GDK_KEY_KP_9:
#endif
      textbox_delete_insert_char(m, '9');
      break;
#ifdef GTK_OLD
      case GDK_KP_Add:
#else
      case GDK_KEY_KP_Add:
#endif
      textbox_delete_insert_char(m, '+');
      break;
#ifdef GTK_OLD
      case GDK_KP_Subtract:
#else
      case GDK_KEY_KP_Subtract:
#endif
      textbox_delete_insert_char(m, '-');
      break;
#ifdef GTK_OLD
      case GDK_KP_Divide:
#else
      case GDK_KEY_KP_Divide:
#endif
      textbox_delete_insert_char(m, '/');
      break;
#ifdef GTK_OLD
      case GDK_KP_Multiply:
#else
      case GDK_KEY_KP_Multiply:
#endif
      textbox_delete_insert_char(m, '*');
      break;
#ifdef GTK_OLD
      case GDK_KP_Decimal:
#else
      case GDK_KEY_KP_Decimal:
#endif
      textbox_delete_insert_char(m, '.');
      break;
      default:
      if (event->keyval > 31 && event->keyval < 127)
        textbox_delete_insert_char(m, event->keyval);
    }
    draw_pixdata_textbox(m->ntextboxes-1);
    redraw();
    return TRUE;
  }

#ifdef GTK_OLD
  if (event->keyval == GDK_Control_L || event->keyval == GDK_Control_R)
#else
  if (event->keyval == GDK_KEY_Control_L || event->keyval == GDK_KEY_Control_R)
#endif
  {
    control_pressed = 1;
    return TRUE;
  }
  if (control_pressed)
  {
    if (event->keyval == 'o') callback_open_file(NULL, NULL);
    else if (event->keyval == 'q') Kill_Gui();
    else if (event->keyval == 'h') luscus_gtk_callback_help(NULL, NULL);
    return TRUE;
  }

  if (!accept_keys) return FALSE;
#ifdef GTK2
  gdk_window_get_pointer(event->window, &x, &y, &mask);
#endif
#ifdef GTK3
  mask=event->state;
/*  gdk_window_get_device_position(event->window, gtk_get_current_event_device(), &x, &y, &mask);*/
#endif

#ifdef GTK_OLD
  if (m->natom+m->nsphere+m->nvector+m->ntriangle+m->nsurf+m->ncells == 0 && (event->keyval != GDK_F1 && event->keyval != GDK_F9)) return TRUE;
#else
  if (m->natom+m->nsphere+m->nvector+m->ntriangle+m->nsurf+m->ncells == 0 && (event->keyval != GDK_KEY_F1 && event->keyval != GDK_KEY_F9))
    return TRUE;
#endif

  switch(event->keyval)
  {
#ifdef GTK_OLD
    case GDK_F1:
#else
    case GDK_KEY_F1:
#endif
    luscus_gtk_callback_help(NULL, NULL);
    break;
#ifdef GTK_OLD
    case GDK_F2:
#else
    case GDK_KEY_F2:
#endif
    callback_save_file(NULL, GINT_TO_POINTER(1));
    break;

#ifdef GTK_OLD
    case GDK_F3:
#else
    case GDK_KEY_F3:
#endif
    if (m->ngrids) luscus_gtk_show_multiorb_window();
    break;

#ifdef GTK_OLD
    case GDK_F5:
#else
    case GDK_KEY_F5:
#endif
      if (mask & GDK_SHIFT_MASK)
        do_screenshot(1);
      else
        do_screenshot(0);
/*    callback_screenshot(NULL, NULL);*/
    break;

#ifdef GTK_OLD
    case GDK_F6:
#else
    case GDK_KEY_F6:
#endif
      if (mask & GDK_SHIFT_MASK)
        remove_watched_coords();
      else
        add_watched_coord();
      break;

#ifdef GTK_OLD
    case GDK_F7:
#else
    case GDK_KEY_F7:
#endif
      mark_neighbor();
      rerender_3d();
      break;

#ifdef GTK_OLD
    case GDK_F8:
#else
    case GDK_KEY_F8:
#endif
      if(m->editable)
      {
        if (mask & GDK_SHIFT_MASK)
          luscus_gtk_set_vec_rot(2);
        else
          luscus_gtk_set_vec_rot(0);
        do_symmetry();
      }
    break;

#ifdef GTK_OLD
    case GDK_F9:
#else
    case GDK_KEY_F9:
#endif
    if (mask & GDK_SHIFT_MASK) callback_change_background(NULL, NULL); 
    else save_settings_data();
    break;

#ifdef GTK_OLD
    case GDK_F10:
#else
    case GDK_KEY_F10:
#endif
    callback_save_file(NULL, 0);
    exit(0);
    break;

#ifdef GTK_OLD
    case GDK_Right:
#else
    case GDK_KEY_Right:
#endif
/*      if(m->editable)
      {*/
        if (mask & GDK_SHIFT_MASK) Do_key_lr(1, 1);
        else Do_key_lr(1, 0);
/*      }*/
      break;

#ifdef GTK_OLD
    case GDK_Left:
#else
    case GDK_KEY_Left:
#endif
/*      if(m->editable)
      {*/
        if (mask & GDK_SHIFT_MASK) Do_key_lr(0, 1);
        else Do_key_lr(0, 0);
/*      }*/
      break;
#ifdef GTK_OLD
    case GDK_Up:
#else
    case GDK_KEY_Up:
#endif
/*      if(m->editable)
      {*/
        if (mask & GDK_SHIFT_MASK) Do_key_ud(1, 1);
        else Do_key_ud(1, 0);
/*      }*/
      break;
#ifdef GTK_OLD
    case GDK_Down:
#else
    case GDK_KEY_Down:
#endif
/*      if(m->editable)
      {*/
        if (mask & GDK_SHIFT_MASK) Do_key_ud(0, 1);
        else Do_key_ud(0, 0);
/*      }*/
    break;
#ifdef GTK_OLD
    case GDK_Page_Up:
#else
    case GDK_KEY_Page_Up:
#endif
      if(m->editable)
      {
        if (mask & GDK_SHIFT_MASK) do_key_page(1, 1);
        else do_key_page(1, 0);
      }
    break;
#ifdef GTK_OLD
    case GDK_Page_Down:
#else
    case GDK_KEY_Page_Down:
#endif
      if(m->editable)
      {
        if (mask & GDK_SHIFT_MASK) do_key_page(0, 1);
        else do_key_page(0, 0);
      }
    break;

#ifdef GTK_OLD
    case GDK_F12:
    case GDK_Insert:
#else
    case GDK_KEY_F12:
    case GDK_KEY_Insert:
#endif
    if (m->editable)
    {
      if (m->n_selected == 0) add_fragment(get_last_fragment());
      else if (m->n_selected == 1)
      {
        if (mask & GDK_SHIFT_MASK) add_atoms(1);
        else /*add_atoms(0); */ add_fragment(get_last_fragment());
      }
      else if (m->n_selected == 2)
      {
        if (mask & GDK_SHIFT_MASK) add_atoms(1);
        else
        {
          i = find_bond(m, m->selected[0], m->selected[1]);
          if (i == -1) change_bond_type(1);
          else change_bond_type((m->bond[i].bond_type+1)%7);
        }
      }
      rerender_3d();
    }
    break;

#ifdef GTK_OLD
    case GDK_Home:
#else
    case GDK_KEY_Home:
#endif
    if (m->editable) set_origin_molecule();
    break;

#ifdef GTK_OLD
    case GDK_End:
#else
    case GDK_KEY_End:
#endif
      if (m->editable)
      {
        add_atoms(1);        
        rerender_3d();
      }
    break;

#ifdef GTK_OLD
    case GDK_BackSpace:
#else
    case GDK_KEY_BackSpace:
#endif
      if (m->editable)
      {
        get_backup();
        rerender_3d();
      }
    break;

#ifdef GTK_OLD
    case GDK_Escape:
#else
    case GDK_KEY_Escape:
#endif
    remove_watched_coords();

      if (m->editable)
      {
        deallocate_vectors(m);
        deallocate_triangles(m);
        deallocate_spheres(m);
        deallocate_surfaces(m);
        deallocate_cells(m);
        luscus_gtk_update_upon_select_or_mark();
        luscus_gtk_update_3Dobject_info();
        rerender_3d();
      }
    break;

#ifdef GTK_OLD
    case GDK_Delete:
#else
    case GDK_KEY_Delete:
#endif

      if (m->editable)
      {
        if (m->n_selected == 1 || m->n_marked)
          delete_coord(0);
/*          delete_atom(m->selected[0]);*/
        else if (m->n_selected == 2)
        {
          i = find_bond(m, m->selected[0], m->selected[1]);
          printf("I = %d\n", i); fflush(stdout);
          if (i >= 0)
            if (m->bond[i].bond_type > 0)
              change_bond_type(m->bond[i].bond_type-1);
        } /*delete_grid should be included*/
        rerender_3d();
      }
    break;

    case 't':
    if (mask & GDK_MOD1_MASK)
      Input_Data.electrostatic_poten_color[0][3] = 
      Input_Data.electrostatic_poten_color[1][3] = 
      Input_Data.neg_pos_color[0][3] = Input_Data.neg_pos_color[1][3] =
      Input_Data.extracolor[3] = 1.0;
    else
    {
      Input_Data.electrostatic_poten_color[0][3] += 0.1;
      Input_Data.electrostatic_poten_color[1][3] += 0.1;
      Input_Data.neg_pos_color[0][3] += 0.1;
      Input_Data.neg_pos_color[1][3] += 0.1;
      Input_Data.extracolor[3] += 0.1;
    }

    if (Input_Data.extracolor[3] > 1.0)
    {
      Input_Data.electrostatic_poten_color[0][3] = 
      Input_Data.electrostatic_poten_color[1][3] = 
      Input_Data.neg_pos_color[0][3] = Input_Data.neg_pos_color[1][3] =
      Input_Data.extracolor[3] = 1.0;
    }
    rerender_3d();
    break;

    case 'T':
      Input_Data.electrostatic_poten_color[0][3] -= 0.1;
      Input_Data.electrostatic_poten_color[1][3] -= 0.1;
      Input_Data.neg_pos_color[0][3] -= 0.1;
      Input_Data.neg_pos_color[1][3] -= 0.1;
      Input_Data.extracolor[3] -= 0.1;

      if (Input_Data.extracolor[3] < 0.0)
      {
        Input_Data.electrostatic_poten_color[0][3] = 0.0;
        Input_Data.electrostatic_poten_color[1][3] = 0.0;
        Input_Data.neg_pos_color[0][3] = 0.0;
        Input_Data.neg_pos_color[1][3] = 0.0;
        Input_Data.extracolor[3] = 0.0;
      }
      rerender_3d();
      break;

    case 'W':
      Input_Data.atomshape++;
      if(Input_Data.atomshape >3) Input_Data.atomshape=0;
      break;

      case 'r':
      case 'g':
      case 'b':
      case 'U':
      case 'G':
      case 'B':

      callback_change_background(NULL, NULL);
      break;
    case 'c':
      if (m->ishow & HAS_MULLIKEN)
        Input_Data.label_atoms &~ (1 << 5);
      else if (m->ishow & HAS_LOPROP)
        Input_Data.label_atoms &~ (1 << 6);
      break;
    case 'l':
      luscus_gtk_move_light_state();
      break;
    case 'm':
      Input_Data.animate=!Input_Data.animate;
      luscus_gtk_start_animation(0);
      break;
    case 'w':
      if (control_pressed) callback_close_file(NULL, NULL);
      else if(m->n_selected==1 && m->editable) add_fragment(0);
      else Input_Data.bw=!Input_Data.bw;
      rerender_3d();
      break;
    case 'x':
      callback_maximize(NULL, NULL);
      break;
    case 'p':
      do_screenshot(3);
      break;
    case 'P':
      do_screenshot(2);
      break;
    case 'A':
      if (m->nvibr)
        Input_Data.frequency_amplitude += 0.1;
      else
      {
        Input_Data.ntexture++;
        if(Input_Data.ntexture > 5 ) Input_Data.ntexture=0;
      }
      break;
    case 'I':
      if (m->editable)
      {
        do_key_insert();
        rerender_3d();
      }
      break;
    case 'a':
      if (m->nvibr) Input_Data.frequency_amplitude -= 0.1;
      else if (m->ngrids) change_orbital_type(ORBITAL_TYPE_2);
      else
      {
        mark_H();
        rerender_3d();
      }
      break;
    case 'f':
      if (m->ngrids) change_orbital_type(ORBITAL_TYPE_F);
      break;
    case 'i':
      if (m->ngrids) change_orbital_type(ORBITAL_TYPE_I);
      break;
#ifdef GTK_OLD
    case GDK_KP_1:
#else
    case GDK_KEY_KP_1:
#endif
    case '1':
      if (m->ngrids) change_orbital_type(ORBITAL_TYPE_1);
      break;
#ifdef GTK_OLD
    case GDK_KP_2:
#else
    case GDK_KEY_KP_2:
#endif
    case '2':
      if (m->ngrids) change_orbital_type(ORBITAL_TYPE_2);
      break;
#ifdef GTK_OLD
    case GDK_KP_3:
#else
    case GDK_KEY_KP_3:
#endif
    case '3':
      if (m->n_selected == 3 && m->editable) 
      {
	add_triangle_selected();
        luscus_gtk_update_upon_select_or_mark();
        luscus_gtk_update_3Dobject_info();
        rerender_3d();
      }
      else if (m->ngrids) change_orbital_type(ORBITAL_TYPE_3);
      break;

#ifdef GTK_OLD
    case GDK_KP_4:
#else
    case GDK_KEY_KP_4:
#endif
    case '4':
      if (m->n_selected == 3 && m->editable)
      {
	add_surface_selected();
        luscus_gtk_update_upon_select_or_mark();
        luscus_gtk_update_3Dobject_info();        
        rerender_3d();
      }
      break;

#ifdef GTK_OLD
    case GDK_KP_6:
#else
    case GDK_KEY_KP_6:
#endif
    case '6':
      if (m->n_selected == 4 && m->editable)
      {
        add_cell_selected();
        luscus_gtk_update_upon_select_or_mark();
        luscus_gtk_update_3Dobject_info();        
        rerender_3d();
      }
      break;

    case 's':
      if (control_pressed) callback_save_file(NULL, (gpointer) 0);
      else if (m->ngrids) change_orbital_type(ORBITAL_TYPE_S);
      break;

    case 'd':
      if (m->ngrids) change_orbital_type(ORBITAL_TYPE_D);
      break;

#ifdef GTK_OLD
    case GDK_plus:
    case GDK_KP_Add:
#else
    case GDK_KEY_plus:
    case GDK_KEY_KP_Add:
#endif
      if (m->nvibr)       Input_Data.frequency_speed += 0.1;
      else if (m->ngrids)
      {
       	Input_Data.lev += 0.01;
        make_surfaces(); 
        rerender_3d();
        redraw();
      }
      else if (m->editable) luscus_gtk_move_coord_step(1);
      break;

#ifdef GTK_OLD
    case GDK_minus:
    case GDK_KP_Subtract:
#else
    case GDK_KEY_minus:
    case GDK_KEY_KP_Subtract:
#endif
      if (m->nvibr)       Input_Data.frequency_speed -= 0.1;
      else if (m->ngrids)
      {
       	Input_Data.lev -= 0.01;
        make_surfaces(); 
        rerender_3d();
        redraw();	
      }
      else if (m->editable) luscus_gtk_move_coord_step(-1);
      break;
    case 'q':
      if(get_number_of_slices() > 4) set_number_of_slices(get_number_of_slices() - 1);
      if(get_number_of_stacks() > 4) set_number_of_stacks(get_number_of_stacks() - 1);
      rerender_3d();
      break;
    case 'Q':
      if(get_number_of_slices() < 30) set_number_of_slices(get_number_of_slices() + 1);
      if(get_number_of_stacks() < 30) set_number_of_stacks(get_number_of_stacks() + 1);
      rerender_3d();
      break;
    case 'Z':
      change_scale(1);
      break;
    case 'z':
      change_scale(-1);
      break;
    case '*':
      reverse_marked();
      rerender_3d();
      break;
    case '#':
      renumber_marked();
      break;
    case ' ':
      if (m->n_selected) luscus_gtk_remove_selection(NULL, NULL);
      else if (m->n_marked) unmark_all();
      rerender_3d();
      break;
    case 'v':
      if(m->n_selected == 2)
      {
        if(get_vec_rot()==0)
        {
          luscus_gtk_set_vec_rot(2);
          luscus_gtk_pop_message_from_statusbar2();
          luscus_gtk_push_message_to_statusbar2("symmetry operation: rotation");
        }
        else
        {
          luscus_gtk_set_vec_rot(0);
          luscus_gtk_pop_message_from_statusbar2();
          luscus_gtk_push_message_to_statusbar2("symmetry operation: translation");
        }
      }
      else  ViewPoint(0);
      break;
     
    case 'V':
      ViewPoint(1);
      break;
    case 'k':
/*      doSlapafDump();*/ printf("doing slapaf dump\n");
      break;
    case 'K':
/*      doIndexDump();*/ printf("doing index dump\n");
      break;
    case '|':
      if (m->n_selected == 2 || m->n_selected == 3)
      {
        mark_one_side();
        rerender_3d();
/*        printf("mark one side!\n"); fflush(stdout);*/
      }
      break;
    default: return TRUE;
  }
  redraw();


  return TRUE;
}

gboolean luscus_gtk_key_release(GtkWidget *da, GdkEventKey *event, gpointer data)
{
  if (!accept_keys) return FALSE;
#ifdef GTK_OLD
  if(event->keyval == GDK_Control_L || event->keyval == GDK_Control_R)
#else
  if(event->keyval == GDK_KEY_Control_L || event->keyval == GDK_KEY_Control_R)
#endif
    control_pressed=0;
  return TRUE;
}
