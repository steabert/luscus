/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*global variables*/
GtkWidget *window;	/*main window*/
GtkWidget *drawingarea; /*main drawingarea*/
GdkPixbuf *pix;		/*icon*/
GtkListStore *atom_list;

enum
{
  COLUMN_NUMBER,
  COLUMN_ATNAME,
  COLUMN_X,
  COLUMN_Y,
  COLUMN_Z,
  COLUMN_NUMER,
  COLUMN_NAME,
  N_COLUMNS
};

#define NONE 0
#define SPHERE 1
#define VECTOR 2
#define TRIANGLE 3
#define SURFACE 4
#define CELL 5
#define TEXTBOX 6

typedef struct gtk_3d_description
{
  GtkWidget *button;
  char i3dtype;
  int inum;
} GTK_3D_DESC;

/*gtk_gui.c*/

void draw_molecule_for_multiview(void);

void set_scale(void);
void gv_gtk_get_screen_size(gint*, gint*);

void callback_maximize(GtkWidget*, gpointer);
gboolean luscus_freq_timer(gpointer);
gboolean check_for_changes(gpointer);
double *get_main_mvm(void);
void redraw(void);

/*gv_notebook.c*/
GtkWidget *make_notebook(void);
void luscus_gtk_diseable_editing(void);
void luscus_gtk_enable_editing(void);
GtkWidget *luscus_gtk_make_treeview(void);
GtkWidget *luscus_make_geometry_objects_treeview(void);
gint luscus_gtk_orbital_compare_func(GtkTreeModel*, GtkTreeIter*, GtkTreeIter*, gpointer);
gboolean luscus_gtk_select_orbital_from_list(GtkTreeSelection*, GtkTreeModel*, GtkTreePath*, gboolean, gpointer);
void luscus_gtk_update_geo_info(void);
void deallocate_vib_buttons(void);

void delete_orbital_from_list(GtkWidget*, gpointer);
void undelete_all_orbitals_from_list(GtkWidget*, gpointer);
void orbital_combo_edited(GtkCellRendererText*, gchar*, gchar*, GtkListStore*);
gint luscus_gtk_change_orbital_type(gint, gint);
void orbital_filter_callback(GtkWidget*, gpointer);
void luscus_gtk_get_filter_str(gchar*);
void luscus_gtk_search_orbital_in_list(gint, gint);
void show_custom_grid(GtkWidget*, gpointer);

/*gv_menubar.c*/
GtkWidget *make_menubar(GtkWidget*);

/*gv_gtk_util.c*/

void luscus_check_writable_pic_files(GdkPixbufFormat*, GSList**);

void luscus_toggle_movement_mode(GtkWidget*, gpointer);
void luscus_gtk_move_up(GtkWidget*, gpointer);
void luscus_gtk_move_down(GtkWidget*, gpointer);
void luscus_gtk_move_forw(GtkWidget*, gpointer);
void luscus_gtk_move_back(GtkWidget*, gpointer);
void luscus_gtk_center(GtkWidget*, gpointer);
void luscus_gtk_zoom(GtkWidget*, gint);
void callback_geo_forw(GtkWidget*, gpointer);
void callback_geo_back(GtkWidget*, gpointer);
void callback_geo_first(GtkWidget*, gpointer);
void callback_geo_last(GtkWidget*, gpointer);
#ifdef GTK3
void callback_adj_geometry(GtkAdjustment*, gpointer);
#endif
#ifdef GTK2
void callback_adj_geometry(GtkObject*, gpointer);
#endif
void callback_geo_play(GtkWidget*, gpointer);
void luscus_gtk_change_transparency_level(GtkWidget*, gpointer);
void luscus_gtk_change_isosurface_level(GtkWidget*, gpointer);
void luscus_gtk_show_all_orbitals(GtkWidget*, gpointer);
void luscus_gtk_show_electrostatic_potential(GtkWidget*, gpointer);

void luscus_gtk_change_vibration_speed(GtkWidget*, gpointer);
void luscus_gtk_change_vibration_amplitude(GtkWidget*, gpointer);

void luscus_gtk_add_fragment(GtkWidget*, gpointer);
void luscus_gtk_remove_fragment(GtkWidget*, gpointer);
void luscus_gtk_change_bond_callback(GtkWidget*, gpointer);
void luscus_gtk_change_angle_callback(GtkWidget*, gpointer);
void callback_change_atom(GtkWidget*, gpointer);
GtkWidget *luscus_gtk_make_periodic_system(void);
gboolean luscus_gtk_select_element(GtkWidget*, gint);
void custom_atom_color_changed_callback(GtkWidget*, gpointer);
void custom_atom_size_changed_callback(GtkWidget*, gpointer);
void custom_bond_size_changed_callback(GtkWidget*, gpointer);
void custom_valency_changed_callback(GtkWidget*, gpointer);
gboolean luscus_gtk_disable_keys(GtkWidget*, GdkEvent*, gpointer);
gboolean luscus_gtk_enable_keys(GtkWidget*, GdkEvent*, gpointer);
void luscus_gtk_change_bond_value(GtkWidget*, gpointer);
void luscus_gtk_change_angle_value(GtkWidget*, gpointer);
void luscus_gtk_change_tors_value(GtkWidget*, gpointer);
/*void luscus_gtk_change_bond_angle_tors_value(GtkWidget*, gpointer);*/
void luscus_gtk_change_xyz_coordinate(GtkWidget*, gpointer);
void luscus_gtk_add_dummy_atoms(GtkWidget*, gpointer);
void luscus_gtk_remove_dummy_atoms(GtkWidget*, gpointer);
void luscus_gtk_remove_selection(GtkWidget*, gpointer);
void luscus_gtk_unmark(GtkWidget*, gpointer);
void luscus_gtk_sort_mark(GtkWidget*, gpointer);
void luscus_gtk_mark_H_atoms(GtkWidget*, gpointer);
void luscus_gtk_mark_reverse(GtkWidget*, gpointer);
void luscus_gtk_mark_element(GtkWidget*, gpointer);
void luscus_gtk_mark_neighbor(GtkWidget*, gpointer);
void luscus_gtk_watch_value(GtkWidget*, gpointer);
void luscus_gtk_unwatch_values(GtkWidget*, gpointer);
void luscus_gtk_select_fragment(GtkWidget*, gpointer);
void luscus_gtk_alpply_symmetry(GtkWidget*, gpointer);
void luscus_gtk_apply_translation(GtkWidget*, gpointer);
void change_translation_magnitude(GtkWidget*, gpointer);
void luscus_gtk_alpply_rot_symmetry(GtkWidget*, gpointer);
void change_rotation_axis_order(GtkWidget*, gpointer);
void luscus_gtk_force_symmetry(GtkWidget*, gpointer);
void luscus_add_arrow(GtkWidget*, gpointer);
void luscus_add_sphere(GtkWidget*, gpointer);
void luscus_add_plain(GtkWidget*, gpointer);
void luscus_add_triangle(GtkWidget*, gpointer);
void luscus_add_cell(GtkWidget*, gpointer);
void luscus_insert_textbox(GtkWidget*, gpointer);
void luscus_clear_drawings(GtkWidget*, gpointer);
void luscus_gtk_modify_3Dobject(GtkWidget*, GTK_3D_DESC*);
void luscus_move_textbox_by_click(GtkWidget*, GTK_3D_DESC*);

void callback_open_file(GtkWidget*, gpointer);
void callback_save_file(GtkWidget*, gpointer);
void set_save_file_filename(GtkWidget*, gint);
void callback_close_file(GtkWidget*, gpointer);

void callback_do_undo(GtkWidget*, gpointer);
void callback_make_editable(GtkWidget*, gpointer);
void callback_change_background(GtkWidget*, gpointer);
void color_t_2_GdkColor(GdkColor *, color_t);
void GdkColor_2_color_t(GdkColor *, color_t);
void callback_change_label(GtkWidget*, gpointer);
void callback_change_pos_orbital(GtkWidget*, gpointer);
void callback_change_neg_orbital(GtkWidget*, gpointer);
void callback_change_pos_epot(GtkWidget*, gpointer);
void callback_change_neg_epot(GtkWidget*, gpointer);
void callback_change_plain(GtkWidget*, gpointer);
void callback_change_font(GtkWidget*, gpointer);
void callback_adjust_atom_properties(GtkWidget*, gpointer);
GtkWidget* luscus_gtk_make_periodic_system1(ELEM_DATA*);
void callback_change_automatic_bonding(GtkWidget*, gpointer);
void luscus_gtk_atom_prop_select_element(GtkWidget*, gpointer);
void atom_prop_atom_size_spin_callback(GtkWidget*, gpointer);
void atom_prop_bond_size_spin_callback(GtkWidget*, gpointer);
void atom_prop_atom_valency_spin_callback(GtkWidget*, gpointer);
void callback_save_settings(GtkWidget*, gpointer);

void callback_move_light(GtkWidget*, gpointer);
void callback_create_viewpoint(GtkWidget*, gpointer);
void callback_restore_viewpoint(GtkWidget*, gpointer);
void callback_grayscale(GtkWidget*, gpointer);
void callback_animate(GtkWidget*, gpointer);
void luscus_gtk_start_animation(int);
gboolean luscus_start_animation1(gpointer);
void callback_teatime(GtkWidget*, gpointer);
void luscus_gtk_show_atoms(GtkWidget*, gpointer);
void luscus_gtk_show_bonds(GtkWidget*, gpointer);
void luscus_gtk_show_axes(GtkWidget*, gpointer);
void luscus_gtk_callback_show_dot_line(GtkWidget*, gpointer);

void luscus_gtk_show_labels(GtkWidget*, gpointer);
void luscus_gtk_show_indices(GtkWidget*, gpointer);
void luscus_gtk_show_numeration(GtkWidget*, gpointer);
void luscus_gtk_show_symbols(GtkWidget*, gpointer);
void luscus_gtk_show_names(GtkWidget*, gpointer);
void luscus_gtk_show_mulliken(GtkWidget*, gpointer);
void luscus_gtk_show_loprop(GtkWidget*, gpointer);
void do_screenshot(int);
void callback_screenshot(GtkWidget*, gpointer);
void luscus_set_filename(GtkWidget*, GtkWidget*);
void luscus_gtk_goto_freq(GtkWidget*, gpointer);

void luscus_start_calculation(GtkWidget*, gint);
void luscus_start_calculation2(GtkWidget*, gint);
gboolean luscus_check_if_calc_runs(char*);
void luscus_calculation_over(GPid, gint, gpointer);

gboolean luscus_gtk_key_press(GtkWidget*, GdkEventKey*, gpointer);
gboolean luscus_gtk_key_release(GtkWidget*, GdkEventKey*, gpointer);

/*gv_gtk_xyz_editor*/
void callback_xyz_edit(GtkWidget*, gpointer);
void luscus_change_x_coordinate(GtkCellRendererText*, const gchar*, const gchar*, gpointer);
void luscus_change_y_coordinate(GtkCellRendererText*, const gchar*, const gchar*, gpointer);
void luscus_change_z_coordinate(GtkCellRendererText*, const gchar*, const gchar*, gpointer);
void luscus_change_atom_numbering(GtkCellRendererText*, const gchar*, const gchar*, gpointer);
void luscus_change_atom_name(GtkCellRendererText*, const gchar*, const gchar*, gpointer);
void luscus_select_atom(GtkTreeView*, GtkTreePath*, GtkTreeViewColumn*, gpointer);
void luscus_xyz_editor_select_row(int);

/*gv_gtk_multiview.c*/
void kill_multi_view_win(GtkWidget*, gpointer);
/*void luscus_gtk_get_cur_color(int, GLfloat*, GLfloat*, GLfloat*);*/
gboolean luscus_gtk_multiview_redraw_list(int, int);
void redraw_list(int);

gboolean luscus_multi_view_key_press(GtkWidget*, GdkEventKey*, gpointer);
/*gboolean luscus_multi_view_win_resize(GtkWidget*, GdkEvent*, gpointer);*/
void arrange_frames_in_layout(void);
void magnify_frames(GtkWidget, gpointer);
void extenuate_frames(GtkWidget, gpointer);
void luscus_gtk_wsub_create(gint);
void luscus_gtk_show_multiorb_window(void);


/*gv_gtk_help.c*/ /*might be more than one function*/
void luscus_gtk_callback_help(GtkWidget*, gpointer);

/*gv_about*/
void luscus_gtk_callback_about(GtkWidget*, gpointer);

/*gv_gtk_make_graph.c*/ /*Only functions that are accessed from other modules are listed!*/
void luscus_gtk_make_graph(GtkWidget*, char*);
void luscus_gtk_show_graphical_energies(GtkWidget*, gpointer);
void luscus_gtk_show_vib_spectrum(GtkWidget*, gpointer);

/*screenshot.c*/
#ifdef GTK_GLEXT
void luscus_save_builtin_graphic_format(gchar*, gchar*, gint, gint);
#else
void luscus_save_builtin_graphic_format(gchar*, gchar*, GdkPixbuf*, gint, gint);
#endif
