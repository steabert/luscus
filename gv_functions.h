/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/

/*function prototypes*/

/*gv_system.c*/
int dir_exists(char*);
void Print_help(void);
void print_plugin_warnings(void);
int check_file_exists(char*, char*);
void init_config_dir_path(void);
void init_variables(void);
void load_settings_data(int, char**);
void init_input_filetypes(void);
void load_default_values(void);
void load_values_from_config_file(void);
void load_values_from_command_line_args(int, char**);
void print_vresion(void);
void open_cmd_line_file(void);
void load_fragment_data(void);
void load_calc_def_data(void);
void save_settings_data(void);
void remove_calc_def_data(void);

 /*gv_geom.c*/
void get_center(double*, double*, double*);
void set_origin_molecule(void);
int find_first_bounded(int, int);
double get_selected_bond_length(void);
double get_selected_angle_value(void);
double get_selected_dihedral_value(void);
double get_bond_length(int, int);
double get_angle_value(int, int, int);
double get_dihedral_value(int, int, int, int);
void accumulate_motion(double);
void luscus_gtk_move_bond(double);
void luscus_gtk_move_angle(double);
void luscus_gtk_move_torsion(double);
void luscus_gtk_move_coord(double);
void luscus_gtk_move_coord_step(int);
int move_angle(double*, double*, double*, double, XYZ, XYZ, XYZ, double [3][3], double, int);
void det_matrix(double*, double [3][3]);
void inv_matrix (double*, double*, double*,
                 double*, double*, double*, 
                 double*, double*, double*,
                 double [3][3]);
void move_d_angle(double*, double*, double*, double*,
                  double,
                  double [3], double [3], double [3], double [3],
                  double, int, double, int);
void translate(XYZ);
void translate_to_origin(XYZ, XYZ, double [3][3]);

void rr_cross(double*, double*, double*,
              double, double, double,
              double, double, double,
              double, double, double);
double rr_vec(XYZ, XYZ, XYZ);
double rr_dist(XYZ, XYZ);
int find_bond(MOL*, int, int);
int find_one_bond_from_atom(int);
void CalcCyllinder(double in[3], double ut[3], double*, double*,  double*, double*);
void Calc2Cyllinder( double, double, double,
                     double, double, double,
                     double, double, double, 
                     double*, double*, double*);
double Calc_Diameter(void);
void pl3to4(double [3], double [3], double [3], double [3],
            double [3], double [3], double [3]);
void pl4to6(double [3], double [3], double [3], double [3],
            double [3], double [3], double [3], double [3]);
void np3(double [3], double [3], double [3], double [3]);

void swap_atoms(int, int);
void add_atoms(int);
void build_rotation_matrix(XYZ, XYZ, double [3][3]);
void change_bond_type(int);
int get_bond_type(int, int);
void delete_bonds_with_atom(int);
void delete_coord(int);
void delete_atom(int);
void delete_bond(int);
int find_next_dummy(void);
void delete_dummy_atoms(void);
void do_inversion(void);

void transform_n_fold(int, int, XYZ, XYZ, XYZ);
void copy_atom_data(int, int);
void add_watched_coord(void);
void remove_watched_coord(int);
void remove_watched_coords(void);
void remove_all_watched_coords(void);
void add_triangle_selected(void);
void add_surface_selected(void);
void add_cell_selected(void);
void rebond(void);

/*gtk_gui.c*/
void get_sizes_data(double*, double*, double*, double*, double*, double*);
unsigned int get_list_3d(void);
void rerender_3d(void);
void set_go_selected(int, int);
void unset_go_selected(void);
int get_number_of_slices(void);
int get_number_of_stacks(void);
void set_number_of_slices(int);
void set_number_of_stacks(int);

void SetAtomColor(int);
void set_sizes(void);
void set_scale(void);
void change_scale(int);
void Do_center(void);
void Do_key_ud(int, int);
void Do_key_lr(int, int);


double grayscale(double, double, double);
void unselect_all(void);
void unmark_all(void);
int is_atom_marked(int);
void mark_atom(int);
void mark_H(void);
void mark_element_as_selected(void);
void mark_neighbor(void);
void mark_nei(int);
void mark_one_side(void);
void unmark_atom(int);
void renumber_marked(void);
void reverse_marked(void);
void change_element_by_one(int);
/*watched data*/
/*void draw_textbox_rectangle(void);*/
void draw_vibrations(void);
void draw_atoms_xyz(XYZ*);
void draw_bonds_xyz(XYZ*);
void draw_atoms(void);
void draw_bonds(void);
///*void compute_screen_coords(void);*/
void draw_labels(int, int);
void draw_axes(void);
void draw_dipole(void);
void draw_vectors(void);
void draw_triangles(void);
void draw_spheres(void);
void draw_surfaces(void);
void draw_cells(void);
void draw_grid(void);
void draw_selected_go(void);
void luscus_draw_3d(void);
void luscus_draw_2d(void);
void draw_rectangle(int, int);
void draw_labels(int, int);
void draw_textboxes(int, int);


///*void draw_textboxes(void);*/
//void draw_epot(void);
//void luscus_gtk_draw_text(char*, double, double);
void Init_Gui(int, char**);
void Start_Gui(void);
void Kill_Gui(void);
void make_warning(char*);
void print_reference(void);

void do_key_insert(void);

void redraw(void);
void print_bond_eps(int, int, double[3]);
void luscus_gtk_push_message_to_statusbar1(char*);
void luscus_gtk_push_message_to_statusbar2(char*);
void luscus_gtk_pop_message_from_statusbar1(void);
void luscus_gtk_pop_message_from_statusbar2(void);

void add_surface(MOL*, XYZ, XYZ, XYZ);
void add_vector(MOL*, XYZ, XYZ);
void add_sphere(MOL*, XYZ, double);
void add_triangle(MOL*, XYZ, XYZ, XYZ);
void add_cell(MOL*, XYZ, XYZ, XYZ, XYZ);
void textbox_delete_last_char(MOL*);
void textbox_delete_insert_char(MOL*, char);
void luscus_set_cursor_text(void);
void luscus_set_cursor_default(void);
void luscus_set_textbox_number(int);
void add_fragment(int);
//int find_connected(int);
int get_vec_rot(void);
void luscus_gtk_set_vec_rot(int);
void do_rotation(void);
void do_translation(void);
void do_expand_or_retract(void);
void transform_mirror(XYZ, XYZ, XYZ, XYZ, XYZ);
void do_mirroring(void);
void do_symmetry(void);
void luscus_gtk_move_light_state(void);
void ViewPoint(int);
void animate(int);
void luscus_gtk_freq_timer(void);
void start_check_for_changes(void);
void luscus_gtk_resort_surfaces(void);
void deallocate_msrf(void);
void allocate_msrf(int);
void make_surfaces(void);
void do_remove_grid(void);
//void do_change_iso(void);
//int tgaSave(char*, int, int);
void initPS(int, char*);
void initPovray(char*);
void do_key_page(int, int);
int get_last_fragment(void);

/*gv_fragment.c*/
void initialize_fragment_data(MOL*);
void free_fragment_data(MOL);
void load_fragment_from_file(char*, MOL*);
MOL load_fragment(int);
void fragment_set_element_data(MOL*, int);
void read_fragment_bond_block(MOL*, FILE*);
void fragment_add_bond(MOL*, int, int, int);

/*gv_notebook.c*/
void luscus_gtk_show_or_hide_widgets(void);
void luscus_gtk_update_upon_select_or_mark(void);
void change_watched_data(void);
void load_grid_data(void);
void delete_all_orbitals_from_list(void);
void change_orbital_type(int);
void luscus_show_substaces(void);
void luscus_gtk_select_orbital_up(void);
void luscus_gtk_select_orbital_down(void);
int get_orbital(int, int);
void luscus_gtk_update_3Dobject_info(void);
void luscus_gtk_update_geo_info(void);
void geo_play_button_play(int);

/*gv_menubar.c*/
void luscus_gtk_menubar_show_or_hide_widgets(void);

/*mystring.c*/
char *my_read_str_value(char*);
int my_read_int_value(char*);
double my_read_dble_value(char*);
char *strdiet(char*);
char *get_one_word(char*);
char *get_ptr_ith_word(char*, int);
int myisnan(double);
int double_ne(double, double);
int line_is_empty(char*);

/*gv_gtk_handle.c*/
void luscus_gtk_draw_sphere(char, double, int, int);
void luscus_gtk_draw_torus(char, double, double, int, int);
void luscus_gtk_draw_teapot(char, double);
void luscus_gtk_draw_cube(char, double);


/*read_file.c*/
void deallocate_mol_data(void);
void deallocate_data(MOL*);
void deallocate_atoms(MOL*);
void deallocate_bonds(MOL*);
void deallocate_dipole(MOL*);
void deallocate_geoms(MOL*);
void deallocate_grids(MOL*);
void deallocate_vibrations(MOL*);
void deallocate_vectors(MOL*);
void deallocate_triangles(MOL*);
void deallocate_spheres(MOL*);
void deallocate_surfaces(MOL*);
void deallocate_cells(MOL*);
void deallocate_textboxes(MOL*);

MOL *new_mol(int);
void allocate_atoms(MOL*, int);
void allocate_bonds(MOL*, int);
/*void allocate_geoms(MOL*, int);*/
void allocate_grids(MOL*, int);
void allocate_vibs(MOL*, int);
void allocate_vectors(MOL*, int);
void allocate_triangles(MOL*, int);
void allocate_spheres(MOL*, int);
void allocate_surfaces(MOL*, int);
void allocate_cells(MOL*, int);
void allocate_textbox(MOL*);
void delete_last_textbox(MOL*);
void delete_ith_textbox(MOL*, int);

char *read_line(FILE*);
char get_file_exist(char*);
void open_file(char*, char*);
FILE *open_current_luscus_file(void);
void open_gv_file(char*);
void parse_gv_file(FILE *in);
int read_geo_positions(FILE *in);
void read_all_sections(FILE *in);
void read_all_grids_from_all_sections(int*);
void read_section(FILE*, int, MOL*);
void read_coord_block(FILE*, MOL*);
void read_additional_atom_block(FILE*, MOL*);
void read_energy_block(FILE*, MOL*);
void read_rms_grad_block(FILE*, MOL*);
void read_max_grad_block(FILE*, MOL*);
void read_bond_block(FILE*, MOL*);
void determine_bonding(MOL*);
int get_number_of_bonds_on_atom(int, MOL*);
/*int get_number_of_bond_types_on_atom(MOL*, int, int);*/
void adapt_bonding(MOL*);
void read_dipole_block(FILE*, MOL*);
void read_vector_block(FILE*, MOL*);
void vector_set_default_values(MOL*, int);
void read_triangle_block(FILE*, MOL*);
void triangle_set_defeault_values(MOL*, int);
void read_sphere_block(FILE*, MOL*);
void sphere_set_defeault_values(MOL*, int);
void read_surface_block(FILE*, MOL*);
void surface_set_default_values(MOL *, int);
void read_cell_block(FILE*, MOL*);
void cell_set_default_values(MOL*, int);

void read_textbox_block(FILE*, MOL*);
void textbox_set_dfault_values(MOL*, int);

void read_grid_block(FILE*, MOL*);
void read_epoten_grid(FILE*, MOL*, int);
void read_vibration_block(FILE*, MOL*);
void read_editable_block(FILE*, MOL*);
void read_sleep_block(FILE*);
void read_write_block(FILE*);

void close_file(void);
void check_if_file_changed(void);
char *get_input_filename(void);
char *get_current_directory(void);
char *nextfname(void);
char *make_luscus_file_name(char*);
void get_new_section(void);

void set_element_data(MOL*, int);
double atom_distance(MOL*, int, int);
void add_bond(MOL*, int, int, int);

void do_load_grid(void);
void read_grid_from_file(void);
int read_bin_block(void*, int, FILE*);
int abfgets(int, char*, int, FILE*);
int read_bin_stream(char*, int, FILE*);
void current_grid_limits(double*, double*, double*);
int current_grid_data(double**);
int get_fmarker(void);
void remove_backup(void);

/*pixdata.c*/

void draw_all_pixdata(void);
void draw_pixdata_indices(void);
void draw_pixdata_symbols(void);
void draw_pixdata_atomnums(void);
void draw_pixdata_atomnames(void);
void draw_pixdata_mulliken(void);
void draw_pixdata_loprop(void);
void draw_pixdata_textbox(int);
void draw_textbox_rectangle(int, int);

/*gv_gtk_atom_list*/
void init_atom_list(void);
void deallocate_atom_list(void);
void insert_atom_into_list(int);
void insert_all_atoms_into_list(void);
void remove_atom_from_list(int);
void change_atom_parameters_in_list(int);

/*screenshot.c*/

int tgaSave(char*, int, int);

/*write_file.c*/

void luscus_gtk_do_save_file(char*, int);
int get_filetype(char*);
void save_gv_file(char*);
void convert_to_target_file(int, char*);
void write_all_data_to_output(FILE*, int);
void write_mol_section(FILE*, MOL*, int);
void write_coord_block(FILE*, MOL*);
void write_additional_atom_info_block(FILE*, MOL*);
void write_bond_block(FILE*, MOL*);
void write_dipole_block(FILE*, MOL*);
void write_vector_block(FILE*, MOL*);
void write_triangle_block(FILE*, MOL*);
void write_sphere_block(FILE*, MOL*);
void write_surface_block(FILE*, MOL*);
void write_cell_block(FILE*, MOL*);
void write_textbox_block(FILE*, MOL*);
void write_grid_block(FILE*, MOL*);
void write_inporb_block(FILE*, FILE*, MOL*);
void write_vibration_block(FILE*, MOL*);
void write_bin_block(FILE*, char*);

/*backup_file.c*/
void close_backup(void);
void open_backup(int);
char *get_backup_filename(void);
void append_backup(void);
void write_backup(void);
char *get_extension_pointer(char*);
char *get_filename_pointer_from_path(char*);
char *get_path_without_file(char*);
char *get_filename_without_extension(char*);
void get_backup(void);
void backup_search_last_backup_point(void);

/*gv_atom.c*/

void load_default_element_data(void);
void load_atom_data(void);
void Write_Atom_data(char*);
void save_atom_data(void);

/*findsym.c*/
int findsym(int, int, char [][4], XYZ*, char*, int,
            char*, int, double*, int, char*, int*, int*);

/*gvgrp.c*/
int gvgrp(char*, SymmetryElement X[]);

/*gveps.c*/
/*#include "gveps.h"*/

/*gv_gtk_make_graph.c*/
void set_current_graph_data(void);

