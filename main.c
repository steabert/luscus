/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<gtk/gtk.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

#ifdef EBUG

char *rcdir;
int number_of_elements;
INPUT_DATA Input_Data;
MOL *mol;
MOL *m;
int n_geometries;
long *geo_file_position;
int igeo;
int accept_keys = 1; 

ELEM_DATA *e;
int n_input_types;
INPUT_FORMAT *input_filetypes;
int n_fragments;
FRAG_DATA *frag;
int n_watched = 0;
WATCH *watch = NULL;

void print_debug_data(void)
{
  int i;
  printf("Input_Data:\n");
  printf("Hide atoms %d \n",             Input_Data.hide_atoms);
  printf("Hide bonds %d \n",             Input_Data.hide_bonds);
  printf("Transparency %f \n",           Input_Data.extracolor[3]);
  printf("Background color Red %f \n",   Input_Data.background_color[0]);
  printf("Background color Green %f \n", Input_Data.background_color[1]);
  printf("Background color Blue %f \n",  Input_Data.background_color[2]);
 
  printf("Extra color Red %f \n",        Input_Data.extracolor[0]);
  printf("Extra color Green %f \n",      Input_Data.extracolor[1]);
  printf("Extra color Blue %f \n",       Input_Data.extracolor[2]);

  printf("Label color Red %f \n",        Input_Data.label_color[0]);
  printf("Label color Green %f \n",      Input_Data.label_color[1]);
  printf("Label color Blue %f \n",       Input_Data.label_color[2]);

  printf("Plus  color Red %f \n",        Input_Data.neg_pos_color[0][0]);
  printf("Plus  color Green %f \n",      Input_Data.neg_pos_color[0][1]);
  printf("Plus  color Blue %f \n",       Input_Data.neg_pos_color[0][2]);
  printf("Minus color Red %f \n",        Input_Data.neg_pos_color[1][0]);
  printf("Minus color Green %f \n",      Input_Data.neg_pos_color[1][1]);
  printf("Minus color Blue %f \n",       Input_Data.neg_pos_color[1][2]);

  printf("Electrostatic plus  color Red %f \n",   Input_Data.electrostatic_poten_color[0][0]);
  printf("Electrostatic plus  color Green %f \n", Input_Data.electrostatic_poten_color[0][1]);
  printf("Electrostatic plus  color Blue %f \n",  Input_Data.electrostatic_poten_color[0][2]);
  printf("Electrostatic minus color Red %f \n",   Input_Data.electrostatic_poten_color[1][0]);
  printf("Electrostatic minus color Green %f \n", Input_Data.electrostatic_poten_color[1][1]);
  printf("Electrostatic minus color Blue %f \n",  Input_Data.electrostatic_poten_color[1][2]);

  printf("Fatness atoms %f \n",          Input_Data.fatness_a);
  printf("Fatness bonds %f \n",          Input_Data.fatness_b);

  printf("Animate %d \n",                Input_Data.animate);
  printf("Initial Screen Size %d \n",    Input_Data.init_screen_size);
  printf("Level %f \n",                  Input_Data.lev);
  
  printf("FILETYPES:\n");

  for(i = 0; i < n_input_types; i++)
  {
    printf("ext = |%s|  desc = |%s|  forw = |%s|  back = |%s|  libpath = |%s|\n", input_filetypes[i].extension, input_filetypes[i].description, input_filetypes[i].forward, input_filetypes[i].backward, input_filetypes[i].libpath);
  }

}
#endif

int main(int argc, char **argv)
{
  init_variables();

  init_input_filetypes();

  load_atom_data();

  load_fragment_data();

  load_calc_def_data();

  load_settings_data(argc, argv);

#ifdef EBUG
  print_debug_data();
#endif

  Init_Gui(argc, argv);

  print_reference();

  open_cmd_line_file();

  start_check_for_changes();

  Start_Gui();

  remove_calc_def_data();

  return EXIT_SUCCESS;
}

