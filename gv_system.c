/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#ifdef WINDOWS
#include<windows.h>
#else
#include<unistd.h>
#endif
#include<dirent.h>
#include<sys/stat.h>
#include<glib.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

char *rcdir;
INPUT_DATA Input_Data;
int n_input_types;
INPUT_FORMAT *input_filetypes;
MOL *m;
MOL *mol;
char *init_filename;
int n_fragments;
FRAG_DATA *frag;
int n_calc_defs;
DEF_CALC_T *calc_defs;
int accept_keys;
int move_camera;
int automatic_rebonding;
int error_config_does_not_exist;

int dir_exists(char* dirname)
{
  struct stat file_info;
  int rc = stat(dirname, &file_info);

  if (!S_ISDIR(file_info.st_mode) || rc!=0) return 0;
  else return 1;
}

void Print_help(void)
{
  printf("%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
  " Visualization of geometry, molecular orbitals and densities\n",
  "Usage: ", PROGRAM_NAME, " [-flag value] [filename]\n",
  "Available flags: \n",
  " -h                     : this help \n",
  " -?                     : this help \n",
  " -m                     : start in animation mode \n",
  " -s integer             : set initial screen size\n",
  " -l real                : isovalue level \n",
  " -v, --version          : version \n",
  " -t real (0.0-1.0)      : transparency level \n",
  " -b white|grey|black    : background \n",
  " --single               : determine single bonds only \n");
}

void print_plugin_warnings(void)
{
  if (error_config_does_not_exist == 1)
  {
#ifdef WINDOWS
    make_warning("Luscus did not find plugin directory. Please locate luscusrc directory in the luscus build directory and copy its contents to the one of the places:\n\n- path defined with the LUSCUS_DIR environment variable\n- \%LOCALAPPDATA\%\\.luscus\n- \%ALLUSERSPROFILE\%\\luscus\n - directory containing luscus.exe");
#else
    make_warning("Luscus can not find plugin directory. Please locate luscusrc directory in the luscus build directory and copy its contents to the one of the places:\n\n- path defined with the LUSCUS_DIR environment variable\n- $HOME/.luscus;\n- /etc/luscus\n- directory containing luscus.exe");
#endif
  }
/*  else if (error_config_does_not_exist == 2)
  {
    make_warning("luscus can not locate atoms.rc plugin file.");
  }*/
  else if (error_config_does_not_exist == 4)
  {
#ifdef WINDOWS
    make_warning("luscus can not locate calculation.rc plugin file. If you want to run calculations from luscus, please make sure that the plugin.rc file is placed in the luscus configuration directory. It should be in one of the following places:\n\n- path defined with the \%LUSCUS_DIR\% environment variable\n- \%LOCALAPPDATA\%\\.luscus\n- \%ALLUSERSPROFILE\%\\luscus\n - directory containing luscus.exe");
#else
    make_warning("luscus can not locate calculation.rc plugin file. If you want to run calculations from luscus, please make sure that the plugin.rc file is placed in the luscus configuration directory. It should be in one of the following places:\n\n- path defined with the LUSCUS_DIR environment variable\n- $HOME/.luscus\n- /etc/luscus\n- directory containing luscus.exe");
#endif
  }
  else if (error_config_does_not_exist == 8)
  {
#ifdef WINDOWS
    make_warning("luscus can not locate fragments.rc plugin file. If you want to use molecular fragments, please locate the fragments.rc file and all files listed in it and place them in the luscus configuration directory. It should be in one of the following places:\n\n- path defined with the \%LUSCUS_DIR\% environment variable\n- \%LOCALAPPDATA\%\\.luscus\n- \%ALLUSERSPROFILE\%\\luscus\n - directory containing luscus.exe");
#else
    make_warning("luscus can not locate fragments.rc plugin file. If you want to use molecular fragments, please locate fragments.rc file and all files listed in it and place them in the luscus configuration directory. It should be in one of the following places:\n\n- path defined with the LUSCUS_DIR environment variable\n- $HOME/.luscus\n- /etc/luscus\n- directory containing luscus.exe");
#endif
  }
  else if (error_config_does_not_exist == 16)
  {
#ifdef WINDOWS
    make_warning("luscus can not locate plugin.rc plugin file. This prevents luscus to open different file formats other than the lus format. Please make sure that the plugin.rc file is placed in the luscus configuration directory. It should be in one of the following places:\n\n- path defined with the \%LUSCUS_DIR\% environment variable\n- \%LOCALAPPDATA\%\\.luscus\n- \%ALLUSERSPROFILE\%\\luscus\n - directory containing luscus.exe");
#else
    make_warning("luscus can not locate plugin.rc plugin file. This prevents luscus to open different file formats other than the lus format. Please make sure that the plugin.rc file is placed in the luscus configuration directory. It should be in one of the following places:\n\n- path defined with the LUSCUS_DIR environment variable\n- $HOME/.luscus\n- /etc/luscus\n- directory containing luscus.exe");
#endif
  }
  else if (error_config_does_not_exist > 2048)
  {
    char tmpstr[1024];
    snprintf(tmpstr, 1024, "Plugin file %s can not be found. Please ensure that the file plugin.rc contains correct path to it.", input_filetypes[n_input_types-1].backward);
    make_warning(tmpstr);
  }
  else if (error_config_does_not_exist > 1024)
  {
    char tmpstr[1024];
    snprintf(tmpstr, 1024, "Plugin file %s can not be found. Please ensure that the file plugin.rc contains correct path to it.", input_filetypes[n_input_types-1].forward);
    make_warning(tmpstr);
  }
}

int check_file_exists(char *path, char *filename)
{
  char *fullpath = malloc(sizeof(char) * (strlen(path) + strlen(filename) + 2));
  int iacc;
  fullpath[0] = 0;
  strcat(fullpath, path);
#ifdef WINDOWS
  strcat(fullpath, "\\");
#else
  strcat(fullpath, "/");
#endif
  strcat(fullpath, filename);

  iacc = access(fullpath, F_OK | R_OK);

#ifdef EBUG
  printf("checking the existence of the file: |%s|: %d (0 means exists and readable)\n", fullpath, iacc);
#endif

  free(fullpath);
  return iacc;
}

void init_config_dir_path(void)
{
  char *ptr;

  error_config_does_not_exist = 0;
#ifdef WINDOWS

/*  1. check $LUSCUS_DIR environment  */
  ptr = getenv("LUSCUS_DIR");
  if (ptr)
  {
    rcdir = strdup(ptr);
    if (dir_exists(rcdir) && !check_file_exists(rcdir, RC_GV)) return;
  }

/*  2. check LOCALAPPDATA directory */  
  ptr = getenv("LOCALAPPDATA");
  if (ptr)
  {
    rcdir=malloc(sizeof(char) * (strlen(ptr) + strlen(RC_DIR) + 2));
    rcdir[0] = 0;
    strcat(rcdir, ptr);
    strcat(rcdir, "\\");
    strcat(rcdir, RC_DIR);

    if (dir_exists(rcdir) && !check_file_exists(rcdir, RC_GV)) return;
    else free(rcdir);
  }

/*  3. check ALLUSERSPROFILE directory */
  ptr = getenv("ALLUSERSPROFILE");
  if (ptr)
  {
    rcdir=malloc(sizeof(char) * (strlen(ptr) + strlen(RC_DIR) + 2));
    rcdir[0] = 0;
    strcat(rcdir, ptr);
    strcat(rcdir, "\\");
    strcat(rcdir, RC_DIR);

    if (dir_exists(rcdir) && !check_file_exists(rcdir, RC_GV)) return;
    else free(rcdir);
  }

/*  4. check luscus exe directory */
  rcdir = malloc(sizeof(char) * 512);
  GetCurrentDirectory(512, rcdir);               /*this should be tested if it works under windows!*/
  if (!check_file_exists(rcdir, RC_GV)) return;

/* If program comes to this point, it means that it can not locate config file */
  error_config_does_not_exist = 1;
#else
/*  1. check $LUSCUS_DIR environment  */
  ptr = getenv("LUSCUS_DIR");

#ifdef EBUG
  printf("checking environment LUSCUS_DIR: |%s|: %d\n", ptr, ptr == NULL ? 0 : 1);
#endif

  if (ptr)
  {
    rcdir = strdup(ptr);
#ifdef EBUG
    printf("checking path: %s\n", rcdir);
#endif
    if (!check_file_exists(rcdir, RC_GV)) return;
  }

/*  2. check $HOME/.luscus directory */
  ptr = getenv("HOME");
#ifdef EBUG
  printf("checking environment HOME: |%s|: %d\n", ptr, ptr == NULL ? 0 : 1);
#endif
  if (ptr)
  {
    rcdir=malloc(sizeof(char) * (strlen(ptr) + strlen(RC_DIR) + 2));
    rcdir[0] = 0;
    strcat(rcdir, ptr);
    strcat(rcdir, "/");
    strcat(rcdir, RC_DIR);

#ifdef EBUG
    printf("checking path: %s\n", rcdir);
#endif

    if (dir_exists(rcdir) && !check_file_exists(rcdir, RC_GV)) return;
    else free(rcdir);
  }

/*  3. check /etc/luscus directory */
#ifdef EBUG
  printf("checking directory: /etc/luscus\n");
#endif

  if (dir_exists("/etc/luscus"))
  {
    rcdir = strdup("/etc/luscus");

#ifdef EBUG
    printf("checking path: %s\n", rcdir);
#endif

    if (!check_file_exists(rcdir, RC_GV)) return;
    return;
  }

/*  4. check luscus exe directory */
  rcdir=malloc(sizeof(char) * 1024);
  getcwd(rcdir, 1024);
/*  rcdir = get_current_dir_name(); */
#ifdef EBUG
    printf("checking path: %s\n", rcdir);
#endif
  if (!check_file_exists(rcdir, RC_GV)) return;

/* If program comes to this point, it means that it can not locate config file */
  error_config_does_not_exist = 1;

#endif
}

void check_plugin_files(void)
{
  int iex;
#ifdef EBUG
  printf("CHECKING PLUGIN FILES\n");
#endif

  /*This is non-critical error*/
/*  iex = check_file_exists(rcdir, RC_ATOM);
  if (iex) error_config_does_not_exist = 2;*/

  iex = check_file_exists(rcdir, RC_CALC);
  if (iex && !error_config_does_not_exist) error_config_does_not_exist = 4;

  iex = check_file_exists(rcdir, RC_FRAG);
  if (iex && !error_config_does_not_exist) error_config_does_not_exist = 8;

  iex = check_file_exists(rcdir, RC_PLUG);
  if (iex && !error_config_does_not_exist) error_config_does_not_exist = 16;
}

void init_variables(void)
{
/*  char *ptr;
  int rc;*/

  accept_keys = 1;
  move_camera = 1;

  control_pressed = 0;

  init_config_dir_path();

  /*check existence of other plugin files*/

  if (!error_config_does_not_exist) check_plugin_files();
/*
  ptr=getenv("LUSCUS_DIR");
  if (ptr == NULL)
  {
#ifdef WINDOWS
    ptr=getenv("LOCALAPPDATA");
#else
    ptr=getenv("HOME");
#endif
#ifdef EBUG
  printf("HOME environment = |%s|\n", ptr);
#endif
    if (ptr == NULL) ptr = strdup(".");

    rcdir = malloc(sizeof(char) * (strlen(ptr) + strlen(RC_DIR) + 2));
    rcdir[0] = 0;
    strcat(rcdir, ptr);
#ifdef WINDOWS
    strcat(rcdir, "\\");
#else
    strcat(rcdir, "/");
#endif
    strcat(rcdir, RC_DIR);
  }
  else rcdir = strdup(ptr);

#ifdef EBUG
  printf("rcdir = |%s|\n", rcdir);
#endif

  if(!dir_exists(rcdir))
  {
#ifdef WINDOWS
    rc= mkdir(rcdir);
#else
    rc= mkdir(rcdir, S_IRWXU | S_IRWXG | S_IRWXO);
#endif
    if (rc) fprintf(stderr, "ERROR: Can't create resource directory %s\n",rcdir);
  }
*/
}

void load_default_values(void)
{
  Input_Data.hide_atoms = 0;
  Input_Data.hide_bonds = 0;
  Input_Data.show_axis = 0;
  Input_Data.show_epot = 0;
  Input_Data.style_atoms = 0;
  Input_Data.background_color[0] = 1.0;
  Input_Data.background_color[1] = 1.0;
  Input_Data.background_color[2] = 1.0;
  Input_Data.label_color[0] = 0.0;
  Input_Data.label_color[1] = 0.0;
  Input_Data.label_color[2] = 0.0;
  Input_Data.extracolor[0] = 0.3;
  Input_Data.extracolor[1] = 0.3;
  Input_Data.extracolor[2] = 0.3;
  Input_Data.extracolor[3] = 0.3;
  Input_Data.neg_pos_color[0][0]=0.953;
  Input_Data.neg_pos_color[0][1]=0.082;
  Input_Data.neg_pos_color[0][2]=0.020;
  Input_Data.neg_pos_color[0][3]=0.7;
  /*transp.*/
  Input_Data.neg_pos_color[1][0]=0.114;
  Input_Data.neg_pos_color[1][1]=0.686;
  Input_Data.neg_pos_color[1][2]=0.286;
  Input_Data.neg_pos_color[1][3]=0.7;

  Input_Data.electrostatic_poten_color[0][0]=0.953;
  Input_Data.electrostatic_poten_color[0][1]=0.082;
  Input_Data.electrostatic_poten_color[0][2]=0.020;
  Input_Data.electrostatic_poten_color[0][3]=0.7;
  /*transp.*/
  Input_Data.electrostatic_poten_color[1][0]=0.005;
  Input_Data.electrostatic_poten_color[1][1]=0.005;
  Input_Data.electrostatic_poten_color[1][2]=0.957;
  Input_Data.electrostatic_poten_color[1][3]=0.7;
  /*font*/
#ifdef WINDOWS
  g_strlcpy(Input_Data.font, "Arial 18", 100);
#else
  g_strlcpy(Input_Data.font, "Helvetica 18", 100);
#endif
  /*transp.*/
  Input_Data.fatness_a = 1.6;
  Input_Data.fatness_b = 1.6;
  Input_Data.animate = 0;
  Input_Data.move_light = 0;
  Input_Data.frequency_amplitude = 0.5;
  Input_Data.frequency_speed = 0.2;
  Input_Data.snapshot_delay = 0.05;
  Input_Data.init_screen_size = 800;
  Input_Data.translation_value = 0.2;
  Input_Data.angle_change_value = 2.0;
  Input_Data.torsion_change_value = 2.0;

  Input_Data.lev = 0.04;
  Input_Data.automatic_rebonding = 0;
  Input_Data.bond_order_one = 0;
/*  Input_Data.transp = 0.3;*/

  Input_Data.bw = 0;
  Input_Data.ntexture = 0;
  Input_Data.atomshape = 0;
  Input_Data.label_atoms = 0;
  Input_Data.povray = 0;
  Input_Data.ps = 0;
  Input_Data.move_camera = 0;
  Input_Data.force_symmetry = 0;

  Input_Data.lev_step = 0.0;
  Input_Data.lev_min = Input_Data.lev_max = 0.0;
  Input_Data.n_lev_steps = 100;
  Input_Data.auto_lev = 1;


}


void load_values_from_config_file(void)
{
  char *ptr;
  FILE *fp;
  char line[256];
  ptr = malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_GV) + 2));
  ptr[0] = 0;
  strcat(ptr, rcdir);
#ifdef WINDOWS
  strcat(ptr, "\\");
#else
  strcat(ptr, "/");
#endif
  strcat(ptr, RC_GV);

  fp = fopen(ptr, "r");
  free(ptr);

  if (fp == NULL) return;

  while(!feof(fp))
  {
    if (fgets(line, sizeof(line), fp) != NULL)
    {
           if (strstr(line, "Hide atoms"))             Input_Data.hide_atoms = my_read_int_value(line);
      else if (strstr(line, "Hide bonds"))             Input_Data.hide_bonds = my_read_int_value(line);
      else if (strstr(line, "Transparency"))           Input_Data.neg_pos_color[0][3] = Input_Data.neg_pos_color[1][3] =
                                                       Input_Data.electrostatic_poten_color[0][3] = 
                                                       Input_Data.electrostatic_poten_color[1][3] =
                                                       Input_Data.extracolor[3] = my_read_dble_value(line);

      else if (strstr(line, "Background color Red"))   Input_Data.background_color[0] = my_read_dble_value(line); 
      else if (strstr(line, "Background color Green")) Input_Data.background_color[1] = my_read_dble_value(line); 
      else if (strstr(line, "Background color Blue"))  Input_Data.background_color[2] = my_read_dble_value(line); 

      else if (strstr(line, "Extra color Red"))        Input_Data.extracolor[0] = my_read_dble_value(line); 
      else if (strstr(line, "Extra color Green"))      Input_Data.extracolor[1] = my_read_dble_value(line); 
      else if (strstr(line, "Extra color Blue"))       Input_Data.extracolor[2] = my_read_dble_value(line); 

      else if (strstr(line, "Label color Red"))        Input_Data.label_color[0] = my_read_dble_value(line); 
      else if (strstr(line, "Label color Green"))      Input_Data.label_color[1] = my_read_dble_value(line); 
      else if (strstr(line, "Label color Blue"))       Input_Data.label_color[2] = my_read_dble_value(line); 

      else if (strstr(line, "Plus  color Red"))        Input_Data.neg_pos_color[0][0] = my_read_dble_value(line); 
      else if (strstr(line, "Plus  color Green"))      Input_Data.neg_pos_color[0][1] = my_read_dble_value(line); 
      else if (strstr(line, "Plus  color Blue"))       Input_Data.neg_pos_color[0][2] = my_read_dble_value(line); 
      else if (strstr(line, "Minus color Red"))        Input_Data.neg_pos_color[1][0] = my_read_dble_value(line); 
      else if (strstr(line, "Minus color Green"))      Input_Data.neg_pos_color[1][1] = my_read_dble_value(line); 
      else if (strstr(line, "Minus color Blue"))       Input_Data.neg_pos_color[1][2] = my_read_dble_value(line); 

      else if (strstr(line, "Epot plus  color Red"))        Input_Data.electrostatic_poten_color[0][0] = my_read_dble_value(line); 
      else if (strstr(line, "Epot plus  color Green"))      Input_Data.electrostatic_poten_color[0][1] = my_read_dble_value(line); 
      else if (strstr(line, "Epot plus  color Blue"))       Input_Data.electrostatic_poten_color[0][2] = my_read_dble_value(line); 
      else if (strstr(line, "Epot minus color Red"))        Input_Data.electrostatic_poten_color[1][0] = my_read_dble_value(line); 
      else if (strstr(line, "Epot minus color Green"))      Input_Data.electrostatic_poten_color[1][1] = my_read_dble_value(line); 
      else if (strstr(line, "Epot minus color Blue"))       Input_Data.electrostatic_poten_color[1][2] = my_read_dble_value(line); 

      else if (strstr(line, "Fatness atoms"))          Input_Data.fatness_a = my_read_dble_value(line);
      else if (strstr(line, "Fatness bonds"))          Input_Data.fatness_b = my_read_dble_value(line);

      else if (strstr(line, "Animate"))                Input_Data.animate = my_read_int_value(line);
      else if (strstr(line, "Initial Screen Size"))    Input_Data.init_screen_size = my_read_int_value(line);
      else if (strstr(line, "Frequency amplitude"))    Input_Data.frequency_amplitude = my_read_dble_value(line);
      else if (strstr(line, "Frequency speed"))        Input_Data.frequency_speed = my_read_dble_value(line);
      else if (strstr(line, "Snapshot delay"))         Input_Data.snapshot_delay = my_read_dble_value(line);
      else if (strstr(line, "Isolevel value"))         Input_Data.lev = my_read_dble_value(line);
      else if (strstr(line, "Translation value"))      Input_Data.translation_value = my_read_dble_value(line);
      else if (strstr(line, "Angle change value"))     Input_Data.angle_change_value = my_read_dble_value(line);
      else if (strstr(line, "Torsion change value"))   Input_Data.torsion_change_value = my_read_dble_value(line);
      else if (strstr(line, "Auto select level"))      Input_Data.auto_lev = my_read_int_value(line);
      else if (strstr(line, "Auto rebond"))            Input_Data.automatic_rebonding = my_read_int_value(line);
    }
  }

  fclose(fp);

}

void print_version(void)
{
  printf("%s version: %s\n", PROGRAM_NAME, VERSION);
}

void load_values_from_command_line_args(int argc, char **argv)
{
  int j;
  int keya = 0;
  init_filename = NULL;
  /*omitted cmd line params
c, z, d, f, p, P, i, a, o, --out, --record, --play, --findsym, f, n, g, G
   */

  for(j = 1; j < argc; j++)
  {
    if(argv[j][0] != '-')
    {
      keya++;
      if (keya == 1) init_filename = argv[j];
    }
    else
    {
      switch(argv[j][1])
      {
        case '?':
        case 'h':
          Print_help();
          exit(0);
        case 'm': Input_Data.animate=1; break;
        case 's': sscanf(argv[j+1],"%d",&Input_Data.init_screen_size); break;
        case 'l': sscanf(argv[j+1],"%lf",&Input_Data.lev); break;
        case 'v':
          print_version();
          exit(0);
        case '-':
          if(strcmp(argv[j],"--version")==0)
          {
            print_version();
            exit(0);
          }
          printf("argv = |%s|\n", argv[j]);
          if (strcmp(argv[j],"--single") == 0)
          {
            printf("setting order-one bonds!\n");
            Input_Data.bond_order_one = 1;
            /*store data in order to show only single bonds*/
          }
          break;
        case 't': sscanf(argv[j+1], "%f", &Input_Data.extracolor[3]); break;
        case 'b':
          switch(argv[j+1][0])
          {
            case 'w':
              Input_Data.background_color[0] =1.0;
              Input_Data.background_color[1] =1.0;
              Input_Data.background_color[2] =1.0;
              break;
            case 'b':
              Input_Data.background_color[0] =0.0;
              Input_Data.background_color[1] =0.0;
              Input_Data.background_color[2] =0.0;
              break;
            case 'g':
              Input_Data.background_color[0] =0.5;
              Input_Data.background_color[1] =0.5;
              Input_Data.background_color[2] =0.5;
              break;
            default:
              Input_Data.background_color[0] =0.8;
              Input_Data.background_color[1] =0.8;
              Input_Data.background_color[2] =0.8;
              break;
          }
      }
    }
  }
}

void load_settings_data(int argc, char** argv)
{
  load_default_values();
  load_values_from_config_file();
  load_values_from_command_line_args(argc, argv);
#ifdef EBUG
  printf("determine single bonds = %d\n", Input_Data.bond_order_one);
#endif
}

void init_input_filetypes(void)
{
  int i;
  char line[1024];
  FILE *pc;
  char *ptr; 
  char *path;
  char *plugpath;

  n_input_types = 0;
  input_filetypes = NULL;

  plugpath = getenv("LUSCUS_DIR");
  if (plugpath)
  {
    path = (char*) malloc(sizeof(char) * (strlen(plugpath) + strlen(RC_PLUG) + 2));
    path[0] = 0;
    strcat(path, plugpath);
  }
  else
  {
    path = (char*) malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_PLUG) + 2));
    path[0] = 0;
    strcat(path, rcdir);
  }

/*  path = (char*) malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_PLUG) + 2));
  path[0] = 0;
  strcat(path, rcdir);*/

#ifdef WINDOWS
  strcat(path, "\\");
#else
  strcat(path, "/");
#endif
  strcat(path, RC_PLUG);
#ifdef EBUG
  printf("PATH = |%s|\n", path);
#endif

  /*first add luscus format*/
  n_input_types++;
  input_filetypes = (INPUT_FORMAT*) realloc(input_filetypes, n_input_types * sizeof(INPUT_FORMAT));
  input_filetypes[n_input_types-1].libpath = NULL;
  input_filetypes[n_input_types-1].forward = NULL;
  input_filetypes[n_input_types-1].backward = NULL;
  input_filetypes[n_input_types-1].extension = strdup(LUSCUS_FILE_EXT);
  input_filetypes[n_input_types-1].description = strdup(LUSCUS_FILE_DES);

  pc = fopen(path, "r");

  if (pc != NULL)
  {
    while(!feof(pc))
    {
      fgets(line, 1024, pc);
      strdiet(line);
#ifdef EBUG
      printf("line = |%s|\n", line);
#endif
      if (line[0] != 0)
      {
        n_input_types++;
        input_filetypes = (INPUT_FORMAT*) realloc(input_filetypes, n_input_types * sizeof(INPUT_FORMAT));

        input_filetypes[n_input_types-1].forward = NULL;
        input_filetypes[n_input_types-1].backward = NULL;
        input_filetypes[n_input_types-1].extension = NULL;
        input_filetypes[n_input_types-1].description = NULL;

      /*CHECK IF THIS IS ACCORDING TO WHAT IS SUPPOSED TO BE IN LUSCUS*/
/*        ptr = getenv("LUSCUS_DIR");
        if (ptr == NULL)
        {*/
          ptr = (char*) strcasestr(line,"libpath");
          input_filetypes[n_input_types-1].libpath = my_read_str_value(ptr);
/*        }
        else input_filetypes[n_input_types-1].libpath = strdup(ptr);*/

        if (input_filetypes[n_input_types-1].libpath != NULL)
        {
          ptr = (char*) strcasestr(line,"forward");
          if (ptr != NULL)
          {
            int iex;
            input_filetypes[n_input_types-1].forward = my_read_str_value(ptr);
            iex = check_file_exists(input_filetypes[n_input_types-1].libpath, input_filetypes[n_input_types-1].forward);
            if (iex) error_config_does_not_exist = 1024 + n_input_types;
          }
          ptr = (char*) strcasestr(line,"backward");
          if (ptr != NULL)
          {
            int iex;
            input_filetypes[n_input_types-1].backward = my_read_str_value(ptr);
            iex = check_file_exists(input_filetypes[n_input_types-1].libpath, input_filetypes[n_input_types-1].backward);
            if (iex) error_config_does_not_exist = 2048 + n_input_types;
          }
          ptr = (char*) strcasestr(line,"extension");
          if (ptr != NULL) input_filetypes[n_input_types-1].extension = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"description");
          if (ptr != NULL) input_filetypes[n_input_types-1].description = my_read_str_value(ptr);
        }
      }
      line[0] = 0;
    }
  }
  free(path);

#ifdef EBUG
  printf("n_input_filetypes = %d\n", n_input_types);
  for(i = 0; i < n_input_types; i++)
    printf("type=|%s|\n", input_filetypes[i].description);
#endif
}

void open_cmd_line_file(void)
{
  if (init_filename) open_file(init_filename, NULL); /*FILE OPENED FROM COMMAND LINE WILL BE RECOGNIZED BY THE FIRST MATCHING EXTENSION*/
  else mol = new_mol(1);
  return;
}

void load_fragment_data(void)
{
  int i;
  FILE *fc;
  char *ptr;
  char *line;

/*  printf("loading fragment data\n"); fflush(stdout);*/

  ptr = malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_FRAG) + 2));
  ptr[0] = 0;
  strcat(ptr, rcdir);
#ifdef WIDNOWS
  strcat(ptr, "\\");
#else
  strcat(ptr, "/");
#endif
  strcat(ptr, RC_FRAG);

/*  printf("fragment path = |%s|\n", ptr); fflush(stdout);*/
  fc = fopen(ptr, "r");
  if (!fc) fprintf(stderr, "ERROR: Can't load fragments\n");
  free(ptr);

  n_fragments = 0;
  frag = NULL;

  if (fc)
  {
    while(!feof(fc))
    {
      line = read_line(fc);

      if (!feof(fc))
      {
        n_fragments++;
        frag = (FRAG_DATA*) realloc(frag, sizeof(FRAG_DATA) * n_fragments);
        frag[n_fragments-1].coord_file_name = get_one_word(get_ptr_ith_word(line,1));
        frag[n_fragments-1].image_file_name = get_one_word(get_ptr_ith_word(line,2));
/*        printf("nfragments = %d; frag = |%s| |%s|\n", n_fragments, frag[n_fragments-1].coord_file_name, frag[n_fragments-1].image_file_name); fflush(stdout);*/
      }
      free(line);
    }
    fclose(fc);
  }

#ifdef EBUG
  printf("FRAGMENT CHECK!\n"); fflush(stdout);
  for(i = 0; i < n_fragments; i++)
  {
    printf("file |%s| image |%s|\n", frag[i].coord_file_name, frag[i].image_file_name); fflush(stdout);
  }
#endif
}

void load_calc_def_data(void)
{
  int i;
  FILE *fc;
  char *ptr;
  char *line;
  char buf[512];

  ptr = malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_CALC) + 2));
  ptr[0] = 0;
  strcat(ptr, rcdir);
#ifdef WIDNOWS
  strcat(ptr, "\\");
#else
  strcat(ptr, "/");
#endif
  strcat(ptr, RC_CALC);

  fc = fopen(ptr, "r");
  if (!fc) fprintf(stderr, "WARNING: Can't find data for calculations\n");
  free(ptr);

  n_calc_defs = 0;
  calc_defs = NULL;

  if (fc)
  {
    while(!feof(fc))
    {
      line = read_line(fc);
/*      strdiet(line);*/
      if (line[0] != 0)
      {
        n_calc_defs++;
        calc_defs = (DEF_CALC_T*) realloc(calc_defs, sizeof(DEF_CALC_T) * n_calc_defs);
        calc_defs[n_calc_defs-1].type = NULL;
        calc_defs[n_calc_defs-1].description = NULL;
        calc_defs[n_calc_defs-1].ext_out = NULL;
        calc_defs[n_calc_defs-1].libpath = NULL;
        calc_defs[n_calc_defs-1].plugin_name = NULL;
        calc_defs[n_calc_defs-1].extraargsym = NULL;
        calc_defs[n_calc_defs-1].extraargval = NULL;
        calc_defs[n_calc_defs-1].extraargdes = NULL;
        calc_defs[n_calc_defs-1].extrafile = 0;
        calc_defs[n_calc_defs-1].iinp = 0;
        calc_defs[n_calc_defs-1].iout = 0;
        ptr = getenv("LUSCUS_DIR");
        if (ptr == NULL)
        {
          ptr = (char*) strcasestr(line,"libpath");
          calc_defs[n_calc_defs-1].libpath = my_read_str_value(ptr);
        }
        else calc_defs[n_calc_defs-1].libpath = strdup(ptr);

        if (calc_defs[n_calc_defs-1].libpath != NULL)
        {
          ptr = (char*) strcasestr(line,"type");
          if (ptr != NULL) calc_defs[n_calc_defs-1].type = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"description");
          if (ptr != NULL) calc_defs[n_calc_defs-1].description = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"ext_out");
          if (ptr != NULL) calc_defs[n_calc_defs-1].ext_out = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"plugin_name");
          if (ptr != NULL) calc_defs[n_calc_defs-1].plugin_name = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"extraargsym");
          if (ptr != NULL) calc_defs[n_calc_defs-1].extraargsym = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"extraargval");
          if (ptr != NULL) calc_defs[n_calc_defs-1].extraargval = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"extraargdes");
          if (ptr != NULL) calc_defs[n_calc_defs-1].extraargdes = my_read_str_value(ptr);
          ptr = (char*) strcasestr(line,"extrafile");
          if (ptr != NULL) calc_defs[n_calc_defs-1].extrafile = 1;

	  if (calc_defs[n_calc_defs-1].ext_out == NULL)
	  {
            printf("ERROR!!!!\n");
            calc_defs[n_calc_defs-1].plugin_name = NULL; /*diseable this plugin*/
            fprintf(stderr, "WARNING: Can't find plugin for %s. Optimization with %s will be diseabled!\n", calc_defs[n_calc_defs-1].description, calc_defs[n_calc_defs-1].description);
	  }
	  else
	  {
            for(calc_defs[n_calc_defs-1].iout = 0; calc_defs[n_calc_defs-1].iout < n_input_types; calc_defs[n_calc_defs-1].iout++)
	    {
              if (strcmp(calc_defs[n_calc_defs-1].ext_out, input_filetypes[calc_defs[n_calc_defs-1].iout].extension) == 0) break;
	    }
            if (calc_defs[n_calc_defs-1].iout == n_input_types) /*this input extension is not registered*/
            {
              if (calc_defs[n_calc_defs-1].plugin_name) free(calc_defs[n_calc_defs-1].plugin_name);
              calc_defs[n_calc_defs-1].plugin_name = NULL; /*diseable this plugin*/
/*            snprintf(buf, 512, "WARNING: Can't find plugin for %s.\nOptimization with %s will be diseabled!", calc_defs[n_calc_defs-1].description, calc_defs[n_calc_defs-1].description);
            make_warning(buf);*/
              fprintf(stderr, "WARNING: Can't find plugin for %s. Optimization with %s will be diseabled!\n", calc_defs[n_calc_defs-1].description, calc_defs[n_calc_defs-1].description);
            }
	  }
        }
      }
      free(line);
    }
    fclose(fc);
  }

#ifdef EBUG
  printf("-------calculation data check--------\n");
  for (i = 0; i < n_calc_defs; i++)
  {
    printf("type=|%s| description=|%s| ext_out=|%s| libpath=|%s| plugin_name=|%s| inum=|%d| onum=|%d| extraargsym=|%s| extraargval=|%s| extrafile=|%d|\n", calc_defs[i].type, calc_defs[i].description, calc_defs[i].ext_out, calc_defs[i].libpath, calc_defs[i].plugin_name, calc_defs[i].iinp, calc_defs[i].iout, calc_defs[i].extraargsym, calc_defs[i].extraargval, calc_defs[i].extrafile);
  }
  printf("-------------------------------------\n");
#endif

}

void save_settings_data(void)
{
  char *ptr;
  FILE *fp;
  char line[256];
  ptr = malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_GV) + 2));
  ptr[0] = 0;
  strcat(ptr, rcdir);
#ifdef WINDOWS
  strcat(ptr, "\\");
#else
  strcat(ptr, "/");
#endif
  strcat(ptr, RC_GV);

  printf("saving settings data in |%s|\n", ptr); fflush(stdout);

  fp = fopen(ptr, "w");
  free(ptr);

  if (fp == NULL) return;

  fprintf(fp, "general settings variables\n");
  fprintf(fp, "Hide atoms = %d\n", Input_Data.hide_atoms);
  fprintf(fp, "Hide bonds = %d\n", Input_Data.hide_bonds);
  fprintf(fp, "Transparency = %f\n", Input_Data.extracolor[3]);

  fprintf(fp, "Background color Red = %f\n", Input_Data.background_color[0]);
  fprintf(fp, "Background color Green = %f\n", Input_Data.background_color[1]);
  fprintf(fp, "Background color Blue = %f\n", Input_Data.background_color[2]);

  fprintf(fp, "Extra color Red = %f\n", Input_Data.extracolor[0]);
  fprintf(fp, "Extra color Green = %f\n", Input_Data.extracolor[1]);
  fprintf(fp, "Extra color Blue = %f\n", Input_Data.extracolor[2]);

  fprintf(fp, "Label color Red = %f\n", Input_Data.label_color[0]);
  fprintf(fp, "Label color Green = %f\n", Input_Data.label_color[1]);
  fprintf(fp, "Label color Blue = %f\n", Input_Data.label_color[2]);

  fprintf(fp, "Plus  color Red = %f\n", Input_Data.neg_pos_color[0][0]);
  fprintf(fp, "Plus  color Green = %f\n", Input_Data.neg_pos_color[0][1]);
  fprintf(fp, "Plus  color Blue = %f\n", Input_Data.neg_pos_color[0][2]);
  fprintf(fp, "Minus color Red = %f\n", Input_Data.neg_pos_color[1][0]);
  fprintf(fp, "Minus color Green = %f\n", Input_Data.neg_pos_color[1][1]);
  fprintf(fp, "Minus color Blue = %f\n", Input_Data.neg_pos_color[1][2]);

  fprintf(fp, "Epot plus  color Red = %f\n", Input_Data.electrostatic_poten_color[0][0]);
  fprintf(fp, "Epot plus  color Green = %f\n", Input_Data.electrostatic_poten_color[0][1]);
  fprintf(fp, "Epot plus  color Blue = %f\n", Input_Data.electrostatic_poten_color[0][2]);
  fprintf(fp, "Epot minus color Red = %f\n", Input_Data.electrostatic_poten_color[1][0]);
  fprintf(fp, "Epot minus color Green = %f\n", Input_Data.electrostatic_poten_color[1][1]);
  fprintf(fp, "Epot minus color Blue = %f\n", Input_Data.electrostatic_poten_color[1][2]);

  fprintf(fp, "Fatness atoms = %f\n", Input_Data.fatness_a);
  fprintf(fp, "Fatness bonds = %f\n", Input_Data.fatness_b);

  fprintf(fp, "Animate = %d\n", Input_Data.animate);
  fprintf(fp, "Initial Screen Size = %d\n", Input_Data.init_screen_size);
  fprintf(fp, "Frequency amplitude = %f\n", Input_Data.frequency_amplitude);
  fprintf(fp, "Frequency speed = %f\n", Input_Data.frequency_speed);
  fprintf(fp, "Snapshot delay = %f\n", Input_Data.snapshot_delay);
  fprintf(fp, "Isolevel value = %f\n", Input_Data.lev);
  fprintf(fp, "Translation value = %f\n", Input_Data.translation_value);
  fprintf(fp, "Angle change value = %f\n", Input_Data.angle_change_value);
  fprintf(fp, "Torsion change value = %f\n",  Input_Data.torsion_change_value);
  fprintf(fp, "Auto select level = %d\n", Input_Data.auto_lev);
  fprintf(fp, "Auto rebond = %d\n", Input_Data.automatic_rebonding);

  fclose(fp);
}

void remove_calc_def_data(void)
{
  int i;

  for (i = 0; i < n_calc_defs; i++)
  {
    if (calc_defs[i].type) free(calc_defs[i].type);
    if (calc_defs[i].description) free(calc_defs[i].description);
    if (calc_defs[i].ext_out) free(calc_defs[i].ext_out);
    if (calc_defs[i].libpath) free(calc_defs[i].libpath);
    if (calc_defs[i].plugin_name) free(calc_defs[i].plugin_name);
    if (calc_defs[i].extraargsym) free(calc_defs[i].extraargsym);
    if (calc_defs[i].extraargval) free(calc_defs[i].extraargval);
    if (calc_defs[i].extraargdes) free(calc_defs[i].extraargdes);
  }
  if (rcdir) free(rcdir);
}
