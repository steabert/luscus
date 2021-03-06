/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include "gvgrp.h"

/*general parameters*/
#define PROGRAM_NAME "luscus"
/*#define VERSION "0.8.5 beta"*/
/*config files parameters*/
#ifdef WINDOWS
#define RC_DIR "luscus"
#else
#define RC_DIR ".luscus"
#endif
#define RC_ATOM "atoms.rc"
#define RC_FRAG "fragments.rc"
#define RC_PLUG "plugin.rc"
#define RC_CALC "calculation.rc"
#define RC_GV "luscus.rc"
/*#define UNDO "gv.undo"*/

/*ishow*/
#define HAS_

#define HAS_ATOMNAMES 0x0001
#define HAS_MULLIKEN 0x0002
#define HAS_LOPROP 0x0004
#define HAS_ATOMNUMS 0x0008
#define HAS_ATOMSYMM 0x0010
#define HAS_DIPOLE 0x0020
#define HAS_IRINT 0x0040
#define HAS_RAMAN 0x0080
#define HAS_ENERGY 0x0100
#define HAS_RMS_GRAD 0x0200
#define HAS_MAX_GRAD 0x0400
#define HAS_E_POT 0x0800
#define HAS_COLOUR 0x1000

/*screenshot file formats*/
#define SCSHOT_TGA 0
#define SCSHOT_PS 1
#define SCSHOT_PS2 2
#define SCSHOT_POV 3

#define NUM_OF_BONDS 7
#define NUM_OF_ANGLES 5
/*definition of bonds*/
#define NO_BOND 0
#define SINGLE_BOND 1
#define DOUBLE_BOND 2
#define TRIPLE_BOND 3
#define PARTIAL_BOND 4
#define S_P_BOND 5
#define LINE_BOND 6

/*definition of orbital types*/
 /*orbital types*/
#define NUM_ORB_TYPES 8

#define ORBITAL_TYPE_U 0
#define ORBITAL_TYPE_F 1
#define ORBITAL_TYPE_I 2
#define ORBITAL_TYPE_1 3
#define ORBITAL_TYPE_2 4
#define ORBITAL_TYPE_3 5
#define ORBITAL_TYPE_S 6
#define ORBITAL_TYPE_D 7

/*grid types*/
#define ORBITAL 1
#define DENSITY 2 
#define E_POTEN 3 
#define CUSTOM 4

/*definition of internal coordinate types*/
#define C_BOND 2
#define C_ANGLE 3
#define C_DIHEDRAL 4

/*maximal extension length*/
#define MAX_EXT_LENGTH 7

/*luscus file format*/
#define LUSCUS_FILE_EXT "lus"
#define LUSCUS_FILE_DES "Luscus native format"

#define MAX_ATOMS_FOR_BONDING 300

/*char *bond_description[NUM_OF_BONDS] =
{
  "No bond",
  "Single bond",
  "double bond",
  "triple bond",
  "partial bond",
  "One and half bond",
  "Line"
};*/

/*char *orbital_description[NUM_ORB_TYPES] =
{
  "undefined",
  "frozen",
  "inactive",
  "RAS 1",
  "RAS 2",
  "RAS 3",
  "secondary",
  "deleted"
};*/

/*const char orbital_description_1[NUM_ORB_TYPES][2] =
{
  "u", "f", "i", "1", "2", "3", "s", "d"
};*/

/*column identifiers in treeview*/
enum 
{
  ORBITAL_TYPE,
  ORBITAL_SYMM,
  ORBITAL_NUM,
  ORBITAL_ENERG,
  ORBITAL_OCC,
  ORBITAL_EDITED,
  NUM_COLUMNS
};

enum
{
  GEOMETRY_OBJECT,
  GEOM_OBJ_COLOR,
  GEOM_OBJ_TRANSP,
  NUM_GEO_COLUMS
};

typedef float color_t[4];
typedef double XYZ[3];

typedef struct bond_data
{
  int iat1;
  int iat2;
  int bond_type;
} BOND;

typedef struct element_data
{
  char *name;
  color_t color;
  float vdw_rad;
  float bond_rad;
  int valency;
  int periodic_pos_x;
  int periodic_pos_y;
} ELEM_DATA;

typedef struct default_elem
{
  int n_default_elem;
  int *default_elem;
} DEFAULT_ELEM;

typedef struct grid_t
{
  int ngrid;
  int ngrid2;
  int nmo;
  int nmo2;
  int npoints;
  int ptotal;
  int npt[3]; /*number of points in 3 dimenstions*/
  int block_size;
  int nBlock;
  int isCutOff;
  double CutOff;
  int nPointsCutOff;
  int *ByteCutOff;
  long *pos;
  long *pos2;
  int *dependence;
  double len[3];
  double origin[3];
  double axisvec1[3], axisvec2[3], axisvec3[3];
  int iniIndex[7];
  char *title;
  char *iType;
  char Gtitle[81];
  double minval, maxval, guess;
  double currentvalue;
  double isod[21];
  double *values;
  int sel[4];
  int is_coord;
/* filter */
  char FltSym[16];
  char FltInd[8];
  double FltMinOcc, FltMaxOcc, FltMinEne, FltMaxEne;
  char **titlesArr;
  int *FltTitle;
  int isFltSym,isFltOcc,isFltEne,isFltInd;
} GRID_T;

typedef struct stxt_t
{
  unsigned int width;
  unsigned int height;
  unsigned char *pixels;
} STXT_T;

typedef struct pixdata
{
  STXT_T pix_index;
  STXT_T pix_symm;
  STXT_T pix_addnum;
  STXT_T pix_name;
  STXT_T pix_charge_m;
  STXT_T pix_charge_l;
} PIXDATA;

typedef struct textbox_t
{
  int coord_x;
  int coord_y;
  char *message;
  char *font;
  color_t color;
  STXT_T pixtext;
} TEXTBOX_T;

typedef struct molecule
{
  unsigned int ishow;

  /*atoms*/
  int natom;
  XYZ *xyz;
  ELEM_DATA *elem;/*data determined from element data*/

  /*adiitional info - atoms*/
  double *charge_m;
  double *charge_l;
  int *additional_numeration;
  int *symmetry;
  char **name;

  PIXDATA *pixdata;

  char editable;

  int n_selected;
  int selected[4];
  int n_marked;
  int *marked;

  int nbond;
  BOND *bond;

  XYZ dipole;
  
  int ngrids;
  int *grid_type;
  double *grid_energy;
  double *grid_occ;
  int *grid_symmetry;
  int *grid_index;
  int *orbital_type;
  int *edited;
  long orb_starting_position;
  long orb_ending_position;

  double geo_energy;
  double max_grad;
  double rms_grad;

  int nvibr; /*number of vibrations*/
  double *freq; /*frequency in cm^-1*/
  double *ir_intensity;
  double *raman_intensity;
  XYZ **normal_mode;

  int nvector;
  XYZ *vector1;
  XYZ *vector2;
  double *radius;
  double *sharpness;
  color_t *vector_color;

  int ntriangle;
  XYZ *triangle1;
  XYZ *triangle2;
  XYZ *triangle3;
  color_t *triangle_color;

  int nsphere;
  XYZ *sphere_center;
  double *sphere_radius;
  color_t *sphere_color;

  int nsurf;
  XYZ *surf1;
  XYZ *surf2;
  XYZ *surf3;
  color_t *surf_color;

  int ncells;
  XYZ *cell1;
  XYZ *cell2;
  XYZ *cell3;
  XYZ *cell4;
  color_t *cell_color;

  int ntextboxes;
  TEXTBOX_T *textboxes;

  GRID_T grid;
  double *epot; /*electrostatic potential values!*/
} MOL;

typedef struct frag
{
  char *coord_file_name;
  char *image_file_name;
} FRAG_DATA;

typedef struct input_data
{
  int hide_atoms;
  int hide_bonds;
  int show_axis;
  int show_epot;
  int style_atoms;

  color_t background_color;
  color_t label_color;
  color_t extracolor;   /*this is default*/
  color_t neg_pos_color[2]; /*orbital negative and positive color*/
  color_t electrostatic_poten_color[2];
  char font[100];
  double fatness_a;      /*relative atom size*/
  double fatness_b;      /*relative bond size*/
  int animate;
  int move_light;
  double frequency_amplitude;
  double frequency_speed;
  int init_screen_size;
  double snapshot_delay;
/*  double isolevel_value;*/
  double symmetry_translation;
  double translation_value;
  double angle_change_value;
  double torsion_change_value;
/*  double transp;*/

  int bw;
  int ntexture;
  int atomshape;
  int label_atoms; /*0->no labels;
		     1->atom symbols
		     2->alternative atom labels
		     3->atom numbers
		     4->alternative atom numbers
		     5->Mulliken charges
		     6->Loprop charges*/
  /*default vector radius and sharpness*/
  int povray;
  int ps;
  int move_camera;
  int force_symmetry;

  /*grid stuff*/
  double lev_step;
  double lev_min,lev_max;
  int n_lev_steps;
  int auto_lev; /* calculate the initial isolevel automatically */
  double lev;
  int max_diam_ok;  /* diameter is recalculated and valid */
  int automatic_rebonding;
  int bond_order_one; /*Determine only single bonds, do not wast time on multiple bonds!*/
} INPUT_DATA;

typedef struct input_format_t
{
  char *extension;
  char *description;
  char *forward; /*conversion to luscus*/
  char *backward; /*conversion from luscus*/
  char *libpath;
} INPUT_FORMAT;

typedef struct def_calc_t
{
  char *type;
  char *description;
  char *ext_out;

  char *libpath;
  char *plugin_name;

  /*extra argument*/
  char *extraargsym;
  char *extraargval;
  char *extraargdes;

  /*extra file*/
  char extrafile;

  /*references on INPUT_FORMAT structure*/
  /*not in rc file; these values are obtained by post-processing*/
  int iinp;
  int iout;
} DEF_CALC_T;

typedef struct watched_coordinate
{
  int coord_type;
  int atom[4];
} WATCH;

extern char *rcdir;
extern int number_of_elements;
extern INPUT_DATA Input_Data;
extern MOL *mol; /*allocatable array that holds all geometries*/
extern MOL *m; /*poitnter that is used to access one geometry from the mol*/

int n_geometries;  /*total number of geometries*/
long *geo_file_position;
int igeo;  /*currenty geometry*/
int ivib; /*current vibration mode*/
int iorb; /*current orbital*/

int accept_keys; /*inactivate hotkeys*/
int move_camera; /*camera movement*/
XYZ move_molecule; /*molecule movement*/

extern ELEM_DATA *e;
extern DEFAULT_ELEM def_elem;
extern int n_input_types;
extern INPUT_FORMAT *input_filetypes;
extern int n_fragments;
extern FRAG_DATA *frag;
extern int n_calc_defs;
extern DEF_CALC_T *calc_defs;
extern int n_watched;
extern WATCH *watch;
extern int control_pressed;
extern int textbox_state;
