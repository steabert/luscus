/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#ifdef WINDOWS
#include<strings.h>
#endif
#include<math.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#define MAXCHAR 127
#define MIN(x,y) (x<y? x : y)
#define MAX(x,y) (x>y? x : y)
#define R529 0.52917721067e0

typedef unsigned char byte_t;
MOL *m;
MOL *mol;
int number_of_elements;
ELEM_DATA *e;
INPUT_DATA Input_Data;
char *input_file_name = NULL;
char *input_file_type = NULL; /*19.02.2013 added in order to discriminate file types*/
char *current_directory = NULL;
int should_bonding_be_determined;
static int isBinary;
static int isPacked;
static int Fmarker = 0;
static time_t init_file_time = 0;
static time_t checked_file_time;
char file_open_success = 0;

static int ReadBlock(double*, int, FILE*);

/* Packing parameters */
struct packing_params_t
{
  int range; /* really, half range */
  int nbytes; /* bytes per packed value */
  
  double xlimit[4];
  double xdelta[3];
  int ylimit[4];
  int ydelta[3];
  
  int ileft, iright;
};

/*definition of function*/
static int getPackingParams(const char*, struct packing_params_t*);
static void unpackData(double*, byte_t*, unsigned int, int, const struct packing_params_t*);

void deallocate_mol_data(void)
{
  int i;
#ifdef EBUG
  printf("DEALLOCATE MOL DATA\n"); fflush(stdout);
#endif
  deallocate_msrf();
  for (i = 0; i < n_geometries; i++)
    deallocate_data(&mol[i]);

  free(mol);
  n_geometries = 0;

  if (geo_file_position) free(geo_file_position);
  geo_file_position = NULL;
  m = NULL;

  remove_all_watched_coords();
}

void deallocate_data(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE DATA\n"); fflush(stdout);
#endif
  deallocate_atoms(mp);
  deallocate_bonds(mp);
  deallocate_dipole(mp);
  deallocate_grids(mp);
  deallocate_vibrations(mp);
  deallocate_vectors(mp);
  deallocate_triangles(mp);
  deallocate_spheres(mp);
  deallocate_surfaces(mp);
  deallocate_cells(mp);
  deallocate_textboxes(mp);
  unmark_all();
  unselect_all();
}

void deallocate_atoms(MOL *mp)
{
  int i;
  if (!mp) return;

  if (mp->xyz) free(mp->xyz);
  mp->xyz = NULL;

  if (mp->charge_m) free(mp->charge_m);
  mp->charge_m = NULL;

  if (mp->charge_l) free(mp->charge_l);
  mp->charge_l = NULL;

  if (mp->additional_numeration) free(mp->additional_numeration);
  mp->additional_numeration = NULL;

  if (mp->symmetry) free(mp->symmetry);
  mp->symmetry = NULL;

  if (mp->name)
  {
    for(i = 0; i < mp->natom; i++)
      if (mp->name[i])
       	free(mp->name[i]);
    free(mp->name);
  }
  mp->name = NULL;

  if (mp->pixdata)
  {
    for(i = 0; i < mp->natom; i++)
    {
      if (mp->pixdata[i].pix_index.pixels) free(mp->pixdata[i].pix_index.pixels);
      if (mp->pixdata[i].pix_addnum.pixels) free(mp->pixdata[i].pix_addnum.pixels);
      if (mp->pixdata[i].pix_symm.pixels) free(mp->pixdata[i].pix_symm.pixels);
      if (mp->pixdata[i].pix_name.pixels) free(mp->pixdata[i].pix_name.pixels);
      if (mp->pixdata[i].pix_charge_m.pixels) free(mp->pixdata[i].pix_charge_m.pixels);
      if (mp->pixdata[i].pix_charge_l.pixels) free(mp->pixdata[i].pix_charge_l.pixels);
    }
    free(mp->pixdata);
  }

  if (mp->elem)
  {
    for(i = 0; i < mp->natom; i++) 
      if (mp->elem[i].name)
        free(mp->elem[i].name);
    free(mp->elem);
  }
  mp->elem = NULL;

  mp->n_selected = 0;
  mp->ishow = 0;
  mp->natom = 0;

/*  if (m->ishow & HAS_ATOMNAMES) m->ishow ^= HAS_ATOMNAMES;
  if (m->ishow & HAS_MULLIKEN) m->ishow ^= HAS_MULLIKEN;
  if (m->ishow & HAS_LOPROP) m->ishow ^= HAS_LOPROP;
  if (m->ishow & HAS_ATOMNUMS) m->ishow ^= HAS_ATOMNUMS;
  if (m->ishow & HAS_ATOMSYMM) m->ishow ^= HAS_ATOMSYMM;*/

  /*deallocate atom list*/
#ifdef EBUG
  printf("deallocate atom list\n"); fflush(stdout);
#endif
  deallocate_atom_list();
}

void deallocate_bonds(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE BONDS\n"); fflush(stdout);
#endif
  if (mp->bond) free(mp->bond);
  mp->bond = NULL;
  mp->nbond = 0;
}

void deallocate_dipole(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE DIPOLE\n"); fflush(stdout);
#endif
  /*no dealocations, just unmark dipole for drawing!*/
  if (mp->ishow & HAS_DIPOLE) mp->ishow ^= HAS_DIPOLE;
}

void deallocate_grids(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE GRIDS\n"); fflush(stdout);
#endif
  if (mp->grid_type) free(mp->grid_type); 
  mp->grid_type = NULL;
  if (mp->grid_energy) free(mp->grid_energy);
  mp->grid_energy = NULL;
  if (mp->grid_occ) free(mp->grid_occ);
  mp->grid_occ = NULL;
  if (mp->grid_symmetry) free(mp->grid_symmetry);
  mp->grid_symmetry = NULL;
  if (mp->grid_index) free(mp->grid_index);
  mp->grid_index = NULL;
  if (mp->orbital_type) free(mp->orbital_type);
  mp->orbital_type = NULL;
  if (mp->edited) free(mp->edited);
  mp->edited = NULL;

  if (mp->grid.ByteCutOff) free(mp->grid.ByteCutOff);
  mp->grid.ByteCutOff = NULL;
  if (mp->grid.pos) free(mp->grid.pos);
  mp->grid.pos = NULL;
  if (mp->grid.title) free(mp->grid.title);
  mp->grid.title = NULL;
  if (mp->grid.titlesArr) free(mp->grid.titlesArr);
  mp->grid.titlesArr = NULL;
  if (mp->grid.FltTitle) free(mp->grid.FltTitle);
  mp->grid.FltTitle = NULL;
  if (mp->grid.iType) free(mp->grid.iType);
  mp->grid.iType = NULL;
  if (mp->grid.dependence) free(mp->grid.dependence);
  mp->grid.dependence = NULL;
  if (mp->grid.values) free(mp->grid.values);
  mp->grid.values = NULL;

  mp->ngrids = 0;
  mp->grid.ngrid = 0;
  mp->grid.npoints = 0;
  if (mp->epot) free(mp->epot);
  mp->epot = NULL;
}

void deallocate_vibrations(MOL *mp)
{
  int i;
#ifdef EBUG
  printf("DEALLOCATE VIBRATIONS\n"); fflush(stdout);
#endif
  if (mp->freq) free(mp->freq);
  mp->freq = NULL;
  if (mp->ir_intensity) free(mp->ir_intensity);
  mp->ir_intensity = NULL;
  if (mp->raman_intensity) free(mp->raman_intensity);
  mp->raman_intensity = NULL;
  if (mp->normal_mode)
  {
    for(i = 0; i < mp->nvibr; i++) free(mp->normal_mode[i]);
    free(mp->normal_mode);
  }
  mp->normal_mode = NULL;
  mp->nvibr = 0;
  if (mp->ishow & HAS_IRINT) mp->ishow ^= HAS_IRINT;
  if (mp->ishow & HAS_RAMAN) mp->ishow ^= HAS_RAMAN;
}

void deallocate_vectors(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE VECTORS\n"); fflush(stdout);
#endif
  if (mp->vector1) free(mp->vector1);
  if (mp->vector2) free(mp->vector2);
  if (mp->radius) free(mp->radius);
  if (mp->sharpness) free(mp->sharpness);
  if (mp->vector_color) free(mp->vector_color);
  mp->vector1 = NULL;
  mp->vector2 = NULL;
  mp->radius = NULL;
  mp->sharpness = NULL;
  mp->vector_color = NULL;
  mp->nvector = 0;
}

void deallocate_triangles(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE TRIANGLES\n"); fflush(stdout);
#endif
  if (mp->triangle1) free(mp->triangle1);
  if (mp->triangle2) free(mp->triangle2);
  if (mp->triangle3) free(mp->triangle3);
  if (mp->triangle_color) free(mp->triangle_color);

  mp->triangle1 = NULL;
  mp->triangle2 = NULL;
  mp->triangle3 = NULL;
  mp->triangle_color = NULL;
  mp->ntriangle = 0;
}

void deallocate_spheres(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE SPHERES\n"); fflush(stdout);
#endif
  if (mp->sphere_center) free(mp->sphere_center);
  if (mp->sphere_radius) free(mp->sphere_radius);
  if (mp->sphere_color) free(mp->sphere_color);

  mp->sphere_center = NULL;
  mp->sphere_radius = NULL;
  mp->sphere_color = NULL;
  mp->nsphere = 0;
}

void deallocate_surfaces(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE SURFACES\n"); fflush(stdout);
#endif
  if (mp->surf1) free(mp->surf1);
  if (mp->surf2) free(mp->surf2);
  if (mp->surf3) free(mp->surf3);
  if (mp->surf_color) free(mp->surf_color);

  mp->surf1 = NULL;
  mp->surf2 = NULL;
  mp->surf3 = NULL;
  mp->surf_color = NULL;
  mp->nsurf = 0;
}

void deallocate_cells(MOL *mp)
{
#ifdef EBUG
  printf("DEALLOCATE CELLS\n"); fflush(stdout);
#endif
  if (mp->cell1) free(mp->cell1);
  if (mp->cell2) free(mp->cell2);
  if (mp->cell3) free(mp->cell3);
  if (mp->cell4) free(mp->cell4);
  if (mp->cell_color) free(mp->cell_color);

  mp->cell1 = NULL;
  mp->cell2 = NULL;
  mp->cell3 = NULL;
  mp->cell4 = NULL;
  mp->cell_color = NULL;
  mp->ncells = 0;
}

void deallocate_textboxes(MOL *mp)
{
  int i;
#ifdef EBUG
  printf("DEALLOCATE TEXTBOXES\n"); fflush(stdout);
#endif
  if (!mp) return;
  for(i = 0; i < mp->ntextboxes; i++)
  {
    if (mp->textboxes[i].message) free(mp->textboxes[i].message);
    if (mp->textboxes[i].font) free(mp->textboxes[i].font);
  }
  if (mp->textboxes) free(mp->textboxes);
  mp->ntextboxes = 0;
}

MOL *new_mol(int number_of_geometries)
{
  int i;
  MOL *mp;
#ifdef EBUG
  printf("NEW MOL\n"); fflush(stdout);
#endif
  n_geometries = number_of_geometries;
  igeo = 0;
  ivib = 0;
  iorb = 0;

  if (number_of_geometries <= 0)
  {
    n_geometries = 0;
    return NULL;
  }
  mp = malloc(sizeof(MOL) * number_of_geometries);

  mp->editable = 1;

  for(i = 0; i < number_of_geometries; i++)
  {
    m = mp + i;

    m->ishow = 0;

    m->natom = 0;
    m->xyz = NULL;
    m->elem = NULL;

    m->charge_m = NULL;
    m->charge_l = NULL;
    m->additional_numeration = NULL;
    m->symmetry = NULL;
    m->name = NULL;

    m->pixdata = NULL;

    m->n_selected = 0;
    m->selected[0] = -1;
    m->selected[1] = -1;
    m->selected[2] = -1;
    m->selected[3] = -1;
    m->n_marked = 0;
    m->marked = NULL;

    m->nbond = 0;
    m->bond = NULL;

    m->dipole[0] = 0.0;
    m->dipole[1] = 0.0;
    m->dipole[2] = 0.0;

    m->ngrids = 0;
    m->grid_type = NULL;
    m->grid_energy = NULL;
    m->grid_occ = NULL;
    m->grid_symmetry = NULL;
    m->grid_index = NULL;
    m->orbital_type = NULL;
    m->edited = NULL;
    m->orb_starting_position = 0;

    m->grid.ByteCutOff = NULL;
    m->grid.pos = NULL;
    m->grid.title = NULL;
    m->grid.titlesArr = NULL;
    m->grid.FltTitle = NULL;
    m->grid.iType = NULL;
    m->grid.dependence = NULL;
    m->grid.values = NULL;
    m->epot = NULL;

    m->geo_energy = 0.0;
    m->max_grad = 0.0;
    m->rms_grad = 0.0;

    m->nvibr = 0;
    m->freq = NULL;
    m->ir_intensity = NULL;
    m->raman_intensity = NULL;
    m->normal_mode = NULL;

    m->nvector = 0;
    m->vector1 = NULL;
    m->vector2 = NULL;
    m->radius = NULL;
    m->sharpness = NULL;
    m->vector_color = NULL;

    m->ntriangle = 0;
    m->triangle1 = NULL;
    m->triangle2 = NULL;
    m->triangle3 = NULL;
    m->triangle_color = NULL;

    m->nsphere = 0;
    m->sphere_center = NULL;
    m->sphere_radius = NULL;
    m->sphere_color = NULL;

    m->nsurf = 0;
    m->surf1 = NULL;
    m->surf2 = NULL;
    m->surf3 = NULL;
    m->surf_color = NULL;

    m->ncells = 0;
    m->cell1 = NULL;
    m->cell2 = NULL;
    m->cell3 = NULL;
    m->cell4 = NULL;
    m->cell_color = NULL;

    m->ntextboxes = 0;
    m->textboxes = NULL;
  }
  m = mp;

  /*init atom list*/
  init_atom_list();
  return mp;
}

void allocate_atoms(MOL *mp, int nnatom)
{
  int i;
  int old_natom = mp->natom;
#ifdef EBUG
  printf("allocating atoms\n"); fflush(stdout);
#endif
  mp->natom = nnatom;
  mp->xyz = (XYZ*) realloc(mp->xyz, sizeof(XYZ) * nnatom);
  mp->charge_m = (double*) realloc(mp->charge_m, sizeof(double) * nnatom);
  mp->charge_l = (double*) realloc(mp->charge_l, sizeof(double) * nnatom);
  mp->additional_numeration = (int*) realloc(mp->additional_numeration, sizeof(int) * nnatom);
  mp->symmetry = (int*) realloc(mp->symmetry, sizeof(int) * nnatom);
  mp->name = (char**) realloc(mp->name, sizeof(char*) * nnatom);
  mp->elem = (ELEM_DATA*) realloc(mp->elem, sizeof(ELEM_DATA) * nnatom);
  mp->pixdata = (PIXDATA*) realloc(mp->pixdata, sizeof(PIXDATA) * nnatom);

  for(i = old_natom; i < nnatom; i++)
  {
    mp->elem[i].vdw_rad = 0.0;
    mp->elem[i].bond_rad = 0.0;
    mp->elem[i].valency = 0;
    mp->elem[i].name = NULL;

    mp->charge_m[i] = 0.0;
    mp->charge_l[i] = 0.0;
    mp->additional_numeration[i] = 0;
    mp->symmetry[i] = 0;
    mp->name[i] = NULL;

    mp->pixdata[i].pix_index.width = 0;
    mp->pixdata[i].pix_index.height = 0;
    mp->pixdata[i].pix_index.pixels = NULL;
    mp->pixdata[i].pix_addnum.width = 0;
    mp->pixdata[i].pix_addnum.height = 0;
    mp->pixdata[i].pix_addnum.pixels = NULL;
    mp->pixdata[i].pix_symm.width = 0;
    mp->pixdata[i].pix_symm.height = 0;
    mp->pixdata[i].pix_symm.pixels = NULL;
    mp->pixdata[i].pix_name.width = 0;
    mp->pixdata[i].pix_name.height = 0;
    mp->pixdata[i].pix_name.pixels = NULL;
    mp->pixdata[i].pix_charge_m.width = 0;
    mp->pixdata[i].pix_charge_m.height = 0;
    mp->pixdata[i].pix_charge_m.pixels = NULL;
    mp->pixdata[i].pix_charge_l.width = 0;
    mp->pixdata[i].pix_charge_l.height = 0;
    mp->pixdata[i].pix_charge_l.pixels = NULL;
  }
#ifdef EBUG
  printf("allocating atoms done\n"); fflush(stdout);
#endif
}

void allocate_bonds(MOL *mp, int nnbonds)
{
  mp->nbond = nnbonds;
  mp->bond = (BOND*) realloc(mp->bond, sizeof(BOND) * (nnbonds + 1));
}

void allocate_grids(MOL *mp, int nngrid)
{
  mp->grid_type = (int*) realloc(mp->grid_type, sizeof(int) * nngrid);
  mp->grid_energy = (double*) realloc(mp->grid_energy, sizeof(double) * nngrid);
  mp->grid_occ = (double*) realloc(mp->grid_occ, sizeof(double) * nngrid);
  mp->grid_symmetry = (int*) realloc(mp->grid_symmetry, sizeof(int) * nngrid);
  mp->grid_index = (int*) realloc(mp->grid_index, sizeof(int) * nngrid);
  mp->orbital_type = (int*) realloc(mp->orbital_type, sizeof(int) * nngrid);
  mp->edited = (int*) realloc(mp->edited, sizeof(int) * nngrid);

  mp->grid.ByteCutOff = (int*) malloc(sizeof(int) * mp->grid.npoints);
  mp->grid.pos = (long*) malloc(sizeof(long) * (unsigned int) mp->grid.nBlock * mp->grid.ngrid);
/*  mp->grid.title = NULL;*/
  mp->grid.titlesArr = (char**) malloc(sizeof(char*) * nngrid);
  mp->grid.FltTitle = (int*) malloc(sizeof(int) * nngrid);
  mp->grid.iType = (char*) malloc(sizeof(char) * nngrid);
  mp->grid.dependence = (int*) malloc(sizeof(int) * nngrid);
  mp->grid.values = (double*) malloc(sizeof(double) * mp->grid.npoints);
}

void allocate_vibs(MOL *mp, int nnvibs)
{
  mp->nvibr = nnvibs;
  mp->freq = (double*) realloc(mp->freq, sizeof(double) * nnvibs);
  mp->ir_intensity = (double*) realloc(mp->ir_intensity, sizeof(double) * nnvibs);
  mp->raman_intensity = (double*) realloc(mp->raman_intensity, sizeof(double) * nnvibs);
  mp->normal_mode = (XYZ**) realloc(mp->normal_mode, sizeof(XYZ*) * nnvibs);
  mp->normal_mode[nnvibs-1] = (XYZ*) malloc(sizeof(XYZ) * mp->natom);
}

void allocate_vectors(MOL *mp, int nnvec)
{
  mp->nvector = nnvec;
  mp->vector1 = (XYZ*) realloc(mp->vector1, sizeof(XYZ) * nnvec);
  mp->vector2 = (XYZ*) realloc(mp->vector2, sizeof(XYZ) * nnvec);
  mp->radius = (double*) realloc(mp->radius, sizeof(double) * nnvec);
  mp->sharpness = (double*) realloc(mp->sharpness, sizeof(double) * nnvec);
  mp->vector_color = (color_t*) realloc(mp->vector_color, sizeof(color_t) * nnvec);
}

void allocate_triangles(MOL *mp, int nntri)
{
  mp->ntriangle = nntri;
  mp->triangle1 = (XYZ*) realloc(mp->triangle1, sizeof(XYZ) * nntri);
  mp->triangle2 = (XYZ*) realloc(mp->triangle2, sizeof(XYZ) * nntri);
  mp->triangle3 = (XYZ*) realloc(mp->triangle3, sizeof(XYZ) * nntri);
  mp->triangle_color = (color_t*) realloc(mp->triangle_color, sizeof(color_t) * nntri);
}

void allocate_spheres(MOL *mp, int nnspher)
{
  mp->nsphere = nnspher;
  mp->sphere_center = (XYZ*) realloc(mp->sphere_center, sizeof(XYZ) * nnspher);
  mp->sphere_radius = (double*) realloc(mp->sphere_radius, sizeof(double) * nnspher);
  mp->sphere_color = (color_t*) realloc(mp->sphere_color, sizeof(color_t) * nnspher);
}

void allocate_surfaces(MOL *mp, int nnsurf)
{
  mp->nsurf = nnsurf;
  mp->surf1 = (XYZ*) realloc(mp->surf1, sizeof(XYZ) * nnsurf);
  mp->surf2 = (XYZ*) realloc(mp->surf2, sizeof(XYZ) * nnsurf);
  mp->surf3 = (XYZ*) realloc(mp->surf3, sizeof(XYZ) * nnsurf);
  mp->surf_color = (color_t*) realloc(mp->surf_color, sizeof(color_t) * nnsurf);
}

void allocate_cells(MOL *mp, int nncells)
{
  mp->ncells = nncells;
  mp->cell1 = (XYZ*) realloc(mp->cell1, sizeof(XYZ) * nncells);
  mp->cell2 = (XYZ*) realloc(mp->cell2, sizeof(XYZ) * nncells);
  mp->cell3 = (XYZ*) realloc(mp->cell3, sizeof(XYZ) * nncells);
  mp->cell4 = (XYZ*) realloc(mp->cell4, sizeof(XYZ) * nncells);
  mp->cell_color = (color_t*) realloc(mp->cell_color, sizeof(color_t) * nncells);
}

void allocate_textbox(MOL *mp)
{
  mp->ntextboxes++;
  mp->textboxes = (TEXTBOX_T*) realloc(mp->textboxes, sizeof(TEXTBOX_T) * mp->ntextboxes);
  mp->textboxes[mp->ntextboxes-1].message = NULL;
  mp->textboxes[mp->ntextboxes-1].font = NULL;
  mp->textboxes[mp->ntextboxes-1].pixtext.pixels = NULL;
}

void delete_last_textbox(MOL *mp)
{
  if (!mp) return;
  if (mp->ntextboxes)
  {
    if (mp->textboxes[mp->ntextboxes-1].message) free(mp->textboxes[mp->ntextboxes-1].message);
    if (mp->textboxes[mp->ntextboxes-1].font) free(mp->textboxes[mp->ntextboxes-1].font);
    if (mp->textboxes[mp->ntextboxes-1].pixtext.pixels) free(mp->textboxes[mp->ntextboxes-1].pixtext.pixels);
    mp->ntextboxes--;
    mp->textboxes = (TEXTBOX_T*) realloc(mp->textboxes, sizeof(TEXTBOX_T) * mp->ntextboxes);
  }
}

void delete_ith_textbox(MOL *mp, int it)
{
  int i;
  if (it >= 0 && it < mp->ntextboxes)
  {
    if (mp->textboxes[it].message) free(mp->textboxes[it].message);
    if (mp->textboxes[it].font) free(mp->textboxes[it].font);
    if (mp->textboxes[it].pixtext.pixels) free(mp->textboxes[it].pixtext.pixels);
    for(i = it; i < mp->ntextboxes-1; i++)
      mp->textboxes[i+1] = mp->textboxes[i];
    mp->ntextboxes--;
    mp->textboxes = (TEXTBOX_T*) realloc(mp->textboxes, sizeof(TEXTBOX_T) * mp->ntextboxes);
  }
}

char *read_line(FILE *in)
{
#ifdef LINUX
  size_t nchar = 0;
  char *line = NULL;
  ssize_t res;

  if (feof(in))
  {
    line = (char*) malloc(sizeof(char));
    line[0] = 0;
    return line;
  }
  res = getline(&line, &nchar, in);
  if (res < 0)
  {
    free(line);
    line = (char*) malloc(sizeof(char));
    line[0] = 0;
  }
  return line;
#endif
#ifdef WINDOWS
  int i = 1;
  char c;
  char *str = NULL;
  str = (char*) malloc(sizeof(char) * i);
  str[i-1] = 0;
  do
  {
    c = fgetc(in);
    if (c == 10 || c == 0 || feof(in)) return str;
    else
    {
      str = (char*) realloc(str, sizeof(char) * ++i);
      str[i-2] = c;
      str[i-1] = 0;
    }
  }
  while(c != 10 && c != 0 && !feof(in));

  return str;
#endif
}

char get_file_exist(char* filename)
{
  struct stat file_info;
  int rc;
  rc=stat(filename, &file_info);
  if (rc == 0) return 1;
  else return 0;
}

time_t get_file_time(char* filename)
{
  struct stat file_info;
  int rc;

  rc=stat(filename, &file_info);

  return file_info.st_mtime;
}

void open_file(char *filename, char *file_description)
{
  int i;
  char *ext = strrchr(filename, '.'); /*POSSIBLE BUG!!!*/
  char *command;
  char *newfile;

#ifdef EBUG
  printf("opening file |%s| file description |%s|\n", filename, file_description); fflush(stdout);
#endif

  delete_all_orbitals_from_list();

  if (input_file_name) free(input_file_name);

  /*if (input_file_type) free(input_file_type);*/
  input_file_name = strdup(filename);

  if (!input_file_type && file_description) input_file_type = strdup(file_description);

  if (current_directory) free(current_directory);

  current_directory = strdup(filename);

#ifdef WINDOWS
  for(i = strlen(current_directory)-1; i > 0 && current_directory[i] != '\\'; i--);
#else
  for(i = strlen(current_directory)-1; i > 0 && current_directory[i] != '/'; i--);
#endif
  if (i != 0) current_directory[i] = 0;

  init_file_time = get_file_time(input_file_name);

  if (ext == NULL)
  {
    make_warning("ERROR: Can't open file: unknown filetype!");
    return;
  }

  ext += 1;

  /*19.02.2013; now determine right plugin for converting file*/
  if (file_description) /*determination based on user-selected file type*/
  {
    if (strcmp(file_description, LUSCUS_FILE_DES) == 0)
    {
/*      if (m)
        if (m->ngrids)
          delete_all_orbitals_from_list();*/

      open_gv_file(input_file_name);
      if (!file_open_success)
      {
        init_file_time = 0;
        close_file();
        /*make false m structure*/
        return;
      }
      Do_center();
      luscus_gtk_show_or_hide_widgets();
      luscus_gtk_menubar_show_or_hide_widgets();
      luscus_gtk_pop_message_from_statusbar1();
      luscus_gtk_push_message_to_statusbar1(input_file_name);
      luscus_gtk_update_3Dobject_info();
      file_open_success = 0;
      return;
    }
    else
    {
      for (i = 0; i < n_input_types; i++)
      {
        if (strcmp(input_filetypes[i].description, file_description) == 0) break;
      }
    }
  }
  else if (strcmp(ext, LUSCUS_FILE_EXT) == 0)  /* special option for LUS files !*/
  {
/*    if (m)
      if (m->ngrids)
        delete_all_orbitals_from_list();*/

    open_gv_file(input_file_name);

    if (!file_open_success)
    {
      init_file_time = 0;
      close_file();
      /*make false m structure*/
      return;
    }
    Do_center();
    luscus_gtk_show_or_hide_widgets();
    luscus_gtk_menubar_show_or_hide_widgets();
    luscus_gtk_pop_message_from_statusbar1();
    luscus_gtk_push_message_to_statusbar1(input_file_name);
    luscus_gtk_update_3Dobject_info();
    file_open_success = 0;
    return;
  }
  else /*determination based on file extension*/
  {
    for (i = 0; i < n_input_types; i++)
    {
      if (strcmp(input_filetypes[i].extension, ext) == 0) break;
    }
  }

  if (i < n_input_types && input_filetypes[i].libpath)
  {
    command = (char*) malloc(sizeof(char) * (strlen(input_filetypes[i].libpath) + strlen(input_filetypes[i].forward) + strlen(filename) + 7));
    command[0] = 0;
/*#ifdef WINDOWS
    strcat(command, "\"");
#endif*/
    strcat(command, input_filetypes[i].libpath);
#ifdef WINDOWS
    strcat(command, "\\");
#else
    strcat(command, "/");
#endif
    strcat(command, input_filetypes[i].forward);
/*#ifdef WINDOWS
    strcat(command, "\" \"");
    strcat(command, filename);
    strcat(command, "\"");
#else*/
    strcat(command, " ");
    strcat(command, filename);
/*#endif*/
    if (system(command) < 0)
    {
      free(command);  /*converting file to GV failed*/
      command = (char*) malloc(sizeof(char) * (strlen(filename) + 27));
      command[0] = 0;
      strcat(command,"ERROR: Can't process file: ");
      strcat(command, filename);
      make_warning(command);
      free(command);
      return;
    }
    free(command);
  }
  else
  {
   /* 08.04.2013: this else diseabled This prevents opening a lus. file!*/
   /* 03.01.2014: this else enabled again*/
    make_warning("ERROR: Can't open file: unknown filetype!");
    return;
  }

  newfile = make_luscus_file_name(filename);

/*  if (m)
    if (m->ngrids)
      delete_all_orbitals_from_list();*/

  open_gv_file(newfile);

  free(newfile);

  if (!file_open_success)
  {
    init_file_time = 0;
    close_file();
    /*make false m structure*/
    return;
  }

  Do_center();
  luscus_gtk_show_or_hide_widgets();
  luscus_gtk_menubar_show_or_hide_widgets();
  luscus_gtk_pop_message_from_statusbar1();
  luscus_gtk_push_message_to_statusbar1(input_file_name);
  luscus_gtk_update_3Dobject_info();
  file_open_success = 0;
}

FILE *open_current_luscus_file(void)
{
  FILE *fp;
  char *luscus_filename;
  luscus_filename = make_luscus_file_name(input_file_name);
  fp = fopen(luscus_filename, "rb");
  free(luscus_filename);
  return fp;
}

void open_gv_file(char *filename)
{
  char *msg;
  FILE *in = NULL;

  in = fopen(filename, "rb");
  if (in == NULL)
  {
    msg = (char*) malloc(sizeof(char) * (strlen(filename) + 25));
    msg[0] = 0;
    strcat(msg,"ERROR: Can't open file: ");
    strcat(msg, filename);
    make_warning(msg);
    free(msg);
    return;
  }

  remove_backup();
  parse_gv_file(in);
  append_backup();
  fclose(in);
  if (m->ngrids) do_load_grid();
  file_open_success = 1;
  luscus_gtk_update_geo_info();
}

void parse_gv_file(FILE *in)
{
  int ngeoms;
  deallocate_mol_data();
  ngeoms = read_geo_positions(in);
  mol = new_mol(ngeoms);

  read_all_sections(in);
  draw_all_pixdata();
/*  luscus_gtk_update_geo_info();*/
}

int read_geo_positions(FILE *in) /*In the case of grid; this functin should jump to the end of grid section*/
{				/*IT WILL MAKE READING MUCH FASTER!!!*/
  int i, j;
  long lastpos;
  char *line;
  int ngeo;
  int sr45 = 0;
  ngeo = 1; /*file beginning one geometry is present*/
  geo_file_position = malloc(sizeof(long));
  geo_file_position[0] = ftell(in);

  line = read_line(in);
  while(!feof(in) && strcasestr(line, "<break>") == 0)
  {
    free(line);
    line = read_line(in);
    if (feof(in))
    {
      if (sr45) ngeo--;
      break;
    }
    if (strcasestr(line, "<break>"))
    {
      if (sr45) ngeo--;
      break;
    }
    if (strcasestr(line, "<end>"))
    {
      ngeo++;
      sr45 = 1;
      geo_file_position = realloc(geo_file_position, sizeof(long) * ngeo);
      geo_file_position[ngeo-1] = ftell(in);
    }
    else sr45 = 0;
  }
  if (line) free(line);

  /*This detects the "<end>" at the end of the file*/
  lastpos = ftell(in);
  geo_file_position = realloc(geo_file_position, sizeof(long) * ngeo);

/*  rewind(in);*/ /*Do I need this?*/
  return ngeo;
}

void read_all_sections(FILE *in)
{
  int i;
  for (i = 0; i < n_geometries; i++)
  {
    should_bonding_be_determined = 1;
    read_section(in, i, mol+i);
  }
  insert_all_atoms_into_list();
}

void read_all_grids_from_all_sections(int *m_ngrids)
{
  int m_igrid = 0;
  FILE *in;
  char *line;

  if (!input_file_name)
  {
    luscus_gtk_pop_message_from_statusbar1();
/*    luscus_gtk_push_message_to_statusbar1("ERROR: Can't find file that contains grid data!");*/
    return;
  }
  in = open_current_luscus_file();
  if (in == NULL)
  {
    luscus_gtk_push_message_to_statusbar1("ERROR: Can't open file that contains grid data!");
    return;
  }

  while(!feof(in))
  {
    line = read_line(in);
    if (strcasestr(line, "<end>")) m_igrid++;
    if (strcasestr(line, "<grid>"))
      if (m_ngrids[m_igrid]) read_grid_block(in, mol + m_igrid);
    free(line);
  }

  fclose(in);
}

void read_section(FILE *in, int isec, MOL *mp)
{
  int next = 1;
  char *line;

  if (isec >= n_geometries) return;
  fseek(in, geo_file_position[isec], 0);
  read_coord_block(in, mp);
  while(next)
  {
    if (feof(in)) next = 0;
    if (next)
    {
      line = read_line(in);
           if (strcasestr(line, "<atom>")) read_additional_atom_block(in, mp);
      else if (strcasestr(line, "<energy>")) read_energy_block(in, mp);
      else if (strcasestr(line, "<rms_grad>")) read_rms_grad_block(in, mp);
      else if (strcasestr(line, "<max_grad>")) read_max_grad_block(in, mp);
      else if (strcasestr(line, "<bond>")) read_bond_block(in, mp);
      else if (strcasestr(line, "<dipole>")) read_dipole_block(in, mp);
      else if (strcasestr(line, "<vector>")) read_vector_block(in, mp);
      else if (strcasestr(line, "<triangle>")) read_triangle_block(in, mp);
      else if (strcasestr(line, "<sphere>")) read_sphere_block(in, mp);
      else if (strcasestr(line, "<surface>")) read_surface_block(in, mp);
      else if (strcasestr(line, "<cell>")) read_cell_block(in, mp);
      else if (strcasestr(line, "<textbox>")) read_textbox_block(in, mp);
      else if (strcasestr(line, "<grid>")) read_grid_block(in, mp);
      else if (strcasestr(line, "<vibration>")) read_vibration_block(in, mp);
      else if (strcasestr(line, "<editable>")) read_editable_block(in, mp); /*19.02.2013 new block!*/
      else if (strcasestr(line, "<sleep>")) read_sleep_block(in);
      else if (strcasestr(line, "<write>")) read_write_block(in);
      else if (strcasestr(line, "<end>")) next = 0;
      free(line);
    }
  }
  if (should_bonding_be_determined) determine_bonding(mp);

  set_sizes();
}

void read_coord_block(FILE *in, MOL *mp)
{
  int i;
  char *line;
  int natom;
  int nread = 0;
  char *tmp;

  line = read_line(in);
  if (line)
  {
    natom = atoi(line);
    allocate_atoms(mp, natom);
  }
  else
  {
    deallocate_mol_data();
    make_warning("ERROR: Can't read file!");
    return;
  }
  free(line);
  line = read_line(in); free(line); /*comment line*/
  for(i = 0; i < natom; i++)
  {
    line = read_line(in);
    if (line)
    {
      if (line[0])
      {
        mp->elem[i].name = get_one_word(get_ptr_ith_word(line,1));
        set_element_data(mp, i);

	tmp = get_ptr_ith_word(line, 2);
        if (tmp) mp->xyz[i][0] = atof(tmp);
	else mp->xyz[i][0] = 0.0;

        tmp = get_ptr_ith_word(line, 3);
        if (tmp) mp->xyz[i][1] = atof(tmp);
	else mp->xyz[i][1] = 0.0;

        tmp = get_ptr_ith_word(line, 4);
        if (tmp) mp->xyz[i][2] = atof(tmp);
	else mp->xyz[i][2] = 0.0;
        nread++;
      }
      free(line);
    }
  }
  if (nread != natom)
  {
    make_warning("WARNING: Broken file; some atoms are missing!\n");
    /*In the case of broken file!!!*/
    allocate_atoms(mp, nread);
  }

}

void read_additional_atom_block(FILE *in, MOL *mp)
{
  char *line;
  int next = 1;
  int iatom;

  if (mp->natom <= 0) return;

  iatom = 0;

  while(next)
  {
    line = read_line(in);
    if (feof(in)) next = 0;
    if (iatom == mp->natom) next = 0; /*Check that file doesn't contain exess of atoms*/
    if (next)
    {
      if (strcasestr(line, "</atom>")) next = 0;
      if (strcasestr(line,"name"))
      {
       	mp->name[iatom] = my_read_str_value((char*) strcasestr(line,"name"));
        if (!(mp->ishow & HAS_ATOMNAMES)) mp->ishow ^= HAS_ATOMNAMES;
      }
      if ((char*) strcasestr(line,"mulliken_charge"))
      {
       	mp->charge_m[iatom] = my_read_dble_value((char*) strcasestr(line,"mulliken_charge"));
        if (!(mp->ishow & HAS_MULLIKEN)) mp->ishow ^= HAS_MULLIKEN;
      }
      if (strcasestr(line,"loprop_charge"))
      {
       	mp->charge_l[iatom] = my_read_dble_value((char*) strcasestr(line,"loprop_charge"));
        if (!(mp->ishow & HAS_LOPROP)) mp->ishow ^= HAS_LOPROP;
      }
      if (strcasestr(line,"number"))
      {
       	mp->additional_numeration[iatom] = my_read_int_value((char*) strcasestr(line,"number"));
        if (!(mp->ishow & HAS_ATOMNUMS)) mp->ishow ^= HAS_ATOMNUMS;
      }
      if (strcasestr(line,"symmetry"))
      {
       	mp->symmetry[iatom] = my_read_int_value((char*) strcasestr(line,"symmetry"));
        if (!(mp->ishow & HAS_ATOMSYMM)) mp->ishow ^= HAS_ATOMSYMM;
      }
/*added 08.04.2013 atom colors*/
      if (strcasestr(line,"red"))
      {
       	mp->elem[iatom].color[0] = my_read_dble_value((char*) strcasestr(line,"red"));
        if (!(mp->ishow & HAS_COLOUR)) mp->ishow ^= HAS_COLOUR;
      }
      if (strcasestr(line,"green"))
      {
       	mp->elem[iatom].color[1] = my_read_dble_value((char*) strcasestr(line,"green"));
        if (!(mp->ishow & HAS_COLOUR)) mp->ishow ^= HAS_COLOUR;
      }
      if (strcasestr(line,"blue"))
      {
       	mp->elem[iatom].color[2] = my_read_dble_value((char*) strcasestr(line,"blue"));
        if (!(mp->ishow & HAS_COLOUR)) mp->ishow ^= HAS_COLOUR;
      }
      iatom++;
    }
    free(line);
  }

}

void read_energy_block(FILE *in, MOL *mp)
{
  char *line;
  line = read_line(in);
  mp->geo_energy = atof(line);
  if (!(mp->ishow & HAS_ENERGY)) mp->ishow ^= HAS_ENERGY;
  free(line);
}

void read_rms_grad_block(FILE *in, MOL *mp)
{
  char *line;
  line = read_line(in);
  mp->max_grad = atof(line);
  if (!(mp->ishow & HAS_RMS_GRAD)) mp->ishow ^= HAS_RMS_GRAD;
  free(line);
}

void read_max_grad_block(FILE *in, MOL *mp)
{
  char *line;
  line = read_line(in);
  mp->rms_grad = atof(line);
  if (!(mp->ishow & HAS_MAX_GRAD)) mp->ishow ^= HAS_MAX_GRAD;
  free(line);
}

void read_bond_block(FILE *in, MOL *mp)
{
  char *line;
  int next = 1;
  int at1, at2, bo;

  line = read_line(in);
  if (!line_is_empty(line))
    if (my_read_int_value((char*) strcasestr(line,"automatic")) == 0) should_bonding_be_determined = 0;
  free(line);

  while(next)
  {
    line = read_line(in);
    if (feof(in)) next = 0;

    if (next)
    {
      if (strcasestr(line, "</bond>")) next = 0;
      else
      {
        sscanf(line, "%d %d %d", &at1, &at2, &bo);
        if (at1 > 0 && at1 <= mp->natom)
          if (at2 > 0 && at2 <= mp->natom)
            add_bond(mp, at1 - 1, at2 - 1, bo);
      }
    }
    free(line);
  }

}

void determine_bonding(MOL *mp)
{
  int i, j;
  int nbond = 0;
  int nne1, nne2;
#ifdef EBUG
  printf("entering function: determine bonding!\n"); fflush(stdout);
#endif
  if (mp->natom > MAX_ATOMS_FOR_BONDING) return;
  /*connect atoms with single bonds*/
  for(i = 0; i < mp->natom - 1; i++)
    for(j = i + 1; j < mp->natom; j++)
    {
#ifndef NO_RESTRICT_BOND
      if (fabs(mp->xyz[i][0] - mp->xyz[j][0]) < 3.0 &&
          fabs(mp->xyz[i][1] - mp->xyz[j][1]) < 3.0 &&
          fabs(mp->xyz[i][2] - mp->xyz[j][2]) < 3.0)
      {
#endif
      /* should be filter that does not check atoms that are some distance apart!*/      
        if (find_bond(mp, i, j) == -1)
        {
          if (atom_distance(mp, i, j) < 1.15 *(mp->elem[i].bond_rad + mp->elem[j].bond_rad))
          {
            nne1 = get_number_of_bonds_on_atom(i, mp);
            nne2 = get_number_of_bonds_on_atom(j, mp);
            if (nne1 < mp->elem[i].valency && nne2 < mp->elem[j].valency)
            {
              nbond++;
              add_bond(mp, i, j, SINGLE_BOND);
            }
          }
        }
#ifndef NO_RESTRICT_BOND
      }
#endif
    }

#ifdef EBUG
  printf("bond order one = %d\n", Input_Data.bond_order_one);
#endif
  if (Input_Data.bond_order_one) return;

  /*find unsatisfied atoms and increase bond order*/
  for(i = 0; i < mp->nbond; i++)
    if (mp->bond[i].bond_type == SINGLE_BOND)
    {
      nne1 = get_number_of_bonds_on_atom(mp->bond[i].iat1, mp);
      nne2 = get_number_of_bonds_on_atom(mp->bond[i].iat2, mp);
      if (nne1 < m->elem[mp->bond[i].iat1].valency && nne2 < m->elem[mp->bond[i].iat2].valency)
        mp->bond[i].bond_type = DOUBLE_BOND;
    }

  /*search for triple bonds*/
  for(i = 0; i < mp->nbond; i++)
    if (mp->bond[i].bond_type == DOUBLE_BOND)
    {
      nne1 = get_number_of_bonds_on_atom(mp->bond[i].iat1, mp);
      nne2 = get_number_of_bonds_on_atom(mp->bond[i].iat2, mp);
      if (nne1 < m->elem[mp->bond[i].iat1].valency && nne2 < m->elem[mp->bond[i].iat2].valency)
        mp->bond[i].bond_type = TRIPLE_BOND;
    }

  /*search for conjugations*/
/*  for(i = 0; i < mp->nbond; i++)
    if (mp->bond[i].bond_type == SINGLE_BOND)
    {
      nne1 = get_number_of_bond_types_on_atom(mp, mp->bond[i].iat1, DOUBLE_BOND);
      nne2 = get_number_of_bond_types_on_atom(mp, mp->bond[i].iat2, DOUBLE_BOND);
      if (nne1 == 1 && nne2 == 1)
      {
        mp->bond[i].bond_type = S_P_BOND;
        for(j = 0; j < mp->nbond; j++)
          if (mp->bond[j].bond_type == DOUBLE_BOND)
          {
            if (mp->bond[j].iat1 == mp->bond[i].iat1 || mp->bond[j].iat2 == mp->bond[i].iat1)
              mp->bond[j].bond_type = S_P_BOND;
            if (mp->bond[j].iat1 == mp->bond[i].iat2 || mp->bond[j].iat2 == mp->bond[i].iat2)
              mp->bond[j].bond_type = S_P_BOND;
          }
      }
    }*/

#ifdef EBUG
  printf("leaving function: determine bonding!\n"); fflush(stdout);
#endif
}

int get_number_of_bonds_on_atom(int iatom, MOL *mp)
{
  int i;
  int ibond = 0;
  for(i = 0; i < mp->nbond; i++)
    if (mp->bond[i].iat1 == iatom || mp->bond[i].iat2 == iatom)
      if (mp->bond[i].bond_type == SINGLE_BOND) ibond++;
      else if (mp->bond[i].bond_type == DOUBLE_BOND) ibond+=2;
      else if (mp->bond[i].bond_type == TRIPLE_BOND) ibond+=3;

  return ibond;
}

/*void adapt_bonding(MOL *mp)
{
  float *bondor;
  float *atomval;
  int i, j, k;

  bondor = (float*) malloc(mp->nbond * sizeof(float));
  atomval = (float*) malloc(mp->natom * sizeof(float));

  /-*initialize values*-/
  for (i = 0; i < mp->nbond; i++)
  {
    switch(mp->bond[i].bond_type)
    {
      case NO_BOND: bondor[i] = 0.F;
      case SINGLE_BOND: bondor[i] = 1.F;
      case DOUBLE_BOND: bondor[i] = 2.F;
      case TRIPLE_BOND: bondor[i] = 3.F;
      case PARTIAL_BOND: bondor[i] = 0.5F;
      case S_P_BOND: bondor[i] = 1.5F;
      case LINE_BOND: bondor[i] = 0.F;
    }
    printf("bondors = %f\n", bondor[i]);
  }

  for(i = 0; i < mp->natom; i++)  /-*number of iterations = number of atoms; could be changed*-/
  {
    /-*1. calculate current valencies*-/
    for(j = 0; j < mp->natom; j++)
      for(k = 0; k < mp->nbond; k++)
        if (j == mp->bond[k].iat1 || j == mp->bond[k].iat2)
          atomval[j] += (float) bondor[k];

    /-*2. correct bond orders*-/
    for(j = 0; j < mp->natom; j++)
      for(k = 0; k < mp->nbond; k++)
      {
        if (j == mp->bond[k].iat1 || j == mp->bond[k].iat2)
          bondor[k] += ((float) mp->elem[j].valency - atomval[j]) / (float) (2* mp->elem[j].valency);
        printf("bondor[%d] = %f\n", k, bondor[k]);
      }
    putchar(10);
  }

  /-*transform corrected bond orders into the bond types*-/
  for(i = 0; i < mp->nbond; i++)
  {
    if (bondor[i] < 0.25F) mp->bond[i].bond_type = NO_BOND;
    else if (bondor[i] < 0.5F) mp->bond[i].bond_type = PARTIAL_BOND;
    else if (bondor[i] < 1.25F) mp->bond[i].bond_type = SINGLE_BOND;
    else if (bondor[i] < 1.75F) mp->bond[i].bond_type = S_P_BOND;
    else if (bondor[i] < 2.5F) mp->bond[i].bond_type = DOUBLE_BOND;
    else mp->bond[i].bond_type = TRIPLE_BOND;
  }

  free(bondor);
  free(atomval);
}*/

/*int get_number_of_bond_types_on_atom(MOL *mp, int iatom, int bond_type)
{
  int i;
  int nt = 0;
  for(i = 0; i < mp->nbond; i++)
    if (mp->bond[i].iat1 == iatom || mp->bond[i].iat2 == iatom)
      if (mp->bond[i].bond_type == bond_type) nt++;
  return nt;
}*/

void read_dipole_block(FILE *in, MOL *mp)
{
  char *line;

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->dipole[0],  &mp->dipole[1],  &mp->dipole[2]);

  if (!(mp->ishow & HAS_DIPOLE)) mp->ishow ^= HAS_DIPOLE;

  free(line);
}

void read_vector_block(FILE *in, MOL *mp)
{
  char *line;
  int nvec = mp->nvector + 1;

  allocate_vectors(mp, nvec);
  vector_set_default_values(mp, nvec-1);
  line = read_line(in);

  if (strcasestr(line,"red"))
    mp->vector_color[nvec-1][0] = my_read_dble_value((char*) strcasestr(line,"red"));
  if (strcasestr(line,"green"))
    mp->vector_color[nvec-1][1] = my_read_dble_value((char*) strcasestr(line,"green"));
  if (strcasestr(line,"blue"))
    mp->vector_color[nvec-1][2] = my_read_dble_value((char*) strcasestr(line,"blue"));
  if (strcasestr(line,"transparency"))
    mp->vector_color[nvec-1][3] = my_read_dble_value((char*) strcasestr(line,"transparency"));
  if (strcasestr(line,"radius"))
    mp->radius[nvec-1] = my_read_dble_value((char*) strcasestr(line,"radius"));
  if (strcasestr(line,"sharpness"))
    mp->sharpness[nvec-1] = my_read_dble_value((char*) strcasestr(line,"sharpness"));

  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->vector1[nvec-1][0], &mp->vector1[nvec-1][1], &mp->vector1[nvec-1][2]);
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->vector2[nvec-1][0], &mp->vector2[nvec-1][1], &mp->vector2[nvec-1][2]);
  free(line);
}

void vector_set_default_values(MOL *mp, int ivec)
{
  int i;
  if (ivec >= mp->nvector) return;
  for (i = 0; i < 4; i++)
    mp->vector_color[ivec][i] = Input_Data.extracolor[i];
  /*default vector radius and sharpness here*/
  mp->radius[ivec] = 0.05F;
  mp->sharpness[ivec] = 0.1F;
}

void read_triangle_block(FILE *in, MOL *mp)
{
  char *line;
  int ntriang = mp->ntriangle + 1;

  allocate_triangles(mp, ntriang);
  triangle_set_defeault_values(mp, ntriang-1);
  line = read_line(in);
  if (strcasestr(line,"red"))
    mp->triangle_color[ntriang-1][0] = my_read_dble_value((char*) strcasestr(line,"red"));
  if (strcasestr(line,"green"))
    mp->triangle_color[ntriang-1][1] = my_read_dble_value((char*) strcasestr(line,"green"));
  if (strcasestr(line,"blue"))
    mp->triangle_color[ntriang-1][2] = my_read_dble_value((char*) strcasestr(line,"blue"));
  if (strcasestr(line,"transparency"))
    mp->triangle_color[ntriang-1][3] = my_read_dble_value((char*) strcasestr(line,"transparency"));
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->triangle1[ntriang-1][0], &mp->triangle1[ntriang-1][1], &mp->triangle1[ntriang-1][2]);
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->triangle2[ntriang-1][0], &mp->triangle2[ntriang-1][1], &mp->triangle2[ntriang-1][2]);
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->triangle3[ntriang-1][0], &mp->triangle3[ntriang-1][1], &mp->triangle3[ntriang-1][2]);
  free(line);
}

void triangle_set_defeault_values(MOL *mp, int itriang)
{
  int i;
  if (itriang >= mp->ntriangle) return;
  for(i = 0; i < 4; i++)
    mp->triangle_color[itriang][i] = Input_Data.extracolor[i];
}

void read_sphere_block(FILE *in, MOL *mp)
{
  char *line;
  int nspheres = mp->nsphere + 1;

  allocate_spheres(mp, nspheres);
  sphere_set_defeault_values(mp, nspheres-1);

  line = read_line(in);
  if (strcasestr(line,"red"))
    mp->sphere_color[nspheres-1][0] = my_read_dble_value((char*) strcasestr(line,"red"));
  if (strcasestr(line,"green"))
    mp->sphere_color[nspheres-1][1] = my_read_dble_value((char*) strcasestr(line,"green"));
  if (strcasestr(line,"blue"))
    mp->sphere_color[nspheres-1][2] = my_read_dble_value((char*) strcasestr(line,"blue"));
  if (strcasestr(line,"transparency"))
    mp->sphere_color[nspheres-1][3] = my_read_dble_value((char*) strcasestr(line,"transparency"));
  if (strcasestr(line,"radius"))
    mp->sphere_radius[nspheres-1] = my_read_dble_value((char*) strcasestr(line,"radius"));

  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->sphere_center[nspheres-1][0], &mp->sphere_center[nspheres-1][1], &mp->sphere_center[nspheres-1][2]);
  free(line);
}

void sphere_set_defeault_values(MOL *mp, int isphere)
{
  int i;
  if (isphere >= mp->nsphere) return;
  for(i = 0; i < 4; i++)
    mp->sphere_color[isphere][i] = Input_Data.extracolor[i];
  mp->sphere_radius[isphere] = 2.0;
}

void read_surface_block(FILE *in, MOL *mp)
{
  char *line;
  int nsurf = mp->nsurf+1;

  allocate_surfaces(mp, nsurf);
  surface_set_default_values(mp, nsurf-1);
  line = read_line(in);

  if (strcasestr(line,"red"))
    mp->surf_color[nsurf-1][0] = my_read_dble_value((char*) strcasestr(line,"red"));
  if (strcasestr(line,"green"))
    mp->surf_color[nsurf-1][1] = my_read_dble_value((char*) strcasestr(line,"green"));
  if (strcasestr(line,"blue"))
    mp->surf_color[nsurf-1][2] = my_read_dble_value((char*) strcasestr(line,"blue"));
  if (strcasestr(line,"transparency"))
    mp->surf_color[nsurf-1][3] = my_read_dble_value((char*) strcasestr(line,"transparency"));


  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->surf1[nsurf-1][0], &mp->surf1[nsurf-1][1], &mp->surf1[nsurf-1][2]); 
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->surf2[nsurf-1][0], &mp->surf2[nsurf-1][1], &mp->surf2[nsurf-1][2]); 
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->surf3[nsurf-1][0], &mp->surf3[nsurf-1][1], &mp->surf3[nsurf-1][2]); 
  free(line);
}

void surface_set_default_values(MOL *mp, int isurf)
{
  int i;
  if (isurf >= mp->nsurf) return;
  for(i = 0; i < 4; i++)
    mp->surf_color[isurf][i] = Input_Data.extracolor[i];
}

void read_cell_block(FILE *in, MOL *mp)
{
  char *line;
  int next = 1;
  int i;
  int ncell;

  ncell = mp->ncells + 1;
  allocate_cells(mp, ncell);
  cell_set_default_values(mp, ncell-1);

  line = read_line(in);

  if (strcasestr(line,"red"))
    mp->cell_color[ncell-1][0] = my_read_dble_value((char*) strcasestr(line,"red"));
  if (strcasestr(line,"green"))
    mp->cell_color[ncell-1][1] = my_read_dble_value((char*) strcasestr(line,"green"));
  if (strcasestr(line,"blue"))
    mp->cell_color[ncell-1][2] = my_read_dble_value((char*) strcasestr(line,"blue"));
  if (strcasestr(line,"transparency"))
    mp->cell_color[ncell-1][3] = my_read_dble_value((char*) strcasestr(line,"transparency"));

  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->cell1[ncell-1][0], &mp->cell1[ncell-1][1], &mp->cell1[ncell-1][2]); 
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->cell2[ncell-1][0], &mp->cell2[ncell-1][1], &mp->cell2[ncell-1][2]); 
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->cell3[ncell-1][0], &mp->cell3[ncell-1][1], &mp->cell3[ncell-1][2]); 
  free(line);

  line = read_line(in);
  sscanf(line, "%lf %lf %lf", &mp->cell4[ncell-1][0], &mp->cell4[ncell-1][1], &mp->cell4[ncell-1][2]); 
  free(line);
}

void cell_set_default_values(MOL *mp, int icell)
{
  int i;
  if (icell >= mp->ncells) return;
  for(i = 0; i < 4; i++)
    mp->cell_color[icell][i] = Input_Data.extracolor[i];
}

void read_textbox_block(FILE *in, MOL *mp)
{
  char *line;
  int itb;

  allocate_textbox(mp);
  itb = mp->ntextboxes-1;
  textbox_set_dfault_values(mp, itb);

  line = read_line(in);

  if (strcasestr(line,"red"))
    mp->textboxes[itb].color[0] = my_read_dble_value((char*) strcasestr(line,"red"));
  if (strcasestr(line,"green"))
    mp->textboxes[itb].color[1] = my_read_dble_value((char*) strcasestr(line,"green"));
  if (strcasestr(line,"blue"))
    mp->textboxes[itb].color[2] = my_read_dble_value((char*) strcasestr(line,"blue"));
  if (strcasestr(line,"x"))
    mp->textboxes[itb].coord_x = my_read_int_value((char*) strcasestr(line,"x"));
  if (strcasestr(line,"y"))
    mp->textboxes[itb].coord_y = my_read_int_value((char*) strcasestr(line,"y"));

  free(line);

  line = read_line(in);
  if (line[strlen(line)-1] == 10) line[strlen(line)-1] = 0;
  mp->textboxes[itb].font = strdup(line);
  free(line);

  line = read_line(in);
  if (line[strlen(line)-1] == 10) line[strlen(line)-1] = 0;
  mp->textboxes[itb].message = strdup(line);
  free(line);
}

void textbox_set_dfault_values(MOL* mp, int itb)
{
  int i;

  if (itb >= mp->ntextboxes) return;
  for(i = 0; i < 3; i++)
    mp->textboxes[itb].color[i] = Input_Data.label_color[i];
  mp->textboxes[itb].coord_x = 1;
  mp->textboxes[itb].coord_x = 2;
  mp->textboxes[itb].font = NULL;
  mp->textboxes[itb].message = NULL;
}

void read_grid_block(FILE *in, MOL *mp)
{
  int i, j, ii;
  char *line;
  char *word;
  char *ptr;
  int next = 1;
  int i0, i1, i2;
  int ie0, ie1, ie2;
  int icount;
  int ishow, iaia;
  double pp[3];
  double myR529 = 1.0;

  long int rec0, rec1;
  char buf[512];
  
  line = read_line(in);
/*failsafe: if the line is empty; return*/
  if (strlen(line) < 2) return;

  deallocate_grids(mp);

#ifdef EBUG
  printf("Start reading grid section\n");
#endif

  if (strcasestr(line,"N_OF_MO"))
    mp->grid.nmo = my_read_int_value((char*) strcasestr(line,"N_OF_MO"));
  if (strcasestr(line,"N_of_Grids"))
    mp->ngrids = mp->grid.ngrid = my_read_int_value((char*) strcasestr(line,"N_of_Grids"));
  if (strcasestr(line,"N_of_Points"))
    mp->grid.npoints = my_read_int_value((char*) strcasestr(line,"N_of_Points"));
  if (strcasestr(line,"Block_Size"))
    mp->grid.block_size = my_read_int_value((char*) strcasestr(line,"Block_Size"));
  if (strcasestr(line,"N_Blocks"))
    mp->grid.nBlock = my_read_int_value((char*) strcasestr(line,"N_Blocks"));
  if (strcasestr(line,"Is_cutoff"))
    mp->grid.isCutOff = my_read_int_value((char*) strcasestr(line,"Is_cutoff"));
  if (strcasestr(line,"CutOff"))
    mp->grid.CutOff = my_read_dble_value((char*) strcasestr(line,"CutOff"));
  if (strcasestr(line,"N_P"))
    mp->grid.nPointsCutOff = my_read_int_value((char*) strcasestr(line,"N_P"));
  if (strcasestr(line,"BOHR")) myR529 = R529;
  free(line);

#ifdef EBUG
  printf("N_OF_MO = %d\n", mp->grid.nmo);
  printf("N_of_Grids = %d\n", mp->ngrids);
  printf("N_of_Points = %d\n", mp->grid.npoints);
  printf("Block_Size = %d\n", mp->grid.block_size);
  printf("N_Blocks = %d\n", mp->grid.nBlock);
  printf("Is_cutoff = %d\n", mp->grid.isCutOff);
  printf("CutOff = %f\n", mp->grid.CutOff);
  printf("N_P = %d\n", mp->grid.nPointsCutOff); fflush(stdout);
  printf("Read header line\n");
#endif

  allocate_grids(mp, mp->grid.ngrid);

  allocate_msrf(mp->grid.ngrid);
  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  sscanf(ptr, "%d %d %d %d %d %d %d", &(mp->grid.iniIndex[0]),
             &(mp->grid.iniIndex[1]), &(mp->grid.iniIndex[2]),
             &(mp->grid.iniIndex[3]), &(mp->grid.iniIndex[4]),
             &(mp->grid.iniIndex[5]), &(mp->grid.iniIndex[6]));
  free(line);

  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%d %d %d", &(mp->grid.npt[0]), &(mp->grid.npt[1]), &(mp->grid.npt[2]));

  free(line);

  if (mp->grid.npoints != mp->grid.npt[0] * mp->grid.npt[1] * mp->grid.npt[2])
  {
    deallocate_grids(mp);
    make_warning("ERROR: Ill-defined number of grids");
    /*read lined untill </GRID>*/
    return;
  }
  if (mp->grid.ngrid == 0)
  {
    deallocate_grids(mp);
    make_warning("ERROR: No data in grid section");
    return;
  }

  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(mp->grid.origin[0]), &(mp->grid.origin[1]), &(mp->grid.origin[2]));
  free(line);

  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(mp->grid.axisvec1[0]),
                              &(mp->grid.axisvec1[1]),
                              &(mp->grid.axisvec1[2]));
  free(line);

  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(mp->grid.axisvec2[0]),
                              &(mp->grid.axisvec2[1]),
                              &(mp->grid.axisvec2[2]));
  free(line);

  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(mp->grid.axisvec3[0]),
                              &(mp->grid.axisvec3[1]),
                              &(mp->grid.axisvec3[2]));
  free(line);

  for(i = 0; i < 3; i++)
  {
    mp->grid.origin[i] *= R529;
    mp->grid.axisvec1[i] *= R529;
    mp->grid.axisvec2[i] *= R529;
    mp->grid.axisvec3[i] *= R529;
  }

  if(mp->grid.isCutOff)
  {
    ie0=mp->grid.npt[0];
    ie1=mp->grid.npt[1];
    ie2=mp->grid.npt[2];

    icount=0;

    for(i0 = 0; i0 < ie0; i0++)
      for(i1 = 0; i1 < ie1; i1++)
        for(i2=0; i2 < ie2; i2++)
        {
          pp[0]=(mp->grid.origin[0]+mp->grid.axisvec1[0]*i0/(ie0-1)+mp->grid.axisvec2[0]*i1/(ie1-1)+mp->grid.axisvec3[0]*i2/(ie2-1));
          pp[1]=(mp->grid.origin[1]+mp->grid.axisvec1[1]*i0/(ie0-1)+mp->grid.axisvec2[1]*i1/(ie1-1)+mp->grid.axisvec3[1]*i2/(ie2-1));
          pp[2]=(mp->grid.origin[2]+mp->grid.axisvec1[2]*i0/(ie0-1)+mp->grid.axisvec2[2]*i1/(ie1-1)+mp->grid.axisvec3[2]*i2/(ie2-1));
          ishow=0;
          for(i = 0; i < mp->natom; i++)
          {
            iaia=0;
            if(fabs(mp->xyz[i][0]-pp[0]) < mp->grid.CutOff*myR529) iaia++;
            if(fabs(mp->xyz[i][1]-pp[1]) < mp->grid.CutOff*myR529) iaia++;
            if(fabs(mp->xyz[i][2]-pp[2]) < mp->grid.CutOff*myR529) iaia++;
            
            if(iaia==3) ishow=1;
          }
          if(ishow==1) mp->grid.ByteCutOff[icount]=1;
          else mp->grid.ByteCutOff[icount]=0;
          icount++;
        }
  }

  /*file positions*/
#ifdef EBUG
  printf("Read up to Axis\n");
#endif  

  ii = 0;
  for (i = 0; i < mp->ngrids; i++)
  {
    for(j = 0; j < mp->grid.nBlock-1; j++)   
    {
      mp->grid.pos[i * mp->grid.nBlock + j] = (mp->ngrids * j + i) * 8 * mp->grid.block_size;
#ifdef EBUG
      printf("orbital position igrid:%d iblock:%d %ld\n", i, j, mp->grid.pos[i * mp->grid.nBlock + j]); fflush(stdout);
#endif
    }
    mp->grid.pos[i * mp->grid.nBlock + mp->grid.nBlock-1] =
        (mp->ngrids * (mp->grid.nBlock-1)) * 8 * mp->grid.block_size + 8 * (mp->grid.npoints - (mp->grid.nBlock - 1) * mp->grid.block_size) * i;
#ifdef EBUG
    printf("orbital position igrid:%d iblock:%d %ld\n", i, mp->grid.nBlock-1, mp->grid.pos[i * mp->grid.nBlock + mp->grid.nBlock-1]); fflush(stdout);
#endif
  }

/*
  {

    {
      mp->grid.pos[i * mp->grid.nBlock + j] = 8 * i * mp->grid.npoints + 8 * j * mp->grid.block_size;
    }
  }*/

/*  mp->orb_ending_position = (long) mp->ngrids * (long) mp->grid.nBloc * (long) 8 * (long) mp->grid.block_size;*/
  /*orbital index position*/
  line = read_line(in);
  ptr = strchr(line, '=') + 1;
  mp->orb_ending_position = atol(ptr);
  free(line);

#ifdef EBUG
  printf("End-of-orbitals location\n");
#endif

  /*orbital data*/
  for(i = 0; i < mp->ngrids; i++)
  {
    line = read_line(in);

    mp->edited[i] = 0;
    if (strcasestr(line,"GridName"))
    {
      mp->grid.titlesArr[i] = my_read_str_value((char*) strcasestr(line,"GridName"));
      if (strcasestr(line, "Orbital"))
      {
        mp->grid_type[i] = ORBITAL;
        if (strcasestr(line, "sym"))
          mp->grid_symmetry[i] = my_read_int_value((char*) strcasestr(line,"sym"));
        else
          mp->grid_symmetry[i] = 0;

        if (strcasestr(line, "index"))
          mp->grid_index[i] = my_read_int_value((char*) strcasestr(line,"index"));
        else
          mp->grid_index[i] = 0; /*some internal numbering should be put here*/

        if (strcasestr(line, "Energ"))
          mp->grid_energy[i] = my_read_dble_value((char*) strcasestr(line,"energ"));
        else
          mp->grid_energy[i] = 0.0;

        if (strcasestr(line, "OCc"))
          mp->grid_occ[i] = my_read_dble_value((char*) strcasestr(line,"occ"));
        else
          mp->grid_occ[i] = 0.0;

        if (strcasestr(line, "type"))
        {
          for(ptr = (char*) strcasestr(line, "type"); ptr[0] != '='; ptr++);
          ptr++;
          for(; ptr[0] == 32; ptr++);
          switch (ptr[0])
          {
            case 'u': mp->orbital_type[i] = ORBITAL_TYPE_U; break;
            case 'f': mp->orbital_type[i] = ORBITAL_TYPE_F; break;
            case 'i': mp->orbital_type[i] = ORBITAL_TYPE_I; break;
            case '1': mp->orbital_type[i] = ORBITAL_TYPE_1; break;
            case '2': mp->orbital_type[i] = ORBITAL_TYPE_2; break;
            case '3': mp->orbital_type[i] = ORBITAL_TYPE_3; break;
            case 's': mp->orbital_type[i] = ORBITAL_TYPE_S; break;
            case 'd': mp->orbital_type[i] = ORBITAL_TYPE_D; break;
            default:  mp->orbital_type[i] = ORBITAL_TYPE_U;
          }
        }
        else          mp->orbital_type[i] = ORBITAL_TYPE_U;
      }
      else if (strcasestr(line, "Electrostatic"))
      {
        mp->ishow ^= HAS_E_POT;
        mp->grid_type[i] = E_POTEN;

        mp->orbital_type[i] = ORBITAL_TYPE_U;
        mp->grid_symmetry[i] = 0;
        mp->grid_index[i] = 0;
        mp->grid_energy[i] = 0.0;
        mp->grid_occ[i] = 0.0;
      }
      else
      {
        mp->grid_type[i] = CUSTOM;

        mp->orbital_type[i] = ORBITAL_TYPE_U;
        mp->grid_type[i] = CUSTOM;
        mp->grid_symmetry[i] = 0;
        mp->grid_index[i] = 0;
        mp->grid_energy[i] = 0.0;
        mp->grid_occ[i] = 0.0;
      }
    }
/*      else if (strcasestr(line, "Density"))     mp->grid_type[i] = DENSITY;
      else if (strcasestr(line, "E_potential")) mp->grid_type[i] = E_POTEN;
      else                                      mp->grid_type[i] = CUSTOM;*/
    free(line);
  }

/*  line = read_line(in);
  free(line);*/
#ifdef EBUG
  printf("Reading <DENSITY>\n");
#endif

  while(next)
  {
    line = read_line(in);
    if (strcasestr(line, "density")) next = 0;
    if (strcasestr(line, "/grid") || feof(in))
    {
      make_warning("Could not find <DENSITY> section in luscus file. You might be attemting to open the old version of luscus file. Blame Valera Veryazov for changing luscus files without backward compatibility! Use old version of luscus or update your molcas.");
      deallocate_grids(mp);
      m->ngrids=0;
      free(line);
      return;
    }
    free(line);
  }
  next = 1;

 /* find where <DENSITY> ends */

  mp->orb_starting_position = ftell(in);

#ifdef EBUG
  printf("Position Set!\n");

  for(i = 0; i < mp->ngrids; i++)
  {
    printf("orb #%d type = %d E = %f occ = %f symm = %d index = %d edit = %d type = %d\n", i,
       	mp->grid_type[i], mp->grid_energy[i], mp->grid_occ[i], mp->grid_symmetry[i],
       	mp->grid_index[i], mp->edited[i], mp->orbital_type[i]);
  }
#endif

/*  if(sizeof(long int) == 4) 
  { 
    if (fread(&rec0, sizeof(int), 1, in)!=1) {fprintf(stderr, "READ ERROR\n"); fclose(in); return;}
    if (fread(buf,(unsigned int)rec0,1,in)!=1) {fprintf(stderr, "READ ERROR\n"); fclose(in); return;}
    if (fread(&rec1, sizeof(int), 1, in)!=1) {fprintf(stderr, "READ ERROR\n"); fclose(in); return;}
    if (rec0!=rec1) Fmarker=1;
  }*/

 /*rewind to the end of the grid section*/

/*  fseek(in, mp->grid.pos[(mp->ngrids - 1) * (mp->grid.nBlock - 1)], SEEK_CUR);*/

  for(i = 0; i < mp->ngrids; i++)
    if (mp->grid_type[i] == E_POTEN)
      read_epoten_grid(in, mp, i);
#ifdef EBUG
  {
    double val1, val2;
    printf("jumping to: %ld\n", mp->grid.pos[1] );
    fseek(in, mp->orb_starting_position+mp->grid.pos[1], SEEK_CUR);
    fread(&val1, (size_t) sizeof(double), 1, in);
    fread(&val2, (size_t) sizeof(double), 1, in);
    printf("GRID 1 BLOCK 2 values: %f %f\n", val1, val2);
  }
#endif


  fseek(in, mp->orb_ending_position, SEEK_CUR);

#ifdef EBUG
  printf("Jump to the Orbitals!\n");
#endif

  while(next)
  {
    line = read_line(in);
    if (strcasestr(line, "/grid")) next = 0;
    free(line);
    if (feof(in)) next = 0;
  }
#ifdef EBUG
  printf("Return from read grid!\n");
#endif

  return;

}

void read_epoten_grid(FILE *in, MOL *mp, int i_grid)
{
  int i;
  int iblock;
  int np;
  double *readbuffer;
  int success;
  int itek = 0;
#ifdef EBUG
  FILE *tst;
#endif

/*
  char *line;
  char word[50];
  int next = 1;
  int i, j;
  int x, y, z;
  long staring_position;
  double bestiso =0.1;
  double g;
 */

  if (in == NULL)
  {
    make_warning("ERROR in reading Epot: Can't open file containig grid data!");
    /*register the existence of epot somewhere if possible*/
    return;
  }

  isBinary = 1;
  isPacked = 0;

  if (i_grid > m->ngrids || i_grid < 0)
  {
    make_warning("ERROR in reading Epot: Grid number out of range!");
    /*register the existence of epot somewhere if possible*/
    return;
  }

  /*allocating data for readbuffer*/
  readbuffer = (double *) calloc ((unsigned int) m->grid.block_size, sizeof (double));
  /*allocating data for epot*/
  if (m->epot == NULL)
    m->epot = (double*) malloc(m->grid.npoints * sizeof(double));

  iblock = 0;
/*  x = y = z = 0;*/
  for (i = 0; i < m->grid.npoints; i++)       /* loop thro' file */
  {
    if (i == m->grid.block_size * iblock)
    {
      fseek(in, m->orb_starting_position + m->grid.pos[i_grid * m->grid.nBlock + iblock], SEEK_SET);
      iblock++;
      np = m->grid.block_size;
      if (iblock == m->grid.nBlock)
      np = m->grid.npoints - m->grid.block_size * (m->grid.nBlock - 1);
      if(m->grid.isCutOff)
      {
        /* calculate number of uncutted data in this block */
        int ic = 0, ii;          
        for(ii = 0; ii < np; ii++)
          if(m->grid.ByteCutOff[ii+i] == 1) ic++;
            success = ReadBlock(readbuffer, ic, in);
      }
      else
      {
        success = ReadBlock(readbuffer, np*sizeof(double), in);
      }
      itek = 0;
    }
    if(m->grid.isCutOff && m->grid.ByteCutOff[i]==0)
    {
      m->epot[i]=0.0;
      itek--;
    }
    else
      m->epot[i] = readbuffer[itek];
    itek++;
    /* For first point, set the min&max values;
       for subsequent points, conditionally set. */
  }  /* endof loop */

  free(readbuffer);

#ifdef EBUG
  tst = fopen("epot.txt", "w");

  for(i = 0; i < m->grid.npoints; i++)
  {
    fprintf(tst, "%d %f\n", i, m->epot[i]);
  }

  fclose(tst);
#endif

}

void read_vibration_block(FILE *in, MOL *mp)
{
  int i;
  char *line;
  int next = 1;
  int nvibs = mp->nvibr + 1;
  allocate_vibs(mp, nvibs);
  line = read_line(in);

  mp->freq[nvibs-1] = 0.0;
  mp->ir_intensity[nvibs-1] = 0.0;
  mp->raman_intensity[nvibs-1] = 0.0;

  if (strcasestr(line,"freq"))
    mp->freq[nvibs-1] = my_read_dble_value((char*) strcasestr(line,"freq"));
  if (strcasestr(line,"ir_int"))
  {
    mp->ir_intensity[nvibs-1] = my_read_dble_value((char*) strcasestr(line,"ir_int"));
    if (!(mp->ishow & HAS_IRINT)) mp->ishow ^= HAS_IRINT;
  }
  if (strcasestr(line,"raman_int"))
  {
    mp->raman_intensity[nvibs-1] = my_read_dble_value((char*) strcasestr(line,"raman_int"));
    if (!(mp->ishow & HAS_RAMAN)) mp->ishow ^= HAS_RAMAN;
  }

  free(line);

  for(i = 0; i < mp->natom; i++)
  {
    line = read_line(in);
    sscanf(line, "%lf %lf %lf", &mp->normal_mode[nvibs-1][i][0], &mp->normal_mode[nvibs-1][i][1], &mp->normal_mode[nvibs-1][i][2]);
    free(line);
  }
}

void read_editable_block(FILE *in, MOL *mp)
{
  char *line;
  line = read_line(in);
  if (strcasestr(line,"yes")) mp->editable=1;
  else if (strcasestr(line,"no")) mp->editable=0;
  free(line);
}

void read_sleep_block(FILE *in)
{
  char *line;
  int seconds;

  line = read_line(in);
  sscanf(line, "%d", &seconds);
  free(line);
/*  redraw();
  do_sleep(seconds);
  redraw();*/ /*BUG reading file should be redesigned if this option should take effect!*/
}

void read_write_block(FILE *in)
{
  char *line = NULL;
  do
  {
    line = read_line(in);
    if (line == NULL) break;
    if (line[strlen(line)-1] == 10) line[strlen(line)-1] = 0;
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2(line);
    free(line);
  } while(!strcasestr(line,"write") && feof(in));
}

void close_file(void)
{
  delete_all_orbitals_from_list();
  deallocate_mol_data();
  mol = new_mol(1); /*Allocate one geometry; user can draw new molecule!*/
  luscus_gtk_show_or_hide_widgets();
  free(input_file_name);
  input_file_name = strdup("new_file.gv");
  remove_backup();
  init_file_time = 0;
}

void check_if_file_changed(void)
{
  char tmppath[512];
  if (!init_file_time) return;
  checked_file_time = get_file_time(input_file_name);
  if (init_file_time != checked_file_time)
  {
    strncpy(tmppath, input_file_name, 511);
    if (input_file_type) open_file(tmppath, input_file_type);
    else open_file(tmppath, NULL);
    rerender_3d();
    redraw();
  }
}

char *get_input_filename(void)
{
  if (input_file_name) return input_file_name;
  else return "newfile.lus";
}

char *get_current_directory(void)
{
  return current_directory;
}

char *nextfname(void)
{
  int n;
  int i;
  char *new_filename;

  if (input_file_name == NULL) return NULL;

  n = strlen(input_file_name);
  for (n = strlen(input_file_name); n > 0 && input_file_name[n] != '.'; n--);
  if (n == 0) n = strlen(input_file_name);

  new_filename = (char*) malloc(sizeof(char) * (n+5));
  new_filename[0] = 0;
  for(i = 0; i < n; i++) new_filename[i] = input_file_name[i];
  new_filename[n] = '_';
  new_filename[n+1] = 'n';
  new_filename[n+2] = 'e';
  new_filename[n+3] = 'w';
  new_filename[n+4] = 0;
  return new_filename;
}

char *make_luscus_file_name(char* original_filename)
{
  int i;
  int extlen;
  char *ext;
  char *luscus_filename;

  for(i = strlen(original_filename); i > 0 && original_filename[i] != '.'; i--, ext = original_filename+i);
  if (strlen(original_filename) - i > MAX_EXT_LENGTH) ext = NULL; /*extensions greater than MAX_EXT_LENGTH do not make sense!*/

  if (ext == NULL)
  {
    luscus_filename = (char*) malloc(sizeof(char) * (strlen(original_filename) + strlen(input_filetypes[n_input_types-1].extension) + 2));
    extlen = 0;

    luscus_filename[0] = 0;
    for(i = 0; i <= strlen(original_filename) - extlen; i++) luscus_filename[i] = original_filename[i];
    luscus_filename[i] = 0;
    strcat(luscus_filename, ".");
    strcat(luscus_filename, input_filetypes[0].extension);
  }
  else
  {
    luscus_filename = (char*) malloc(sizeof(char) * (strlen(original_filename) - strlen(ext) + strlen(input_filetypes[n_input_types-1].extension) + 2));
    extlen = strlen(ext);

    luscus_filename[0] = 0;
    for(i = 0; i <= strlen(original_filename) - extlen; i++) luscus_filename[i] = original_filename[i];
    luscus_filename[i] = 0;
    strcat(luscus_filename, input_filetypes[0].extension);
  }

  return luscus_filename;
}

void get_new_section(void)
{
  m = mol + igeo;

  /*populate atom-list with new data*/
  deallocate_atom_list();
  insert_all_atoms_into_list();
}

void set_element_data(MOL *mp, int iat)
{
  int i;
  for(i = 0; i < number_of_elements; i++)
    if (strcasecmp(mp->elem[iat].name, e[i].name) == 0) break;

  mp->elem[iat].color[0] = e[i].color[0];
  mp->elem[iat].color[1] = e[i].color[1];
  mp->elem[iat].color[2] = e[i].color[2];
  mp->elem[iat].vdw_rad = e[i].vdw_rad;
  mp->elem[iat].bond_rad = e[i].bond_rad;
  mp->elem[iat].valency = e[i].valency;
  mp->elem[iat].periodic_pos_x = e[i].periodic_pos_x;
  mp->elem[iat].periodic_pos_y = e[i].periodic_pos_y;
}

double atom_distance(MOL *mp, int iat1, int iat2)
{
  if (m == NULL) return 0;
  return sqrt((mp->xyz[iat1][0] - mp->xyz[iat2][0]) * (mp->xyz[iat1][0] - mp->xyz[iat2][0]) +
              (mp->xyz[iat1][1] - mp->xyz[iat2][1]) * (mp->xyz[iat1][1] - mp->xyz[iat2][1]) + 
              (mp->xyz[iat1][2] - mp->xyz[iat2][2]) * (mp->xyz[iat1][2] - mp->xyz[iat2][2]));
}

void add_bond(MOL *mp, int iat1, int iat2, int bond_order)
{
  int i;
  int nbond;
  for(i = 0; i < mp->nbond; i++)
    if ((mp->bond[i].iat1 == iat1 && mp->bond[i].iat2 == iat2) ||
        (mp->bond[i].iat1 == iat2 && mp->bond[i].iat2 == iat1))
    {
      mp->bond[i].bond_type = bond_order;
      return;
    }

  /*bond is not recorded yet; make new record*/
  nbond = mp->nbond+1;
  allocate_bonds(mp, nbond);
  mp->bond[mp->nbond-1].iat1 = iat1;
  mp->bond[mp->nbond-1].iat2 = iat2;
  mp->bond[mp->nbond-1].bond_type = bond_order;
}

void do_load_grid(void)
{
  read_grid_from_file();
  make_surfaces();
}

void read_grid_from_file(void)
{
  FILE *in;
  char *line;
  char *luscus_filename;
  char word[50];
  int next = 1;
  int i, j;
  int np;
  int x, y, z;
  int itek = 0;
  int iblock;
  int success;
  long staring_position;
  double *readbuffer;
  double bestiso =0.1;
  double g;

  if (!input_file_name)
  {
    make_warning("ERROR: File containing grid data is undefined!");
    deallocate_grids(m);
    return;
  }

  luscus_filename = make_luscus_file_name(input_file_name);

  in = fopen(luscus_filename, "rb");
  free(luscus_filename);
  if (in == NULL)
  {
    make_warning("ERROR: Can't open file containig grid data!");
    deallocate_grids(m);
    return;
  }

  isBinary = 1;
  isPacked = 0;

  if (iorb > m->ngrids || iorb < 0)
  {
    make_warning("ERROR: Grid number out of range!");
    deallocate_grids(m);
    return;
  }

  /*allocating data for readbuffer*/

  readbuffer = (double *) calloc ((unsigned int) m->grid.block_size, sizeof (double));

  iblock = 0;
  x = y = z = 0;
  for (i = 0; i < m->grid.npoints; i++, x++)       /* loop thro' file */
  {
    if (i == m->grid.block_size * iblock)
    {
/*#ifdef EBUG
      printf("seeking position %ld\n", m->orb_starting_position + m->grid.pos[iorb * m->grid.nBlock + iblock]);
#endif*/
      fseek(in, m->orb_starting_position + m->grid.pos[iorb * m->grid.nBlock + iblock], SEEK_SET);
      iblock++;
      np = m->grid.block_size;
      if (iblock == m->grid.nBlock)
      np = m->grid.npoints - m->grid.block_size * (m->grid.nBlock - 1);
      if(m->grid.isCutOff)
      {
        /* calculate number of uncutted data in this block */
        int ic = 0, ii;          
        for(ii = 0; ii < np; ii++)
          if(m->grid.ByteCutOff[ii+i] == 1) ic++;
            success = ReadBlock(readbuffer, ic, in);
      }
      else
      {
        success = ReadBlock(readbuffer, np*8, in);
#ifdef EBUG
        printf("NEW POSITION: %ld\n", ftell(in));
        printf("address #%d val0 = %f val1 = %f\n", i, readbuffer[0], readbuffer[1]);
#endif
      }
      itek = 0;
    }
    if (x == m->grid.npt[0])
    {
      x = 0;
      y++;
      if (y == m->grid.npt[1])
      {
        y = 0;
        z++;
      }
    }
    if(m->grid.isCutOff && m->grid.ByteCutOff[i]==0)
    {
      m->grid.values[i]=0.0;
      itek--;
    }
    else
      m->grid.values[i] = readbuffer[itek];
    itek++;
    /* For first point, set the min&max values;
       for subsequent points, conditionally set. */
    if (i == 0)
    {
      m->grid.minval = 10000.0;
      m->grid.maxval = -10000.0;
    }
    else
    {
      /* let's cut off huge values */
      if (m->grid.values[i] > 2) m->grid.values[i] = 2;
      m->grid.minval = MIN (m->grid.minval, m->grid.values[i]);
      m->grid.maxval = MAX (m->grid.maxval, m->grid.values[i]);
    }
  }  /* endof loop */

  {
    int iguess;
    double cutoff[21];
    double gguess[21];
    double mysmall;
    double mybig; 
    double searchval;
    int isd=0;
    int ptotal;
    if (m->grid_type[iorb] == DENSITY) isd=1;
/*    if (strcasestr(m->grid.title,"Densit")) isd=1; */
    g = MAX(fabs(m->grid.maxval), fabs(m->grid.minval)) / 20;

    ptotal = m->grid.npt[0] * m->grid.npt[1] * m->grid.npt[2];
    for (iguess = 0; iguess < 21; iguess++)
    {
      gguess[iguess] = 0;
      cutoff[iguess] = (iguess-.001)*g;
      if(iguess==0) cutoff[0]=0;
      for (i = 0; i < ptotal; i++)
        if(fabs(m->grid.values[i])>cutoff[iguess])
        {
          if(isd) gguess[iguess]+=m->grid.values[i];
          else gguess[iguess]+=m->grid.values[i]*m->grid.values[i];
        }
    }

    for(i = 1; i < 21; i++)
    { 
      searchval=gguess[0]*(i/20.);
      for (iguess = 0; iguess < 20; iguess++)
      {
        if(gguess[iguess]>=searchval && searchval>gguess[iguess+1])     
        {
          mybig = gguess[iguess];
          mysmall = gguess[iguess+1];
          if(fabs(mybig-mysmall) < 1e-10) bestiso = .5*(cutoff[iguess]+cutoff[iguess+1]) ;
          else bestiso = ((searchval-mysmall)/(mybig-mysmall)) * (cutoff[iguess+1]-cutoff[iguess]) + cutoff[iguess];
         }
       } 
       if(bestiso > m->grid.maxval) bestiso = -bestiso ;
       m->grid.isod[i]=bestiso;
    }
  }

  m->grid.guess = m->grid.isod[14];
  sprintf(word, "%7.0le", m->grid.guess);
  sscanf(word, "%lf", &m->grid.guess);
      /* close the file and fly away */

  free(readbuffer);
  fclose(in);

#ifdef EBUG
  in = fopen("test.txt", "w");

  for(i = 0; i < m->grid.npoints; i++)
  {
    fprintf(in, "%d %f\n", i, m->grid.values[i]);
  }

  fclose(in);
#endif
}

/***
 Read or skip block. When reading, unpack values if necessary.

 In: fd  -- file descriptor
     buf -- pointer to buffer to read into
     np  -- expected number of values in block (buffer length, in doubles)
 Returns: success flag
 Notes:
  1) If pointer is NULL, does fake read (i.e. skips the block)
  2) After successfull reading/skipping the file pointer is positioned
     just after the block (presumably at the beginning of the next one)
***/ 
static int ReadBlock(double* buf, int np, FILE* fp)
{
  int j;
  char line[MAXCHAR];
  struct packing_params_t pk_prm;
  byte_t* packed_data;

  if (fread(buf, (size_t) np, 1, fp) !=1) return -1;

  return np;
}   

#ifdef NOCOMPILE
/** Read block/record from unformatted Fortran file. 
 Returns:
   -1  : on failure
   number of bytes
   read into buffer : on success
   
   If buf==NULL, fakes read.
***/
int read_bin_block(void *buf, int bufsz, FILE * fp)
{
  int rec_len, rec_len1, rec_8;
 
  if (fread(&rec_len, sizeof(rec_len), 1, fp)!=1) return -1;
 
  if(Fmarker==1)
    if (fread(&rec_8, sizeof(rec_8), 1, fp)!=1) return -1;
   
  if (buf)
  {
    if (rec_len>bufsz) return -1;
    if (fread(buf,(unsigned int)rec_len,1,fp)!=1) return -1;
  }
  else
    fseek(fp,rec_len,SEEK_CUR);
 
  if (fread(&rec_len1, sizeof(rec_len1), 1, fp)!=1) return -1;

  if(Fmarker==1)
    if (fread(&rec_8, sizeof(rec_8), 1, fp)!=1) return -1;
 
  if (rec_len1!=rec_len)
  {
    printf("ERROR reading BLOCK!\n"); fflush(stdout);
    return -1;
  }
  return rec_len;
}


/***
  Parse header, fill a structure with packing params.
***/
static int getPackingParams(const char* line, struct packing_params_t* prm)
{
  int flags;
  int i;

  if (sscanf(line,"BHeader= %d %lf %lf %lf %lf %d %d %d",
             &flags,                 /* not used yet */
             &(prm->xlimit[0]),      /* xmin */ 
             &(prm->xlimit[1]),      /* xleft */ 
             &(prm->xlimit[2]),      /* xright */ 
             &(prm->xlimit[3]),      /* xmax */ 
             &(prm->ydelta[0]),      /* dy1 */ 
             &(prm->ydelta[1]),      /* dy2 */ 
             &(prm->ydelta[2])       /* dy3 */)!=8) return 0;

  prm->ileft=0;  prm->iright=2;
  if (prm->ydelta[0]==0)
  {
    prm->ileft++;
    prm->xdelta[0]=0.;
  }
  if (prm->ydelta[2]==0)
  { 
    prm->iright--;
    prm->xdelta[2]=0.;
  }

  prm->ylimit[0]=0;
  prm->ylimit[1]=prm->ydelta[0];
  prm->ylimit[2]=prm->ylimit[1] + prm->ydelta[1];
  prm->ylimit[3]=prm->ylimit[2] + prm->ydelta[2];
  prm->range=prm->ylimit[3];
 
  prm->nbytes=(prm->range==128)?1:2;
 
  for (i=prm->ileft; i<=prm->iright; i++)
    prm->xdelta[i]=prm->xlimit[i+1] - prm->xlimit[i];

  fprintf(stderr,"LIMITS:\n"
          "range=%d\nxlimit=%f %f %f %f\nxdelta=%f %f %f\n"
          "ylimit=%d %d %d %d\nydelta=%d %d %d\n",
          prm->range,
          prm->xlimit[0],prm->xlimit[1],prm->xlimit[2],prm->xlimit[3],
          prm->xdelta[0],prm->xdelta[1],prm->xdelta[2],
          prm->ylimit[0],prm->ylimit[1],prm->ylimit[2],prm->ylimit[3],
          prm->ydelta[0],prm->ydelta[1],prm->ydelta[2]);
  return 1;
}
#endif

int abfgets (int isBinary, char *str, int max, FILE * fp)
{
  int i = 1;
  if (isBinary) i = read_bin_stream (str, max, fp);
  else if (fgets (str, max, fp) == NULL) i = 0;
  return i;
}

/***
 Unpack data either from array of bytes or from integer value.
 In: dest -- array to unpack into
     packed_data -- packed values, represented as 1-byte or MSB-first 2-bytes
     packed_val  -- alternatively, single integer to unpack
     ndata  -- length (number of values, not size in bytes!) of packed_data,
               or, the same, length (in doubles) of the dest buffer
     prm -- packing parameters

***/
static void unpackData(double* dest,
                       byte_t* packed_data, unsigned int packed_val,
                       int ndata,
                       const struct packing_params_t* prm)
{
  if (packed_data==0) { ndata=1; } /* we're called to unpack single value */
 
  while (ndata--)
  {
    if (packed_data)
    {  /* we're called to unpack array */
      packed_val=*(packed_data++);
      if (prm->nbytes==2)
      { /* 2 bytes, MSB first!!! */
        packed_val<<=8;
        packed_val += *(packed_data++);
      }
    } /* else use supplied single value of packed_val */
   
   /** Unpacking begins. **/
   /** Inlining -- for effectiveness. Separate block -- for clarity **/
    {
      int i;
      double x=0.0;
      int y=packed_val; /* for brevity */
      int minus=0;

      y-=prm->range; /* BEWARE: reperesentable int range and sign! */
      if (y<0) { y=-y; minus=1; }

      if (y<prm->ylimit[0])  { x=prm->xlimit[0]; goto x_done; }
      if (y>=prm->ylimit[3]) { x=prm->xlimit[3]; goto x_done; }

      for (i=prm->ileft; i<=prm->iright; i++)
      {
        if (y>=prm->ylimit[i] && y<prm->ylimit[i+1])
       	{
          x=((double)(y - prm->ylimit[i]))/ prm->ydelta[i] * prm->xdelta[i] + prm->xlimit[i];
          goto x_done;
        }
      }
      x_done:
      *(dest)=(minus?-x:x); 
      if( *(dest) > prm->xlimit[3]) *(dest)=0;
      *(dest++);
    }
     /** Unpacking ends. **/
  }  /* cycle through array ends */
}

/** Read string from unformatted Fortran file
 Returns:
   0  : on failure
   1  : on success
***/
int read_bin_stream(char *str, int max, FILE * fp)
{
  int j;

  j = fread(str, (size_t) max-1, 1, fp);

/*  j=read_bin_block(str,max-1,fp);*/
  if (j==-1) return 0;

  str[j]=0;
  return 1;
}

void current_grid_limits(double* min, double* max, double* guess)
{
  if (m->ngrids)
  {
    *min = m->grid.minval;
    *max = m->grid.maxval;
    *guess = m->grid.guess;
  }
}

int current_grid_data(double **data)
{
  if (m->ngrids)
  {
    *data = m->grid.values;
    return 1;
  }
  return 0;
}

int get_fmarker(void)
{
  return Fmarker;
}

