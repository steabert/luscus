/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*Fragments should be read with same functions as all other molecules!!!*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#ifdef WINDOWS
#include<strings.h>
#endif
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

ELEM_DATA *e;
int number_of_elements;

/*Only coordinate and bond sections of fragments are read*/
/*TODO: read all sections for fragments*/

void initialize_fragment_data(MOL *f)
{
  f->natom = 0;
  f->xyz = NULL;

  f->elem = NULL;
  f->charge_m = NULL;
  f->charge_l = NULL;
  f->additional_numeration = NULL;
  f->symmetry = NULL;
  f->name = NULL;

  /*selected/marked not used*/

  f->nbond = 0;
  f->bond = NULL;

  f->dipole[0] = 0.F;
  f->dipole[1] = 0.F;
  f->dipole[2] = 0.F;

  f->ngrids = 0;
  f->grid_type = NULL;
  f->grid_energy = NULL;
  f->grid_occ = NULL;
  f->grid_symmetry = NULL;
  f->grid_index = NULL;
  f->orbital_type = NULL;
  f->edited = NULL;
  f->orb_starting_position = 0;

  f->grid.ByteCutOff = NULL;
  f->grid.pos = NULL;
  f->grid.title = NULL;
  f->grid.titlesArr = NULL;
  f->grid.FltTitle = NULL;
  f->grid.iType = NULL;
  f->grid.dependence = NULL;
  f->grid.values = NULL;
  f->epot = NULL;

  f->geo_energy = 0.F;
  f->max_grad = 0.F;
  f->rms_grad = 0.F;

  f->nvibr = 0;
  f->freq = NULL;
  f->ir_intensity = NULL;
  f->raman_intensity = NULL;
  f->normal_mode = NULL;

  f->nvector = 0;
  f->vector1 = NULL;
  f->vector2 = NULL;
  f->radius = NULL;
  f->sharpness = NULL;
  f->vector_color = NULL;

  f->ntriangle = 0;
  f->triangle1 = NULL;
  f->triangle2 = NULL;
  f->triangle3 = NULL;
  f->triangle_color = NULL;

  f->nsphere = 0;
  f->sphere_center = NULL;
  f->sphere_radius = NULL;
  f->sphere_color = NULL;

  f->nsurf = 0;
  f->surf1 = NULL;
  f->surf2 = NULL;
  f->surf3 = NULL;
  f->surf_color = NULL;

  f->ncells = 0;
  f->cell1 = NULL;
  f->cell2 = NULL;
  f->cell3 = NULL;
  f->cell4 = NULL;
  f->cell_color = NULL;
}

void free_fragment_data(MOL f)
{
  int i;
  if (f.xyz) free(f.xyz);
  if (f.elem)
  {
    for(i = 0; i < f.natom; i++)
      if (f.elem[i].name)
        free(f.elem[i].name);
    free(f.elem);
  }
  if (f.bond) free(f.bond);
  if (f.charge_m) free(f.charge_m);
  if (f.charge_l) free(f.charge_l);
  if (f.additional_numeration) free(f.additional_numeration);
  if (f.symmetry) free(f.symmetry);
  if (f.name)
  {
    for(i = 0; i < f.natom; i++)
      if (f.name[i])
        free(f.name[i]);
    free(f.name);
  }

  f.natom = 0;
  f.nbond = 0;
}

void load_fragment_from_file(char *path, MOL *f)
{
  FILE *ffc;
  int i;
  int natom;
  char *line;
  int next = 1;

  ffc = fopen(path, "r");
/*coord block*/
  line = read_line(ffc);
  if (line)
  {
    natom = atoi(line);
    /*allocate atoms*/
    f->natom = natom;
    f->xyz = (XYZ*) malloc(sizeof(XYZ) * natom);
    f->elem = (ELEM_DATA*) malloc(sizeof(ELEM_DATA) * natom);
  }
  else
  {
    make_warning("ERROR: Can't load fragment!");
    return;
  }
  free(line);
  line = read_line(ffc); free(line); /*comment line*/
  /*XYZ coordinates*/
  for(i = 0; i < natom; i++)
  {
    line = read_line(ffc);
    f->elem[i].name = get_one_word(get_ptr_ith_word(line,1));
    fragment_set_element_data(f, i);
/*    f->xyz[i][0] = atof(get_ptr_ith_word(line, 2));
    f->xyz[i][1] = atof(get_ptr_ith_word(line, 3));
    f->xyz[i][2] = atof(get_ptr_ith_word(line, 4));*/
    sscanf(get_ptr_ith_word(line, 2), "%lf", f->xyz[i]);
    sscanf(get_ptr_ith_word(line, 3), "%lf", f->xyz[i]+1);
    sscanf(get_ptr_ith_word(line, 4), "%lf", f->xyz[i]+2);/*corrected 2014.01.24*/
    free(line);
  }

  while(next)
  {
    if (feof(ffc)) next = 0;
    if (next)
    {
      line = read_line(ffc);
           if (strcasestr(line, "<bond>")) read_fragment_bond_block(f, ffc);
      else if (strcasestr(line, "<end>")) next = 0;
      free(line);
    }
  }

  fclose(ffc);
}

MOL load_fragment(int ifrag)
{
#ifdef EBUG
  int i;
#endif

  char *ptr;
  MOL fragment;

  initialize_fragment_data(&fragment);

  ptr = malloc(sizeof(char) * (strlen(rcdir) + strlen(frag[ifrag].coord_file_name) + 2));
  ptr[0] = 0;

  strcat(ptr, rcdir);
#ifdef WINDOWS
  strcat(ptr, "\\");
#else
  strcat(ptr, "/");
#endif
  strcat(ptr, frag[ifrag].coord_file_name);

  load_fragment_from_file(ptr, &fragment);

  free(ptr);

#ifdef EBUG
  printf("FRAGMENT: atoms\n");
  for(i = 0; i < fragment.natom; i++)
  {
    printf(" |%s| %f %f %f\n", fragment.elem[i].name, fragment.xyz[i][0], fragment.xyz[i][1], fragment.xyz[i][2]);
  }
  printf("FRAGMENT: bonds\n");
  for(i = 0; i < fragment.nbond; i++)
  {
    printf("%d %d order = %d\n", fragment.bond[i].iat1, fragment.bond[i].iat2, fragment.bond[i].bond_type);
  }
#endif

  return fragment;
}

void fragment_set_element_data(MOL *f, int iat)
{
  int i;
  for(i = 0; i < number_of_elements; i++)
    if (strcasecmp(f->elem[iat].name, e[i].name) == 0) break;

  f->elem[iat].color[0] = e[i].color[0];
  f->elem[iat].color[1] = e[i].color[1];
  f->elem[iat].color[2] = e[i].color[2];
  f->elem[iat].vdw_rad = e[i].vdw_rad;
  f->elem[iat].bond_rad = e[i].bond_rad;
  f->elem[iat].valency = e[i].valency;
  f->elem[iat].periodic_pos_x = e[i].periodic_pos_x;
  f->elem[iat].periodic_pos_y = e[i].periodic_pos_y;
}

void read_fragment_bond_block(MOL* f, FILE *ffc)
{
  char *line;
  int next = 1;
  int at1, at2, bo;

  line = read_line(ffc);
  
 /*read keywords from line*/

  free(line);
  while(next)
  {
    line = read_line(ffc);
    if (feof(ffc)) next = 0;

    if (next)
    {
      if (strcasestr(line, "</bond>")) next = 0;
      else
      {
        sscanf(line, "%d %d %d", &at1, &at2, &bo);
        fragment_add_bond(f, at1 - 1, at2 - 1, bo);
      }
    }
    free(line);
  }
}

void fragment_add_bond(MOL* f, int at1, int at2, int bo)
{
  int i;

  /*search if atoms are allready bonded*/
  for(i = 0; i < f->nbond; i++)
  {
    if ((f->bond[i].iat1 == at1 && f->bond[i].iat2 == at2) ||
       	(f->bond[i].iat1 == at2 && f->bond[i].iat2 == at1))
    {
      f->bond[i].bond_type = bo;
      return;
    }
  }
  /*atoms are not bonded; insert new bond*/
  f->nbond++;
  f->bond = (BOND*) realloc(f->bond, sizeof(BOND) * f->nbond);
  f->bond[f->nbond - 1].iat1 = at1;
  f->bond[f->nbond - 1].iat2 = at2;
  f->bond[f->nbond - 1].bond_type = bo;
}

