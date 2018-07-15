/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"gv_atoms.h"

char *rcdir;
ELEM_DATA *e;
DEFAULT_ELEM def_elem;
int number_of_elements;

void load_default_element_data(void)
{
  int i;
  int default_atoms[]={1,9,10,11,12,19,20};
  e = e_default;

  def_elem.n_default_elem=7;
  def_elem.default_elem = (int*) malloc(7*sizeof(int));
  for(i = 0; i < 7; i++) def_elem.default_elem[i] = default_atoms[i];
}

void Write_Atom_data(char *path)
{
  int i = 0, j;
  int default_atoms[]={1,9,10,11,12,19,20};
  char next;
  FILE *ap;

  ap = fopen(path, "w");
  printf("writing atom data in |%s|\n", path); fflush(stdout);
  fprintf(ap, "%s", "# Configuration file with atom settings\n"
                    "# Atom Name, color(3), Atom size, Bond size, Max valency, periodic system coord x, periodic system coord y, default atoms\n"
                    "# Atom name can NOT be altered!\n#\n#\n");
  while(e[i].name)
  {
    fprintf(ap, "\"%s\" %7.3f %7.3f %7.3f %7.3f %7.3f %4d %4d %4d", e[i].name, e[i].color[0], e[i].color[1], e[i].color[2], e[i].vdw_rad, e[i].bond_rad, e[i].valency, e[i].periodic_pos_x, e[i].periodic_pos_y);
    
    next = 0;
    for(j = 0; j < 7; j++)
      if (i == default_atoms[j])
      {
        next = 1;
        break;
      }
    fprintf(ap, " %4d", next);
    if (e[i++].name) fputc(10, ap);
  }
  number_of_elements = i-1;
  fclose(ap);
}

void load_atom_data(void)
{
  int i;
  char *path;
  FILE *ap;
  char *line;
  int iatom = 0;
  int isdefault;
  e = NULL;
  def_elem.n_default_elem = 0;
  def_elem.default_elem = NULL;

  path = malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_ATOM) + 2));
  path[0] = 0;
  strcat(path, rcdir);
#ifdef WINDOWS
  strcat(path, "\\");
#else
  strcat(path, "/");
#endif
  strcat(path, RC_ATOM);

  ap = fopen(path, "r");

  if (ap == NULL)
  {
    load_default_element_data();
    Write_Atom_data(path);
    free(path);
    return;
  }
  free(path);

  while(!feof(ap))
  {
    line = read_line(ap);
#ifdef EBUG
    printf("LOAD ATOM DATA |%s|\n", line);
#endif
    if (line)
    {
      if (line[0])
      {
        if (line[0] != '#' && strlen(line) > 5)
        {
          e = (ELEM_DATA*) realloc(e, sizeof(ELEM_DATA) * (++iatom));
          e[iatom-1].name = (char*) malloc(sizeof(char) * 3);
          e[iatom-1].name[0] = line[1];
          if (line[2] == 34)
            e[iatom-1].name[1] = 0;
          else
          {
            e[iatom-1].name[1] = line[2];
            e[iatom-1].name[2] = 0;
          }

          sscanf(line+5, "%f %f %f %f %f %d %d %d %d", &e[iatom-1].color[0], &e[iatom-1].color[1], &e[iatom-1].color[2], &e[iatom-1].vdw_rad, &e[iatom-1].bond_rad, &e[iatom-1].valency, &e[iatom-1].periodic_pos_x, &e[iatom-1].periodic_pos_y, &isdefault);
          if (isdefault)
          {
            def_elem.n_default_elem++;
            def_elem.default_elem = (int*) realloc(def_elem.default_elem, def_elem.n_default_elem*sizeof(int));
            def_elem.default_elem[def_elem.n_default_elem-1] = iatom-1;
          }
/*        printf("data = |%s| %f %f %f %f %f %d %d %d\n", e[iatom-1].name, e[iatom-1].color[0], e[iatom-1].color[1], e[iatom-1].color[2], e[iatom-1].vdw_rad, e[iatom-1].bond_rad, e[iatom-1].valency, e[iatom-1].periodic_pos_x, e[iatom-1].periodic_pos_y);*/
        }
      }
    }
    free(line);
  }

#ifdef EBUG
  printf("DEFAULT ATOMS:"); for(i = 0; i < def_elem.n_default_elem; i++) printf(" %d", def_elem.default_elem[i]); putchar(10);
#endif

  fclose(ap);
  number_of_elements = iatom;
}

void save_atom_data(void)
{
  FILE *ap;
  char *path;
  path = malloc(sizeof(char) * (strlen(rcdir) + strlen(RC_ATOM) + 2));
  path[0] = 0;
  strcat(path, rcdir);
  strcat(path, "/"); /*POSSIBLE BUG IN WINDOWS !!!*/
  strcat(path, RC_ATOM);

  Write_Atom_data(path);

  free(path);
}

void free_atom_data(void)
{
  int i;
  if (def_elem.default_elem) free(def_elem.default_elem);
   
  for(i = 0; i < number_of_elements; i++)
    if (e[i].name) free(e[i].name);
  if (e) free(e);
}

