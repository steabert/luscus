/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"
#include"vstring.h"
#define R529 0.52917721067e0

char *orb_desc_short[] = 
{
  "u", "f", "i", "1", "2", "3", "s", "d"
};

MOL *mol;

void luscus_gtk_do_save_file(char *filename, int filetype)
{
  int i;
  char *ext;
  char *gv_filename;
  char *command;

  gv_filename = make_luscus_file_name(filename);
  save_gv_file(gv_filename);

  if (input_filetypes[filetype].backward)
    convert_to_target_file(filetype, gv_filename);

  if (gv_filename) free(gv_filename);
}

int get_filetype(char *filename)
{
  int i;
  char *ext;

  for(i = strlen(filename); i > 0 && filename[i] != '.'; i--);
  if (i == 0) return n_input_types - 1;
  ext = &filename[i+1];

  for(i = 0; i < n_input_types; i++)
    if (strcmp(ext, input_filetypes[i].extension) == 0) break;
  return i;
}

void save_gv_file(char *filename)
{
  int i;
  char c;
  char rndword[35];
  char *command;
  FILE *out;

  for(i = 0; i < 30; i++)
  {
    c = 0;
    while(!isalnum(c))
      c = (char) rand() % 128;
    rndword[i] = c;
  }
  rndword[30] = 0;
  strcat(rndword, ".");
  strcat(rndword, input_filetypes[0].extension);

  out = fopen(rndword, "wb");
  write_all_data_to_output(out, 1);

  fclose(out);

  /*command mv rndword filename*/
#ifdef WINDOWS
  command = (char*) malloc(sizeof(char) * (strlen(rndword) + strlen(filename) + 7));
  command[0] = 0;
  strcat(command, "copy ");
  strcat(command, rndword);
  strcat(command, " ");
  strcat(command, filename);
#else
  command = (char*) malloc(sizeof(char) * (strlen(rndword) + strlen(filename) + 8));
  command[0] = 0;
  strcat(command, "mv -f ");
  strcat(command, rndword);
  strcat(command, " ");
  strcat(command, filename);
#endif

  if (system(command) < 0)
    make_warning("ERROR: system command failed! Can't rename file!");

  free(command);

#ifdef WINDOWS
  command = (char*) malloc(sizeof(char) * (strlen(rndword) + 7));
  command[0] = 0;
  strcat(command, "del ");
  strcat(command, rndword);

  if (system(command) < 0)
    make_warning("ERROR: system command failed! Can't move file!");
  free(command);
#endif
}

/*converts gv file to some other file format*/
void convert_to_target_file(int filetype, char *filename)
{
  int n;
  int i;
  char *command;
/*  char *new_filename; */

/*  new_filename = (char*) malloc(sizeof(char) * (strlen(filename) + strlen(input_filetypes[filetype].extension) - strlen(input_filetypes[n_input_types-1].extension)));
  for(n = strlen(filename); n > 0 && filename[n] != '.'; n--);
  for(i = 0; i <= n; i++) new_filename[i] = filename[i];
  new_filename[i] = 0;
  strcat(new_filename, input_filetypes[filetype].extension);*/

  command = (char*) malloc(sizeof(char) * (strlen(input_filetypes[filetype].libpath) + strlen(input_filetypes[filetype].backward) + strlen(filename) + 5));
  command[0] = 0;
  strcat(command, input_filetypes[filetype].libpath);
#ifdef WINDOWS
  strcat(command, "\\");
#else
  strcat(command, "/");
#endif
  strcat(command, input_filetypes[filetype].backward);
  strcat(command, " ");
  strcat(command, filename);

  if (system(command) < 0) make_warning("ERROR: Can't convert file.\nFile will be saved in gv format.");

/*  free(new_filename);*/
  free(command);
}

void write_all_data_to_output(FILE *out, int write_grid)
{
  int i;
  char *filename;
/*  FILE *tmpin;*/

/*  close_input_file();
  filename = get_input_filename();
  tmpin = fopen(filename, "r");
  */

  for(i = 0; i < n_geometries; i++)
    write_mol_section(out, mol+i, write_grid);

  /*fclose(tmpin);*/
}

void write_mol_section(FILE *out, MOL *mp, int write_grid)
{
#ifdef EBUG
  printf("writing coord block\n"); fflush(stdout);
#endif
  write_coord_block(out, mp);
#ifdef EBUG
  printf("writing atom info block\n"); fflush(stdout);
#endif
  write_additional_atom_info_block(out, mp);
#ifdef EBUG
  printf("writing bond block!\n"); fflush(stdout);
#endif
  write_bond_block(out, mp);
#ifdef EBUG
  printf("writing dipole block\n"); fflush(stdout);
#endif
  write_dipole_block(out, mp);
#ifdef EBUG
  printf("writing vector block\n"); fflush(stdout);
#endif
  write_vector_block(out, mp);
#ifdef EBUG
  printf("writing triangle block\n"); fflush(stdout);
#endif
  write_triangle_block(out, mp);
#ifdef EBUG
  printf("writing sphere block\n"); fflush(stdout);
#endif
  write_sphere_block(out, mp);
#ifdef EBUG
  printf("writing surface block\n"); fflush(stdout);
#endif
  write_surface_block(out, mp);
#ifdef EBUG
  printf("writing cell block\n"); fflush(stdout);
#endif
  write_cell_block(out, mp);
#ifdef EBUG
  printf("writing textbox block");
#endif
  write_textbox_block(out, mp);
#ifdef EBUG
  printf("write_grid = %d\n", write_grid); fflush(stdout);
  if (write_grid) printf("writing grid block\nwriting inporb block\n"); fflush(stdout);
#endif
  if (write_grid) write_grid_block(out, mp);
/*#ifdef EBUG
  printf("writing inporb block\n"); fflush(stdout);
#endif
  if (write_grid) write_inporb_block(out, mp);*/
#ifdef EBUG
  printf("writing vibration block\n"); fflush(stdout);
#endif
  write_vibration_block(out, mp);
  fprintf(out, " <END>\n");
}

void write_coord_block(FILE *out, MOL *mp)
{
  int i;
  fprintf(out, "  %d\n", mp->natom);
  fprintf(out, "File generated by %s version %s\n", PROGRAM_NAME, VERSION);
  for(i = 0; i < mp->natom; i++)
    fprintf(out, " %s  %15.9f  %15.9f  %15.9f\n", mp->elem[i].name, mp->xyz[i][0], mp->xyz[i][1], mp->xyz[i][2]);
}

void write_additional_atom_info_block(FILE *out, MOL *mp)
{
  int i;

  if (!(mp->ishow & HAS_ATOMNAMES) && !(mp->ishow & HAS_MULLIKEN) && !(mp->ishow & HAS_LOPROP) && !(mp->ishow & HAS_ATOMNUMS) && !(mp->ishow & HAS_ATOMSYMM) && !(mp->ishow & HAS_COLOUR))
    return;

  fprintf(out, " <ATOM>\n");
  for(i = 0; i < mp->natom; i++)
  {
    if (mp->ishow & HAS_ATOMNAMES && mp->name[i])
      fprintf(out, "NAME = %s ", mp->name[i]);
    if (mp->ishow & HAS_MULLIKEN && mp->charge_m)
      fprintf(out, "MULLIKEN_CHARGE = %f ", mp->charge_m[i]);
    if (mp->ishow & HAS_LOPROP && mp->charge_l)
      fprintf(out, "LOPROP_CHARGE = %f ", mp->charge_l[i]);
    if (mp->ishow & HAS_ATOMNUMS && mp->additional_numeration)
      fprintf(out, "NUMBER = %d ", mp->additional_numeration[i]);
    if (mp->ishow & HAS_ATOMSYMM && mp->symmetry)
      fprintf(out, "SYMMETRY = %d ", mp->symmetry[i]);
    fprintf(out, "RED=%f GREEN=%f BLUE=%f", m->elem[i].color[0], m->elem[i].color[1], m->elem[i].color[2]);
    fputc(10, out);
  }
  fprintf(out, " </ATOM>\n");
}

void write_bond_block(FILE *out, MOL *mp)
{
  int i;
  if (!mp->nbond) return;
  fprintf(out, " <BOND>\n");
  fprintf(out, " AUTOMATIC = 0\n");
  for(i = 0; i < mp->nbond; i++)
    fprintf(out, " %d  %d  %d\n", mp->bond[i].iat1 + 1, mp->bond[i].iat2 + 1, mp->bond[i].bond_type);
  fprintf(out, " </BOND>\n");
}

void write_dipole_block(FILE *out, MOL *mp)
{
  if (mp->ishow & HAS_DIPOLE)
  {
    fprintf(out, " <DIPOLE>\n");
    fprintf(out, "%.6f  %.6f  %.6f\n", mp->dipole[0], mp->dipole[1], mp->dipole[2]);
    fprintf(out, " </DIPOLE>\n");
  }
}

void write_vector_block(FILE *out, MOL *mp)
{
  int i;
  for(i = 0; i < mp->nvector; i++)
  {
    fprintf(out, "<VECTOR>\n");
    fprintf(out, "RED = %f GREEN = %f BLUE = %f TRANSPARENCY = %f RADIUS = %f SHARPNESS = %f\n",
       	mp->vector_color[i][0], mp->vector_color[i][1], mp->vector_color[i][2], mp->vector_color[i][3],
        mp->radius[i], mp->sharpness[i]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->vector1[i][0],  mp->vector1[i][1],  mp->vector1[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->vector2[i][0],  mp->vector2[i][1],  mp->vector2[i][2]);
    fprintf(out, "</VECTOR>\n");
  }
}

void write_triangle_block(FILE *out, MOL *mp)
{
  int i;
  for(i = 0; i < mp->ntriangle; i++)
  {
    fprintf(out, "<TRIANGLE>\n");
    fprintf(out, "RED = %f GREEN = %f BLUE = %f TRANSPARENCY = %f\n ",
       	mp->triangle_color[i][0], mp->triangle_color[i][1], mp->triangle_color[i][2], mp->triangle_color[i][3]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->triangle1[i][0],  mp->triangle1[i][1],  mp->triangle1[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->triangle2[i][0],  mp->triangle2[i][1],  mp->triangle2[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->triangle3[i][0],  mp->triangle3[i][1],  mp->triangle3[i][2]);
    fprintf(out, "</TRIANGLE>\n");
  }
}

void write_sphere_block(FILE *out, MOL *mp)
{
  int i;
  for(i = 0; i < mp->nsphere; i++)
  {
    fprintf(out, "<SPHERE>\n");
    fprintf(out, "RED = %f GREEN = %f BLUE = %f TRANSPARENCY = %f RADIUS = %f \n",
       	mp->sphere_color[i][0], mp->sphere_color[i][1], mp->sphere_color[i][2], mp->sphere_color[i][3], mp->sphere_radius[i]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->sphere_center[i][0],  mp->sphere_center[i][1],  mp->sphere_center[i][2]);
    fprintf(out, "</SPHERE>\n");

  }
}

void write_surface_block(FILE *out, MOL *mp)
{
  int i;
  for(i = 0; i < mp->nsurf; i++)
  {
    fprintf(out, "<SURFACE>\n");
    fprintf(out, "RED = %f GREEN = %f BLUE = %f TRANSPARENCY = %f \n", mp->surf_color[i][0], mp->surf_color[i][1], mp->surf_color[i][2], mp->surf_color[i][3]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->surf1[i][0],  mp->surf1[i][1],  mp->surf1[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->surf2[i][0],  mp->surf2[i][1],  mp->surf2[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->surf3[i][0],  mp->surf3[i][1],  mp->surf3[i][2]);
    fprintf(out, "</SURFACE>\n");
  }
}

void write_cell_block(FILE *out, MOL *mp)
{
  int i;
  for(i = 0; i < mp->ncells; i++)
  {
    fprintf(out, "<CELL>\n");
    fprintf(out, "RED = %f GREEN = %f BLUE = %f TRANSPARENCY = %f \n", mp->cell_color[i][0], mp->cell_color[i][1], mp->cell_color[i][2], mp->cell_color[i][3]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->cell1[i][0],  mp->cell1[i][1],  mp->cell1[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->cell2[i][0],  mp->cell2[i][1],  mp->cell2[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->cell3[i][0],  mp->cell3[i][1],  mp->cell3[i][2]);
    fprintf(out, " %.6f  %.6f  %.6f\n", mp->cell4[i][0],  mp->cell4[i][1],  mp->cell4[i][2]);
    fprintf(out, "</CELL>\n");
  }
}

void write_textbox_block(FILE *out, MOL *mp)
{
  int i;
  for(i = 0; i < mp->ntextboxes; i++)
  {
    fprintf(out, "<TEXTBOX>\n");
    fprintf(out, "RED = %f GREEN = %f BLUE = %f X=%d Y=%d\n", mp->textboxes[i].color[0], mp->textboxes[i].color[1], mp->textboxes[i].color[2], mp->textboxes[i].coord_x, mp->textboxes[i].coord_y);
    fprintf(out, "%s\n", mp->textboxes[i].font);
    fprintf(out, "%s\n", mp->textboxes[i].message);
    fprintf(out, "</TEXTBOX>\n");
  }
}

void write_grid_block(FILE *out, MOL *mp) /*WRONG!!! MUST BE UPDATED!!!*/
{
  int i, j;
  FILE *gridin;
  char c;
  int n;
  double l;
  int size;

#ifdef EBUG
  printf("SAVING GRID BLOCK!\n"); fflush(stdout);
#endif

  if (!mp->ngrids) return;
  fprintf(out, " <GRID>");
  fputc(10, out); 
  fprintf(out, " N_of_MO=%d", mp->grid.nmo);
  fprintf(out, " N_of_Grids=%d", mp->ngrids);
  fprintf(out, " N_of_Points=%d", mp->grid.npoints);/*Not needed!???*/
  fprintf(out, " Block_Size=%d", mp->grid.block_size);
  fprintf(out, " N_Blocks=%d", mp->grid.nBlock); /*Not needed!???*/
  fprintf(out, " Is_cutoff=%d", mp->grid.isCutOff);
  fprintf(out, " CutOff=%8.4f", mp->grid.CutOff);
  fprintf(out, " N_P=%d", mp->grid.nPointsCutOff);
  fputc(10, out);

  fprintf(out, " N_INDEX= %d %d %d %d %d %d %d", mp->grid.iniIndex[0],
                   mp->grid.iniIndex[1], mp->grid.iniIndex[2],
                   mp->grid.iniIndex[3], mp->grid.iniIndex[4],
                   mp->grid.iniIndex[5], mp->grid.iniIndex[6]);
  fputc(10, out);
  fprintf(out, " Net= %d %d %d", mp->grid.npt[0],
                                 mp->grid.npt[1],
                                 mp->grid.npt[2]);
  fputc(10, out);
  fprintf(out, " Origin= %lf %lf %lf", mp->grid.origin[0] / R529,
                                       mp->grid.origin[1] / R529,
                                       mp->grid.origin[2] / R529);
  fputc(10, out);

  fprintf(out, " Axis_1= %lf %lf %lf", mp->grid.axisvec1[0] / R529,
                                       mp->grid.axisvec1[1] / R529,
                                       mp->grid.axisvec1[2] / R529);
  fputc(10, out);
  fprintf(out, " Axis_2= %lf %lf %lf", mp->grid.axisvec2[0] / R529,
                                       mp->grid.axisvec2[1] / R529,
                                       mp->grid.axisvec2[2] / R529);
  fputc(10, out);
  fprintf(out, " Axis_3= %lf %lf %lf", mp->grid.axisvec3[0] / R529,
                                       mp->grid.axisvec3[1] / R529,
                                       mp->grid.axisvec3[2] / R529);
  fputc(10, out);

/*  for(i = 0; i < mp->ngrids; i++)
  {
    fprintf(out, " File_pointers = ");
    for(j = 0; j < mp->grid.nBlock; j++) 
      fprintf(out, " %ld", mp->grid.pos[i * mp->grid.nBlock + j] - mp->grid.pos[0]);
    fputc(10, out);
  }*/
  fprintf(out," ORBOFF = %ld", mp->orb_ending_position);
  fputc(10, out);

  for(i = 0; i < mp->ngrids; i++)
  {
    if (strstr(mp->grid.titlesArr[i], "Orbital"))
    {
      fprintf(out, " GridName= Orbital  sym= %d index=%d Energ=%f occ=%f type=%s",
              mp->grid_symmetry[i], mp->grid_index[i], mp->grid_energy[i],
              mp->grid_occ[i], orb_desc_short[mp->orbital_type[i]]);
              fputc(10, out);
/*      fputc(mp->grid.titlesArr[i][1], out);
      fputc(mp->grid.titlesArr[i][2], out); 
      fprintf(out, " index=");
      fputc(mp->grid.titlesArr[i][5], out);
      fputc(mp->grid.titlesArr[i][6], out); 
      fputc(mp->grid.titlesArr[i][7], out); 
      fprintf(out, "  Energ=");
      for(j = 11; j < 20; j++)
        fputc(mp->grid.titlesArr[i][j], out);
      fprintf(out, " occ= ");
      for(j = 22; j < 26; j++)
        fputc(mp->grid.titlesArr[i][j], out);
      fprintf(out, " type= ");
      fputc(mp->grid.titlesArr[i][28], out);*/
    }
    else
    {
      fprintf(out, " GridName=%s", mp->grid.titlesArr[i]);
      fputc(10, out);
    }
  }
  gridin = open_current_luscus_file();
  fseek(gridin, mp->orb_starting_position, SEEK_SET);
  /*loop to almost to the orbital ending position*/
/*  for(ii = 0; ii < mp->orb_ending_position - mp->orb_starting_position; ii++)
  {
    fread(&c, sizeof(char), 1, gridin);
    fwrite(&c, sizeof(char), 1, out);
  }*/
  /*loop to the real orbital ending position; this is to ensure that output is written up to the right bite!*/
  fprintf(out, " <DENSITY>\n");
  while(ftell(gridin) < mp->orb_ending_position + mp->orb_starting_position)
  {
    size = fread(&c, sizeof(char), 1, gridin);
    if (size) /*failsafe! Shit happens!*/
      fwrite(&c, sizeof(char), 1, out);
    else
      break;
  }
  fprintf(out, "\n</DENSITY>\n");
  write_inporb_block(gridin, out, mp);

  fclose(gridin);
}

void write_inporb_block(FILE *gridin, FILE *out, MOL *mp)
{
  int i;
  int next = 1;
  int nirr;
  int iBas[8];
  char line[1024];
  char **iTypeIndex;
  /*search for the INPORB*/
  while(next)
  {
    abfgets(0, line, 1024, gridin);
    chomp(line);
    if (feof(gridin)) return; /*No INPORB; return*/
    if (strstr(line, "INPORB") != NULL) next = 0;
  }
  chomp(line); 
  fprintf(out, "%s\n", line);
  abfgets(0, line, 1024, gridin);
  nirr = atoi(line);
  if (nirr > 8) nirr = 8;
  iTypeIndex = (char**) malloc(nirr * sizeof(char*));
  for(i = 0; i < nirr; i++)
  {
    iBas[i] = atoi(line+8+i*8);
    iTypeIndex[i]=(char*) malloc(sizeof(char) * (1+iBas[i]));
  }
  chomp(line); 
  fprintf(out, "%s\n", line);
  /*search for the #INDEX and print every line*/
  next = 1;
  while(next)
  {
    abfgets(0, line, 1024, gridin);
    chomp(line); 
    fprintf(out, "%s\n", line);
    if (strstr(line, "#INDEX") != NULL) next = 0;
  }
  /*read the body of the #INDEX section*/

  for(i = 0; i < nirr; i++)
  {
    int j;
    int iline = -1, ich = 0;
    abfgets(0, line, 1024, gridin); /*  * 1234567890  */
    chomp(line);
    for(j = 0; j < iBas[i]; j++)
    {
      if (!(j%10))
      {
        abfgets(0, line, 1024, gridin);
        chomp(line);
        iline++;
        ich=0;
      }
      iTypeIndex[i][j] = line[ich+++2];
    }
  }

  /*put new data in the iTypeIndex*/

  for(i = 0; i < mp->ngrids; i++)
    if (mp->edited[i])
      iTypeIndex[mp->grid_symmetry[i]-1][mp->grid_index[i]-1] = orb_desc_short[mp->orbital_type[i]][0];

  /*print the body of the #INDEX*/

  for(i = 0; i < nirr; i++)
  {
    int j;
    int iline = -1;
    fprintf(out, "* 1234567890");
    for(j = 0; j < iBas[i]; j++)
    {
      if (!(j%10))
      {
        iline++;
        fprintf(out, "\n%d ", iline%10);
      }
      fprintf(out, "%c", iTypeIndex[i][j]);
    }
    fprintf(out, "\n");
  }
  fprintf(out, " </INPORB>\n </GRID>\n");
  for(i = 0; i < nirr; i++)
    free(iTypeIndex[i]);
  free(iTypeIndex);
}

#ifdef NOCOMPILE
void write_inporb_block(FILE *out, MOL *mp)
{
  int i, j;
  int k, k1, ik, kk, kkk;
  int iLab;
  char Lab[11]="0123456789";
  char c;
  char line[1024]; /*temporeraly*/
  char **iTypeIndex;
  FILE *gridin;
  int nirr, is, it, n;
  double l;
  int isIndex=0;
  char i18[19];
  int iBas[8];
  int errnum;

#ifdef EBUG
  printf("writing INPORB BLOCK!\n");
#endif
  if (!mp->ngrids) return;

  gridin = open_current_luscus_file();
  errnum = fseek(gridin, mp->orb_ending_position + mp->orb_starting_position, SEEK_SET);

  /*search for the INPORB!*/
  if (errnum == 0)
  {
    int next = 1;
    char *lin;
    while(next)
    {
      lin = read_line(gridin);
#ifdef EBUG
      printf("READ LINE |%s|\n", lin); fflush(stdout);
#endif
      if (strcasestr(lin, "inporb")) next = 0;
      if (strcasestr(lin, "/grid") || feof(gridin))
      {
        make_warning("Could not find INPORB section in luscus input. Therefore, it will not be saved.");

        free(lin);
        return;
      }
      free(lin);
    }
  }

  fprintf(out, "<INPORB>\n");

  /*write the data footer*/

  abfgets(0, line, 1024, gridin);
  is = 8;
  sscanf(line, "%d", &nirr);
  fprintf(out, "%s", line);

  for(i=0; i<nirr; i++)
  {
    sscanf(line+is, "%d", &j);
    iBas[i]=j;
    is+=8;
  }

  abfgets(0, line, 1024, gridin);
  chomp(line); 
  fprintf(out, "%s\n", line);

  abfgets(0, line, 1024, gridin);
  chomp(line); 
  fprintf(out, "%s\n", line);

  abfgets(0, line, 1024, gridin);
  chomp(line); 
  fprintf(out, "%s\n", "* orbitals generated by luscus");

  abfgets(0, line, 1024, gridin);
  chomp(line); 
  fprintf(out, "%s\n", line);

  abfgets(0, line, 1024, gridin);
  chomp(line); 
  fprintf(out, "%s\n", line);

  abfgets(0, line, 1024, gridin);
  chomp(line); 
  fprintf(out,"%s\n",line);

  iTypeIndex = (char**) malloc(sizeof(char*) * nirr);
  for(i = 0; i < nirr; i++)
    iTypeIndex[i]=(char*) malloc(sizeof(char) * (1+iBas[i]));

  while(!feof(gridin))
  {
    if(!abfgets(0, line, 1024, gridin)) continue;
    if (strcasestr(line, "/inporb")) break;
    chomp(line);
    fprintf(out, "%s\n", line);

    if(strstr(line, "* OCCUPATION NUM") != NULL)
    {
      for(i = 0; i < nirr; i++)
      {
        it=0;
        if(iBas[i] == 0) continue;
        while(it < iBas[i])
        {
          abfgets(0, line, 1024, gridin);
          chomp(line);
          fprintf(out,"%s\n",line);

          n = ((int) strlen(line)) / 18;
          for(j = 0; j < n; j++)
          {
            strncpy(i18, line+j*18, 18);
            i18[18] = 0;
            sscanf(i18+16, "%d", &k);

	    kkk = 1;
            if (i18[15]=='-') kkk=-1;
            i18[14]=0;
            sscanf(i18, "%lf", &l);
            kk=1;
            for (ik = 0; ik<k; ik++) kk=kk*10;
            if (kkk == 1) l=kk*l;
            else l = (1.0f/kk)*l;  
            c='2';
            if(l > 1.99) c='I';
            if(l < 0.01) c='S';
            iTypeIndex[i][it]=c;
            it++;
          }
        }
        iTypeIndex[i][it]=0;
      }
    }

    if (strstr(line,"#INDEX") != NULL)
    {
      isIndex=1;

      for(i = 0; i < nirr; i++)
      {
        abfgets(0, line, 1024, gridin);
	printf("line: %s\n", line); fflush(stdout);
        it = 0;
        if(iBas[i] == 0) continue;
        while(it < iBas[i])
        {
          abfgets(0, line, 1024, gridin);
          chomp(line);
          n=(int) strlen(line)-2;
          if(n < 0)
          {
            fprintf(stderr, "InpOrb is in Old format\n");
            exit(1);
          }
          for(j=0; j<n; j++)
          {
            iTypeIndex[i][it]=line[j+2];
            it++;
          }
        }
        iTypeIndex[i][it]=0;
      }
    }
  }

  /*put new data in the iTypeIndex*/

  for(i = 0; i < mp->ngrids; i++)
    if (mp->edited[i])
      iTypeIndex[mp->grid_symmetry[i]-1][mp->grid_index[i]-1] = orb_desc_short[mp->orbital_type[i]][0];

  line[0] = 0;
  if (!isIndex) fprintf(out,"#INDEX\n");

  for(i = 0; i < nirr; i++)
  {
    fprintf(out, "* 1234567890\n");
    if (iBas[i] == 0) continue;

    k=0;
    k1=0;
    iLab=0;
    line[k1]=Lab[iLab];
    line[k1+1]=' ';
    k1=k1+2;
    iLab=iLab+1;

    for(j = 0; j < iBas[i]; j++)
    {
      line[k1]=iTypeIndex[i][j];
      k++;
      k1++;
      line[k1]=0;
      if(k == 10)
      {
        fprintf(out, "%s\n", line);
        k1 = 0;
        k = 0;
        line[k1] = Lab[iLab];
        line[k1 + 1] = ' ';
        k1 = k1 + 2;
        iLab = iLab + 1;
        if (iLab > 9) iLab = 0;
      }
    }

    fprintf(out,"%s\n",line);
  }

  fprintf(out, "</INPORB>\n</GRID>\n");
  fputc(10, out);


  fclose(gridin);

  for(i = 0; i < nirr; i++)
    free(iTypeIndex[i]);
  free(iTypeIndex);

}
#endif

void write_vibration_block(FILE *out, MOL *mp)
{
  int i, j;
  for(i = 0; i < mp->nvibr; i++)
  {
    fprintf(out, "<VIBRATION>\n");
    fprintf(out, " FREQ = %.6f ", mp->freq[i]);
    if (mp->ishow & HAS_IRINT)
      fprintf(out, " IR_INT = %.6f ", mp->ir_intensity[i]);
    if (mp->ishow & HAS_RAMAN)
      fprintf(out, " RAMAN_INT = %.6f ", mp->raman_intensity[i]);
    fputc(10, out);
    for(j = 0; j < mp->natom; j++)
      fprintf(out, " %.6f  %.6f  %.6f\n", mp->normal_mode[i][j][0],
                mp->normal_mode[i][j][1], mp->normal_mode[i][j][2]);
    fprintf(out, "</VIBRATION>\n");
  }
}

/*reverse of read_bin_block*/
void write_bin_block(FILE *fo, char *line)
{
  int rec_len, rec_8;
  int fmarker = get_fmarker();
  rec_8 = 0;

  rec_len = strlen(line);

  fwrite(&rec_len, sizeof(rec_len), 1, fo);
  if (fmarker) fwrite(&rec_8, sizeof(rec_8), 1, fo);

  fwrite(line, (unsigned int) rec_len, 1, fo);

  fwrite(&rec_len, sizeof(rec_len), 1, fo);
  if (fmarker) fwrite(&rec_8, sizeof(rec_8), 1, fo);
}

