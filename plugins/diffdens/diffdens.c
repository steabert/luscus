#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"diffdens.h"

int is_num(char *str)
{
  int i;
  for(i = 0; i < strlen(str); i++)
    if (!isdigit(str[i]) && str[i] != 45 && str[i] != 46) return 0;
  return 1;
}

void skip_grid(FILE *fp)
{
  char *line;
  int inipos, endpos = 0;
  int next = 1;

  while(next)
  {
    char *ptr;
    line = read_line(fp);
    if (find_substring_case(line, "orboff"))
    {
      ptr = strchr(line, '=') + 1;
      endpos = atol(ptr);
      next = 0;
    }
    if (find_substring_case(line, "</grid>") || feof(fp))
    {
      free(line);
      return;
    }
 
    free(line);
  }
  next = 1;

  while(next)
  {
    char *ptr;
    line = read_line(fp); 
    if (find_substring_case(line, "<density>")) next = 0;
    if (find_substring_case(line, "</grid>") || feof(fp))
    {
      free(line);
      return;
    }
    free(line);
  }
  next = 1;
  inipos = ftell(fp);
  fseek(fp, inipos+endpos, SEEK_SET);
  while(next)
  {
    char *ptr;
    line = read_line(fp);

    if (find_substring_case(line, "</grid>") || feof(fp))
    {
      free(line);
      return;
    }

    free(line);
  }

  return;
}

int read_luscus_file(char *fname, char **buffer)
{
  FILE *fp;
  char *buf;
  ssize_t chartot = 1;

  buf = (char*) malloc(sizeof(char));
  buf[0] = 0;
  fp = fopen(fname, "r");
  if (!fp) return 1;
  while(!feof(fp))
  {
    char *line = NULL;
    line = read_line(fp);
    if (find_substring_case(line, "<grid>"))
    {
      skip_grid(fp);
    }
    else if (find_substring_case(line, "<end>"))
    {
      free(line);
      break;
    }
    else
    {
      ssize_t nchar;
      nchar = strlen(line)+1;
      buf = (char*) realloc(buf, (chartot+nchar) *sizeof(char));
      strncat(buf, line, nchar-1);
      buf[chartot+nchar-2] = 10;
      buf[chartot+nchar-1] = 0;
      chartot += nchar;
    }
    free(line);
  }

  *buffer = buf;
  return 0;
}

void write_output_file(char *filename, char *buffer)
{
  FILE *of;
  of = fopen(filename, "w");

  fputs(buffer, of);

  fclose(of);
  return;
}

void write_output_grid(char *filename, GRID_T resgrid)
{
  FILE *of;
  char buffer[25];

  of = fopen(filename, "a");

  fputs("<GRID>", of);
  fputc(10, of);
  fputs(" N_of_MO=1", of);
  fputs(" N_of_Grids=1", of);
  fputs(" N_of_Points=", of);
  snprintf(buffer, 25, "%d", resgrid.npoints);
  fputs(buffer, of);
  fputs(" Block_Size=", of);
  fputs(buffer, of); 
  fputs(" N_Blocks=1", of);
  fputs(" Is_cutoff=0", of);
  fputs(" CutOff=0.0", of);
  fputs(" N_P=", of);
  fputs(buffer, of);
  fputc(10, of);
  fputs(" N_INDEX= 0 0 0 0 0 0 0", of);
  fputc(10, of);
  fputs(" Net=", of);
  snprintf(buffer, 25, " %d", resgrid.npt[0]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %d", resgrid.npt[1]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %d", resgrid.npt[2]);
  fputs(buffer, of);
  fputc(10, of);
  fputs(" Origin=", of);
  snprintf(buffer, 25, " %lf", resgrid.origin[0]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.origin[1]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.origin[2]);
  fputs(buffer, of);
  fputc(10, of);
  fputs(" Axis_1=", of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec1[0]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec1[1]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec1[2]);
  fputs(buffer, of);
  fputc(10, of);
  fputs(" Axis_2=", of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec2[0]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec2[1]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec2[2]);
  fputs(buffer, of);
  fputc(10, of);
  fputs(" Axis_3=", of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec3[0]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec3[1]);
  fputs(buffer, of);
  snprintf(buffer, 25, " %lf", resgrid.axisvec3[2]);
  fputs(buffer, of);
  fputc(10, of);
  fputs(" ORBOFF=", of);
  snprintf(buffer, 25, " %d", /*ftell(of)+*/resgrid.npoints*8);
  fputs(buffer, of);
  fputc(10, of);
  fputs( " GridName=  Density (difference)", of);
  fputc(10, of);
  fputs(" <DENSITY>", of);
  fputc(10, of);
  fclose(of);
  of = fopen(filename, "ab");
  fwrite(resgrid.values, 8, resgrid.npoints, of);
  fputc(10, of);
  fputs(" </DENSITY>", of);
  fputc(10, of);
  fputs("</GRID>", of);
  fputc(10, of);
  fclose(of);
}

int main(int argc, char **argv)
{
  int i;
  double weight;
  char *buffer;
  GRID_T grid1, grid2, resgrid;

  if (!strcmp(argv[1], "--version"))
  {
    printf("diffdens version %s\n", VERSION);
    return EXIT_SUCCESS;
  }

  if (!is_num(argv[3]))
  {
    fprintf(stderr, "ERROR: ARGUMENT 3 (WEIGHT) SHOULD BE A NUMERICAL VALUE!\n");
    return EXIT_FAILURE;
  }

  weight=atof(argv[3]);

  if (read_luscus_file(argv[1], &buffer))
  {
    fprintf(stderr, "ERROR reading %s\n", argv[1]);
    return EXIT_FAILURE;
  }

  /*open luscus file 1*/
  open_luscus_file(argv[1], &grid1);
  if (grid1.npoints == 0)
  {
    fprintf(stderr, "ERROR READING GRID 1\n");
    return EXIT_FAILURE;
  }
#ifdef EBUG
  printf("grid 1:\n");
  printf("ngrid = %d\n", grid1.ngrid);
  printf("nmo = %d\n", grid1.nmo);
  printf("npoints = %d\n", grid1.npoints);
  printf("npt = %d %d %d\n", grid1.npt[0], grid1.npt[1], grid1.npt[2]);
  printf("nBlock = %d\n", grid1.nBlock);
/*  printf("len = %f %f %f\n", grid1.len[0], grid1.len[1], grid1.len[2]);*/
  printf("origin = %f %f %f\n", grid1.origin[0], grid1.origin[1], grid1.origin[2]);
#endif
  /*open luscus file 2*/
  open_luscus_file(argv[2], &grid2);
  if (grid2.npoints == 0)
  {
    deallocate_grids(&grid1);
    fprintf(stderr, "ERROR READING GRID 2\n");
    return EXIT_FAILURE;
  }

#ifdef EBUG
  printf("grid 2:\n");
  printf("ngrid = %d\n", grid2.ngrid);
  printf("nmo = %d\n", grid2.nmo);
  printf("npoints = %d\n", grid2.npoints);
  printf("npt = %d %d %d\n", grid2.npt[0], grid2.npt[1], grid2.npt[2]);
  printf("nBlock = %d\n", grid2.nBlock);
/*  printf("len = %f %f %f\n", grid2.len[0], grid2.len[1], grid2.len[2]);*/
  printf("origin = %f %f %f\n", grid2.origin[0], grid2.origin[1], grid2.origin[2]);
#endif

  /*test if grid sizes match*/
  if (grid1.npt[0] != grid2.npt[0] &&
      grid1.npt[1] != grid2.npt[1] && 
      grid1.npt[2] != grid2.npt[2])
  {
    fprintf(stderr, "ERROR: grids are not equal!\n");
    return EXIT_FAILURE;
  }

  /*initialize the resulting grid*/
  resgrid.nmo = 1;
  resgrid.ngrid = 1;
  resgrid.npoints = grid1.npoints;
  resgrid.block_size = grid1.npoints;
  resgrid.nBlock = 1;
  resgrid.isCutOff = grid1.isCutOff;
  resgrid.CutOff = grid1.CutOff;
  resgrid.nPointsCutOff = grid1.nPointsCutOff;
  resgrid.npt[0] = grid1.npt[0];
  resgrid.npt[1] = grid1.npt[1];
  resgrid.npt[2] = grid1.npt[2];
  resgrid.origin[0] = grid1.origin[0];
  resgrid.origin[1] = grid1.origin[1];
  resgrid.origin[2] = grid1.origin[2];
  resgrid.axisvec1[0] = grid1.axisvec1[0];
  resgrid.axisvec1[1] = grid1.axisvec1[1];
  resgrid.axisvec1[2] = grid1.axisvec1[2];
  resgrid.axisvec2[0] = grid1.axisvec2[0];
  resgrid.axisvec2[1] = grid1.axisvec2[1];
  resgrid.axisvec2[2] = grid1.axisvec2[2];
  resgrid.axisvec3[0] = grid1.axisvec3[0];
  resgrid.axisvec3[1] = grid1.axisvec3[1];
  resgrid.axisvec3[2] = grid1.axisvec3[2];
  allocate_grids(&resgrid);

  /*make density difference*/
  for(i = 0; i < grid1.npoints; i++)
  {
    resgrid.values[i] = grid1.values[i] + weight * grid2.values[i];
#ifdef EBUG
/*    printf("%d %f %f %f\n", i, grid1.values[i], grid2.values[i], resgrid.values[i]); fflush(stdout);*/
#endif
  }


  write_output_file(argv[1], buffer);
  
  free(buffer);

  write_output_grid(argv[1], resgrid);

  deallocate_grids(&grid1);
  deallocate_grids(&grid2);
  deallocate_grids(&resgrid);

  return EXIT_SUCCESS;
}

