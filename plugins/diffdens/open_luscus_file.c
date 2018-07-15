#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"diffdens.h"
#define MIN(x,y) (x<y? x : y)
#define MAX(x,y) (x>y? x : y)


static long starting_position;

/* ---- helper functions ---- */

char *read_line(FILE* fp)
{
  int i = 1;
  char c;
  char *str = NULL;
  str = (char*) malloc(sizeof(char) * i);
  str[i-1] = 0;
  do
  {
    c = fgetc(fp);
    if (c == 10 || c == 0 || feof(fp)) return str;
    else
    {
      str = (char*) realloc(str, sizeof(char) * ++i);
      str[i-2] = c;
      str[i-1] = 0;
    }
  }
  while(c != 10 && c != 0 && !feof(fp));

  return str;
}

char *find_substring_case(char *haystack, const char *needle)
{
  int i, j;
  int len_hay, len_nee;
  int nmatch;

  len_hay = strlen(haystack);
  len_nee = strlen(needle);

  for(i = 0; i <= len_hay - len_nee; i++)
  {
    nmatch = 0;
    for(j = 0; j < len_nee; j++)
    {
      if (toupper(haystack[i+j]) == toupper(needle[j])) nmatch++;
      else break;
    }
    if (nmatch == len_nee) return haystack + i;
  }
  return NULL;
}

char *my_read_str_value(char *str)
{
  int i, j;
  char *ptr;
  char *val = NULL;
  int next, quoted = 0;

  if (str == NULL) return NULL;
  if (str[0] == 0) return NULL;
  ptr = strchr(str, '=') + 1;
  if (ptr == NULL) return NULL;
  if (ptr[0] == 0) return NULL;

  while(ptr[0] == 32 || ptr[0] == 9)
  {
    if (ptr[0] == 0) return NULL;
    ptr++;
  }
  if (ptr[0] == 0) return NULL;

  i = j = 0;
  val = (char*) malloc(sizeof(char) * (i+1));
  val[i] = 0;
  next = 1;
  if (ptr[i] < 33 || ptr[i] > 126) next = 0;
  while(next)
  {
    if (ptr[j] == 0) return val;
    else if (ptr[j] == 34)
    {
      if (quoted) quoted = 0;
      else quoted = 1;
    }
    else
    {
      val = (char*) realloc(val, sizeof(char) * (++i+1));
      val[i-1] = ptr[j];
      val[i] = 0;
    }
    j++;
    if (quoted && (ptr[j] < 9 || (ptr[j] > 13 && ptr[j] < 32) || ptr[j] > 126)) next = 0;
    else if (!quoted && (ptr[j] < 33 || ptr[j] > 126)) next = 0;
  }

  return val;
}

int my_read_int_value(char *str)
{
  char *ptr;
  ptr = strchr(str, '=') + 1;
  if (ptr) return atoi(ptr); 
  else return 0;
}

double my_read_dble_value(char *str)
{
  char *ptr;
  double tmp;
  ptr = strchr(str, '=') + 1;
  if (ptr)
  {
    sscanf(ptr, "%lf", &tmp);
    return tmp; /*corrected 2014.01.24*/
    /*return atof(ptr);*/
  }
  else return 0.F;
}

/*char *strdiet(char* in)
{
  int i, j, k;
  j = strlen(in) - 1;
  while(in[j] == 32 || in[j] == 9) j--;
  in[j] = 0;
  for(i = 0; in[i] == 32 || in[i] == 9; i++);
  for(k = 0; k < j-i; k++) in[k] = in[k+i];
  in[j-1] = 0;
  return in;
}*/

/*char *get_one_word(char*);*/
/*char *get_ptr_ith_word(char*, int);*/
/*int line_is_empty(char*);*/

/* ---- luscus dealing functions ---- */

void allocate_grids(GRID_T *grid)
{
  grid->ByteCutOff = (int*) malloc(sizeof(int) * grid->npoints);
  grid->pos = (long*) malloc(sizeof(long) * (unsigned int) grid->nBlock * grid->ngrid);
  grid->titlesArr = (char**) malloc(sizeof(char*) * grid->ngrid);
  grid->FltTitle = (int*) malloc(sizeof(int) * grid->ngrid);
  grid->iType = (char*) malloc(sizeof(char) * grid->ngrid);
  grid->dependence = (int*) malloc(sizeof(int) * grid->ngrid);
  grid->values = (double*) malloc(sizeof(double) * grid->npoints);
}

void deallocate_grids(GRID_T *grid)
{
/*  printf("deallocating ByteCutOff\n"); fflush(stdout);*/
  if (grid->ByteCutOff) free(grid->ByteCutOff);
  grid->ByteCutOff = NULL;
/*  printf("deallocating pos\n"); fflush(stdout);*/
  if (grid->pos) free(grid->pos);
  grid->pos = NULL;
/*  printf("deallocating title\n"); fflush(stdout);
  if (grid->title) free(grid->title);
  grid->title = NULL;*/
/*  printf("deallocating titlesArr\n"); fflush(stdout);*/
  if (grid->titlesArr) free(grid->titlesArr);
  grid->titlesArr = NULL;
/*  printf("deallocating FltTitle\n"); fflush(stdout);*/
  if (grid->FltTitle) free(grid->FltTitle);
  grid->FltTitle = NULL;
/*  printf("deallocating itype\n"); fflush(stdout);*/
  if (grid->iType) free(grid->iType);
  grid->iType = NULL;
/*  printf("deallocating dependence\n"); fflush(stdout);*/
  if (grid->dependence) free(grid->dependence);
  grid->dependence = NULL;
/*  printf("deallocating values\n"); fflush(stdout);*/
  if (grid->values) free(grid->values);
  grid->values = NULL;

  grid->ngrid = 0;
  grid->npoints = 0;
}

static int ReadBlock(double* buf, int np, FILE* fp)
{
  if (fread(buf, (size_t) np, 1, fp) !=1) return -1;
  return np;
}

void read_grid_block(FILE *fp, GRID_T *grid)
{
  char *line;
  char *ptr;
  long unix2dos;
  long orb_ending_position, pos1, pos2, ip;
  int ii, i, j;
  int next = 1;

  line = read_line(fp);
/*  printf("line = |%s|\n", line);*/
  if (find_substring_case(line,"N_OF_MO"))
    grid->nmo = my_read_int_value((char*) find_substring_case(line,"N_OF_MO"));
  if (find_substring_case(line,"N_of_Grids"))
    grid->ngrid = my_read_int_value((char*) find_substring_case(line,"N_of_Grids"));
  if (find_substring_case(line,"N_of_Points"))
    grid->npoints = my_read_int_value((char*) find_substring_case(line,"N_of_Points"));
  if (find_substring_case(line,"Block_Size"))
    grid->block_size = my_read_int_value((char*) find_substring_case(line,"Block_Size"));
  if (find_substring_case(line,"N_Blocks"))
    grid->nBlock = my_read_int_value((char*) find_substring_case(line,"N_Blocks"));
  if (find_substring_case(line,"Is_cutoff"))
    grid->isCutOff = my_read_int_value((char*) find_substring_case(line,"Is_cutoff"));
  if (find_substring_case(line,"CutOff"))
    grid->CutOff = my_read_dble_value((char*) find_substring_case(line,"CutOff"));
  if (find_substring_case(line,"N_P"))
    grid->nPointsCutOff = my_read_int_value((char*) find_substring_case(line,"N_P"));
  free(line);

  allocate_grids(grid);

  line = read_line(fp);
/*  printf("line = |%s|\n", line);*/
  ptr = strchr(line, '=') + 1;
  sscanf(ptr, "%d %d %d %d %d %d %d", &(grid->iniIndex[0]),
                &(grid->iniIndex[1]), &(grid->iniIndex[2]),
                &(grid->iniIndex[3]), &(grid->iniIndex[4]),
                &(grid->iniIndex[5]), &(grid->iniIndex[6]));
  free(line);

  line = read_line(fp);
/*  printf("line = |%s|\n", line);*/
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%d %d %d", &(grid->npt[0]), &(grid->npt[1]), &(grid->npt[2]));
  free(line);

  if (grid->npoints != grid->npt[0] * grid->npt[1] * grid->npt[2])
  {
    deallocate_grids(grid);
    fprintf(stderr, "ERROR: Ill-defined number of grids");
    grid->npoints = 0;
    return;
  }
  if (grid->ngrid == 0)
  {
    deallocate_grids(grid);
    fprintf(stderr, "ERROR: No data in grid section");
    grid->npoints = 0;
    return;
  }

  line = read_line(fp);
/*  printf("line = |%s|\n", line);*/
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(grid->origin[0]), &(grid->origin[1]), &(grid->origin[2]));
  free(line);

  line = read_line(fp);
/*  printf("line = |%s|\n", line);*/
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(grid->axisvec1[0]),
                              &(grid->axisvec1[1]),
                              &(grid->axisvec1[2]));
  free(line);

  line = read_line(fp);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(grid->axisvec2[0]),
                              &(grid->axisvec2[1]),
                              &(grid->axisvec2[2]));
  free(line);

  line = read_line(fp);
  ptr = strchr(line, '=') + 1;
  sscanf (ptr, "%lf %lf %lf", &(grid->axisvec3[0]),
                              &(grid->axisvec3[1]),
                              &(grid->axisvec3[2]));
  free(line);

  for(i = 0; i < 3; i++)
  {
    grid->origin[i];
    grid->axisvec1[i];
    grid->axisvec2[i];
    grid->axisvec3[i];
  }

  ii = 0;
  for (i = 0; i < grid->ngrid; i++)
  {
    for(j = 0; j < grid->nBlock-1; j++)   
    {
      grid->pos[i * grid->nBlock + j] = (grid->ngrid * j + i) * 8 * grid->block_size;
    }
    grid->pos[i * grid->nBlock + grid->nBlock-1] =
        (grid->ngrid * (grid->nBlock-1)) * 8 * grid->block_size + 8 * (grid->npoints - (grid->nBlock - 1) * grid->block_size) * i;
  }

  /*orbital index position*/
  line = read_line(fp);
  ptr = strchr(line, '=') + 1;
  orb_ending_position = atol(ptr);
  free(line);

  /*orbital data*/
  for(i = 0; i < grid->ngrid; i++)
  {
    line = read_line(fp);
    if (find_substring_case(line,"GridName"))
      grid->titlesArr[i] = my_read_str_value((char*) find_substring_case(line,"GridName"));
    free(line);
  }

  while(next)
  {
    pos1=ftell(fp);
    line = read_line(fp);
    if (find_substring_case(line, "density")) next = 0;
    if (find_substring_case(line, "/grid") || feof(fp))
    {
      fprintf(stderr, "ERROR can't find <density> section\n");
      deallocate_grids(grid);
      grid->npoints=0;
      free(line);
      return;
    }
    free(line);
  }
  next = 1;

  pos2=ftell(fp);
  ip=pos2-pos1;
  switch(ip)
  {
    case 11: unix2dos=0; break;
    case 12: unix2dos=1; break;
    default:puts("I can't determine EOL characters\n");
    deallocate_grids(grid);
    grid->npoints=0;
/*    free(line);*/
    return;
  }

  starting_position = ftell(fp) + unix2dos;
}

void load_grid_data(GRID_T *grid, FILE *in)
{
  double *readbuffer;
  int i;
  int found = 0;
  int igrid = 0;
  int iblock;
  int x, y, z;
  int itek, success;


  for(i = 0; i < grid->ngrid; i++)
    if (find_substring_case(grid->titlesArr[i], "Density"))
    {
      igrid = i;
      found = 1;
    }

  if (!found)
  {
    fprintf(stderr, "grid differece error: Can't find density!\n");
    grid->npoints = 0;
    return;
  }

  readbuffer = (double *) calloc ((unsigned int) grid->block_size, sizeof (double));

  iblock = 0;
  x = y = z = 0;

  for (i = 0; i < grid->npoints; i++, x++)
  {
    if (i == grid->block_size * iblock)
    {
      int np;
      fseek(in, starting_position + grid->pos[igrid * grid->nBlock + iblock], SEEK_SET);
      iblock++;
      np = grid->block_size;
      if (iblock == grid->nBlock)
      np = grid->npoints - grid->block_size * (grid->nBlock - 1);
      if(grid->isCutOff)
      {
        /* calculate number of uncutted data in this block */
        int ic = 0, ii;          
        for(ii = 0; ii < np; ii++)
          if(grid->ByteCutOff[ii+i] == 1) ic++;
            success = ReadBlock(readbuffer, ic, in);
      }
      else
      {
        success = ReadBlock(readbuffer, np*8, in);
      }
#ifdef EBUG
      printf("address %d values: %f %f\n", i, readbuffer[0], readbuffer[1]);
#endif
      itek = 0;
    }
    if (x == grid->npt[0])
    {
      x = 0;
      y++;
      if (y == grid->npt[1])
      {
        y = 0;
        z++;
      }
    }
    if(grid->isCutOff && grid->ByteCutOff[i]==0) /*NOT NEEDED;*/
    {                                            /*NOT NEEDED;*/
      grid->values[i]=0.0;                       /*NOT NEEDED;*/
      itek--;                                    /*NOT NEEDED;*/
    }                                            /*NOT NEEDED;*/
    else                                         /*NOT NEEDED;*/
      grid->values[i] = readbuffer[itek];
    itek++;

    if (i == 0)
    {
      grid->minval = 10000.0;
      grid->maxval = -10000.0;
    }
    else
    {
      /* let's cut off huge values */
      if (grid->values[i] > 2) grid->values[i] = 2;
      grid->minval = MIN(grid->minval, grid->values[i]);
      grid->maxval = MAX(grid->maxval, grid->values[i]);
    }
  }
  free(readbuffer);

}

void open_luscus_file(char *filename, GRID_T *grid)
{
  FILE *fp;
  int next = 1;
  char *line;

  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    grid->npoints = 0;
    fprintf(stderr, "ERROR OPENING FILE |%s|\n", filename);
    return;
  }

  while(next)
  {
    if (feof(fp)) next = 0;
    if (next)
    {
      line = read_line(fp);
      if (find_substring_case(line, "<grid>"))
      {
        read_grid_block(fp, grid);
        if (grid->npoints == 0) return;
      }
      free(line);
    }
  }

  load_grid_data(grid, fp);

  fclose(fp);
}

