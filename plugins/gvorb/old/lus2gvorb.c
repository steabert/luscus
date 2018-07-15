#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#define MAXCHAR 1024

int Fmarker = 0;

char *chomp(char *str)
{
  int i;
  i = (int)strlen (str) - 1;
  if (str[i] == '\n') str[i] = 0;
  return str;
}

/*This is the strcasestr function; Some opeerating systems don't have it included!*/
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

static int read_bin_block(void *buf, int bufsz, FILE * fp)
{
  int rec_len,rec_len1, rec_8;

  if (fread(&rec_len, sizeof(rec_len), 1, fp)!=1) return -1;

  if(Fmarker==1)  
    if (fread(&rec_8, sizeof(rec_8), 1, fp)!=1) return -1;
  
  if (buf)
  {
    if (rec_len>bufsz) return -1;
    if (fread(buf,(unsigned int)rec_len,1,fp)!=1) return -1;
  }
  else
  {
    fseek(fp,rec_len,SEEK_CUR);
  }

  if (fread(&rec_len1, sizeof(rec_len1), 1, fp)!=1) return -1;
  if(Fmarker==1)  
    if (fread(&rec_8, sizeof(rec_8), 1, fp)!=1) return -1;

  if (rec_len1!=rec_len) return -1;
  return rec_len;
}

int a1(FILE *fp, FILE *ouf)
{
  int i, j;
  int ngrids;
  int next = 1;
  char line[MAXCHAR];
  long orbend;
  char *ptr;
  int nindexlin; /*number of lines in #INDEX*/
  int nstar; /*number of lines in #INDEX that start with '*'*/
  int iBas[8];
  int nIrr;
  int is;
  
  while(next)
  {
    fgets(line, MAXCHAR, fp);
    if (feof(fp)) return 0;
    if (strstr(line, "<GRID>")) next = 0;
  }
  fgets(line, MAXCHAR, fp);
  if (feof(fp)) return 0;
  ptr = find_substring_case(line, "N_of_Grids");

  ptr = strchr(ptr, '=') + 1;
  ngrids = atoi(ptr);

  for(i = 0; i < 7 + ngrids; i++)
  {
    fgets(line, MAXCHAR, fp);
    if (feof(fp)) return 0;
  }

  ptr = strchr(line, '=');
  if (!ptr) return 0;
  ptr++;
  orbend = atol(ptr);

  for(i = 0; i < ngrids; i++)
  {
    fgets(line, MAXCHAR, fp);
    if (feof(fp)) return 0;
  }

  fseek(fp, orbend, SEEK_CUR);

  for(j = 0; j < MAXCHAR; j++) line[j] = 0;
  if (read_bin_block(line, MAXCHAR, fp) == -1) return 0;
  chomp(line);
  is = 8;
  sscanf(line, "%d", &nIrr);
  for(i = 0; i < nIrr; i++)
  {
    sscanf(line+is, "%d", &j);
    iBas[i]=j;
    is+=8;
  }

/*  printf("NIRR = %d\n", nIrr); fflush(stdout);*/

  for(i = 0; i < 7; i++)
  {
    for(j = 0; j < MAXCHAR; j++) line[j] = 0;
    if (read_bin_block(line, MAXCHAR, fp) == -1) return 0;
/*    printf("LINE = |%s|\n", line);*/
    if (feof(fp)) return 0;
    fprintf(ouf,"%s\n", line);
  }

  next = 1;

  while(!feof(fp) && next)
  {
    for(j = 0; j < MAXCHAR; j++) line[j] = 0;
    if (!read_bin_block(line, MAXCHAR, fp)) continue;
    chomp(line);
    fprintf(ouf,"%s\n",line);
    if (strstr(line, "#INDEX")) next = 0;
  }

  /*calculate the number of lines to be read*/

  nindexlin = 0;
  nstar = 0;
  for(i = 0; i < nIrr; i++)
  {
    nstar += (1 + (iBas[i] - 1)/100);
    if (iBas[i]) nindexlin += (1 + (iBas[i] - 1)/10);
  }
  nindexlin += nstar;

  /*read next nindexlin lines*/
  for(i = 0; i < nindexlin; i++)
  {
    for(j = 0; j < MAXCHAR; j++) line[j] = 0;
    if (read_bin_block(line, MAXCHAR, fp) == -1) return 0;
    if (feof(fp)) return 0;
    fprintf(ouf,"%s\n", line);
  }

  return 1;
}

int main(int argc, char **argv)
{
  FILE *fp;
  FILE *ouf;
  char outputfile[MAXCHAR];
  char *ptr;

  if (strstr(argv[argc-1], ".lus") == NULL)
  {
    fprintf(stderr, "ERROR: Input must be a luscus file!\n");
    return EXIT_FAILURE;
  }

  strcpy(outputfile, argv[argc-1]);
  ptr = strrchr(outputfile, '.');
  if (strstr(ptr, "lus") == NULL)
    ptr = argv[argc-1] + strlen(argv[argc-1]);
  else
    ptr[0] = 0;
  strcpy(ptr, ".GvOrb");

  fp = fopen(argv[argc - 1], "rb");

  if (fp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open luscus file: %s\n", argv[argc-1]);
    return EXIT_FAILURE;
  }

  ouf = fopen(outputfile, "w");

  if (!a1(fp, ouf))
  {
    fprintf(stderr, "ERROR: Can't read grid data in luscus file!\n");
    return EXIT_FAILURE;
  }

  fclose(fp);
  fclose(ouf);
}

