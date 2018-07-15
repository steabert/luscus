#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#define MAXCHAR 1024
#define VERSION "1.0.0"

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
/*    if (strcasestr(line, "<GRID>")) next = 0;*/
    if (find_substring_case(line, "<GRID>")) next = 0;
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

  next = 1;
  while(next)
  {
    fgets(line, MAXCHAR, fp);
    if (feof(fp)) return 0;
/*    if (strcasestr(line, "<DENSITY>")) next = 0;*/
    if (find_substring_case(line, "<DENSITY>")) next = 0;
  }

  fseek(fp, orbend, SEEK_CUR);

  next = 1;

  while(next)
  {
    fgets(line, MAXCHAR, fp);
    if (feof(fp)) return 0;
/*    if (strcasestr(line, "<INPORB>")) next = 0;*/
    if (find_substring_case(line, "<INPORB>")) next = 0;
  }
  fgets(line, MAXCHAR, fp);/*read one lines and ignore them!*/

  next = 1;

  while(next)
  {
    fgets(line, MAXCHAR, fp);
    if (find_substring_case(line, "</INPORB>")) next = 0;
    else
    {
       fprintf(ouf, "%s", line);
    }
    if (feof(fp)) return 0;
/*    if (strcasestr(line, "</INPORB>")) next = 0;*/
  }

  return 1;
}

int main(int argc, char **argv)
{
  FILE *fp;
  FILE *ouf;
  char outputfile[MAXCHAR];
  char *ptr;

  if (argc < 2) return EXIT_FAILURE;

  if (!strcmp(argv[1], "--version"))
  {
    printf("lus2gvorb version %s\n", VERSION);
    return EXIT_SUCCESS;
  }

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

