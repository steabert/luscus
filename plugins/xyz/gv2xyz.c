#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define VERSION "1.0.0"

void len_trim(char *str)
{
  int i = strlen(str) - 1;
  while(str[i] == 32 || str[i] == 9) i--;
  str[i] = 0;
}

/*
int line_contains_non_numerical_characters(char *line)
{
  int i = 0;
  if (!line) return 1;
  while(line[i] != 0)
  {
    printf("i = %d\n", i);
    if (line[i] > 20 && line[i] < 48 && line[i] > 57 && line[i] < 127) return 1;
    i++;
  }
  return 0;
}*/

int gv_xyz(char *path)
{
  int i, j;
  FILE *in, *out;
  char line[2048];
  char *ptr;
  int rdgeo, next;
  int natom;

  in = fopen(path, "r");
  if (in == NULL) return 0;

 if (path[strlen(path) - 4] != '.') return 0;

/*  ptr = (char*) malloc(sizeof(char) * (strlen(path) + 2));
  strcpy(ptr, path);
  j = strlen(ptr);

  for(i = 0; i < 3; i++) ptr[j-2+i] = 120+i;
  ptr[j+1] = 0;*/

  ptr = strdup(path);
  for(i = 0; i < 3; i++) ptr[strlen(ptr)-3+i] = 120+i;

  out = fopen(ptr, "w");

  if (out == NULL) return 0;

  rdgeo = 1;
  while(!feof(in))
  {
    fgets(line, sizeof(line), in);
    if (feof(in)) break;

    len_trim(line);
    if (line[0] == 0) next = 0;
    else next = 1;
    if (rdgeo)
    {
      if (next) natom = atoi(line);
      if (next)
      {
        fprintf(out, "%s\n", line);
        for(i = 0; i <= natom; i++)
        {
          fgets(line, sizeof(line), in);
          len_trim(line);
          fprintf(out, "%s\n", line);
        }
      }

      rdgeo = 0;
    }

    if (strstr(line, "<END>")) rdgeo = 1;
  }

  return 1;
}

int main(int argc, char **argv)
{
  if (argc < 2) return EXIT_FAILURE;
  if (!strcmp(argv[1], "--version"))
  {
    printf("gv2xyz version %s\n", VERSION);
    return EXIT_SUCCESS;
  }
  if (gv_xyz(argv[1])) return EXIT_SUCCESS;
  else return EXIT_FAILURE;
}

