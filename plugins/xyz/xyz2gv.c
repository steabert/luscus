#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define VERSION "1.0.0"

void len_trim(char *str)
{
  int i = strlen(str) - 1;
  while(i >= 0)
  {
    if (str[i] == 32 || str[i] == 9 || str[i] == 10) str[i] = 0;
    else return;
    i--;
  }
}

int xyz2gv(char *path)
{
  FILE *in, *out;
  char *ptr;
  char line[2048];
  int natom;
  int i;
  int next = 1;
  int notfirst = 0;

  if (path == NULL) return 0;
  ptr = strrchr(path, '.');
  if (ptr == NULL) return 0;
  
  if (strcmp(ptr+1, "xyz") != 0) return 0;

  in = fopen(path, "r");
  if (in == NULL) return 0;

  i = strlen(path);
  ptr = (char*) malloc(sizeof(char)*(i+1));
  strncpy(ptr, path, i);
  ptr[i-3] = 'l';
  ptr[i-2] = 'u';
  ptr[i-1] = 's';
  ptr[i] = 0;

  out = fopen(ptr, "w");
  if (out == NULL) return 0;

  while(!feof(in) && next)
  {
    fgets(line, sizeof(line), in);
    len_trim(line);
    if (line[0] == 0) next = 0;
    if (next) next = sscanf(line, "%d", &natom);
    if (next)
    {
      if (notfirst) fprintf(out, "<EDITABLE>\n NO\n</EDITABLE>\n<END>\n");
      fprintf(out, "  %d\n", natom);
      fgets(line, sizeof(line), in);
      len_trim(line);
      fprintf(out, "%s\n", line);
      for(i = 0; i < natom; i++)
      {
        fgets(line, sizeof(line), in);
        len_trim(line);
        fprintf(out, "%s\n", line);
      }
    }
    notfirst = 1;
  }

  return 1;
}

int main(int argc, char **argv)
{
  if (argc < 2) return EXIT_FAILURE;
  if (!strcmp(argv[1],"--version"))
  {
    printf("xyz2gv version %s\n", VERSION);
    return EXIT_SUCCESS;
  }
  if (xyz2gv(argv[1])) return EXIT_SUCCESS;
  else return EXIT_FAILURE;
}

