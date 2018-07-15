#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"mystring.h"
#include"atom.h"

#define BOHR 0.5291772083

#define ENERGY 0
#define OPTIMIZE 1
#define IRC 2
#define HESSIAN 3
#define SADPOINT 4
#define GRADIENT 5

#define min(x,y) x < y ? x : y

int calc;
int natom;

char *read_line(FILE *in)
{
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
}

char *atsym_from_num(int atnum)
{
  if (atnum < 0 || atnum >= NATOMS) return "  ";
  else return atom_symbols[atnum - 1];
}

void read_geom(FILE *in, FILE *out)
{
  int i;
  char *line;
  int nat = 0;
  double x, y, z;
  int atnum;
  long pos;

  for(i = 0; i < 4; i++)
  {
    line = read_line(in);
    free(line);
    line = NULL;
  }

  pos = ftell(in);

  do
  {
    if (line) free(line);
    line = read_line(in);
    if (!strstr(line,"---------------------------------------")) nat++;
  }
  while(!strstr(line,"---------------------------------------"));
  if (line) free(line);

  fseek(in, pos, SEEK_SET);

  natom = nat;
  fprintf(out,"  %d\n File created by luscus\n", nat);

  for(i = 0; i < nat; i++)
  {
    line = read_line(in);
    atnum = atoi(get_ptr_ith_word(line, 2));
    x = atof(get_ptr_ith_word(line, 4));
    y = atof(get_ptr_ith_word(line, 5));
    z = atof(get_ptr_ith_word(line, 6));
    free(line);
    fprintf(out, " %s  %f %f %f\n", atsym_from_num(atnum), x, y, z);
  }


}

void read_energy(char *line, FILE *out)
{
  double energy;

  energy = atof(strstr(line,"=")+1);
  fprintf(out, "  <ENERGY>\n");
  fprintf(out, "  %f\n", energy);
  fprintf(out, "  </ENERGY>\n");

}

void read_freq(FILE *in, FILE *out)
{
  int i, j, k;
  double *x, *y, *z;
  double freq[3];
  double irint[3], raman_int[3];
 /*missing raman1, raman2, raman3;*/
  int ifreq = 0;
  int nfreq;
  char *line;
  char imaraman = 0;

  x = (double*) malloc(3*natom*sizeof(double));
  y = (double*) malloc(3*natom*sizeof(double));
  z = (double*) malloc(3*natom*sizeof(double));

  for(i = 0; i < 3; i++)
  {
    line = read_line(in);
    free(line);
  }

  for(i = 0; i < natom; i++)
  {
    printf("reading group #%d\n", i);
    line = read_line(in);
    if (line_is_empty(line))
    {
      free(line);
      return;
    }
    nfreq = my_string_get_num_words(line);
    printf("nfreq = %d\n", nfreq);

    while(!strstr(line,"Frequencies"))
    {
      if (line) free(line);
      line = NULL;
      line = read_line(in);
      printf("searching for Frequencies: |%s|\n", line);
    }

    for(j = 0; j < nfreq; j++)
      freq[j] = atof(get_ptr_ith_word(line, 3+j));

    while(!strstr(line,"IR Inten"))
    {
      printf("TU SAM0\n");  fflush(stdout);
      if (line) free(line);
      printf("TU SAM1\n");  fflush(stdout);
      line = NULL;
      printf("TU SAM2\n");  fflush(stdout);
      line = read_line(in);
      printf("searching for IR INT |%s|\n", line); fflush(stdout);
    }

    for(j = 0; j < nfreq; j++)
      irint[j] = atof(get_ptr_ith_word(line, 3+j));

    while(!strstr(line,"Atom  AN"))
    {
      if (line) free(line);
      line = NULL;
      line = read_line(in);
      printf("searching for Atom  AN |%s|\n", line); fflush(stdout);
    }

    if (line) free(line);
    line = NULL;

    for(j = 0; j < natom; j++)
    {
      line = read_line(in);
      for(k = 0; k < nfreq; k++)
      {
        x[j+natom*k] = atof(get_ptr_ith_word(line, 3+3*k));
        y[j+natom*k] = atof(get_ptr_ith_word(line, 4+3*k));
        z[j+natom*k] = atof(get_ptr_ith_word(line, 5+3*k));
      }

      free(line);
    }

    for(j = 0; j < nfreq; j++)
    {
      fprintf(out, "  <VIBRATION>\n");
      fprintf(out, "freq=%f ir_int=%f ", freq[j], irint[j]);
      if (imaraman) fprintf(out, "raman_int=%f", raman_int[j]);
      fputc(10, out);
      for(k = 0; k < natom; k++)
        fprintf(out," %f %f %f\n", x[j*natom+k], y[j*natom+k], z[j*natom+k]);

      fprintf(out, "  </VIBRATION>\n");

    }

  }

  free(x);
  free(y);
  free(z);

}

void analyse_gauss_output(FILE *in, FILE *out)
{
  char *line;
  int natom;
  int notfirstgeo = 0;

  while(!feof(in))
  {
    line = read_line(in);

    if (strstr(line,"Standard orientation:"))
    {
      if (notfirstgeo) fprintf(out, "  <EDITABLE>\n  NO\n  </EDITABLE>\n  <END>\n");
      read_geom(in, out);
      notfirstgeo = 1;
    }
    if (strstr(line,"SCF Done")) read_energy(line, out);
    if (strstr(line,"Harmonic frequencies")) read_freq(in, out);

    free(line);
  }
}

int main(int argc, char **argv)
{  
  char *output;
  char *ptr;
  FILE *in, *out;
  calc = ENERGY;

  if (argc < 2) return EXIT_FAILURE;
  output = (char*) malloc(sizeof(char) * (strlen(argv[argc-1]) + 1));
  strncpy(output, argv[argc-1], strlen(argv[argc-1]));
  ptr = strrchr(output, '.');
  if (ptr == NULL) ptr = output + strlen(argv[argc-1]);
  output = (char*) realloc(output, sizeof(char) * (ptr - output + 5));
  ptr[0] = '.';
  ptr[1] = 'l';
  ptr[2] = 'u';
  ptr[3] = 's';
  ptr[4] = 0;

  in = fopen(argv[argc-1], "r");
  if (in == NULL) return EXIT_FAILURE;

  out = fopen(output, "w");
  if (out == NULL) return EXIT_FAILURE;

  analyse_gauss_output(in, out);

  free(output);
  return EXIT_SUCCESS;
}

