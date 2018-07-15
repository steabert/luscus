#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"mystring.h"
#include"atom.h"

#define VERSION "1.0.0"

#define BOHR 0.5291772083

#define ENERGY 0
#define OPTIMIZE 1
#define IRC 2
#define HESSIAN 3
#define SADPOINT 4
#define GRADIENT 5

#define U_ANGSTR 0
#define U_BOHR 1

int units;
double *x=NULL, *y=NULL, *z=NULL;
int natom = 0;

char *atsym_from_num(int atnum)
{
  if (atnum < 0 || atnum >= NATOMS) return "  ";
  else return atom_symbols[atnum - 1];
}

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

void rd_0th_geom(FILE *in, FILE *out)
{
  int i;
/*  int natom = 0;*/
  char *line;
  char **atom_name = NULL;
  int *atom_number = NULL;
/*  double *x=NULL, *y=NULL, *z=NULL;*/
  char *tmp;

  rewind(in);

  line = read_line(in);

  while(strstr(line,"COORDINATES (BOHR)") == NULL)
  {
    free(line);
    line = read_line(in);
  }

  free(line);
  line = read_line(in);

  while(strlen(line) > 5)
  {
    free(line);
    line = read_line(in);
    if (strlen(line) > 5)
    {
      natom++;
      atom_name = (char**) realloc(atom_name, sizeof(char*) * natom);
      atom_number = (int*) realloc(atom_number, sizeof(int) * natom);
      x = (double*) realloc(x, sizeof(double) * natom);
      y = (double*) realloc(y, sizeof(double) * natom);
      z = (double*) realloc(z, sizeof(double) * natom);

      atom_name[natom-1] = get_one_word(get_ptr_ith_word(line, 1));
      tmp = get_one_word(get_ptr_ith_word(line, 2));
      atom_number[natom-1] = atoi(tmp);
      free(tmp);
      tmp = get_one_word(get_ptr_ith_word(line, 3));
      x[natom-1] = atof(tmp) * BOHR;
      free(tmp);
      tmp = get_one_word(get_ptr_ith_word(line, 4));
      y[natom-1] = atof(tmp) * BOHR;
      free(tmp);
      tmp = get_one_word(get_ptr_ith_word(line, 5));
      z[natom-1] = atof(tmp) * BOHR;
      free(tmp);
    }
  }

  fprintf(out, "  %d\n\n", natom);
  for(i = 0; i < natom; i++)
    fprintf(out, " %s  %f  %f  %f\n", atsym_from_num(atom_number[i]), x[i], y[i], z[i]);

  fprintf(out, " <ATOM>\n");

  for(i = 0; i < natom; i++)
    fprintf(out, "name = %s\n", atom_name[i]);

  fprintf(out, " </ATOM>\n");

  for(i = 0; i < natom; i++) free(atom_name[i]);
  free(atom_name);
  free(atom_number);
/*  free(x);
  free(y);
  free(z);*/


  free(line);
}

void rd_dipoles(FILE *in, FILE *out)
{
  char *line;
  double d;

  line = read_line(in);
  free(line);

  line = read_line(in);
  free(line);

  line = read_line(in);
  free(line);

  line = read_line(in);
  free(line);

  line = read_line(in);
  free(line);

  line = read_line(in);
  free(line);

 /*!!! should check if dipole is > 0*/
  fprintf(out, " <DIPOLE>\n");

  d = atof(get_ptr_ith_word(line, 1));
  fprintf(out, "  %f", d);

  d = atof(get_ptr_ith_word(line, 1));
  fprintf(out, "  %f", d);

  d = atof(get_ptr_ith_word(line, 1));
  fprintf(out, "  %f\n", d);

  fprintf(out, " </DIPOLE>\n");
}

void read_gradient(FILE *in, FILE *out, int natom)
{
  char *line;
  int next = 1;
  int i, j;
  double xyz[3];
/*  double len;*/

  if (!x) return;
  if (!y) return;
  if (!z) return;

  while(next)
  {
    line = read_line(in);
    if (strstr(line, "UNITS ARE HARTREE/BOHR")) next = 0;
    free(line);
  }

  for (i = 0; i < natom; i++)
  {
    line = read_line(in);
    fprintf(out, " <VECTOR>\n RED = 0.700  GREEN = 0.100 BLUE = 0.100 RADIUS = 0.020 SHARPNESS = 0.200 TRANSPARENCY=1.000\n");

    for(j = 0; j < 3; j++)
    {
      xyz[j] = atof(get_ptr_ith_word(line, 3+j));
      if (units == U_BOHR) xyz[j] /= BOHR;
    }
/*    len = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);*/
    /*asuming 0.3 as smallest atomic radius*/
    /*scale factor 5 applied!*/
    fprintf(out, " %f %f %f\n %f %f %f\n </VECTOR>\n", x[i], y[i], z[i], x[i]-xyz[0]*65.0, y[i]-xyz[1]*65.0, z[i]-xyz[2]*65.0);
  }
}

void rd_statpt(FILE *in, FILE *out, int natom)
{
  char *line;
  int next = 1;
  int i, j;
  char **atname;
  int atbr;
  double xyz[3];

  while(next)
  {
    line = read_line(in);
    if (strstr(line, "COORDINATES OF ALL ATOMS ARE")) next = 0;
    free(line);
  }

  atname = (char**) malloc(sizeof(char*) * natom);

  fprintf(out,"  %d\n\n", natom);

  line = read_line(in);
  free(line);

  line = read_line(in);
  free(line);

  for(i = 0; i < natom; i++)
  {
    line = read_line(in);
    atname[i] = get_one_word(get_ptr_ith_word(line, 1));
    atbr = atoi(get_ptr_ith_word(line, 2));
    for(j = 0; j < 3; j++)
    {
      xyz[j] = atof(get_ptr_ith_word(line, 3+j));
      if (units == U_BOHR) xyz[j] *= BOHR;
    }
    free(line);
    fprintf(out, " %s  %f %f %f\n", atsym_from_num(atbr), xyz[0], xyz[1], xyz[2]);
  }

  fprintf(out,"\n  <ATOM>\n");
  for(i = 0; i < natom; i++)
  {
    fprintf(out, "  name=%s\n", atname[i]);
  }
  fprintf(out,"  </ATOM>\n");

  free(atname);


}

void rd_irc(FILE *in, FILE *out, int natom)
{
  int i, j;
  double distance, energy;
  char *line;

  char **atname;
  int atbr;
  double xyz[3];

  atname = (char**) malloc(sizeof(char*) * natom);

  fprintf(out, "  %d\n\n", natom);
  line = read_line(in);
  free(line);
  distance = atof(line+37);
  line = read_line(in);
  free(line);
  line = read_line(in);
  energy = atof(line+37);
  free(line);

  for(i = 0; i < 5; i++)
  {
    line = read_line(in);
    free(line);
  }
  for(i = 0; i < natom; i++)
  {
    line = read_line(in);

    atname[i] = get_one_word(get_ptr_ith_word(line, 1));
    atbr = atoi(get_ptr_ith_word(line, 2));
    for(j = 0; j < 3; j++)
    {
      xyz[j] = atof(get_ptr_ith_word(line, 3+j));
      if (units == U_BOHR) xyz[j] *= BOHR;
    }
    fprintf(out, " %s  %f %f %f\n", atsym_from_num(atbr), xyz[0], xyz[1], xyz[2]);


    free(line);
  }

  fprintf(out, "  <ENERGY>\n");
  fprintf(out, "  %f\n", energy);
  fprintf(out, "  </ENERGY>\n");

  fprintf(out, "<DISTANCE>\n");
  fprintf(out, "  %f\n", distance);
  fprintf(out, "</DISTANCE>\n");


  fprintf(out, " <ATOM>\n");

  for(i = 0; i < natom; i++)
    fprintf(out, "name = %s\n", atname[i]);

  fprintf(out, " </ATOM>\n  <EDITABLE>\n  NO\n  </EDITABLE> <END>\n");



  free(atname);

}

void rd_freqs(FILE *in, FILE *out, int natom)
{
  int i, j, k;
  char *line = NULL;
  double freq[5];
  double **normod_x, **normod_y, **normod_z;
  double ir_int[5], raman_int[5];
  int nfields;
  char imaraman = 0;

  normod_x = (double**) malloc(natom * sizeof(double*));
  normod_y = (double**) malloc(natom * sizeof(double*));
  normod_z = (double**) malloc(natom * sizeof(double*));

  for(i = 0; i < natom; i++)
  {
    normod_x[i] = (double*) malloc(5 * sizeof(double));
    normod_y[i] = (double*) malloc(5 * sizeof(double));
    normod_z[i] = (double*) malloc(5 * sizeof(double));
  }

  do
  {
    line = read_line(in);
    free(line);
  }
  while(strstr(line,"FREQUENCIES IN CM**-1") == NULL);

  line = read_line(in);
  free(line);
  line = read_line(in);
  free(line);

  for(k = 0; k < 3*natom/5; k++)
  {
    line = read_line(in);
    nfields = 0;
    if (line == NULL) return;
    if (line[0] != 0)
      for(i = 1; i < strlen(line); i++)
        if (line[i-1] == 32 && line[i] != 32) nfields++;
    free(line);
 
    line = read_line(in);
    if (strstr(line,"FREQUENCY"))
      for(i = 0; i < nfields; i++)
        freq[i] = atof(line+22+12*i);
    free(line);
 
    line = read_line(in);
    free(line);
    line = read_line(in);
    free(line);
 
    line = read_line(in);
    if (strstr(line,"IR INTENSITY"))
      for(i = 0; i < nfields; i++)
        ir_int[i] = atof(line+22+12*i);
    free(line);
 
    line = read_line(in);
    if (strstr(line,"RAMAN INTENSITY"))
    {
      imaraman = 1;
      for(i = 0; i < nfields; i++)
        ir_int[i] = atof(line+22+12*i);
      free(line); /*one extra line*/
      read_line(in);
    }
    free(line);
 
    for(i = 0; i < natom; i++)
    {
      line = read_line(in);
      for(j = 0; j < nfields; j++)
        normod_x[i][j] = atof(line+22+12*j);
      free(line);
      line = read_line(in);
      for(j = 0; j < nfields; j++)
        normod_y[i][j] = atof(line+22+12*j);
      free(line);
      line = read_line(in);
      for(j = 0; j < nfields; j++)
        normod_z[i][j] = atof(line+22+12*j);
      free(line);
    }
 
    for(i = 0; i < 11; i++)
    {
      line = read_line(in);
      free(line);
    }

    for(j = 0; j < nfields; j++)
    {
      fprintf(out, "  <VIBRATION>\n");
      fprintf(out, "freq=%f ir_int=%f ", freq[j], ir_int[j]);
      if (imaraman) fprintf(out, "raman_int=%f", raman_int[j]);
      fputc(10, out);
      for(i = 0; i < natom; i++)
	fprintf(out," %f %f %f\n", normod_x[i][j], normod_y[i][j], normod_z[i][j]);

      fprintf(out, "  </VIBRATION>\n");
    }
  }

  for(i = 0; i < natom; i++)
  {
    free(normod_x[i]);
    free(normod_y[i]);
    free(normod_z[i]);
  }
  free(normod_x);
  free(normod_y);
  free(normod_z);
  /*end block*/
}

void skip_initial_hessian(FILE *in)
{
  char *line;
  int next = 1;

  while(next)
  {
    line = read_line(in);
    if (feof(in)) next = 0;
    if (strstr(line, "END OF NUMERICAL HESSIAN CALCULATION")) next = 0;
    if (strstr(line,"DONE WITH CPHF CONTRIBUTIONS")) next = 0;

    free(line);
  }
}

void analyse_gamess_output(FILE *in, FILE *out)
{
  int units;
  int calc;

  char *line;

  int natom = 0;
  double energy;
  double xyz[3];
  char notfirstgeo = 0;
  long pos;

  units = U_ANGSTR;

  while(!feof(in))
  {
    line = read_line(in);
/*    if (strstr(line,"COORDINATES (BOHR)")) rd_0th_geom(in, out);*/
    if (strstr(line,"TOTAL NUMBER OF ATOMS")) natom = my_read_int_value(line);
    if (strstr(line,"    RUNTYP="))
    {
      if (strstr(line,"RUNTYP=ENERGY"))
      {
        pos = ftell(in);
        rd_0th_geom(in, out);
        calc = ENERGY;
        fseek(in, pos, SEEK_SET);
      }
      if (strstr(line,"RUNTYP=GRADIENT"))
      {
        pos = ftell(in);
        rd_0th_geom(in, out);
        calc = GRADIENT;
        fseek(in, pos, SEEK_SET);
      }
      if (strstr(line,"RUNTYP=OPTIMIZE")) calc = OPTIMIZE;
      if (strstr(line,"RUNTYP=IRC")) calc = IRC;
      if (strstr(line,"RUNTYP=HESSIAN"))
      {
        pos = ftell(in);
        rd_0th_geom(in, out);
        calc = HESSIAN;
        fseek(in, pos, SEEK_SET);
      }
      if (strstr(line,"RUNTYP=SADPOINT")) calc = SADPOINT;
    }
    if (strstr(line, "THERE ARE ATOMS LESS THAN"))
    {
      pos = ftell(in);
      rd_0th_geom(in, out);
      calc = ENERGY;
      fseek(in, pos, SEEK_SET);
/*----------*/
    }
    if (strstr(line, "UNITS ="))
    {
      if (strstr(line, "ANGS")) units=U_ANGSTR;
      if (strstr(line, "BOHR")) units=U_BOHR;
    }
    if (calc != IRC)
      if (strstr(line, "FINAL") && strstr(line,"ENERGY IS"))
      {
        energy = atof(line + 24);
        fprintf(out, "  <ENERGY>\n");
        fprintf(out, "  %f\n", energy);
        fprintf(out, "  </ENERGY>\n");
      }
    if (strstr(line,"ELECTROSTATIC MOMENTS")) rd_dipoles(in, out);
    if (calc == OPTIMIZE || calc == SADPOINT)
    {
      if (strstr(line,"BEGINNING GEOMETRY SEARCH POINT"))
      {
        if (notfirstgeo) fprintf(out, "<EDITABLE>\n  NO\n  </EDITABLE> <END>\n");
        rd_statpt(in, out, natom);
        notfirstgeo = 1;
      }
      if (strstr(line,"PROCEEDING DIRECTLY TO THE GEOMETRY SEARCH"))
      {
        if (notfirstgeo) fprintf(out, "<EDITABLE>\n  NO\n  </EDITABLE> <END>\n");
        rd_statpt(in, out, natom);
        notfirstgeo = 1;
      }

    }
    if (calc == GRADIENT)
      if (strstr(line,"GRADIENT OF THE ENERGY"))
	read_gradient(in, out, natom);

    if (strstr(line,"MAXIMUM GRADIENT ="))
      fprintf(out, "\n  <MAX_GRAD>\n  %f\n  </MAX_GRAD>\n", my_read_dble_value(strstr(line,"MAXIMUM GRADIENT =")));
    if (strstr(line,"RMS GRADIENT ="))
      fprintf(out, "\n  <RMS_GRAD>\n  %f\n  </RMS_GRAD>\n", my_read_dble_value(strstr(line,"RMS GRADIENT =")));

    if (strstr(line,"OBTAINING INITIAL HESSIAN, HESS=CALC"))
      skip_initial_hessian(in);

    if (calc == HESSIAN)
      if (strstr(line, "ANALYZING SYMMETRY OF NORMAL MODES"))
       	rd_freqs(in, out, natom);

    if (calc == IRC)
    {
      if (strstr(line, "ON THE REACTION PATH"))
	rd_irc(in, out, natom);
    }
    free(line);
  }
  if (x) free(x);
  if (y) free(y);
  if (z) free(z);

}

int main(int argc, char **argv)
{  
  char *output;
  char *ptr;
  FILE *in, *out;

  if (argc < 2) return EXIT_FAILURE;
  if (!strcmp(argv[1],"--version"))
  {
    printf("gamess2lus version %s\n", VERSION);
    return EXIT_SUCCESS;
  }


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

  analyse_gamess_output(in, out);

  free(output);
  return EXIT_SUCCESS;
}

