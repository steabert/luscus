#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

double cx=0.0, cy, cz; /*unit cell parameters*/
double ca, cb, cg;

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

void read_alat_coordiantes(int natom, double *x, double *y, double *z, char *atsym0, char *atsym1, FILE *in)
{
  int i;
  char *line;
  line = read_line(in); free(line);
  line = read_line(in); free(line);
  for(i = 0; i < natom; i++)
  {
    printf("READING LINE\n"); fflush(stdout);
    line = read_line(in);

    printf("TEST read cartesian |%s|\n", line); fflush(stdout);

    atsym0[i] = line[21];
    atsym1[i] = line[22];
    printf("i = %d line_x = |%s|\n", i, line+38); fflush(stdout);   
    x[i] = atof(line+38);
    y[i] = atof(line+50);
    z[i] = atof(line+62);
    printf("TU SAM\n"); fflush(stdout);
    free(line);
  }
  printf("x = %f %f %f\n", x[0], x[1], x[2]);
  printf("y = %f %f %f\n", y[0], y[1], y[2]);
  printf("z = %f %f %f\n", z[0], z[1], z[2]);

}

void print_xyz(int natom, double *x, double *y, double *z, char *atsym0, char *atsym1, FILE *out)
{
  int i;
  fprintf(out, "  %d\nfile created by luscus\n", natom);
  for(i = 0; i < natom; i++)
    fprintf(out, " %c%c  %f  %f  %f\n", atsym0[i], atsym1[i], x[i], y[i], z[i]);
}

void convert_alat_2_xyz(int natom, double *x, double *y, double *z)
{
  int i;
  double tmpx, tmpy, tmpz;
  double /*alpha, beta,*/ gamma;
  double /*sa, sb,*/ sg;
  double v;
/*  alpha = acos(ca);
  beta = acos(cb);*/
  gamma = acos(cg);
/*  sa = sin(alpha);
  sb = sin(beta);*/
  sg = sin(gamma);
  v = sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg);

  for(i = 0; i < natom; i++)
  {
    tmpx = ca * x[i] + cy * y[i] * cg + cz * z[i] * cb;
    tmpy = cy * y[i] * sg + cz * (ca - cb*cg) / sg;
    tmpz = cz * z[i] * v / sg;

    x[i] = tmpx;
    y[i] = tmpy;
    z[i] = tmpz;
  }
}

void analyse_QE_output(FILE *in, FILE *out)
{
  int i;
  char *line;

  int natom = 0;
  double energy;
  double *x, *y, *z;
  double tot_force;
  double force[3];
  int firstgeo = 1;
  char *atsym0, *atsym1;

  while(!feof(in))
  {
    line = read_line(in);
    printf("TEST checking line |%s|\n", line); fflush(stdout);
/*    if (strstr(line,"lattice parameter")) cx = my_read_dble_value(line);*/
    if (strstr(line,"number of atoms/cell")) natom = my_read_int_value(line);
    if (strstr(line,"celldm(1)="))
    {
      cx = my_read_dble_value(line);
      cx *= 0.5291772083;
      cy = cx * my_read_dble_value(strstr(line,"celldm(2)="));
      cz = cx * my_read_dble_value(strstr(line,"celldm(3)="));
    }
    if (strstr(line,"celldm(4)="))
    {
      ca = my_read_dble_value(strstr(line,"celldm(4)="));
      cb = my_read_dble_value(strstr(line,"celldm(5)="));
      cg = my_read_dble_value(strstr(line,"celldm(6)="));
    }
    if (strstr(line,"total energy              =")) energy = my_read_dble_value(line);
    if (strstr(line,"Cartesian axes"))
    {
      printf("ALLOCATING MEMORY FOR THE GEOMETRY: NATOM = %d\n", natom); fflush(stdout);
      if (firstgeo)
      {
        x = (double*) malloc(natom * sizeof(double));
        y = (double*) malloc(natom * sizeof(double));
        z = (double*) malloc(natom * sizeof(double));
        atsym0 = (char*) malloc(natom * sizeof(char));
        atsym1 = (char*) malloc(natom * sizeof(char));
        firstgeo = 0;
      }
      read_alat_coordiantes(natom, x, y, z, atsym0, atsym1, in);
      printf("x = %f %f %f\n", x[0], x[1], x[2]);
      convert_alat_2_xyz(natom, x, y, z);
      printf("x = %f %f %f\n", x[0], x[1], x[2]);
      print_xyz(natom, x, y, z, atsym0, atsym1, out);
    }
    free(line);
    
    
  }
}

int main(int argc, char **argv)
{
  char *output;
  char *ptr;
  FILE *in, *out;

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

  analyse_QE_output(in, out);

  fclose(in);
  fclose(out);
  free(output);
  return EXIT_SUCCESS;
}

