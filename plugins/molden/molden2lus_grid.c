#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define CARTESIAN 1
#define BOHR 2

#define BOHR2CART 0.529177249
#define MAXLEN 1024
#define MAXELEM 104
char *element[] = 
{
" ",  
"H",                                                            "He", 
"Li", "Be",                                                     "B",  "C",  "N",  "O",  "F",  "Ne", 
"Na", "Mg",                                                     "Al", "Si", "P",  "S",  "Cl", "Ar", 
"K",  "Ca", "Sc","Ti","V","Cr", "Mn", "Fe","Co","Ni","Cu","Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
"Rb", "Sr", "Y", "Zr","Nb","Mo","Tc","Ru", "Rh","Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I","Xe",
"Cs", "Ba",
   "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
           "Hf","Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
"Fr", "Ra", 
   "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
};

double minx[3], maxx[3];
double *x_ = NULL, *y_ = NULL, *z_ = NULL;

int has_geoconv = 0;

void calculate_epot_grid(int natom, double *charge, FILE *pout)
{
  int i, j, k, l;
  double mminx[3], mmaxx[3], dv[3], dx[3];
  double cutoff = 1.0e-5;
  int npt[3];
  long file_pointers, orboff;
  FILE *tmp;
  char c;
  char tmpfilename[65];
  time_t ttm;
  double psi; /*some stupid wavefunction*/
  double r2;
  
  ttm = time(NULL);
  srand(ttm);
  
  for(i = 0; i < 64; i++)
  {
    do
    {
      c = rand()%74+48;
    }
    while((c < 65 && c > 57) || (c < 97 && c > 90));
    tmpfilename[i] = c;
  }

  tmpfilename[i] = 0;
  printf("tmpfilename = |%s|\n", tmpfilename);

  tmp=fopen(tmpfilename, "rb");

#ifdef EBUG
  return;
#endif

  for (i = 0; i < 3; i++)
  {
    mminx[i] = minx[i] - 5.0e-1 * (maxx[i] - minx[i]);
    mmaxx[i] = maxx[i] + 5.0e-1 * (maxx[i] - minx[i]);
    printf("i = %d mminx=%f mmaxx=%f\n", i, mminx[i], mmaxx[i]);
    npt[i] = (int) (mmaxx[i] - mminx[i]) / 6.0;
    if (npt[i] < 30) npt[i] = 30;
    dv[i] = (mmaxx[i] - mminx[i]) / (double) npt[i];
  }

  /*print potential to the temporary file*/
  for(i = 0; i < npt[0]+1; i++)
  {
    dx[0] = mminx[0] + dv[0] * (double) i;
    for(j = 0; j < npt[1]+1; j++)
    {
      dy[1] = mminx[1] + dv[1] * (double) j;
      for(k = 0; k < npt[2]+1; k++)
      {
        dy[2] = mminx[2] + dv[2] * (double) z;
        psi = 0.0;
        for(l = 0; l < natom; l++)
        {
          r2 =((dx[0] - x_[i])*(dx[0] - x_[i]) + 
               (dx[1] - y_[i])*(dx[1] - y_[i]) +
               (dx[2] - z_[i])*(dx[2] - z_[i]));
          psi += exp(-r2);
        }
        printf("%f", );
      }
    }
  }




  /*print the grid section of the lus file*/
  fprintf(pout, "<GRID>\n");

  fprintf(pout, " N_of_MO=1");
  fprintf(pout, " N_of_Grids=1");
  fprintf(pout, " N_of_Points=27000");/*Not needed!???*/
  fprintf(pout, " Block_Size=%d", 1);/*???*/
  fprintf(pout, " N_Blocks=%d", 1);/*???*/ /*Not needed!???*/
  fprintf(pout, " Is_cutoff=0");
  fprintf(pout, " CutOff=%8.4f", cutoff);
  fprintf(pout, " N_P=%d", 50);
  fputc(10, pout);
  fprintf(pout, " N_INDEX= %d %d %d %d %d %d %d\n", 0, 0, 0, 0, 0, 0, 0);
  fprintf(pout, " Net= %d %d %d\n", npt[0], npt[1], npt[2]);
  fprintf(pout, " Origin= %lf %lf %lf\n", mminx[0], mminx[1], mminx[2]);
  fprintf(pout, " Axis_1= %lf %lf %lf\n", (float) mmaxx[0]-mminx[0], 0.0, 0.0);
  fprintf(pout, " Axis_2= %lf %lf %lf\n", 0.0, (float) mmaxx[1]-mminx[1], 0.0);
  fprintf(pout, " Axis_3= %lf %lf %lf\n", 0.0, 0.0, (float) mmaxx[2]-mminx[2]);

  fprintf(pout, " File_pointers = ");
  for(i = 0; i < 3; i++)
  {
    fprintf(pout, " %ld", (long) 1); /*...*/
  }
  fputc(10, pout);
  fprintf(pout, " ORBOFF = %ld", orboff);
  fprintf(pout, " Grid Name=   Elecrostatic potential\n");

  fprintf(pout, "</GRID>\n");

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

void read_array(FILE *pin, int ngeom, double *array)
{
  int i;
  char *line;

  for(i = 0; i < ngeom; i++)
  {
    line = read_line(pin);
    array[i] = atof(line);
    free(line);
  }
  
}

void read_geoconv(FILE *pin, int ngeom, double *energy, double *rms_grad, double *max_grad)
{
  int i;
  char *line;
  has_geoconv = 1;
  for(i = 0; i < 3; i++)
  {
    line = read_line(pin);
    if (strstr(line, "energy")) read_array(pin, ngeom, energy);
    else if (strstr(line, "max-force")) read_array(pin, ngeom, rms_grad);
    else if (strstr(line, "rms-force")) read_array(pin, ngeom, max_grad);
    free(line);
  }

}

void read_geometries(FILE *pin, FILE *pout, int ngeom, double *energy, double *rms_grad, double *max_grad, int iunit)
{
  int i, j;
  int natom;
  char *line;
  double c_unit = 1.0;
  char atsym[6];
  double x, y, z;
  if (iunit == BOHR) c_unit = BOHR2CART;

  minx[0] = minx[1] = minx[2] = 9.99e+99;
  maxx[0] = maxx[1] = maxx[2] =-9.99e+99;

  for(i = 0; i < ngeom; i++)
  {
    line = read_line(pin);
    printf("first_line = |%s|\n", line);
    natom = atoi(line);
    fprintf(pout, "  %d\n", natom - 1);
    free(line);
    line = read_line(pin);
    fprintf(pout, "%s\n", line);
    free(line);
    for(j = 0; j < natom - 1; j++)
    {
      line = read_line(pin);
      sscanf(line, "%s %lf %lf %lf", atsym, &x, &y, &z);
      fprintf(pout, " %s   %f   %f   %f \n", atsym, x * BOHR2CART, y * BOHR2CART, z * BOHR2CART);
      if (x * BOHR2CART < minx[0]) minx[0] = x * BOHR2CART;
      if (x * BOHR2CART > maxx[0]) maxx[0] = x * BOHR2CART;
      if (y * BOHR2CART < minx[1]) minx[1] = y * BOHR2CART;
      if (y * BOHR2CART > maxx[1]) maxx[1] = y * BOHR2CART;
      if (z * BOHR2CART < minx[2]) minx[2] = z * BOHR2CART;
      if (z * BOHR2CART > maxx[2]) maxx[2] = z * BOHR2CART;
/*      fprintf(pout, "%s\n", line);*/
      free(line);
    }
    if (has_geoconv)
    {
      fprintf(pout, " <ENERGY>\n %f\n </ENERGY>\n", energy[i]);
      fprintf(pout, " <RMS_GRAD>\n %f\n </RMS_GRAD>\n", rms_grad[i]);
      fprintf(pout, " <MAX_GRAD>\n %f\n </MAX_GRAD>\n", max_grad[i]);
    }
    fprintf(pout, "<END>\n");
    line = read_line(pin);
    free(line);
  }
}

void read_geometry(FILE *pin, FILE *pout, int natom, int cunit)
{
  int i;
  char *line;
  char **name;
  int *nnum, atbr;
  double x, y, z;
  double unitc = 1.0;

  if (cunit == BOHR) unitc = BOHR2CART;

  name = (char**) malloc(sizeof(char*) * natom);
  for(i = 0; i < natom; i++)
    name[i] = (char*) malloc(sizeof(char) * natom);
  nnum = (int*) malloc(sizeof(int) * natom);
  if (x_ == NULL) x_ = (double) malloc(natom * sizeof(double));
  if (y_ == NULL) y_ = (double) malloc(natom * sizeof(double));
  if (z_ == NULL) z_ = (double) malloc(natom * sizeof(double));

  fprintf(pout, "  %d\n\n", natom);
  for(i = 0; i < natom; i++)
  {RDGEOM
    line = read_line(pin);
    sscanf(line, "%s %d %d %lf %lf %lf", name[i], &nnum[i], &atbr, &x, &y, &z);
    x_[i] = x;
    y_[i] = y;
    z_[i] = z;
    if (atbr > MAXELEM || atbr < 0)
      fprintf(pout, " %s  %f  %f  %f\n", "Q ", x * unitc, y * unitc, z * unitc);
    else
      fprintf(pout, " %s  %f  %f  %f\n", element[atbr], x * unitc, y * unitc, z * unitc);
    free(line);

    if (x * unitc < minx[0]) minx[0] = x * unitc;
    if (x * unitc > maxx[0]) maxx[0] = x * unitc;
    if (y * unitc < minx[1]) minx[1] = y * unitc;
    if (y * unitc > maxx[1]) maxx[1] = y * unitc;
    if (z * unitc < minx[2]) minx[2] = z * unitc;
    if (z * unitc > maxx[2]) maxx[2] = z * unitc;
  }
  fprintf(pout, "<ATOM>\n");
  for(i = 0; i < natom; i++)
  {
    fprintf(pout, " NAME=%s  NUMBER=%d\n", name[i], nnum[i]);
  }
  fprintf(pout, "</ATOM>\n");


  for(i = 0; i < natom; i++) free(name[i]);
  free(name);
  free(nnum);
}

void read_geom1(FILE *pin, FILE *pout, int natom)
{
  int i;
  char *line;
  char atsym[6];
  double x, y, z;

  fprintf(pout, "  %d\n\n", natom);
  for(i = 0; i < natom; i++)
  {
    line = read_line(pin);
    sscanf(line, "%s %lf %lf %lf", atsym, &x, &y, &z);
    fprintf(pout, " %s   %f   %f   %f \n", atsym, x * BOHR2CART, y * BOHR2CART, z * BOHR2CART);
    free(line);
  }
}

void read_charge(FILE *pin, FILE *pout, int natom, int is_mulliken)
{
  int i;
  double *charge;
  char *line;

  charge = (double*) malloc(natom * sizeof(double));

  fprintf(pout, "<ATOM>\n");
  for(i = 0; i < natom; i++)
  {
    line = read_line(pin);
    charge[i] = atof(line);
    if (is_mulliken) fprintf(pout, "  MULLIKEN_CHARGE = ");
    else fprintf(pout, "  LOPROP_CHARGE = ");
    fprintf(pout, "%f\n", charge[i]);
    free(line);
  }
  fprintf(pout, "</ATOM>\n");

  calculate_epot_grid(natom, charge, pout);

  free(charge);
}

void read_normod(FILE *pin, FILE *pout, int natom, int nfreq, double *freq, double *freq_int)
{
  int i, j;
  char *line;
  double x, y, z;
  for(i = 0; i < nfreq; i++)
  {
    line = read_line(pin);
    free(line);
    fprintf(pout, "<VIBRATION>\n");
    fprintf(pout," FREQ=%f  IR_INT=%f\n", freq[i], freq_int[i]);
    for(j = 0; j < natom; j++)
    {
      line = read_line(pin);
      sscanf(line, " %lf %lf %lf", &x, &y, &z);
      fprintf(pout, "  %f  %f  %f\n", x, y, z);
      free(line);
    }
    fprintf(pout, "</VIBRATION>\n");
  }
}

void convert_molden_to_luscus(FILE *pin, FILE *pout)
{
  char *line;
  int ngeom, nfreq, natom;
  int next = 1;
  double *energy = NULL;
  double *rms_grad = NULL;
  double *max_grad = NULL;
  double *freq = NULL;
  double *freq_int = NULL;

  line = read_line(pin);
  if (strstr(line, "MOLDEN FORMAT") == NULL) next = 0;
  free(line);

  if (!next) return;

  while(!feof(pin))
  {
    line = read_line(pin);
    if (strstr(line, "N_GEO"))
    {
      free(line);
      line = read_line(pin);
      ngeom = atoi(line);
      if (ngeom)
      {
        energy = (double*) malloc(sizeof(double) * ngeom);
        rms_grad = (double*) malloc(sizeof(double) * ngeom);
        max_grad = (double*) malloc(sizeof(double) * ngeom);
      }
    }
    else if (strstr(line, "N_ATOMS") || strstr(line, "NATOM"))
    {
      free(line);
      line = read_line(pin);
      natom = atoi(line);
    }
    else if (strstr(line, "N_FREQ"))
    {
      free(line);
      line = read_line(pin);
      nfreq = atoi(line);
      if (nfreq)
      {
        freq = (double*) malloc(sizeof(double) * nfreq);
        freq_int = (double*) malloc(sizeof(double) * nfreq);
      }
    }
    else if (strstr(line, "GEOCONV"))
      read_geoconv(pin, ngeom, energy, rms_grad, max_grad);
    else if (strstr(line, "GEOMETRIES"))
      if (strstr(line, "XYZ"))
        read_geometries(pin, pout, ngeom, energy, rms_grad, max_grad, CARTESIAN);
      else if (strstr(line, "AU"))
        read_geometries(pin, pout, ngeom, energy, rms_grad, max_grad, BOHR);
      else
        read_geometries(pin, pout, ngeom, energy, rms_grad, max_grad, CARTESIAN);
    else if (strstr(line, "ATOMS"))
      if (strstr(line, "XYZ")) read_geometry(pin, pout, natom, CARTESIAN);
      else if (strstr(line, "AU")) read_geometry(pin, pout, natom, BOHR);
      else read_geometry(pin, pout, natom, CARTESIAN);
    else if (strstr(line, "FR-COORD"))
      read_geom1(pin, pout, natom);
    else if (strstr(line, "CHARGE"))
    {
      if (strstr(line, "MULLIKEN")) read_charge(pin, pout, natom, 1);
      else read_charge(pin, pout, natom, 0);
    }
    else if (strstr(line, "FREQ"))
      read_array(pin, nfreq, freq);
    else if (strstr(line, "INT"))
      read_array(pin, nfreq, freq_int);
    else if (strstr(line, "FR-NORM-COORD"))
      read_normod(pin, pout, natom, nfreq, freq, freq_int);

    free(line);
  }
}

int main(int argc, char **argv)
{
  int i;
  char outname[MAXLEN];
  FILE *pin;
  FILE *pout;
  if (argc < 2) return;

  strncpy(outname, argv[argc-1], MAXLEN-4);
  for(i = strlen(outname); i > 0 && outname[i] != 46; i--);
  if (i == 0) i = strlen(outname);
  strncpy(outname+i, ".lus", 4);
  outname[i+4] = 0;

  pin = fopen(argv[argc-1], "r");
  if (pin == NULL)
  {
    fprintf(stderr, "ERROR: Can't open input file: %s\n", argv[argc-1]);
    return EXIT_FAILURE;
  }

  pout = fopen(outname, "w");
  if (pout == NULL)
  {
    fprintf(stderr, "ERROR: Can't open output file: %s\n", outname);
    return EXIT_FAILURE;
  }
 
  convert_molden_to_luscus(pin, pout);

  return EXIT_SUCCESS;
}

