#define VERSION "1.0.0"
char *read_line(FILE*);
char *find_substring_case(char*, const char*);
char *my_read_str_value(char*);
int my_read_int_value(char*);
double my_read_dble_value(char*);
/*char *strdiet(char*);*/
/*char *get_one_word(char*);*/
/*char *get_ptr_ith_word(char*, int);*/
/*int line_is_empty(char*);*/


typedef struct grid_t
{
  int ngrid;
  int ngrid2;
  int nmo;
  int nmo2;
  int npoints;
  int ptotal;
  int npt[3]; /*number of points in 3 dimenstions*/
  int block_size;
  int nBlock;
  int isCutOff;
  double CutOff;
  int nPointsCutOff;
  int *ByteCutOff;
  long *pos;
  long *pos2;
  int *dependence;
  double len[3];
  double origin[3];
  double axisvec1[3], axisvec2[3], axisvec3[3];
  int iniIndex[7];
  char *title;
  char *iType;
  char Gtitle[81];
  double minval, maxval, guess;
  double currentvalue;
  double isod[21];
  double *values;
  int sel[4];
  int is_coord;
/* filter */
  char FltSym[16];
  char FltInd[8];
  double FltMinOcc, FltMaxOcc, FltMinEne, FltMaxEne;
  char **titlesArr;
  int *FltTitle;
  int isFltSym,isFltOcc,isFltEne,isFltInd;
} GRID_T;

void open_luscus_file(char*, GRID_T*);
void allocate_grids(GRID_T*);
void deallocate_grids(GRID_T*);


