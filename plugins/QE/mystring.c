#include<stdio.h>
#include<stdlib.h>
#include<string.h>

char *my_read_str_value(char *str)
{
  int i;
  char *ptr;
  char *val = NULL;

  ptr = strchr(str, '=') + 1;
  while(ptr[0] == 32 || ptr[0] == 9) ptr++;

  if (ptr == NULL) return NULL;
  i = 0;
  val = (char*) malloc(sizeof(char) * (i+1));
  val[i] = 0;
  while(ptr[i] > 32 && ptr[i] < 127)
  {
    val = (char*) realloc(val, sizeof(char) * (++i+1));
    val[i-1] = ptr[i-1];
    val[i] = 0;
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
  ptr = strchr(str, '=') + 1;
  if (ptr) return atof(ptr);
  else return 0.F;
}

/*This function is rewritten from mydiet; it is simpler and handles string in place*/
char *strdiet(char* in)
{
  int i, j, k;
  j = strlen(in) - 1;
  while(in[j] == 32 || in[j] == 9) j--;
  in[j] = 0;
  for(i = 0; in[i] == 32 || in[i] == 9; i++);
  for(k = 0; k < j-i; k++) in[k] = in[k+i];
  in[j-1] = 0;
  return in;
}

char *get_one_word(char *str)
{
  int i = 0;
  char *n;
  n = (char*) malloc(sizeof(char));
  n[0] = 0;
  while(str[i] != 32 && str[i] != 9 && str[i] != 0)
  {
    n = (char*) realloc(n, sizeof(char) * ++i + 1);
    n[i-1] = str[i-1];
    n[i] = 0;
  }
  return n;
}

char *get_ptr_ith_word(char *str, int i)
{
  int j = 0;
  while(str[j] == 32 || str[j] == 9) j++;
  while(i-1)
  {
    while(str[j] != 32 && str[j] != 9 && str[j] != 0)
    {
      if (str[j] == 0) return NULL;
      j++;
    }
    while(str[j] == 32 || str[j] == 9) j++;
    i--;
  }
  return str+j;
}

int line_is_empty(char *line)
{
  int i = 0;
  if (line == NULL) return 1;
  while(line[i] != 0)
  {
    if (line[i] > 32) return 0;
    i++;
  }
  return 1;
}

