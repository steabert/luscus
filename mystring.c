/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

/*char *my_read_str_value(char *str)
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
}*/

/*05.07.2013. This function can read quoted string as one word*/
char *my_read_str_value(char *str)
{
  int i, j;
  char *ptr;
  char *val = NULL;
  int next, quoted = 0;

  if (str == NULL) return NULL;
  if (str[0] == 0) return NULL;
  ptr = strchr(str, '=') + 1;
  if (ptr == NULL) return NULL;
  if (ptr[0] == 0) return NULL;

  while(ptr[0] == 32 || ptr[0] == 9)
  {
    if (ptr[0] == 0) return NULL;
    ptr++;
  }
  if (ptr[0] == 0) return NULL;

  i = j = 0;
  val = (char*) malloc(sizeof(char) * (i+1));
  val[i] = 0;
  next = 1;
  if (ptr[i] < 33 || ptr[i] > 126) next = 0;
  while(next)
  {
    if (ptr[j] == 0) return val;
    else if (ptr[j] == 34)
    {
      if (quoted) quoted = 0;
      else quoted = 1;
    }
    else
    {
      val = (char*) realloc(val, sizeof(char) * (++i+1));
      val[i-1] = ptr[j];
      val[i] = 0;
    }
    j++;
    if (quoted && (ptr[j] < 9 || (ptr[j] > 13 && ptr[j] < 32) || ptr[j] > 126)) next = 0;
    else if (!quoted && (ptr[j] < 33 || ptr[j] > 126)) next = 0;
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
  double tmp;
  ptr = strchr(str, '=') + 1;
  if (ptr)
  {
    sscanf(ptr, "%lf", &tmp);
    return tmp; /*corrected 2014.01.24*/
    /*return atof(ptr);*/
  }
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

/*char *get_one_word(char *str)
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
}*/

 /*05.07.2013 this function is upgraded; now it can accept quoted strings*/
char *get_one_word(char *str)
{
  int i = 0, j = 0;
  char *n;
  int next = 1, quoted = 0;
  n = (char*) malloc(sizeof(char));
  n[0] = 0;

  while(next)
  {
    if (str[j] == 0) return n;
    else if (str[j] == 34)
    {
      if (quoted) quoted = 0;
      else quoted = 1;
    }
    else
    {
      n = (char*) realloc(n, sizeof(char) * ++i + 1);
      n[i-1] = str[j];
      n[i] = 0;
    }
    j++;
    if (quoted && (str[j] < 8 || str[j] > 126)) next = 0;
    else if (!quoted && (str[j] < 33 || str[j] > 126)) next = 0;
  }

  return n;
}

/*char *get_ptr_ith_word(char *str, int i)
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
}*/

/*05.07.2013. this function is upgraded so it can recognize quoted string as one word*/
char *get_ptr_ith_word(char *str, int i)
{
  int j = 0;
  int quoted = 0;
  int next;
  while(str[j] == 32 || str[j] == 9) j++;
  while(i-1)
  {
    next = 1;

    while(next)
    {
      if (str[j] == 0) return NULL;
      else if (str[j] == 34)
      {
	if (quoted) quoted = 0;
       	else quoted = 1;
      }
      j++;
      if (quoted && str[j] < 8) next = 0;
      else if (!quoted && (str[j] == 32 || str[j] == 9)) next = 0;
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

