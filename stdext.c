/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<string.h>
#include<ctype.h>

/*this is extension to string.h for oerating systems that do not support this extension*/

char *strcasestr(char *haystack, const char *needle)
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

