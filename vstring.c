/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck mystring.c $Revision: 7.7 $ */


/*****************************************************************************
*                                                                            *
* Any modifications are prohibited unless explicitly permitted by the luscus *
* authors.                                                                   *
*                                                                            *
*****************************************************************************/           


/****************************************************************************/
/*                                                                          */
/* General string operations                                                */
/*                                                                          */
/****************************************************************************/

#include "vstring.h"
/* See desctription in header file */

int
 myNext (char *in, char c, int *shift, int *len)
{
  char *ptr;
  int tmp;
  tmp=(int)strlen(in);
  if (*shift >= tmp)
    return -1;
  ptr = strchr (in + *shift, c);

  if (ptr == NULL)
    ptr = strchr (in + *shift, '\n');
  if (ptr == NULL)
    {
      *len = (int)strlen (in) - *shift;
    }
  else
    {
      *len = (int)(ptr - in - *shift);
    }
  if (*len < 1)
    return -1;
  return 0;
}

int
mytoken (char *in, unsigned int c, char *out, char *innew)
{
  int i;
  int flag = 0, il = 0;
  for (i = 0; in[i] != 0; i++)
    {
      if (flag == 0)
	{
	  out[i] = in[i];
	  if (in[i] == c)
	    {
	      flag = 1;
	      out[i] = 0;
	      il = i;
	    }
	}
      else
	{
	  innew[i - il - 1] = in[i];
	}
    }
  innew[i - il - 1] = 0;
  return flag;
}


int
mycut (char *in, char c, char *out)
{
  int i = 0;
  while (in[i] != 0 && in[i] != c)
    {
      out[i] = in[i];
      i++;
    }
  out[i] = 0;
  return 1;
}


int
myterminate (char *in, char *term)
{
  int i, j;
  int l = (int)strlen (term);
  for (i = 0; in[i] != 0; i++)
    for (j = 0; j < l; j++)
      if (in[i] == term[j])
	{
	  in[i] = 0;
	  return 1;
	}
  return 1;
}


int
my7cut (char *in, char c)
{
  int i = 0;
  while (in[i] != 0 && in[i] != c)
    {
      i++;
    }
  in[i] = 0;
  return 1;
}

char *
myexchange (char *in, unsigned int a, unsigned int b)
{
  int i;
  for (i = 0; in[i] != 0; i++)
    {
      if (in[i] == a)
	in[i] = b;
    }
  return in;
}

int
myfgets (char *str, int len, FILE * inp)
{
  char temp[80];
  char *ptr;
  if (fgets (str, len - 1, inp) == NULL)
    return 0;
  if ((ptr = strchr (str, '\n')) != NULL)
    {
      *ptr = 0;
      return 1;
    }

  str[len - 1] = 0;
  if (fgets (temp, 80, inp) == NULL)
    return 1;
  while ((strchr (temp, '\n')) == NULL)
    if (fgets (temp, 80, inp) == NULL)
      break;

  return 1;
}

int
mysubdelete (char *str, char *token, char c)
{
  char *ptr, *ptr1;
  ptr = strstr (str, token);
  if (ptr == NULL)
    return 0;
  ptr1 = strchr (ptr, c);
  if (ptr1 == NULL)
    {
      *ptr = 0;
      return 1;
    }
  strcpy (ptr, ptr1 + 1);
  return 1;
}

int
mysubstitute (char *in, char *out, char *token, char *ins)
{
  char *ptr;
  ptr = strstr (in, token);
  if (ptr == NULL)
    return 0;
  strncpy (out, in, (unsigned int) (ptr - in));
  out[ptr - in] = 0;
  strcat (out, ins);
  strcat (out, ptr + strlen (token));
  return 1;
}

int
mylastchar (char *str, char c)
{
  if (str[strlen (str) - 1] == c)
    return 1;
  else
    return 0;
}

char *
strcatc (char *str, unsigned int c)
{
  int i;
  i = (int)strlen (str);
  str[i] = c;
  str[i + 1] = 0;
  return str;
}

char *mystrcatl (char *str, unsigned int c)
{
  int i;
  i = (int)strlen (str);
  if (str[i-1] == c)
    return str;
  str[i] = c;
  str[i + 1] = 0;
  return str;
}

int
mylastvchar (char *str, char c)
{
  int i;
  int j;
  j = (int)strlen (str) - 1;
  for (i = j; i > 0; i--)
    {
      if (str[i] != '\n' && str[i] != ' ')
	break;
    }
  if (str[i] == c)
    return 1;
  else
    return 0;
}


int
mylastdel (char *str, unsigned int c)
{
  int i;
  i = (int)strlen (str) - 1;
  if (str[i] == c)
    str[i] = 0;
  return 0;
}

char *
chop (char *str)
{
  int i;
  i = (int)strlen (str) - 1;
  str[i] = 0;
  return str;
}

char *
chomp (char *str)
{
  int i;
  i = (int)strlen (str) - 1;
  if (str[i] == '\n')
    str[i] = 0;
  return str;
}


char *
mydiet (char *in)
{
  int i, ii, j = 0, l;
  char *temp, *ptr;
  if (!strlen (in))
    return in;

  temp = (char *) malloc ((strlen (in) + 2) * sizeof (char));
  if (temp == NULL)
    {
      fprintf(stderr,"Memory is over\n");
      return NULL;
    }
  if ((ptr = strchr (in, '\n')) != NULL)
    *ptr = ' ';
  l = (int)strlen (in) - 1;
  for (i = 0; i < l; i++)
    {
      if (in[i] == ' ' && in[i + 1] == ' ')
	continue;
      temp[j++] = in[i];
    }
  temp[j++] = in[i];

  temp[j] = 0;
  if (temp[0] == ' ')
    {
      for (ii=0; ii<j; ii++)
       temp[i]=temp[i+1];
/*      strcpy (temp, temp + 1);  */
      j--;
    }

  if (j > 0)
    if (temp[j - 1] == ' ')
      temp[j - 1] = 0;

  strcpy (in, temp);

  free (temp);
  return in;
}

char *
mydietspace (char *in)
{
  int i, j = 0, l, k;
  char *temp;
  if (!strlen (in))
    return in;
  temp = (char *) malloc ((strlen (in) + 2) * sizeof (char));
  if (temp == NULL)
    {
      fprintf(stderr,"Memory is over\n");
      return NULL;
    }
  l = (int)strlen (in) - 1;
  for (i = 0; i < l; i++)
    {
      if (in[i] == ' ' && in[i + 1] == ' ')
	continue;
      temp[j++] = in[i];
    }
  temp[j++] = in[i];

  temp[j] = 0;
  if (temp[0] == ' ')
    {
      for(k=0; k<j; k++)
      {
      temp[k]=temp[k+1];
     /* 
      strcpy (temp, temp + 1);
      */
      }
      j--;
    }

  if (j > 0)
    if (temp[j - 1] == ' ')
      temp[j - 1] = 0;
  if (j > 1)
    if (temp[j - 1] == '\n' && temp[j - 2] == ' ')
      temp[j - 1] = '\n';
  strcpy (in, temp);

  free (temp);
  return in;
}


char *
mystrupr (char *str)
{
/* simple substitution of strupr function */
  int i;
  unsigned int c;

  for (i = 0; str[i] != 0; i++)
    if (islower (c = *(str + i)))
      *(str + i) = toupper ((int)c);
  return str;

}

char *
mystrlwr (char *str)
{
/* simple substitution of strlwr function */
  int i;
  unsigned int c;

  for (i = 0; str[i] != 0; i++)
    if (isupper (c = *(str + i)))
      *(str + i) = tolower ((int)c);
  return str;

}

int
mysubtranc (char *in, char *token, char *ins)
{
  char *ptr;
  int l_ins;
  int l_token;
  int l = 0;
  unsigned int i;
  l_ins = (int)strlen (ins);
  l_token = (int)strlen (token);
  if (l_token < l_ins)
    {
      fprintf(stderr,"Wrong usage of mysubtranc\n");
      return 0;
    }
  while ((ptr = strstr (in + l, token)) != NULL)
    {
    
      strncpy (ptr, ins,(unsigned int)l_ins);
      for(i=0; i<strlen (ptr + l_token); i++)
       {
         *(ptr+l_ins+i)=*(ptr+l_token+i);
       }
/*       
      strncpy (ptr + l_ins, ptr + l_token, strlen (ptr + l_token));
*/
      ptr[l_ins + strlen (ptr + l_token)] = 0;
      l = l_token + 1;
    }
  return 1;
}


/*
int itoa(int i, char* s)
{
char c;
int j,k=0;
  if(i<=0) return -1;
while(i>0)
 {
  j=i-i/10*10;
  i=i/10;
  c=char(j+int('0'));
  s[k++]=c;
if(k>sizeof(s)) return -1;
  }
s[k]=0;
 for(j=0; j<k/2; j++)
 {
 c=s[j];
 s[j]=s[k-1-j];
 s[k-1-j]=c;
 }
 return 0;
}
*/


int
mycount (char *in, unsigned int c)
{
  int i = 0, j = 0;
  while (in[i])
    {
      if (in[i++] == c)
	j++;
    }
  return j;
}

int
mycount2 (char *in, unsigned int c, int n)
{
  int i = 0, j = 0;
  while (in[i])
    {
      if (i == n)
	break;
      if (in[i++] == c)
	j++;

    }
  return j;
}

int
mycounts (char *in, char *sub, char c)
{
  int i = 0, j = 0, l;
  char *ptr;
  ptr=strstr(in,sub);
  if(ptr==NULL) return 0;
  l=(int)(ptr-in);
  while (i<l)
    {
      if (in[i++] == c)
	j++;

    }
  return j;
}

int
mycounts_next (char *in, char *sub, char c, char *next, int len)
{
  int i = 0, j = 0, l;
  char *ptr;
  ptr=strstr(in,sub);
  if(ptr==NULL) return 0;
  
  l=(int)(ptr-in);

  strncpy(next,ptr,(unsigned int)len);
  next[len]=0;
  ptr=strchr(next,c);
  if(ptr) *ptr=0;
  while (i<l)
    {
      if (in[i++] == c)
	j++;

    }
  return j;
}

#include <stdlib.h>

void *
xmalloc (size_t size)
{

  register void *value = malloc (size);
  if (size == 0)
    {
      fprintf(stderr,"malloc call with zero size!\n");
      exit (-1);
    }
  if (value == 0)
    {
      fprintf(stderr, "Virtual memory exhausted\n");
      exit (-1);
    }
  return value;
}

char *
myextradiet (char *s, char c)
{
  int i, j = 0, l;
  char *s1;
  if (s[0] == 0)
    return s;
  if ((s1 = (char *) xmalloc ((strlen (s) + 1) * sizeof (char))) == NULL)
    return NULL;
  strcpy (s1, s);
  l = (int)strlen (s);
  for (i = 0; i < l; i++)
    {
      if (s1[i] != c)
	{
	  s[j] = s1[i];
	  j++;
	}
    }
  s[j] = 0;
  free (s1);

  return s;
}

static int double_ne(double, double);

int
mycutsub (char *in, char c)
{
  char *ptr;
  ptr = strchr (in, c);
  if (ptr != NULL)
    strcpy (in, ptr + 1);
  return 0;
}

int 
myisnan(double n)
{
 return double_ne(n,n);

}
static int double_ne(double n1,double n2)
{
 return n1!=n2;
}

