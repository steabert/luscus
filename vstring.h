/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/* $Revision: 7.7 $ */

/*****************************************************************************
*                                                                            *
* This software is copyrighted. You are granted the use of this software and *
* can redistribute it freely provided that the software is redistributed in  *
* its original form without omissions.                                       *
*                                                                            *
* Any modifications are prohibited unless explicitly permitted by the MOLCAS *
* development group.                                                         *
* URL: http://www.teokem.lu.se/molcas                                        *
*                                                                            *
* Copyright: Theoretical Chemistry                                           *
*           Lund University                                                  *
*           Sweden                                                           *
*                                                                            *
*****************************************************************************/           


#ifndef _mystring_H_
#define _mystring_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void *xmalloc (size_t size);

int mytoken (char *in, unsigned int c, char *out, char *innew);
   /* cut string _in_ at the position _c_
      in == out // c // innew 
      return 1 if OK.
    */

int myNext (char *, char, int *, int *);
   /*
      NOTE : input sring must be mydiet(); 
      while (myNext(in,' ', &shift, &len)==0)
      { 
      strncpy(out,in+shift,len);
      out[len]=0;
      shift=shift+len+1;
      }

    */

int mycut (char *in, char c, char *out);
   /* cut string _in_ at the position _c_
      in == out // c // ... 
      return 1
    */

char *strcatc (char *str, unsigned int c);
  /* strcat for char */
char *mystrcatl (char *str, unsigned int c);
  /* strcatc if last char not c */

char *chomp (char *str);
char *chop (char *str);

int myterminate (char *in, char *term);
    /* cut string, if any char from terminator are exist */

int my7cut (char *in, char c);
   /* cut string _in_ at the position _c_
      out == out // c // ... 
      return 1
    */

char *myexchange (char *in, unsigned int a, unsigned int b);
   /*  replace in _in_ char _a_ to _b_
      return in.
    */


int myfgets (char *str, int len, FILE * inp);
   /*  like fgets, but
    *    -remove \n
    *    -skip chars before next line 
    *   return 0 at eof, else 1.
    */

int mylastchar (char *str, char c);
   /*  is last char in str _c_ ?
    *   return 0 if not, else 1.
    */

int mylastvchar (char *str, char c);
   /*  is last visible char in str _c_ ?
    *  return 0 is not, else 1.
    */

int mylastdel (char *str, unsigned int c);
   /*  is last char in str is _c_ delete it
    */

int mysubdelete (char *str, char *token, char c);
   /*  if token exist in str, cut token before char c
    */

int mysubstitute (char *in, char *out, char *token, char *ins);
   /*  if token exist in in, substitute it to ins
    */
int mysubtranc (char *in, char *token, char *ins);
   /*  substitute all token to ins on place. strlen(token)>=strlen(ins)
    */

char *mydiet (char *str);
   /*  remove all extra space and \n from str
    */

char *myextradiet (char *str, char c);
   /*  remove all c from str;
    */

char *mydietspace (char *str);
   /*  remove all extra space  from str
    */

char *mystrupr (char *str);
   /*  strupr
    */
char *mystrlwr (char *str);
   /*  strlwr
    */

int mycount (char *in, unsigned int c);
    /* return number of chars c in in */
int mycount2 (char *in, unsigned int c, int n);
int mycounts (char *in, char *sub, char c);

int mycounts_next (char *in, char *sub, char c, char *next, int len);

int mycutsub (char *in, char c);
   /* remove beginning of in */
int myisnan(double);

void printlog(char*);
#endif
