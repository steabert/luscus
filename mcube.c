/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck mcube.c $Revision: 7.7 $ */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"luscus.h"
#include"gv.h"
#include"mcube.h"
#include"gveps.h"
#include"surface.h"


static int calcgr(int ind,int jnd,int knd,double *fvals,
           int griddim[3], 
           double cubgr[3][8]);

int mcubes(surface_t* surf,
           int nx,int ny,int nz, double fvals[], double evals[], double origin[3],
           double cntval,double rad[3], int idominus)
{

	   
int i,j,k;
int  griddim[3], yzdim;
int casemask[]={1,2,4,8,16,32,64,128};
int tptr,caseptr,findex,posneg;
int  cneg;
int lines[][12]={
     {1,2,4,1,5,6,8,5,1,2,4,3},
     {2,3,3,4,6,7,7,8,5,6,8,7}
    };

      
int linepts[2];
double factor,xtmp,ytmp,ztmp,contval;
int l,m;
double scale[3];
double txyz[3][3],tgrd[3][3],epval[3];
double cubval[8], cubxyz[8][3], cubgr[3][8];
double cubeval[8];
double norm;

      griddim[0] = nx;
      griddim[1] = ny;
      griddim[2] = nz;

       idominus=!!idominus; 

   for (i=0;i<3;i++)
          scale[i] = rad[i]/(griddim[i]-1.0);

         yzdim = griddim[2] * griddim[1];
	 for (posneg=0; posneg<idominus+1; posneg++)
	 {
            contval = cntval*(1-2*posneg);
            cneg = 0;
            if (contval<0.0) cneg = 1;
            for( k=0; k<griddim[2]-1; k++)
	      {
               cubxyz[0][2] = origin[2] + k*scale[2];
               ztmp = origin[2] + (k+1)*scale[2];
               for (j=0; j<griddim[1]-1; j++)
	        {
                  cubxyz[0][1] = origin[1] + j*scale[1];
                  ytmp = origin[1] + (j+1)*scale[1];

                  for( i=0; i<griddim[0]-1; i++)
		  {
                     findex = k + j*griddim[2] + i*yzdim;
                     cubval[0] = fvals[findex];
                     cubval[1] = fvals[findex +yzdim];
                     cubval[2] = fvals[findex +yzdim+griddim[2]];
                     cubval[3] = fvals[findex +griddim[2]];
                     cubval[4] = fvals[findex +1];
                     cubval[5] = fvals[findex +1+yzdim];
                     cubval[6] = fvals[findex +1+griddim[2]+yzdim];
                     cubval[7] = fvals[findex +1+griddim[2]];
                     if (evals)
                     {
                       cubeval[0] = evals[findex];
                       cubeval[1] = evals[findex +yzdim];
                       cubeval[2] = evals[findex +yzdim+griddim[2]];
                       cubeval[3] = evals[findex +griddim[2]];
                       cubeval[4] = evals[findex +1];
                       cubeval[5] = evals[findex +1+yzdim];
                       cubeval[6] = evals[findex +1+griddim[2]+yzdim];
                       cubeval[7] = evals[findex +1+griddim[2]];
                     }
                     caseptr = 0;
                     for( l=0; l<8; l++)
		     {
                         if (cneg) {
                            if (cubval[l]>=contval) caseptr = caseptr + casemask[l];
                          }
                         else {
                            if (cubval[l]<=contval) caseptr = caseptr + casemask[l];
                          }
                     }

                     if (! (caseptr==0 || caseptr==255 ))
                    {		      
/*                    create Cube points */
                     cubxyz[0][0] = origin[0] + i*scale[0];
                     xtmp = origin[0] + (i+1)*scale[0];
                     cubxyz[1][0] = xtmp;           cubxyz[1][1] = cubxyz[0][1];   cubxyz[1][2] = cubxyz[0][2];
                     cubxyz[2][0] = xtmp;           cubxyz[2][1] = ytmp;	   cubxyz[2][2] = cubxyz[0][2];
                     cubxyz[3][0] = cubxyz[0][0];   cubxyz[3][1] = ytmp;	   cubxyz[3][2] = cubxyz[0][2];
                     cubxyz[4][0] = cubxyz[0][0];   cubxyz[4][1] = cubxyz[0][1];   cubxyz[4][2] = ztmp;
                     cubxyz[5][0] = xtmp;           cubxyz[5][1] = cubxyz[0][1];   cubxyz[5][2] = ztmp;
                     cubxyz[6][0] = xtmp;           cubxyz[6][1] = ytmp;	   cubxyz[6][2] = ztmp;
                     cubxyz[7][0] = cubxyz[0][0];   cubxyz[7][1] = ytmp;	   cubxyz[7][2] = ztmp;
/*                    create Cube gradients */
                     calcgr(i,  j,  k,  fvals,griddim,cubgr);
                     for(l=0; l<8; l++)
		      {
		         cubgr[0][l]=cubgr[0][l]/scale[0];
		         cubgr[1][l]=cubgr[1][l]/scale[1];
		         cubgr[2][l]=cubgr[2][l]/scale[2];
		      }
                     tptr = 0;
                     while (cases[caseptr][tptr]>-1)
		     {
                        for( l=0; l<3; l++)
			{
                           linepts[0] = lines[0][cases[caseptr][tptr+l]]-1;
                           linepts[1] = lines[1][cases[caseptr][tptr+l]]-1;
                           factor = (contval - cubval[linepts[0]])/(cubval[linepts[1]] - cubval[linepts[0]]);

                           for( m=0; m<3; m++)
			    {
                              txyz[m][l] = cubxyz[linepts[0]][m] + factor*(cubxyz[linepts[1]][m] - cubxyz[linepts[0]][m]);
                              tgrd[m][l] = cubgr[m][linepts[0]] + factor*(cubgr[m][linepts[1]] - cubgr[m][linepts[0]]);
                            }
                           norm = sqrt(tgrd[0][l]*tgrd[0][l] + tgrd[1][l]*tgrd[1][l] + tgrd[2][l]*tgrd[2][l]);
                           if (norm>0.0) { tgrd[0][l]=tgrd[0][l]/norm; tgrd[1][l]=tgrd[1][l]/norm; tgrd[2][l]=tgrd[2][l]/norm; }
                           if (contval>0.0) 
			     {
			      tgrd[0][l]=-tgrd[0][l];      tgrd[1][l]=-tgrd[1][l];      tgrd[2][l]=-tgrd[2][l];
			     }

                           if (evals)
                             epval[l] = cubeval[linepts[0]] + factor*(cubeval[linepts[1]] - cubeval[linepts[0]]);
                           else
                             epval[l] = 0.0;
                        }
                         if ( ! (txyz[0][0]==txyz[0][1] && txyz[1][0]==txyz[1][1] && txyz[2][0]==txyz[2][1]) 
			     || (txyz[0][0]==txyz[0][2] && txyz[1][0]==txyz[1][2] && txyz[2][0]==txyz[2][2]) 
			     || (txyz[0][1]==txyz[0][2] && txyz[1][1]==txyz[1][2] && txyz[2][1]==txyz[2][2]))
			   { 
/*			   printf("txz %lf %lf %lf\n",txyz[0][0],txyz[0][1],txyz[0][2]); */
                            if (surf) srf_Add_Triangle(surf,txyz,tgrd,epval,posneg);
                           }
                             
                        tptr = tptr + 3;
                     }
                    }
                  }
               }
            }
         }

      
  return 0;
}


int calcgr(int ind,int jnd,int knd,double *fvals,
int griddim[3], 
double cubgr[3][8])
{   
int dir[3], dir1[3], dir2[3];
int i,j;
double f1,f2, fact;
int yzdim;
int it;
int shiftI[]={0,1,1,0,0,1,1,0};
int shiftJ[]={0,0,1,1,0,0,1,1};
int shiftK[]={0,0,0,0,1,1,1,1};

     yzdim = griddim[2] * griddim[1];

     for(it=0; it<8; it++)
     {
       dir[0] = ind+shiftI[it];
       dir[1] = jnd+shiftJ[it];
       dir[2] = knd+shiftK[it];
       for (i=0; i<3; i++)
        {
          fact = 1;
           for (j=0; j<3; j++)
	   {
            dir1[j] = dir[j];
            dir2[j] = dir[j];
	   }   
        
           if (dir[i]==0) { dir2[i] ++; }
	    
           else
	     {
	       if (dir[i]==(griddim[i]-1)) {  dir1[i] --; }
               else
	        {
                 dir1[i] --;
                 dir2[i] ++;
                 fact = 0.5;
	        }
	      }
        f2 = fvals[dir2[2] + dir2[1]*griddim[2] + dir2[0]*yzdim];
        f1 = fvals[dir1[2] + dir1[1]*griddim[2] + dir1[0]*yzdim];
        cubgr[i][it] = fact * (f2 - f1) ;
       }
     }	 
      return 0;
}

