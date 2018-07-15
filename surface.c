/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*cdeck surface.c $Revision: 7.7 $ */
/*****
Project: Iso-surface drawing

Multi-surface support

*******/

/**** Implementation ****/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"luscus.h"
#include"gv.h"
#include"gveps.h"
#include"surface.h"

void make_warning(char*);
/* cleanup a surface */
static void srf_Clean(surface_t*);

/* Calculate the max. diameter of the surfaces */
/*GLdouble msrf_Diameter(multisrf_t*);*//*No such functions -Goran*/

/* Internal functions */

/* comparison function for depth sorting */
static int compare_triangles(const void*, const void*);
/* sort all triangles */
static void Sort_by_Depth(multisrf_t*);

static void Compute_Distances(multisrf_t*);
/* Construct "index" array of pointers to all triangles */
static void Make_Index(multisrf_t*);


/* Initializer/constructor */
multisrf_t* msrf_Init(multisrf_t* self, int n_srf)
{
 if (self==0) {
   self=(multisrf_t*)calloc(sizeof(multisrf_t),1);
 }
 
 self->index=0;
 
 self->n_surfaces=n_srf;
 self->child=(surface_t**)calloc(sizeof(surface_t*),(unsigned int)n_srf);
 
 return self;
}

/* Cleanup */
void msrf_Clean(multisrf_t* self)
{
 int i;
 
 for (i=0; i<self->n_surfaces; i++) {
   surface_t* child=self->child[i];
   if (child) {
     srf_Clean(child);
     free(child);
   }
 }
 
 if (self->index) free(self->index);
 self->index=0;
 self->is_valid=0;
}

/* Construct & init new surface */
surface_t* msrf_New_Surface(multisrf_t* self)
{
 int child_no;
 surface_t* new_child=0;
 int i;
 
 for (child_no=0; child_no<self->n_surfaces; child_no++)
 {
   if (self->child[child_no]==0) goto place; /* place for a new surface */
 }
 /* if we're here, there's no place for a new surface */
 return 0;
 
 place:
#ifdef TRACE_MALLOC
 puts("allocate new_child");
#endif
 new_child=(surface_t*)malloc(sizeof(surface_t));
 new_child->head=0;
 new_child->n_triangles=0;

 new_child->parent=self;
 new_child->offset=child_no;

 for (i=0; i<3; i++) {
   new_child->max_coor[i]=new_child->min_coor[i]=0;
 }
 
 self->is_valid=0; /* invalidate triangle array */
 self->child[child_no]=new_child;
 return new_child;
}
 

/* Delete a surface */
void msrf_Delete_Surface(surface_t* self)
{
 multisrf_t* parent=self->parent;
 
 /* assertion */
 if (parent->child[self->offset]!=self) { /* can't happen */
/*   SayError("msrf_Delete_Surface(): structure corrupted");*/
   make_warning("msrf_Delete_Surface(): structure corrupted");
   exit(1);
 }
 
 parent->is_valid=0; /* invalidate triangle array */
 parent->child[self->offset]=0; /* free the slot */
 srf_Clean(self);
 free(self);
}

/* add triangle with normals */
triangle_t* srf_Add_Triangle(surface_t* self,
                             double txyz[3][3],
                             double tgrd[3][3],
                             double val[3],
                             int clr)
{

 triangle_t* tri=(triangle_t*)malloc(sizeof(triangle_t));
 int vi;
 
 if (self->n_triangles==0) { /* the 1st triangle */
   int k;
   for (k=0; k<3; k++) {
     self->min_coor[k]=self->max_coor[k]=txyz[k][0]; /* init min/max vaues */
   }
 }
   
 for (vi=0; vi<3; vi++) { /* for each vertex */
   int k;
   for (k=0; k<3; k++) { /* for each coor */
     tri->vertex[vi][k]=txyz[k][vi];
     tri->normal[vi][k]=tgrd[k][vi];

     /* update bounding box coors */
     if (self->max_coor[k] < txyz[k][vi]) { self->max_coor[k]=txyz[k][vi]; }
     if (self->min_coor[k] > txyz[k][vi]) { self->min_coor[k]=txyz[k][vi]; }
   }
   tri->value[vi] = val[vi];
 } 
 tri->distance=0;
 tri->clr_idx=clr;
 
 /* add this node-triangle to list */ 
 tri->lnk_prev=self->head;
 self->head=tri;
 
 self->n_triangles++;
 
 self->parent->is_valid=0; /* invalidate array */
 return tri;
}



/* draw all surfaces in current OpenGL context */
void msrf_Draw(multisrf_t* self, color_t color[],int povray, FILE *fl_povray, int ps, FILE *fl_ps, gvepsdef *ctl)
{
 int i;

#ifdef TRACE
  puts("Trace: msrf_Draw");
#endif
 if (!self->is_valid) {
   if (self->index) free(self->index);
   Make_Index(self);
 }
 
 if (!self->is_sorted) {
   Sort_by_Depth(self);
 }
 
 if(povray)
 {
   fprintf(fl_povray,"\n#declare colorA =\n texture {\n pigment { rgbt <%7.4f, %7.4f, %7.4f, %7.4f>}\n",
        color[0][0],color[0][1],color[0][2], 1.0-color[0][3]);
   fprintf(fl_povray,"finish {ambient 0.2 diffuse 0.8 specular 0.8}\n}");

   fprintf(fl_povray,"\n#declare colorB =\n texture {\n pigment { rgbt <%7.4f, %7.4f, %7.4f, %7.4f>}\n",
        color[1][0],color[1][1],color[1][2], 1.0-color[1][3]);
   fprintf(fl_povray,"finish {ambient 0.2 diffuse 0.8 specular 0.8}\n}");
 
 }
 
 glBegin(GL_TRIANGLES);
 for (i=0; i<self->n_triangles; i++) {
   triangle_t* triangle=self->index[i];
   int j;
   double bw;
   if(Input_Data.bw==0)
   {
   glColor4fv(color[triangle->clr_idx]);
   }
   else
   {
   bw=grayscale(color[triangle->clr_idx][0],color[triangle->clr_idx][1],color[triangle->clr_idx][2]);
   glColor4d(bw,bw,bw,color[triangle->clr_idx][3]);
   }
/*   printf("triangle #%d\n %f %f %f\n %f %f %f\n %f %f %f\n", i,
       triangle->vertex[0][0], triangle->vertex[0][1], triangle->vertex[0][2],
       triangle->vertex[1][0], triangle->vertex[1][1], triangle->vertex[1][2],
       triangle->vertex[2][0], triangle->vertex[2][1], triangle->vertex[2][2] );*/
   for (j=0; j<3; j++) {
/*   printf("tria %lf %lf %lf\n",triangle->normal[j][0],triangle->normal[j][1],triangle->normal[j][2]); */
     glNormal3dv(triangle->normal[j]);
     glVertex3dv(triangle->vertex[j]);
   }
/* povray starts here */

  if(povray) 
   {
   	
      fprintf(fl_povray, "\nsmooth_triangle {\n"); 
       for (j=0; j<3; j++) {
      fprintf(fl_povray, "<%10.6f,%10.6f,%10.6f>, <%10.6f, %10.6f, %10.6f>\n", 
       triangle->vertex[j][0],
       triangle->vertex[j][1],
       triangle->vertex[j][2],
       triangle->normal[j][0],
       triangle->normal[j][1],
       triangle->normal[j][2]);
       }
      if ( triangle->clr_idx ==0) fprintf(fl_povray,"texture { colorA }\n}\n");
      else                        fprintf(fl_povray,"texture { colorB }\n}\n");
/*
      fprintf(fl_povray, "  pigment { color rgb <%7.4f, %7.4f, %7.4f>\n  }\n}\n", 
      color[triangle->clr_idx][0],
      color[triangle->clr_idx][1],
      color[triangle->clr_idx][2]);
*/      
    }  
    
/* povray ends here */

/* gveps starts here */

  if(ps) 
   {
       gveps_triangle(fl_ps, 0,  
                      triangle->vertex[0][0],triangle->vertex[0][1],triangle->vertex[0][2],
                      triangle->vertex[1][0],triangle->vertex[1][1],triangle->vertex[1][2],
                      triangle->vertex[2][0],triangle->vertex[2][1],triangle->vertex[2][2],

                      triangle->normal[0][0],triangle->normal[0][1],triangle->normal[0][2],
                      triangle->normal[1][0],triangle->normal[1][1],triangle->normal[1][2],
                      triangle->normal[2][0],triangle->normal[2][1],triangle->normal[2][2],

                      color[triangle->clr_idx][0],
                      color[triangle->clr_idx][1],
                      color[triangle->clr_idx][2],
		      ctl);
		      
   }  
    
/* gveps ends here */
   
 }
 glEnd();

} 
 
/* draw all surfaces in current OpenGL context */
void msrf_Draw1(multisrf_t* self, color_t color[],int povray, FILE *fl_povray, int ps, FILE *fl_ps, gvepsdef *ctl)
{
 int i;
 double tp_color[4];
 double white_region = 0.025;
 double max_value = 0.10;
 tp_color[3] = color[0][3];

#ifdef TRACE
  puts("Trace: msrf_Draw");
#endif
 if (!self->is_valid) {
   if (self->index) free(self->index);
   Make_Index(self);
 }
 
 if (!self->is_sorted) {
   Sort_by_Depth(self);
 }
 
 if(povray) /*this must be changed!*/
 {
   fprintf(fl_povray,"\n#declare colorA =\n texture {\n pigment { rgb <%7.4f, %7.4f, %7.4f>}\n",
        color[0][0],color[0][1],color[0][2]);
   fprintf(fl_povray,"finish {ambient 0.2 diffuse 0.8 specular 0.8}\n}");

   fprintf(fl_povray,"\n#declare colorB =\n texture {\n pigment { rgb <%7.4f, %7.4f, %7.4f>}\n",
        color[1][0],color[1][1],color[1][2]);
   fprintf(fl_povray,"finish {ambient 0.2 diffuse 0.8 specular 0.8}\n}");
 
 }
 
 glBegin(GL_TRIANGLES);
 for (i=0; i<self->n_triangles; i++) {
   triangle_t* triangle=self->index[i];
   int j;
   double bw;
/*   if(Input_Data.bw==0)
   {
   glColor4fv(color[triangle->clr_idx]);
   }
   else
   {
   bw=grayscale(color[triangle->clr_idx][0],color[triangle->clr_idx][1],color[triangle->clr_idx][2]);
   glColor4d(bw,bw,bw,color[triangle->clr_idx][3]);
   }*/
/*   printf("triangle #%d %f %f %f\n", i, triangle->value[0], triangle->value[1], triangle->value[2]);*/
/*   printf("triangle #%d\n %f %f %f\n %f %f %f\n %f %f %f\n", i,
       triangle->vertex[0][0], triangle->vertex[0][1], triangle->vertex[0][2],
       triangle->vertex[1][0], triangle->vertex[1][1], triangle->vertex[1][2],
       triangle->vertex[2][0], triangle->vertex[2][1], triangle->vertex[2][2] );*/
   for (j=0; j<3; j++)
   {
     int k;
     if (Input_Data.bw==1)
     {
       if (triangle->value[j] > 1.0) bw = 1.0;
       else if (triangle->value[j] < -1.0) bw = -1.0;
       else if (triangle->value[j] > -white_region && triangle->value[j] < white_region) bw = 0.5;
       else if (triangle->value[j] > white_region) bw = 0.5 + 0.5 *(triangle->value[j] - white_region) / (1.0 - white_region);
       else if (triangle->value[j] < -white_region) bw = 0.5 - 0.5 *(white_region + triangle->value[j]) / (1.0 - white_region);
       glColor4d(bw,bw,bw,tp_color[3]);
     }
     else
     {
       for(k = 0; k < 3; k++)
       {
         if (triangle->value[j] > max_value) tp_color[k] = color[0][k];
         else if (triangle->value[j] < -max_value) tp_color[k] = color[1][k];
         else if (triangle->value[j] < white_region && triangle->value[j] > -white_region) tp_color[k] = 1.0;
         else if (triangle->value[j] > white_region)
           tp_color[k] = 1.0 - (1.0 - color[0][k]) * (triangle->value[j] - white_region)/(max_value-white_region);
         else if (triangle->value[j] < -white_region)
           tp_color[k] = 1.0 + (1.0 - color[1][k]) * (white_region + triangle->value[j])/(max_value-white_region);
       }
       glColor4dv(tp_color);
     }

/*   printf("tria %lf %lf %lf\n",triangle->normal[j][0],triangle->normal[j][1],triangle->normal[j][2]); */
     glNormal3dv(triangle->normal[j]);
     glVertex3dv(triangle->vertex[j]);
   }
/* povray starts here */

  if(povray) 
   {
   	
      fprintf(fl_povray, "\nsmooth_triangle {\n"); 
       for (j=0; j<3; j++) {
      fprintf(fl_povray, "<%10.6f,%10.6f,%10.6f>, <%10.6f, %10.6f, %10.6f>\n", 
       triangle->vertex[j][0],
       triangle->vertex[j][1],
       triangle->vertex[j][2],
       triangle->normal[j][0],
       triangle->normal[j][1],
       triangle->normal[j][2]);
       }
      if ( triangle->clr_idx ==0) fprintf(fl_povray,"texture { colorA }\n}\n");
      else                        fprintf(fl_povray,"texture { colorB }\n}\n");
/*
      fprintf(fl_povray, "  pigment { color rgb <%7.4f, %7.4f, %7.4f>\n  }\n}\n", 
      color[triangle->clr_idx][0],
      color[triangle->clr_idx][1],
      color[triangle->clr_idx][2]);
*/      
    }  
    
/* povray ends here */

/* gveps starts here */

  if(ps) 
   {
       gveps_triangle(fl_ps, 0,  
                      triangle->vertex[0][0],triangle->vertex[0][1],triangle->vertex[0][2],
                      triangle->vertex[1][0],triangle->vertex[1][1],triangle->vertex[1][2],
                      triangle->vertex[2][0],triangle->vertex[2][1],triangle->vertex[2][2],

                      triangle->normal[0][0],triangle->normal[0][1],triangle->normal[0][2],
                      triangle->normal[1][0],triangle->normal[1][1],triangle->normal[1][2],
                      triangle->normal[2][0],triangle->normal[2][1],triangle->normal[2][2],

                      color[triangle->clr_idx][0],
                      color[triangle->clr_idx][1],
                      color[triangle->clr_idx][2],
		      ctl);
		      
   }  
    
/* gveps ends here */
   
 }
 glEnd();

} 



/* Internal functions */

/* cleanup a surface */
static void srf_Clean(surface_t* self)
{
 triangle_t* tp=self->head;

#ifdef TRACE
  puts("Trace: srf_Clean");
#endif

 while (tp) { /* go through the linked list of triangles */
   triangle_t* link=tp->lnk_prev;
   free(tp);
   tp=link;
 }
}

/* comparison function for depth sorting */
static int compare_triangles(const void* tptr1, const void* tptr2)
{
 const triangle_t* tre1_ptr=*(const triangle_t **)tptr1;
 const triangle_t* tre2_ptr=*(const triangle_t **)tptr2;
 
 if (tre1_ptr->distance > tre2_ptr->distance) return -1;
 if (tre1_ptr->distance < tre2_ptr->distance) return  1;
 /* Debug("Triangles are eq: %g %g",tre1_ptr->distance,tre2_ptr->distance); */
 return 0;
}

/* sort all triangles */
static void Sort_by_Depth(multisrf_t* self)
{
 /* recalculate distances */
 Compute_Distances(self);
 /* reorder triangle indices, sorting by depth */
 qsort(self->index,self->n_triangles,sizeof(triangle_t*),compare_triangles);
 
 /* mark as sorted */
 self->is_sorted=1;
}

static void Compute_Distances(multisrf_t* self)
{
 GLdouble mvm[16];
 int idx;

/*
 if (!Surface.index) {
   Make_Surface_Index();
 }
*/
 glGetDoublev(GL_MODELVIEW_MATRIX,mvm); /* get model-view matrix */
 for (idx=0; idx<self->n_triangles; idx++) {
   GLdouble avg_z=0;
   triangle_t* tri=self->index[idx];
   int ivert;
   for (ivert=0; ivert<3; ivert++) {
     GLdouble z_dist=0;
     int k;
     for (k=0; k<3; k++) {
       z_dist+=tri->vertex[ivert][k]*mvm[2+4*k];
     }
     z_dist+=mvm[2+4*3]; /* add 4th component of scalar product */
     
     avg_z+=z_dist;
   }
   tri->distance=-avg_z/3;
 }
}

/* Construct "index" array of pointers to all triangles */
static void Make_Index(multisrf_t* self)
{
 int ichild;
 triangle_t** tptr;

 /* compute total number of triangles */
 self->n_triangles=0;
 for (ichild=0; ichild<self->n_surfaces; ichild++) {
   surface_t* child=self->child[ichild];
   if (child) {
     self->n_triangles+=child->n_triangles;
   }
 }
 
 /* allocate memory for array */
 self->index=tptr=(triangle_t**)calloc((unsigned)self->n_triangles,sizeof(triangle_t*));
 
 /* fill the array: */
 /* for each child surface */
 for (ichild=0; ichild<self->n_surfaces; ichild++) {
   surface_t* child=self->child[ichild];
   triangle_t* link;
   int i;
   
   if (!child) continue;

   /* for each triangle within surface */
   for (i=0, link=child->head; i<child->n_triangles; i++) {
     if (!link) { /* can't happen */
       /*SayError("Make_Index(): inconsistency, too few triangles!");*/
       make_warning("Make_Index(): inconsistency, too few triangles!");
       exit(1);
     }
     *(tptr++)=link; /* store pointer to triangle */
     link=link->lnk_prev;
   } /* endfor: by triangles */
   
   if (link) { /* can't happen */
/*     SayError("Make_Index(): inconsistency, too many triangles!");*/
       make_warning("Make_Index(): inconsistency, too many triangles!");
     exit(1);
   }
 } /* endfor: by surfaces */
 
 /* mark as valid, but unsorted */
 self->is_valid=1;
 self->is_sorted=0;
}
