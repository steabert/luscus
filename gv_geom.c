/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

MOL *m;
MOL *mol;
int n_watched;
WATCH *watch;
INPUT_DATA Input_Data;

double cumulative_move = 0.0;

void get_center(double *x, double *y, double *z)
{
  int i;
  int ntotal = 0;
  *x = 0.F;
  *y = 0.F;
  *z = 0.F;

  if (!m->natom) /*rare case: if there are no atoms, determine the center of geometrical objects*/
  {
    ntotal = m->nsphere + m->nvector + m->ntriangle + m->ncells;
    if (ntotal)
    {
      for(i = 0; i < m->nsphere; i++)
      {
        *x += m->sphere_center[i][0];
        *y += m->sphere_center[i][0];
        *z += m->sphere_center[i][0];
      }
      for(i = 0; i < m->nvector; i++)
      {
        *x += 0.5*(m->vector1[i][0] + m->vector2[i][0]);
        *y += 0.5*(m->vector1[i][1] + m->vector2[i][1]);
        *z += 0.5*(m->vector1[i][2] + m->vector2[i][2]);
      }
      for(i = 0; i < m->ntriangle; i++)
      {
        *x += 0.3333333*(m->triangle1[i][0] + m->triangle2[i][0] + m->triangle3[i][0]);
        *y += 0.3333333*(m->triangle1[i][1] + m->triangle2[i][1] + m->triangle3[i][1]);
        *z += 0.3333333*(m->triangle1[i][2] + m->triangle2[i][2] + m->triangle3[i][2]);
      }
      for(i = 0; i < m->ncells; i++)
      {
        *x += 0.25*(m->cell1[i][0] + m->cell2[i][0] + m->cell3[i][0] + m->cell3[i][0]);
        *y += 0.25*(m->cell1[i][1] + m->cell2[i][1] + m->cell3[i][1] + m->cell3[i][1]);
        *z += 0.25*(m->cell1[i][2] + m->cell2[i][2] + m->cell3[i][2] + m->cell3[i][2]);
      }
      *x /= (double) ntotal;
      *y /= (double) ntotal;
      *z /= (double) ntotal;
    }
  }
  else if (m->n_marked)
  {
    for(i = 0; i < m->n_marked; i++)
    {
      *x += m->xyz[m->marked[i]][0];
      *y += m->xyz[m->marked[i]][1];
      *z += m->xyz[m->marked[i]][2];
    }

    *x /= (double) m->n_marked;
    *y /= (double) m->n_marked;
    *z /= (double) m->n_marked;
  }
  else
  {
    for(i = 0; i < m->natom; i++)
    {
      *x += m->xyz[i][0];
      *y += m->xyz[i][1];
      *z += m->xyz[i][2];
    }

    *x /= (double) m->natom;
    *y /= (double) m->natom;
    *z /= (double) m->natom;
  }
  return;
}

void set_origin_molecule(void)
{
  int i, j;
  XYZ cm;

  if (m->n_selected == 0) get_center(&cm[0], &cm[1], &cm[2]);
  else if (m->n_selected == 1)
    for(i = 0; i < 3; i++)
      cm[i] = m->xyz[m->selected[0]][i];
  else if (m->n_selected == 2)
    for(i = 0; i < 3; i++)
      cm[i] = 0.5F * (m->xyz[m->selected[0]][i] + m->xyz[m->selected[1]][i]);

  for(i = 0; i < m->natom; i++)
  {
    for(j = 0; j < 3; j++)
      m->xyz[i][j] -= cm[j];
    change_atom_parameters_in_list(i);
  }

  for(i = 0; i < m->nvector; i++)
    for(j = 0; j < 3; j++)
    {
      m->vector1[i][j] -= cm[j];
      m->vector2[i][j] -= cm[j];
    }

  for(i = 0; i < m->ntriangle; i++)
    for(j = 0; j < 3; j++)
    {
      m->triangle1[i][j] -= cm[j];
      m->triangle2[i][j] -= cm[j];
      m->triangle3[i][j] -= cm[j];
    }

  for(i = 0; i < m->nsphere; i++)
    for(j = 0; j < 3; j++)
      m->sphere_center[i][j] -= cm[j];

  for(i = 0; i < m->nsurf; i++)
    for(j = 0; j < 3; j++)
    {
      m->surf1[i][j] -= cm[j];
      m->surf2[i][j] -= cm[j];
      m->surf3[i][j] -= cm[j];
    }

  for(i = 0; i < m->ncells; i++)
    for(j = 0; j < 3; j++)
    {
      m->cell1[i][j] -= cm[j];
      m->cell2[i][j] -= cm[j];
      m->cell3[i][j] -= cm[j];
      m->cell4[i][j] -= cm[j];
    }

  if (m->ngrids)
  {
    for(j = 0; j < 3; j++)
      m->grid.origin[j] -= cm[j];
    make_surfaces();
  }

  rerender_3d();
}

int find_first_bounded(int ia1, int ia2)
{
  int i;
  int ii=-1, type;

  /* find bounded atom */
  for (i = 0; i < m->nbond; i++)
  {
    int ib1=m->bond[i].iat1;
    int ib2=m->bond[i].iat2;
    type=m->bond[i].bond_type;

    if(type>0)
    {
      if((ib1==ia1 && ib2!=ia2) || (ib1==ia2 && ib2!=ia1) ||
         (ib2==ia1 && ib1!=ia2) || (ib2==ia2 && ib1!=ia1))
      {
        if(ib1==ia1 || ib1==ia2)
        {
          ii=ib2;
          break;
        }
        if(ib2==ia1 || ib2==ia2)
        {
          ii=ib1;
          break;
        }
      }
    }
  }
  return ii;
}

double get_selected_bond_length(void)
{
  return sqrt((m->xyz[m->selected[1]][0] - m->xyz[m->selected[0]][0]) *
              (m->xyz[m->selected[1]][0] - m->xyz[m->selected[0]][0]) +
              (m->xyz[m->selected[1]][1] - m->xyz[m->selected[0]][1]) *
              (m->xyz[m->selected[1]][1] - m->xyz[m->selected[0]][1]) + 
              (m->xyz[m->selected[1]][2] - m->xyz[m->selected[0]][2]) *
              (m->xyz[m->selected[1]][2] - m->xyz[m->selected[0]][2]));
}

double get_selected_angle_value(void)
{
  double sp, r1, r2;

  sp = rr_vec(m->xyz[m->selected[0]], m->xyz[m->selected[1]], m->xyz[m->selected[2]]);
  r1 = rr_dist(m->xyz[m->selected[0]], m->xyz[m->selected[1]]);
  r2 = rr_dist(m->xyz[m->selected[2]], m->xyz[m->selected[1]]);
  return acos(sp / (r1 * r2)) * 180.F / M_PI;
}

double get_selected_dihedral_value(void)
{
  double xi0[3], xi1[3], xi2[3], xi3[3];
  double r12[3], r23[3], r34[3];
  double l, cosphiNy, sinphiNy, o[3], p[3], q[3];

  xi0[0] = m->xyz[m->selected[0]][0];  xi0[1] = m->xyz[m->selected[0]][1]; xi0[2]=m->xyz[m->selected[0]][2];
  xi1[0] = m->xyz[m->selected[1]][0];  xi1[1] = m->xyz[m->selected[1]][1]; xi1[2]=m->xyz[m->selected[1]][2];
  xi2[0] = m->xyz[m->selected[2]][0];  xi2[1] = m->xyz[m->selected[2]][1]; xi2[2]=m->xyz[m->selected[2]][2];
  xi3[0] = m->xyz[m->selected[3]][0];  xi3[1] = m->xyz[m->selected[3]][1]; xi3[2]=m->xyz[m->selected[3]][2];

  r12[0] = xi0[0] - xi1[0];  r12[1] = xi0[1] - xi1[1];  r12[2] = xi0[2] - xi1[2];
  r23[0] = xi1[0] - xi2[0];  r23[1] = xi1[1] - xi2[1];  r23[2] = xi1[2] - xi2[2];
  r34[0] = xi2[0] - xi3[0];  r34[1] = xi2[1] - xi3[1];  r34[2] = xi2[2] - xi3[2];

  o[0] = r12[1]*r23[2]-r12[2]*r23[1];
  o[1] = r12[2]*r23[0]-r12[0]*r23[2];
  o[2] = r12[0]*r23[1]-r12[1]*r23[0];

  l = sqrt(o[0]*o[0]+o[1]*o[1]+o[2]*o[2]);
  o[0] = o[0]/l; o[1] = o[1]/l; o[2] = o[2]/l;

  p[0] = r23[1]*r34[2]-r23[2]*r34[1]; 
  p[1] = r23[2]*r34[0]-r23[0]*r34[2]; 
  p[2] = r23[0]*r34[1]-r23[1]*r34[0];

  l = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  p[0] = p[0]/l; p[1] = p[1]/l; p[2] = p[2]/l;

  q[0] = r23[1]*o[2]-r23[2]*o[1]; 
  q[1] = r23[2]*o[0]-r23[0]*o[2]; 
  q[2] = r23[0]*o[1]-r23[1]*o[0];

  l = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
  q[0] = q[0]/l; q[1] = q[1]/l; q[2] = q[2]/l;

  cosphiNy = (o[0]*p[0]+o[1]*p[1]+o[2]*p[2]);
  sinphiNy = (q[0]*p[0]+q[1]*p[1]+q[2]*p[2]);

  return atan2(sinphiNy, cosphiNy) * 180.F / M_PI;
}

double get_bond_length(int at1, int at2)
{
  return sqrt((m->xyz[at2][0] - m->xyz[at1][0]) *
              (m->xyz[at2][0] - m->xyz[at1][0]) +
              (m->xyz[at2][1] - m->xyz[at1][1]) *
              (m->xyz[at2][1] - m->xyz[at1][1]) + 
              (m->xyz[at2][2] - m->xyz[at1][2]) *
              (m->xyz[at2][2] - m->xyz[at1][2]));
}

double get_angle_value(int at1, int at2, int at3)
{
  double sp, r1, r2;

  sp = rr_vec(m->xyz[at1], m->xyz[at2], m->xyz[at3]);
  r1 = rr_dist(m->xyz[at1], m->xyz[at2]);
  r2 = rr_dist(m->xyz[at3], m->xyz[at2]);
  return acos(sp / (r1 * r2)) * 180.F / M_PI;
}

double get_dihedral_value(int at0, int at1, int at2, int at3)
{
  double xi0[3], xi1[3], xi2[3], xi3[3];
  double r12[3], r23[3], r34[3];
  double l, cosphiNy, sinphiNy, o[3], p[3], q[3];

  xi0[0] = m->xyz[at0][0];  xi0[1] = m->xyz[at0][1]; xi0[2]=m->xyz[at0][2];
  xi1[0] = m->xyz[at1][0];  xi1[1] = m->xyz[at1][1]; xi1[2]=m->xyz[at1][2];
  xi2[0] = m->xyz[at2][0];  xi2[1] = m->xyz[at2][1]; xi2[2]=m->xyz[at2][2];
  xi3[0] = m->xyz[at3][0];  xi3[1] = m->xyz[at3][1]; xi3[2]=m->xyz[at3][2];

  r12[0] = xi0[0] - xi1[0];  r12[1] = xi0[1] - xi1[1];  r12[2] = xi0[2] - xi1[2];
  r23[0] = xi1[0] - xi2[0];  r23[1] = xi1[1] - xi2[1];  r23[2] = xi1[2] - xi2[2];
  r34[0] = xi2[0] - xi3[0];  r34[1] = xi2[1] - xi3[1];  r34[2] = xi2[2] - xi3[2];

  o[0] = r12[1]*r23[2]-r12[2]*r23[1];
  o[1] = r12[2]*r23[0]-r12[0]*r23[2];
  o[2] = r12[0]*r23[1]-r12[1]*r23[0];

  l = sqrt(o[0]*o[0]+o[1]*o[1]+o[2]*o[2]);
  o[0] = o[0]/l; o[1] = o[1]/l; o[2] = o[2]/l;

  p[0] = r23[1]*r34[2]-r23[2]*r34[1]; 
  p[1] = r23[2]*r34[0]-r23[0]*r34[2]; 
  p[2] = r23[0]*r34[1]-r23[1]*r34[0];

  l = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  p[0] = p[0]/l; p[1] = p[1]/l; p[2] = p[2]/l;

  q[0] = r23[1]*o[2]-r23[2]*o[1]; 
  q[1] = r23[2]*o[0]-r23[0]*o[2]; 
  q[2] = r23[0]*o[1]-r23[1]*o[0];

  l = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
  q[0] = q[0]/l; q[1] = q[1]/l; q[2] = q[2]/l;

  cosphiNy = (o[0]*p[0]+o[1]*p[1]+o[2]*p[2]);
  sinphiNy = (q[0]*p[0]+q[1]*p[1]+q[2]*p[2]);

  return atan2(sinphiNy, cosphiNy) * 180.F / M_PI;
}

void accumulate_motion(double value)
{
  cumulative_move += fabs(value);

  if (cumulative_move > 0.1)
  {
    append_backup();
    cumulative_move = 0.0;
  }

  change_watched_data();
  if (Input_Data.automatic_rebonding) rebond();
}

void luscus_gtk_move_bond(double value)
{
  int i, j;
  int i0, i1, ii0;
  double r01[3];
  double tempr;
  if (m->n_selected < 2 || m->n_selected > 4) return;

  i0=m->selected[0];
  i1=m->selected[1];

  for(i = 0; i < 3; i++) r01[i] = m->xyz[i0][i] - m->xyz[i1][i];
  tempr=sqrt(r01[0]*r01[0]+r01[1]*r01[1]+r01[2]*r01[2]);

  for(i = 0; i < 3; i++) r01[i]=r01[i] * value/tempr;

  if (m->n_marked == 0)
  {
    for(i = 0; i < 3; i++)
      m->xyz[i0][i] = m->xyz[i1][i] + r01[i];
      change_atom_parameters_in_list(i0);
  }
  else
  {
    for(i = 0; i < 3; i++)
      r01[i] = - m->xyz[i0][i] + m->xyz[i1][i] + r01[i];

    for(j = 0; j < m->n_marked; j++)
    {
      ii0 = m->marked[j];
      for(i = 0; i < 3; i++) m->xyz[ii0][i] += r01[i];
      change_atom_parameters_in_list(ii0);
    }
  }
  cumulative_move += fabs(value - tempr);

  if (cumulative_move > 0.1)
  {
    append_backup();
    cumulative_move = 0.0;
  }

  change_watched_data();
  if (Input_Data.automatic_rebonding) rebond();
  return;
}

void luscus_gtk_move_angle(double value)
{
  int i, j;
  int i0, i1, i2, ii0;
  double xi0[3], xi1[3], xi2[3], xi3[3];
  double newx, newy, newz;
  double A[3][3];
  double scale;
  int angle_return;

  int fri;
/**/

  if (m->n_selected != 3) return;

  fri = 0;

  i0 = m->selected[0];
  i1 = m->selected[1];
  i2 = m->selected[2];

  for(i = 0; i < 3; i++) 
  {
    xi0[i] = m->xyz[i0][i];
    xi1[i] = m->xyz[i1][i];
    xi2[i] = m->xyz[i2][i];
  }

  scale = 1.F;
  angle_return = move_angle(&newx, &newy, &newz, scale, xi0, xi1, xi2, A, value, 1);

  if (angle_return)
  {
    for(i = 0; i < 3; i++) /*move atom for small, random value*/
      xi0[i] = m->xyz[i0][i] = m->xyz[i0][i] + 0.001 * (double) rand() / (double) RAND_MAX;

    angle_return = move_angle(&newx, &newy, &newz, scale, xi0, xi1, xi2, A, value, 1);
  }

  if (m->n_marked == 0)
  {
    translate_to_origin(m->xyz[i0], m->xyz[i1], A);
    change_atom_parameters_in_list(i0);
  }
  else
  {
    for(j = 0; j < m->n_marked; j++)
    {
      ii0 = m->marked[j];
      for(i = 0; i < 3; i++) xi3[i] = m->xyz[ii0][i];
      translate_to_origin(xi3, xi1, A);
      for(i = 0; i < 3; i++) m->xyz[ii0][i] = xi3[i];

      change_atom_parameters_in_list(ii0);
    }
  }

  if (cumulative_move > 0.1)
  {
    append_backup();
    cumulative_move = 0.0;
  }

  change_watched_data();
  if (Input_Data.automatic_rebonding) rebond();
  return;
}

void luscus_gtk_move_torsion(double value)
{
  int i, j;
  int i0, i1, i2, i3, ii0;
  double xi0[3], xi1[3], xi2[3], xi3[3];
  double newx, newy, newz;
  double a, ang;
  double scale;
  
  if (m->n_selected != 4) return;

  i0=m->selected[0];
  i1=m->selected[1];
  i2=m->selected[2];
  i3=m->selected[3];

  for (i=0; i<3; i++)
  {
    xi0[i]=m->xyz[i0][i];
    xi1[i]=m->xyz[i1][i];
    xi2[i]=m->xyz[i2][i];
    xi3[i]=m->xyz[i3][i];
  }

  move_d_angle(&newx, &newy, &newz, &ang, scale,
               xi0, xi1, xi2, xi3, value, 1, 0.0, 0);
  m->xyz[i0][0] = newx;
  m->xyz[i0][1] = newy;
  m->xyz[i0][2] = newz;
  change_atom_parameters_in_list(i0);

  if (m->n_marked > 0)
  {
    a = ang;
    for(j = 0; j < m->n_marked; j++)
    {
      ii0 = m->marked[j];
      if ((ii0 != i1) && (ii0 != i2) && (ii0 != i3))
      {
        for(i = 0; i < 3; i++) xi0[i] = m->xyz[ii0][i];

        if (ii0 != i0)
        {
          move_d_angle(&newx, &newy, &newz, &ang, scale,
                       xi0, xi1, xi2, xi3, 0.0, 1, a, 1);
          m->xyz[ii0][0] = newx;
          m->xyz[ii0][1] = newy;
          m->xyz[ii0][2] = newz;
          change_atom_parameters_in_list(ii0);
        }
      }
    }
  }

  if (cumulative_move > 0.1)
  {
    append_backup();
    cumulative_move = 0.0;
  }

  change_watched_data();
  if (Input_Data.automatic_rebonding) rebond();
  return;
}

void luscus_gtk_move_coord(double value)
{
  double scale;
  double r01[3];
  int i0, i1, i2, i3, i,j, ii0;
  double xi0[3], xi1[3], xi2[3], xi3[3];
  double xyz[4][3];
  double newx, newy, newz;
  double tempr;
  double a, ang;
  double A[3][3];
  int angle_return;
  double dlt;
  int fri;

  if (m->n_selected > 4) return;
  else if (m->n_selected == 2)
  {
     i0=m->selected[0];
     i1=m->selected[1];

     for(i = 0; i < 3; i++) r01[i] = m->xyz[i0][i] - m->xyz[i1][i];
     tempr=sqrt(r01[0]*r01[0]+r01[1]*r01[1]+r01[2]*r01[2]);

     for(i = 0; i < 3; i++) r01[i]=r01[i] * value/tempr;

     if (m->n_marked == 0)
     {
       for(i = 0; i < 3; i++)
         m->xyz[i0][i] = m->xyz[i1][i] + r01[i];
       change_atom_parameters_in_list(i0);
     }
     else
     {
       for(i = 0; i < 3; i++)
         r01[i] = - m->xyz[i0][i] + m->xyz[i1][i] + r01[i];

       for(j = 0; j < m->n_marked; j++)
       {
         ii0 = m->marked[j];
         for(i = 0; i < 3; i++) m->xyz[ii0][i] += r01[i];
       }
       change_atom_parameters_in_list(ii0);
     }
     cumulative_move += fabs(value - tempr);
/*     tempr=sqrt(r01[0]*r01[0]+r01[1]*r01[1]+r01[2]*r01[2]);*/
  }
  else if (m->n_selected == 3)
  {
    fri = 0;

    i0 = m->selected[0];
    i1 = m->selected[1];
    i2 = m->selected[2];

    for(i = 0; i < 3; i++)
    {
      xi0[i] = m->xyz[i0][i];
      xi1[i] = m->xyz[i1][i];
      xi2[i] = m->xyz[i2][i];
    }

    scale = 1.F;
    angle_return = move_angle(&newx, &newy, &newz, scale, xi0, xi1, xi2, A, value, 1);

    if (angle_return)
    {
      for(i = 0; i < 3; i++) /*move atom for small, random value*/
        xi0[i] = m->xyz[i0][i] = m->xyz[i0][i] + 0.001 * (double) rand() / (double) RAND_MAX;

      angle_return = move_angle(&newx, &newy, &newz, scale, xi0, xi1, xi2, A, value, 1);
    }

/*    if (angle_return==1)
    {
      dlt = (xi0[1]-xi1[1])-(xi0[0]-xi1[0]);
      xi0[1] -= dlt;
      xi0[0] += dlt;

      if (ang==0)
      {
        ang = 1.0;
        if (scale<1) ang = -1.0;
      }
      angle_return = move_angle(&newx, &newy, &newz, scale, xi0, xi1, xi2, A, ang, 1);
    }

    if (angle_return == 2)
    {
      dlt = (xi0[1]-xi1[1])-(xi0[0]-xi1[0]);

      xi0[1] -= dlt;
      xi0[0] += dlt;

      if (ang==0)
      {
        ang = 181.0;
        if (scale < 1) ang = 179.0;
      }
      angle_return = move_angle(&newx, &newy, &newz, scale, xi0, xi1, xi2, A, ang, 1);
    }*/

    if (m->n_marked == 0)
    {
      translate_to_origin(m->xyz[i0], m->xyz[i1], A);
/*      for(i = 0; i < 3; i++) m->xyz[i0][i] = xi0[i];*/
      change_atom_parameters_in_list(i0);
    }
    else
    {
      for(j = 0; j < m->n_marked; j++)
      {
        ii0 = m->marked[j];
        for(i = 0; i < 3; i++) xi3[i] = m->xyz[ii0][i];
        translate_to_origin(xi3, xi1, A);
        for(i = 0; i < 3; i++) m->xyz[ii0][i] = xi3[i];
      }
      change_atom_parameters_in_list(ii0);
    }

  }
  else if (m->n_selected == 4)
  {
    i0=m->selected[0];
    i1=m->selected[1];
    i2=m->selected[2];
    i3=m->selected[3];

    for (i=0; i<3; i++)
    {
      xi0[i]=m->xyz[i0][i];
      xi1[i]=m->xyz[i1][i];
      xi2[i]=m->xyz[i2][i];
      xi3[i]=m->xyz[i3][i];
    }

    move_d_angle(&newx, &newy, &newz, &ang, scale,
                 xi0, xi1, xi2, xi3, value, 1, 0.0, 0);
    m->xyz[i0][0] = newx;
    m->xyz[i0][1] = newy;
    m->xyz[i0][2] = newz;
    change_atom_parameters_in_list(i0);

    if (m->n_marked > 0)
    {
      a = ang;
      for(j = 0; j < m->n_marked; j++)
      {
        ii0 = m->marked[j];
        if ((ii0 != i1) && (ii0 != i2) && (ii0 != i3))
        {
          for(i = 0; i < 3; i++) xi0[i] = m->xyz[ii0][i];

          if (ii0 != i0)
          {
            move_d_angle(&newx, &newy, &newz, &ang, scale,
                         xi0, xi1, xi2, xi3, 0.0, 1, a, 1);
            m->xyz[ii0][0] = newx;
            m->xyz[ii0][1] = newy;
            m->xyz[ii0][2] = newz;
          }
        }
        change_atom_parameters_in_list(ii0);
      }
    }
  }

  if (cumulative_move > 0.1)
  {
    append_backup();
    cumulative_move = 0.0;
  }

  change_watched_data();
  if (Input_Data.automatic_rebonding) rebond();
  return;
}

void luscus_gtk_move_coord_step(int isign)
{
  if (abs(isign) != 1) return;
  if (m->n_selected == 2)
  {
    luscus_gtk_move_coord(isign * Input_Data.translation_value + get_selected_bond_length());
  }
  else if (m->n_selected == 3)
    luscus_gtk_move_coord(isign * Input_Data.angle_change_value + get_selected_angle_value());
  else if (m->n_selected == 4)
    luscus_gtk_move_coord(isign * Input_Data.torsion_change_value + get_selected_dihedral_value());
  return;
}

int move_angle(double *newx, double *newy, double *newz,
               double scale,
               XYZ xi0, XYZ xi1, XYZ xi2,
               double A[3][3], double ang, int fri)
{
  int i, j, k;
  double delta, theta;
  double deter, al;
  double matx[3];
  double A_Matrix[3][3], A_Inverted[3][3], tmpM[3][3];
  double xyz[4][3];
  double rBA[3], rBC[3];
  double rax[3][3], ray[3][3], raz[3][3];

  double angz, angx;

  /* Set to old point in case of failure */

  (*newx) = xi0[0];
  (*newy) = xi0[1];
  (*newz) = xi0[2];

  /* compute the current angle */

  for(i=0; i<3; i++)
  {
    xyz[0][i]=xi0[i];
    xyz[1][i]=xi1[i];
    xyz[2][i]=xi2[i];
  }

  theta = get_selected_angle_value();  /* old angle */
  if (!finite(theta)) return 1;

  if (ang==0)
  {
    if (fri==1) delta=-(theta-ang) * M_PI/180.0e0;
    else delta = M_PI/180.0e0;
  }
  else
  {
    delta = -(theta-ang) * M_PI/180.0e0;
  }

  if (ang==0 && fri==0)
  {
    if(scale>1) delta=delta;    /* angle increase/decrease, 1 degree */
    else delta=-delta;
  }

  /* Translate the system */
  /* Save the old value of atom B */

  /* Normal vector */
  rBA[0] = xi0[0]-xi1[0]; rBA[1] = xi0[1]-xi1[1]; rBA[2] = xi0[2]-xi1[2];
  rBC[0] = xi2[0]-xi1[0]; rBC[1] = xi2[1]-xi1[1]; rBC[2] = xi2[2]-xi1[2];

  matx[0] = rBA[1]*rBC[2]-rBA[2]*rBC[1];
  matx[1] = rBA[2]*rBC[0]-rBA[0]*rBC[2];
  matx[2] = rBA[0]*rBC[1]-rBA[1]*rBC[0];

  /* lets normalize n */
  al = sqrt(matx[0]*matx[0]+matx[1]*matx[1]+matx[2]*matx[2]);

  matx[0] = matx[0]/al;
  matx[1] = matx[1]/al;
  matx[2] = matx[2]/al;

  /* Angle 1 */
  if(sqrt(matx[0]*matx[0]+matx[1]*matx[1]) < 1e-4)
    angz = 90.0e0 * M_PI/180.0e0;
  else
    angz = acos(matx[1]/sqrt(matx[0]*matx[0]+matx[1]*matx[1]));

  /* Angle 2 */
  angx = acos(sqrt(matx[0]*matx[0]+matx[1]*matx[1])/
              sqrt(matx[0]*matx[0]+matx[1]*matx[1]+matx[2]*matx[2]));

  if (matx[2]>0.0) angx= angx;
  else angx = -angx;

  if (matx[0]>0.0) angz = -angz;
  else angz = angz;

  /* Matrix for rotation around x axis */
  rax[0][0] = 1.0;  rax[0][1] =        0.0;  rax[0][2] =         0.0;
  rax[1][0] = 0.0;  rax[1][1] =  cos(angx);  rax[1][2] =   sin(angx);
  rax[2][0] = 0.0;  rax[2][1] = -sin(angx);  rax[2][2] =   cos(angx);

  /* Matrix for rotation around y axis */
  ray[0][0] = cos(delta);  ray[0][1] = 0.0;  ray[0][2] = -sin(delta);
  ray[1][0] =        0.0;  ray[1][1] = 1.0;  ray[1][2] =         0.0;
  ray[2][0] = sin(delta);  ray[2][1] = 0.0;  ray[2][2] =  cos(delta);

  /* Matrix for rotation around z axis */
  raz[0][0] =  cos(angz);  raz[0][1] = sin(angz);  raz[0][2] = 0.0;
  raz[1][0] = -sin(angz);  raz[1][1] = cos(angz);  raz[1][2] = 0.0;
  raz[2][0] =        0.0;  raz[2][1] =       0.0;  raz[2][2] = 1.0;

  /* make sure that there is no info in A or tmpM */
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
    {
      A[i][j] = 0.0;
      tmpM[i][j] = 0.0;
      A_Matrix[i][j] = 0.0;
    }

   /* First raz*rax */
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++)
        A_Matrix[i][j] = A_Matrix[i][j]+rax[i][k]*raz[k][j];

  det_matrix(&deter, A_Matrix);

  if(!finite(deter))
  {
/*    luscus_gtk_push_message_to_statusbar1("NaN in determinant, will try to change anyway");*/
    /* We will return identity matrix */
    A[0][0] = 1.0; A[0][1] = 0.0;  A[0][2] = 0.0;
    A[1][0] = 0.0; A[1][1] = 1.0;  A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0;  A[2][2] = 1.0;

    if (theta==0) return 1;
    if (theta==180) return 2;

    return 3;
  }

  if(fabs(deter)<1.0e-4)
  {
/*    luscus_gtk_push_message_to_statusbar1("Cannot move it due to determinant equal zero");*/
   /* We will return identity matrix */
    A[0][0] = 1.0; A[0][1] = 0.0;  A[0][2] = 0.0;
    A[1][0] = 0.0; A[1][1] = 1.0;  A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0;  A[2][2] = 1.0;
    return 3;
  }    /* to avoid division by 0 */

  inv_matrix(&A_Inverted[0][0], &A_Inverted[0][1], &A_Inverted[0][2],
             &A_Inverted[1][0], &A_Inverted[1][1], &A_Inverted[1][2],
             &A_Inverted[2][0], &A_Inverted[2][1], &A_Inverted[2][2],
             A_Matrix);

  /* First matrix product */
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++)
        tmpM[i][j] = tmpM[i][j]+A_Inverted[i][k]*ray[k][j];
  
  /* Second matrix product */
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++)
        A[i][j] = A[i][j]+tmpM[i][k]*A_Matrix[k][j];
  
  (*newx) = xi0[0];
  (*newy) = xi0[1];
  (*newz) = xi0[2];

  cumulative_move += fabs(theta-ang) * M_PI / 90.0; 

  return 0;
}

void det_matrix (double *d, double mm[3][3])
{
 double mm0, mm1, mm2;
 
 mm0 = mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1];
 mm1 = mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0];
 mm2 = mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0];
 
 *d = mm[0][0] * mm0 - mm[0][1] * mm1 + mm[0][2] * mm2;
 
 return;
}


void inv_matrix (double *inv00, double *inv01, double *inv02,
                 double *inv10, double *inv11, double *inv12, 
                 double *inv20, double *inv21, double *inv22,
                 double mm[3][3])
{
  double mm00, mm10, mm20,
         mm01, mm11, mm21,
         mm02, mm12, mm22,
         d;

  mm00 = mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1];
  mm10 = mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0];
  mm20 = mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0];

  mm01 = mm[0][1] * mm[2][2] - mm[0][2] * mm[2][1];
  mm11 = mm[0][0] * mm[2][2] - mm[0][2] * mm[2][0];
  mm21 = mm[0][0] * mm[2][1] - mm[0][1] * mm[2][0];

  mm02 = mm[0][1] * mm[1][2] - mm[0][2] * mm[1][1];
  mm12 = mm[0][0] * mm[1][2] - mm[0][2] * mm[1][0];
  mm22 = mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0];

  d = mm[0][0] * mm00 - mm[0][1] * mm10 + mm[0][2] * mm20;

  d = 1/d;

  *inv00 = mm00 * d;
  *inv02 = mm02 * d;
  *inv11 = mm11 * d;
  *inv20 = mm20 * d;
  *inv22 = mm22 * d;

  d = -1 * d;

  *inv01 = mm01 * d;
  *inv10 = mm10 * d;
  *inv12 = mm12 * d;
  *inv21 = mm21 * d;
 
  return;
}

void move_d_angle(double *newx, double *newy, double *newz, double *ang,
                  double scale,
                  double xi0[3], double xi1[3], double xi2[3], double xi3[3],
                  double setangle, int ifsetangle, double fixedmove, int ifixedmove)
{
  double theta, alpha;
  double rr, t1, t2;

  double a,b,c,d,e,f;
  double d1,d2;
  double ux1, ux2, ux3, uy1, uy2, uy3, uz1, uz2, uz3;
  double x2, y2, z2, radius, cosa, sinphi, cosphi, cdx, cdy, cdz;
  double delta;

  double r12[3], r23[3], r34[3], o[3], p[3], q[3];
  double cosphiNy, sinphiNy, l;
 
  /*in the case of failure, coordinates should be unchanged!*/
  *newx=xi0[0];
  *newy=xi0[1];
  *newz=xi0[2];  

  if(scale>1) delta = 1;           /* angle increase/decrease */
  else delta = -1;
  delta = delta * M_PI/180.0e0;

/* calculate old theta */ 
  r12[0] = xi0[0] - xi1[0];  r12[1] = xi0[1] - xi1[1];  r12[2] = xi0[2] - xi1[2];
  r23[0] = xi1[0] - xi2[0];  r23[1] = xi1[1] - xi2[1];  r23[2] = xi1[2] - xi2[2];
  r34[0] = xi2[0] - xi3[0];  r34[1] = xi2[1] - xi3[1];  r34[2] = xi2[2] - xi3[2];

  o[0] = r12[1]*r23[2]-r12[2]*r23[1]; 
  o[1] = r12[2]*r23[0]-r12[0]*r23[2]; 
  o[2] = r12[0]*r23[1]-r12[1]*r23[0];

  l = sqrt(o[0]*o[0]+o[1]*o[1]+o[2]*o[2]);
  if (l < 0.01)
  {
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2("ERROR: dihedral undefined");
    printf("ERROR: dihedral undefined\n");
    return;
  }
  o[0] = o[0]/l; o[1] = o[1]/l; o[2] = o[2]/l;

  p[0] = r23[1]*r34[2]-r23[2]*r34[1]; 
  p[1] = r23[2]*r34[0]-r23[0]*r34[2]; 
  p[2] = r23[0]*r34[1]-r23[1]*r34[0];

  l = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  if (l < 0.01)
  {
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2("ERROR: dihedral undefined");
    printf("ERROR: dihedral undefined\n");
    return;
  }
  p[0] = p[0]/l; p[1] = p[1]/l; p[2] = p[2]/l;

  q[0] = r23[1]*o[2]-r23[2]*o[1]; 
  q[1] = r23[2]*o[0]-r23[0]*o[2]; 
  q[2] = r23[0]*o[1]-r23[1]*o[0];

  l = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
  if (l < 0.01)
  {
    luscus_gtk_pop_message_from_statusbar2();
    luscus_gtk_push_message_to_statusbar2("ERROR: dihedral undefined");
    printf("ERROR: dihedral undefined\n");
    return;
  }
  q[0] = q[0]/l; q[1] = q[1]/l; q[2] = q[2]/l; 

  cosphiNy = (o[0]*p[0]+o[1]*p[1]+o[2]*p[2]);
  sinphiNy = (q[0]*p[0]+q[1]*p[1]+q[2]*p[2]);

  theta = -atan2(sinphiNy, cosphiNy);

/* calculate alpha */  
  rr = rr_vec(xi0,  xi1,  xi2);

  t1 = rr_dist(xi0, xi1);
  t2 = rr_dist(xi2, xi1);

  /*avoid numerical errors!*/
  if (rr/(t1*t2) >= 1.0) alpha = 0.0;
  else if (rr/(t1*t2) <= -1.0) alpha = M_PI;
  else alpha=acos(rr/(t1*t2));

/* calculate dist */
  
  rr = rr_dist(xi0, xi1);
/*  rr=sqrt(t1);*/

/*  printf("mes %f %f %f \n", rr, alpha*180/M_PI, theta*180/M_PI); */

  if (ifsetangle==1)
  {
    if (ifixedmove!=0)
    {
      theta=theta+fixedmove;
    }
    else
    {
      *ang = -theta-setangle * M_PI/180.0e0;      
      theta = -setangle * M_PI/180.0e0;
    }
  }
  else
  {
    theta=theta+delta;
    *ang = delta;
  }

/* hope that connectivity is ok */

  a=xi1[0]-xi2[0];
  b=xi1[1]-xi2[1];
  c=xi1[2]-xi2[2];
  d1=sqrt(a*a+b*b+c*c);

  if(fabs(d1) < 0.01)
  {
    ux1=0; uy1=0; uz1=0;
  }
  else
  {
    ux1=a/d1; uy1=b/d1; uz1=c/d1;
  }

  d=xi3[0]-xi2[0];
  e=xi3[1]-xi2[1];
  f=xi3[2]-xi2[2];

  x2 = b*f-c*e;
  y2 = c*d-a*f;
  z2 = a*e-b*d;

  d2 = sqrt(x2*x2+y2*y2+z2*z2);

  if(fabs(d2) < 0.01)
  {
    ux2=0; uy2=0; uz2=0;
  }
  else
  {
    ux2=x2/d2; uy2=y2/d2; uz2=z2/d2;
  }

  ux3 = uy2*uz1-uz2*uy1;
  uy3 = uz2*ux1-ux2*uz1;
  uz3 = ux2*uy1-uy2*ux1;

  radius = rr*sin(M_PI-alpha);
  cosa   = cos(M_PI-alpha);
  sinphi = sin(theta);
  cosphi = cos(theta);
  cdx = rr*cosa*ux1+radius*(ux2*sinphi+ux3*cosphi);
  cdy = rr*cosa*uy1+radius*(uy2*sinphi+uy3*cosphi);
  cdz = rr*cosa*uz1+radius*(uz2*sinphi+uz3*cosphi);

  *newx=xi1[0]+cdx;
  *newy=xi1[1]+cdy;
  *newz=xi1[2]+cdz;

  if (*ang < -M_PI) *ang += 2.0 * M_PI;
  if (*ang > M_PI) *ang -= 2.0 * M_PI;
  cumulative_move += fabs(*ang * 3.0); 

  return;
}

void translate(XYZ trxyz)
{
  int i, j;
  int ii0;
  if (m->n_selected)
  {
    for(j = 0; j < 3; j++)
      m->xyz[m->selected[0]][j] += trxyz[j];
  }
  else if (m->n_marked)
  {
    for(i = 0; i < m->n_marked; i++)
    {
      ii0 = m->marked[i];
      for(j = 0; j < 3; j++) m->xyz[ii0][j] += trxyz[j];
    }
  }
  else
  {
    for(i = 0; i < m->natom; i++)
      for(j = 0; j < 3; j++)
        m->xyz[i][j] += trxyz[j];

    for(i = 0; i < m->nvector; i++)
      for(j = 0; j < 3; j++)
      {
        m->vector1[i][j] += trxyz[j];
        m->vector2[i][j] += trxyz[j];
      }

  for(i = 0; i < m->ntriangle; i++)
    for(j = 0; j < 3; j++)
    {
      m->triangle1[i][j] += trxyz[j];
      m->triangle2[i][j] += trxyz[j];
      m->triangle3[i][j] += trxyz[j];
    }

  for(i = 0; i < m->nsphere; i++)
    for(j = 0; j < 3; j++)
      m->sphere_center[i][j] += trxyz[j];

  for(i = 0; i < m->nsurf; i++)
    for(j = 0; j < 3; j++)
    {
      m->surf1[i][j] += trxyz[j];
      m->surf2[i][j] += trxyz[j];
      m->surf3[i][j] += trxyz[j];
    }

  for(i = 0; i < m->ncells; i++)
    for(j = 0; j < 3; j++)
    {
      m->cell1[i][j] += trxyz[j];
      m->cell2[i][j] += trxyz[j];
      m->cell3[i][j] += trxyz[j];
      m->cell4[i][j] += trxyz[j];
    }

    if (m->ngrids)
    {
      for(j = 0; j < 3; j++)
        m->grid.origin[j] += trxyz[j];
      make_surfaces();
    }
  }
  rerender_3d();
}

void translate_to_origin(XYZ xi0, XYZ xi1, double A[3][3])
{
  double a, b, c;

  /* Translate all atoms selected to the computed origin */
  xi0[0] = xi0[0] - xi1[0];
  xi0[1] = xi0[1] - xi1[1];
  xi0[2] = xi0[2] - xi1[2];

  /* Rotate all atoms selected for rotation */
  a = A[0][0] * xi0[0] + A[0][1] * xi0[1] + A[0][2] * xi0[2];
  b = A[1][0] * xi0[0] + A[1][1] * xi0[1] + A[1][2] * xi0[2];
  c = A[2][0] * xi0[0] + A[2][1] * xi0[1] + A[2][2] * xi0[2];
  xi0[0] = a; xi0[1] = b; xi0[2] = c;

  /* Translate back to correct origin */
  xi0[0] = xi0[0]+xi1[0];
  xi0[1] = xi0[1]+xi1[1];
  xi0[2] = xi0[2]+xi1[2];
  return;

}

void rr_cross (double *rr1, double *rr2, double *rr3,
               double x0, double x1, double x2,
               double y0, double y1, double y2,
               double z0, double z1, double z2)
{
 double rx1, ry1, rz1, rx2, ry2, rz2;
 rx1=x0-y0;
 ry1=x1-y1;
 rz1=x2-y2;
 
 rx2=z0-y0;
 ry2=z1-y1;
 rz2=z2-y2;

 *rr1=ry1*rz2-rz1*ry2;
 *rr2=rz1*rx2-rx1*rz2;
 *rr3=rx1*ry2-ry1*rx2;
 
 return;
}

double rr_vec(XYZ x0, XYZ x1, XYZ x2)
{
  return (x0[0] - x1[0]) * (x2[0] - x1[0]) +
         (x0[1] - x1[1]) * (x2[1] - x1[1]) +
	 (x0[2] - x1[2]) * (x2[2] - x1[2]);
}

double rr_dist(XYZ a1, XYZ a0)
{
  return sqrt((a1[0] - a0[0]) *
              (a1[0] - a0[0]) +
              (a1[1] - a0[1]) *
              (a1[1] - a0[1]) + 
              (a1[2] - a0[2]) *
              (a1[2] - a0[2]));
}

int find_bond(MOL *mp, int iat1, int iat2)
{
  int i;

  for(i = 0; i < mp->nbond; i++)
    if (iat1 == mp->bond[i].iat1 && iat2 == mp->bond[i].iat2) return i;
    else if (iat2 == mp->bond[i].iat1 && iat1 == mp->bond[i].iat2) return i;
  return -1;
}

int find_one_bond_from_atom(int iat)
{
  int i;

  for (i = 0; i < m->nbond; i++)
    if (iat == m->bond[i].iat1) return i;
    else if (iat == m->bond[i].iat2) return i;
  return -1;
}

void CalcCyllinder(double in[3], double ut[3], double *angle, double *vn0,  double *vn1, double *vl)
{
  double cosa, r;
  double v[3];
  v[0] = ut[0] - in[0];
  v[1] = ut[1] - in[1];
  v[2] = ut[2] - in[2];
   
  *vl=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  cosa=(v[2])/(*vl);

  if (fabs(cosa) == 1.0e0)
  {
    *vn0 = 0.0;
    *vn1 = 1.0;
    *angle = 0.0;
    if (cosa < 0) *angle = 180.0e0;

  }
  else
  {
    *angle = (double) acos(cosa)*180.0e0/M_PI;
    r=sqrt( v[0]*v[0]+v[1]*v[1]);
    *vn0 = -v[1]/r;
    *vn1 =  v[0]/r;
  }
}

void Calc2Cyllinder(double in0, double in1,  double in2, 
                    double ut0, double ut1, double ut2, 
                    double c0, double c1, double c2, 
                    double *vn0,  double *vn1, double *vn2)
{
  double v[3];
  double vl;
  v[0] = ut0 - in0;
  v[1] = ut1 - in1;
  v[2] = ut2 - in2;
  
  rr_cross(vn0, vn1, vn2, in0, in1, in2, ut0, ut1, ut2, c0, c1, c2);
  
  if (*vn0 + *vn1 + *vn2 < 0.0001) *vn1 = 1.0;
  
  vl=sqrt((*vn0)*(*vn0)+(*vn1)*(*vn1)+(*vn2)*(*vn2));
  if(vl<0.001) vl=0.01;
  *vn0=*vn0/(vl*15);
  *vn1=*vn1/(vl*15);
  *vn2=*vn2/(vl*15);
      
  return;     

}

double Calc_Diameter(void)
{
  int i, j, np;
  double d, dmax = 0.0, diam = 0.0;
  double cent[3] = {0.0, 0.0, 0.0};

  if (m->ngrids)
  {
    d = sqrt((m->grid.axisvec1[0]-m->grid.origin[0])*(m->grid.axisvec1[0]-m->grid.origin[0]) +
             (m->grid.axisvec1[1]-m->grid.origin[1])*(m->grid.axisvec1[1]-m->grid.origin[1]) +
             (m->grid.axisvec1[2]-m->grid.origin[2])*(m->grid.axisvec1[2]-m->grid.origin[2]));
    if (d > dmax) dmax = d;
    d = sqrt((m->grid.axisvec2[0]-m->grid.origin[0])*(m->grid.axisvec2[0]-m->grid.origin[0]) +
             (m->grid.axisvec2[1]-m->grid.origin[1])*(m->grid.axisvec2[1]-m->grid.origin[1]) +
             (m->grid.axisvec2[2]-m->grid.origin[2])*(m->grid.axisvec2[2]-m->grid.origin[2]));
    if (d > dmax) dmax = d;
    d = sqrt((m->grid.axisvec3[0]-m->grid.origin[0])*(m->grid.axisvec3[0]-m->grid.origin[0]) +
             (m->grid.axisvec3[1]-m->grid.origin[1])*(m->grid.axisvec3[1]-m->grid.origin[1]) +
             (m->grid.axisvec3[2]-m->grid.origin[2])*(m->grid.axisvec3[2]-m->grid.origin[2]));
    if (d > dmax) dmax = d;
  }

  for(i = 0; i < m->natom; i++)
  {
    d = sqrt(m->xyz[i][0]*m->xyz[i][0] +
             m->xyz[i][1]*m->xyz[i][1] +
             m->xyz[i][2]*m->xyz[i][2]);
    if (d > dmax) dmax = d;
  }

  for(i = 0; i < m->nsphere; i++)
  {
    d = sqrt(m->sphere_center[i][0]*m->sphere_center[i][0] +
             m->sphere_center[i][1]*m->sphere_center[i][1] +
             m->sphere_center[i][2]*m->sphere_center[i][2]) + m->sphere_radius[i];
    if (d > dmax) dmax = d;
  }

  for(i = 0; i < m->nvector; i++)
  {
    d = sqrt(m->vector1[i][0]*m->vector1[i][0] +
             m->vector1[i][1]*m->vector1[i][1] +
             m->vector1[i][2]*m->vector1[i][2]);
    if (d > dmax) dmax = d;
    d = sqrt(m->vector2[i][0]*m->vector2[i][0] +
             m->vector2[i][1]*m->vector2[i][1] +
             m->vector2[i][2]*m->vector2[i][2]);
    if (d > dmax) dmax = d;
  }

  for(i = 0; i < m->ntriangle; i++)
  {
    d = sqrt(m->triangle1[i][0]*m->triangle1[i][0] +
             m->triangle1[i][1]*m->triangle1[i][1] +
             m->triangle1[i][2]*m->triangle1[i][2]);
    if (d > dmax) dmax = d;
    d = sqrt(m->triangle2[i][0]*m->triangle2[i][0] +
             m->triangle2[i][1]*m->triangle2[i][1] +
             m->triangle2[i][2]*m->triangle2[i][2]);
    if (d > dmax) dmax = d;
    d = sqrt(m->triangle3[i][0]*m->triangle3[i][0] +
             m->triangle3[i][1]*m->triangle3[i][1] +
             m->triangle3[i][2]*m->triangle3[i][2]);
    if (d > dmax) dmax = d;
  }

  for(i = 0; i < m->ncells; i++)
  {
    d = sqrt(m->cell1[i][0]*m->cell1[i][0] +
             m->cell1[i][1]*m->cell1[i][1] +
             m->cell1[i][2]*m->cell1[i][2]);
    if (d > dmax) dmax = d;
    d = sqrt(m->cell2[i][0]*m->cell2[i][0] +
             m->cell2[i][1]*m->cell2[i][1] +
             m->cell2[i][2]*m->cell2[i][2]);
    if (d > dmax) dmax = d;
    d = sqrt(m->cell3[i][0]*m->cell3[i][0] +
             m->cell3[i][1]*m->cell3[i][1] +
             m->cell3[i][2]*m->cell3[i][2]);
    if (d > dmax) dmax = d;
    d = sqrt(m->cell4[i][0]*m->cell4[i][0] +
             m->cell4[i][1]*m->cell4[i][1] +
             m->cell4[i][2]*m->cell4[i][2]);
    if (d > dmax) dmax = d;
  }

  if (dmax < 0.1) diam = 1.0; /*arbitrary value*/
  else diam = 1.5 * dmax;

  return diam;
}

void pl3to4(double in1[3],double in2[3],double in3[3], double ut1[3], double ut2[3], double ut3[3], double ut4[3])
{
  double in4[3], cen[3]/*, r[3], diam, d*/, v1[3], v2[3], d1, d2, d, dot[3];
  int i;
/*  in4[0] = in2[0] + in3[0] - in1[0];
  in4[1] = in2[1] + in3[1] - in1[1];
  in4[2] = in2[2] + in3[2] - in1[2];

  r[0] = (in2[0] + in3[0]) / 2.0;
  r[1] = (in2[1] + in3[1]) / 2.0;
  r[2] = (in2[2] + in3[2]) / 2.0;

  diam=(int)(Calc_Diameter()/2 + 1);
  d=(sqrt((r[0]-in1[0]) * (r[0]-in1[0]) + (r[1]-in1[1]) * (r[1]-in1[1]) + (r[2]-in1[2]) * (r[2]-in1[2]))+
     sqrt((r[0]-in2[0]) * (r[0]-in2[0]) + (r[1]-in2[1]) * (r[1]-in2[1]) + (r[2]-in3[2]) * (r[2]-in3[2])))/2;
  d=diam/d*1.5;
  for(i = 0; i < 3; i++)
  {
    ut1[i] = r[i] + (in1[i] - r[i]) * d; 
    ut2[i] = r[i] + (in2[i] - r[i]) * d; 
    ut3[i] = r[i] + (in3[i] - r[i]) * d; 
    ut4[i] = r[i] + (in4[i] - r[i]) * d; 
  }*/
  d = 0.0;
  d1 = 0.0;
  d2 = 0.0;
  for (i = 0; i < 3; i++)
  {
    cen[i] = 0.0;/*2016.09.15*/
    v1[i] = in2[i] - in1[i];
    v2[i] = in3[i] - in1[i];
    d  += v1[i]*v2[i];
    d1 += v1[i]*v1[i];
    d2 += v2[i]*v2[i];
  }
  d /= (d1);
  d1 = sqrt(d1);
  d2 = sqrt(d2);
  d2 = 0.0;
  for (i = 0; i < 3; i++)
  {
    cen[i] += in1[i];/*2016.09.15*/
    cen[i] += in2[i];/*2016.09.15*/
    cen[i] += in3[i];/*2016.09.15*/
    v2[i] -= d * v1[i];
    v1[i] /= d1;
    d2 += v2[i]*v2[i];
/*    v2[i] *= d / d2;*/
  }
  d = 0.5 * (d1 + d2);
  d2 = sqrt(d2);
  for(i = 0; i < 3; i++)
  {
    v1[i] *= d;
    v2[i] *= (d/d2);
    cen[i] /= 3.0;/*2016.09.15*/
  }
  for (i = 0; i < 3; i++) /*2014.01.23*//*2016.09.15*/
  {
/*    ut1[i] = in1[i] - v1[i] - v2[i];
    ut2[i] = in1[i] - v1[i] + v2[i];
    ut3[i] = in1[i] + v1[i] + v2[i];
    ut4[i] = in1[i] + v1[i] - v2[i];*/
    ut1[i] = cen[i] - v1[i] - v2[i];
    ut2[i] = cen[i] - v1[i] + v2[i];
    ut3[i] = cen[i] + v1[i] + v2[i];
    ut4[i] = cen[i] + v1[i] - v2[i];
  }
  return;
}

void pl4to6(double in0[3],  double in1[3],  double in2[3],  double in3[3],
            double in12[3], double in13[3], double in23[3], double in123[3])
{
  int i; 
  for(i = 0; i < 3; i++)
  {
    in12[i]  = in1[i] + in2[i] - in0[i];
    in13[i]  = in1[i] + in3[i] - in0[i];
    in23[i]  = in2[i] + in3[i] - in0[i];
    in123[i] = in1[i] + in2[i] + in3[i] - 2 * in0[i]; 
  }
  return;
}

void np3(double p1[3], double p2[3], double p3[3], double n[3])
{
  /*corrected cross product 2014.01.24*/
  n[0] = p1[0] + (p2[1] - p1[1])*(p3[2] - p1[2]) - (p2[2] - p1[2])*(p3[1] - p1[1]);
  n[1] = p1[1] + (p2[2] - p1[2])*(p3[0] - p1[0]) - (p2[0] - p1[0])*(p3[2] - p1[2]);
  n[2] = p1[2] + (p2[0] - p1[0])*(p3[1] - p1[1]) - (p2[1] - p1[1])*(p3[0] - p1[0]);
/*  n[0] = p1[1]*(p2[2] - p3[2]) + p2[1]*(p3[2] - p1[2]) + p3[1]*(p1[2] - p2[2]);
  n[1] = p1[2]*(p2[0] - p3[0]) + p2[2]*(p3[0] - p1[0]) + p3[2]*(p1[0] - p2[0]);
  n[2] = p1[0]*(p2[1] - p3[1]) + p2[0]*(p3[1] - p1[1]) + p3[0]*(p1[1] - p2[1]);*/
}

void swap_atoms(int ati, int atj)
{
  int i;
  XYZ tmpxyz;
  ELEM_DATA tmpe;
  double tmp_charg_m, tmp_charg_l;
  int tmp_add_numer, tmp_symm;
  char *tmp_name;

  if (ati >= m->natom || atj >= m->natom) return;

  /*atj -> tmp*/
  for(i = 0; i < 3; i++) tmpxyz[i] = m->xyz[atj][i];
  tmpe.name = m->elem[atj].name;
  tmpe.vdw_rad = m->elem[atj].vdw_rad;
  tmpe.bond_rad = m->elem[atj].bond_rad;
  tmpe.valency = m->elem[atj].valency;
  for(i = 0; i < 4; i++) tmpe.color[i] = m->elem[atj].color[i];
  tmp_charg_m = m->charge_m[atj];
  tmp_charg_l = m->charge_l[atj];
  tmp_add_numer = m->additional_numeration[atj];
  tmp_symm = m->symmetry[atj];
  tmp_name = m->name[atj];

  /*ati -> atj*/
  for(i = 0; i < 3; i++) m->xyz[atj][i] = m->xyz[ati][i];
  m->elem[atj].name = m->elem[ati].name;
  m->elem[atj].vdw_rad = m->elem[ati].vdw_rad;
  m->elem[atj].bond_rad = m->elem[ati].bond_rad;
  m->elem[atj].valency = m->elem[ati].valency;
  for(i = 0; i < 4; i++) m->elem[atj].color[i] = m->elem[ati].color[i];
  m->charge_m[atj] = m->charge_m[ati];
  m->charge_l[atj] = m->charge_l[ati];
  m->additional_numeration[atj] = m->additional_numeration[ati];
  m->symmetry[atj] = m->symmetry[ati];
  m->name[atj] = m->name[ati];

  /*tmp->ati*/
  for(i = 0; i < 3; i++) m->xyz[ati][i] = tmpxyz[i];
  m->elem[ati].name = tmpe.name;
  m->elem[ati].vdw_rad = tmpe.vdw_rad;
  m->elem[ati].bond_rad = tmpe.bond_rad;
  m->elem[ati].valency = tmpe.valency;
  for(i = 0; i < 4; i++) m->elem[ati].color[i] = tmpe.color[i];
  m->charge_m[ati] = tmp_charg_m;
  m->charge_l[ati] = tmp_charg_l;
  m->additional_numeration[ati] = tmp_add_numer;
  m->symmetry[ati] = tmp_symm;
  m->name[ati] = tmp_name;

  /*correct bonding*/
  if (m->nbond)
    for (i = 0; i < m->nbond; i++)
      if (m->bond[i].iat1 == ati) m->bond[i].iat1 = atj;
      else if (m->bond[i].iat1 == atj) m->bond[i].iat1 = ati;
}

void add_atoms(int is_dummy_sym) /*this function does not symmetry operations!*/
{
  int i, j;
  int i0;
  int old_natom = m->natom;
  int atom_in_zero = 0;
  double diam;
  XYZ cm;
  XYZ rm;
  double d_rm = 0.0;
  int nbonded;
  int *ibonded;

  if (!m)
  {
    if (!mol) mol = new_mol(1);
    else m = mol;
  }

  if (is_dummy_sym == 1)
  {
    if (!m->n_selected && m->n_marked)
    {
      double center[3] = {0.0, 0.0, 0.0};
      for(i = 0; i < m->n_marked; i++)
        for(j = 0; j < 3; j++)
          center[j] += m->xyz[m->marked[i]][j];

      for(i = 0; i < 3; i++)
        center[i] /= (double) m->n_marked;

      allocate_atoms(m, m->natom+1);
      for(i = 0; i < 3; i++)
        m->xyz[m->natom-1][i] = center[i];
      m->elem[m->natom-1].name = strdup("Q");

      set_element_data(m, m->natom-1);
      insert_atom_into_list(m->natom-1);
    }
    else if (m->n_selected == 0)
    {
      diam = Calc_Diameter();
      for(i = 0; i < m->natom; i++)
      {
        if (fabs(m->xyz[i][0]) < 0.01 &&
            fabs(m->xyz[i][1]) < 0.01 &&
            fabs(m->xyz[i][2]) < 0.01) atom_in_zero = 1;
      }

      if (atom_in_zero) allocate_atoms(m, m->natom + 6);
      else allocate_atoms(m, m->natom + 7);

      for(i = old_natom; i < m->natom; i++)
      {
        for(j = 0; j < 3; j++)
          m->xyz[i][j] = 0.F;
        m->elem[i].name = strdup("Q");
        set_element_data(m, i);
      }
#ifdef EBUG
      printf("DIAM = %f\n", diam);
      printf("atom in zero = %d\n", atom_in_zero);
#endif
      m->xyz[old_natom  ][0] = diam;
      m->xyz[old_natom+1][0] =-diam;
      m->xyz[old_natom+2][1] = diam;
      m->xyz[old_natom+3][1] =-diam;
      m->xyz[old_natom+4][2] = diam;
      m->xyz[old_natom+5][2] =-diam;

      /*add atoms to atom_list*/
      for(i = old_natom; i < m->natom; i++) insert_atom_into_list(i);

    }
    else if (m->n_selected == 2)
    {
      allocate_atoms(m, m->natom+1);
      for(i = 0; i < 3; i++)
        m->xyz[m->natom-1][i] = 0.5F * (m->xyz[m->selected[0]][i] + m->xyz[m->selected[1]][i]);
      m->elem[m->natom-1].name = strdup("Q");

      set_element_data(m, m->natom-1);
      /*add atom to atom_list*/
      insert_atom_into_list(m->natom-1);
    }
  }
  else if (is_dummy_sym == 0)
  {
    if (m->n_selected == 0)
    {
      if (m->natom == 0)
      {
        allocate_atoms(m, m->natom+1);
        for(i = 0; i < 3; i++)
          m->xyz[0][i] = 0.F;
      }
      else
      {
        /*1. calculate geometrical center*/
        get_center(&cm[0], &cm[1], &cm[2]);
        /*2. find the most distant atom from the center!*/
        i0 = 0;
        for(j = 0; j < 3; j++) rm[j] = 0.F;
        for(i = 0; i < m->natom; i++)
        {
          diam = rr_dist(cm, m->xyz[i]);
          if (diam > d_rm)
          {
            d_rm = diam;
            i0 = i;
          }
        }
        if (d_rm < 0.1)
        {
          for(i = 0; i < 3; i++) rm[i] = 0.5;
          d_rm = 0.866;
        }
        else
          for(i = 0; i < 3; i++)
            rm[i] = m->xyz[i0][i] - cm[i];

        allocate_atoms(m, m->natom+1);

        for(i = 0; i < 3; i++)
          m->xyz[m->natom-1][i] = m->xyz[i0][i] + rm[i] * (m->elem[i0].bond_rad + e[1].bond_rad) / d_rm;

      }

      m->elem[m->natom-1].name = strdup("H");
      set_element_data(m, m->natom-1);
      /*add atom to atom_list*/
      insert_atom_into_list(m->natom-1);
    }
    else if (m->n_selected == 1)
    {
      i0 = m->selected[0];
      nbonded = 0;

      /*1. determine the number of neighboring atoms*/
      for(i = 0; i < m->nbond; i++)
        if (m->bond[i].iat1 == i0 || m->bond[i].iat2 == i0) 
          nbonded++;

      ibonded = (int*) malloc(sizeof(int) * nbonded);

      /*2 determine neighboring atoms*/
      j = 0;
      for(i = 0; i < m->nbond; i++)
      {
        if (m->bond[i].iat1 == i0)
          ibonded[j++] = m->bond[i].iat2;
        if (m->bond[i].iat2 == i0)
          ibonded[j++] = m->bond[i].iat1;
      }

      /*3. find geometrical center of neighboring atoms*/
      if (nbonded)
      {
        for(i = 0; i < 3; i++) cm[i] = 0.0;
        for(i = 0; i < nbonded; i++)
	{
          for(j = 0; j < 3; j++)
            cm[j] += m->xyz[ibonded[i]][j];
	}

        for(i = 0; i < 3; i++)
          cm[i] /= (float) nbonded;

        free(ibonded);

      /*4. determine vector */

        for(i = 0; i < 3; i++)
          rm[i] = cm[i] - m->xyz[i0][i];

      /*5. calculate vector length*/
        d_rm = sqrt(rm[0] * rm[0] + rm[1] * rm[1] + rm[2] * rm[2]);
      }

      /*6. if cm and i0 are at the same place, use the geometry center of entire molecule*/
      if (d_rm < 0.1)
      {
        for(i = 0; i < 3; i++) cm[i] = 0.0;
        for(i = 0; i < m->natom; i++)
          for(j = 0; j < 3; j++)
            cm[j] += m->xyz[i][j];
        for(i = 0; i < 3; i++)
          cm[j] /= (float) m->natom;

        for(i = 0; i < 3; i++)
          rm[i] = cm[i] - m->xyz[i0][i];
        d_rm = sqrt(rm[0] * rm[0] + rm[1] * rm[1] + rm[2] * rm[2]);
      }
      /*7. if cm and i0 are at the same place
       	(there is a larger probability that goat
         can calmly withstend pouring with cold water...)
         use arbitary vector*/
      if (d_rm < 0.1)
      {
        rm[0] = rm[1] = rm[2] = 0.5;
        d_rm = 0.866;
      }

      /*8. put new atom in the direction specified by the vector*/
      allocate_atoms(m, m->natom+1);

      for(i = 0; i < 3; i++)
        m->xyz[m->natom-1][i] = m->xyz[i0][i] - (m->elem[i0].bond_rad + e[1].bond_rad) * rm[i] / d_rm;

      m->elem[m->natom-1].name = strdup("H"); /*could be another atom...*/
      set_element_data(m, m->natom-1);
      /*add atom to atom_list*/
      insert_atom_into_list(m->natom-1);
    }

    /*check if new atom coincides with some of the old atoms*/

/*    nbonded = 0;
    for(i = 0; i < m->natom - 1; i++)
    {
      d_rm = rr_dist(m->xyz[i], m->xyz[m->natom-1]);
      if (d_rm < 0.5)
      {
        i0 = i;
	nbonded = 1;
	break;
      }
    }
    if (nbonded)
    {
    }*/

  }
  append_backup();
  if (Input_Data.automatic_rebonding) rebond();
}

void build_rotation_matrix(XYZ v1, XYZ v2, double A[3][3])
{
  int s;
  int i, j;
  XYZ ax;
  double ang;
  double norm;
  double cosa, sina;

  norm = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  for(i = 0; i < 3; i++) v1[i] /= norm;

  norm = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
  for(i = 0; i < 3; i++) v2[i] /= norm;

  A[0][0] = 1.0; A[0][1] = 0.0; A[0][2] = 0.0;
  A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 0.0;
  A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 1.0;

  if (fabs(v1[0]-v2[0]) + fabs(v1[1]-v2[1]) + fabs(v1[2]-v2[2]) < 1.0e-4) /*same vectors; no rotation*/
  {
/*    printf("fabss = %f\n", fabs(v1[0]-v2[0]) + fabs(v1[1]-v2[1]) + fabs(v1[2]-v2[2]));
    printf("UNIMATRIX!!\n");*/
    return;
  }

  for(i = 0; i < 3; i++)
  {
    s = 0;
    for(j = 0; j < 3; j++)
      if (i == j)
      {
        if (fabs(v1[j] + v2[j]) < 1.0e-4) s++;
      }
      else
      {
	if (fabs(v1[j] - v2[j]) < 1.0e-4) s++;
      }
    if (s == 3)
    {
      A[i][i] = -1.0; /*180 deg. rotation; rot. axis undefined*/
      return;
    }
  }

  ang = acos(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);

  if (!finite(ang))
  {
    /*angle is undefined -> do not rotate; return unit matrix*/
    A[0][0] = 1.0; A[0][1] = 0.0; A[0][2] = 0.0;
    A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 1.0;
    return;
  }

  ax[0] = v1[1]*v2[2] - v1[2]*v2[1];
  ax[1] = v1[2]*v2[0] - v1[0]*v2[2];
  ax[2] = v1[0]*v2[1] - v1[1]*v2[0];

  norm = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
  for(i = 0; i < 3; i++) ax[i] /= norm;
 
  cosa = cos(ang);
  sina = sin(ang);

  A[0][0] = ax[0]*ax[0] + (1.0 - ax[0]*ax[0])*cosa;  A[0][1] = ax[0]*ax[1]*(1.0-cosa) + ax[2]*sina;     A[0][2] = ax[0]*ax[2]*(1.0-cosa) - ax[1]*sina;
  A[1][0] = ax[0]*ax[1]*(1.0-cosa) - ax[2]*sina;     A[1][1] = ax[1]*ax[1] + (1.0 - ax[1]*ax[1])*cosa;  A[1][2] = ax[1]*ax[2]*(1.0-cosa) + ax[0]*sina;
  A[2][0] = ax[0]*ax[2]*(1.0-cosa) + ax[1]*sina;     A[2][1] = ax[1]*ax[2]*(1.0-cosa) - ax[0]*sina;     A[2][2] = ax[2]*ax[2] + (1.0 - ax[2]*ax[2])*cosa;
  return;
}

void change_bond_type(int new_bond_type)
{
  int ibond = find_bond(m, m->selected[0], m->selected[1]);
  if (m->selected[0] < 0 || m->selected[1] < 0) return;
  if (ibond >= 0 && m->bond[ibond].bond_type != new_bond_type)
  {
    m->bond[ibond].bond_type = new_bond_type;
    append_backup();
  }
  else if (ibond == -1 && new_bond_type != 0)
  {
    add_bond(m, m->selected[0], m->selected[1], new_bond_type);
    append_backup();
  }
}

int get_bond_type(int at1, int at2)
{
  int ibond = find_bond(m, at1, at2);
  if (ibond < 0 || ibond >= m->nbond) return NO_BOND;
  return m->bond[ibond].bond_type;
}

void delete_bonds_with_atom(int iatom)
{
  int ibond = 0;

  while(ibond >= 0)
  {
    ibond = find_one_bond_from_atom(iatom);
    if (ibond >= 0) delete_bond(ibond);
    /**/
  }
}

void delete_coord(int type)
{
  int i, j, k;

  if (m->n_selected)
  {
    if (m->selected[0] < 0) return;
    delete_atom(m->selected[0]);

  }
  else if (m->n_marked)
  {
    for(i = 0; i < m->n_marked; i++)
    {
      delete_bonds_with_atom(m->marked[i]);
      if (m->pixdata[m->marked[i]].pix_index.pixels)
      {
        free(m->pixdata[m->marked[i]].pix_index.pixels);
        m->pixdata[m->marked[i]].pix_index.pixels = NULL;
      }
      if (m->pixdata[m->marked[i]].pix_addnum.pixels)
      {
        free(m->pixdata[m->marked[i]].pix_addnum.pixels);
        m->pixdata[m->marked[i]].pix_addnum.pixels = NULL;
      }
      if (m->pixdata[m->marked[i]].pix_symm.pixels)
      {
        free(m->pixdata[m->marked[i]].pix_symm.pixels);
        m->pixdata[m->marked[i]].pix_symm.pixels = NULL;
      }
      if (m->pixdata[m->marked[i]].pix_name.pixels)
      {
        free(m->pixdata[m->marked[i]].pix_name.pixels);
        m->pixdata[m->marked[i]].pix_name.pixels = NULL;
      }
      if (m->pixdata[m->marked[i]].pix_charge_m.pixels)
      {
        free(m->pixdata[m->marked[i]].pix_charge_m.pixels);
        m->pixdata[m->marked[i]].pix_charge_m.pixels = NULL;
      }
      if (m->pixdata[m->marked[i]].pix_charge_l.pixels)
      {
        free(m->pixdata[m->marked[i]].pix_charge_l.pixels);
        m->pixdata[m->marked[i]].pix_charge_l.pixels = NULL;
      }

      for(j = m->marked[i] + 1; j < m->natom; j++)
      {
        for(k = 0; k < 3; k++)
          m->xyz[j-1][k] = m->xyz[j][k];
        m->elem[j-1] = m->elem[j];
        m->charge_m[j-1] =  m->charge_m[j];
        m->charge_l[j-1] =  m->charge_l[j];
        m->additional_numeration[j-1] =  m->additional_numeration[j];
        m->symmetry[j-1] =  m->symmetry[j];
        m->name[j-1] =  m->name[j];
        m->pixdata[j-1] = m->pixdata[j];
      }

      for(j = 0; j < m->nbond; j++)
      {
        if (m->bond[j].iat1 >= m->marked[i]) m->bond[j].iat1--;
        if (m->bond[j].iat2 >= m->marked[i]) m->bond[j].iat2--;
      }

      for(j = i+1; j < m->n_marked; j++)
        if (m->marked[j] > m->marked[i]) m->marked[j]--;

      /*remove atom from the list*/
      remove_atom_from_list(m->marked[i]);
    }
    allocate_atoms(m, m->natom - m->n_marked);
    unmark_all();
  }
  append_backup();
}

void delete_atom(int iatom)
{
  int i, j;

  /*delete bonds near the atom*/
  delete_bonds_with_atom(iatom);

  /*renumerate selected atoms*/
  for(i = 0; i < m->n_selected; i++)
    if (m->selected[i] > iatom) m->selected[i]--;

  /*delete atom name*/
  if (m->name[iatom])
    free(m->name[iatom]);

  /*delete pixdata associated with atom iatom*/
  if (m->pixdata[iatom].pix_index.pixels) free(m->pixdata[iatom].pix_index.pixels);
  if (m->pixdata[iatom].pix_addnum.pixels) free(m->pixdata[iatom].pix_addnum.pixels);
  if (m->pixdata[iatom].pix_symm.pixels) free(m->pixdata[iatom].pix_symm.pixels);
  if (m->pixdata[iatom].pix_name.pixels) free(m->pixdata[iatom].pix_name.pixels);
  if (m->pixdata[iatom].pix_charge_m.pixels) free(m->pixdata[iatom].pix_charge_m.pixels);
  if (m->pixdata[iatom].pix_charge_l.pixels) free(m->pixdata[iatom].pix_charge_l.pixels);

  /*renumerate atoms*/
  for (i = iatom + 1; i < m->natom; i++)
  {
    for(j = 0; j < 3; j++)
      m->xyz[i-1][j] = m->xyz[i][j];
    m->elem[i-1] = m->elem[i];
    m->charge_m[i-1] =  m->charge_m[i];
    m->charge_l[i-1] =  m->charge_l[i];
    m->additional_numeration[i-1] =  m->additional_numeration[i];
    m->symmetry[i-1] =  m->symmetry[i];
    m->name[i-1] =  m->name[i];
    m->pixdata[i-1] = m->pixdata[i];
  }

  /*renumerate bond definitions*/
  for(i = 0; i < m->nbond; i++)
  {
    if (m->bond[i].iat1 >= iatom) m->bond[i].iat1--;
    if (m->bond[i].iat2 >= iatom) m->bond[i].iat2--;
  }

  unselect_all();

  /*reallocate atoms*/
  allocate_atoms(m, m->natom-1);
  /*remove atom from the list*/
  remove_atom_from_list(iatom);
}

void delete_bond(int ibond)
{
  int i;
  if (ibond < 0) return;
  if (ibond > m->nbond) return;
  for(i = ibond + 1; i < m->nbond; i++)
    m->bond[i-1] = m->bond[i];
  allocate_bonds(m, m->nbond - 1);
}

int find_next_dummy(void)
{
  int i;
  if (m->natom <= 0) return -1;
  for(i = 0; i < m->natom; i++)
    if (m->elem[i].name[0] == 'Q' && m->elem[i].name[1] == 0)
      return i;
  return -1; 
}

void delete_dummy_atoms(void)
{
  int next = 1;
  while(next)
  {
    next = find_next_dummy();
    if (next > -1)
    {
      m->selected[0] = next;
      m->n_selected = 1; 
      delete_coord(0);
    }
    next++;
  }
  append_backup();
}

void do_inversion(void) /*MARKED ATOMS !!!!*/
{
  int i, j, k;
  int natom = m->natom;
  int nbond = m->nbond;
  int iat1, iat2;
  int nb;

  if (m->n_selected <= 0) return;

  if (m->n_marked)
  {
    if (is_atom_marked(m->selected[0])) allocate_atoms(m, m->natom + m->n_marked - 1);
    else allocate_atoms(m, m->natom + m->n_marked);

    /*copy atoms*/
    for(i = 0, j = 0; i < m->n_marked; i++)
      if (m->marked[i] != m->selected[0])
      {
        copy_atom_data(m->marked[i], natom + j);
        for(k = 0; k < 3; k++)
          m->xyz[natom+j][k] = 2.0 * m->xyz[m->selected[0]][k] - m->xyz[m->marked[i]][k];
        j++;
      }

    /*copy bonds*/
    for(i = 0; i < nbond; i++)
    {
      nb = 0;
      for(k = 0; k < m->n_marked; k++)
      {
        if (m->marked[k] == m->bond[i].iat1)
	{
          nb++;
          iat1 = k;
	}
	if (m->marked[k] == m->bond[i].iat2)
	{
          nb++;
          iat2 = k;
	}
      }
      if (nb == 2)
        add_bond(m, iat1 + natom,
                    iat2 + natom, m->bond[i].bond_type);
    }
  }
  else
  {
    allocate_atoms(m, 2*m->natom - 1);
    /*copy atoms*/
    for(i = 0, j = 0; i < natom; i++)
      if (i != m->selected[0])
      {
        copy_atom_data(i, natom+j);
        for(k = 0; k < 3; k++)
          m->xyz[natom+j][k] = 2.0 * m->xyz[m->selected[0]][k] - m->xyz[i][k];
        j++;
      }
    /*copy bonds*/

    for(i = 0; i < nbond; i++)
    {
      iat1 = -1;
      iat2 = -1;
      if (m->bond[i].iat1 > m->selected[0]) iat1 = natom + m->bond[i].iat1 - 1;
      else if (m->bond[i].iat1 < m->selected[0]) iat1 = natom + m->bond[i].iat1;
      if (m->bond[i].iat2 > m->selected[0]) iat2 = natom + m->bond[i].iat2 - 1;
      else if (m->bond[i].iat2 < m->selected[0]) iat2 = natom + m->bond[i].iat2;
      if (iat1 >= 0 && iat2 >= 0) add_bond(m, iat1, iat2, m->bond[i].bond_type);
    }
  }
}

void transform_n_fold(int n, int k, XYZ p0, XYZ p1, XYZ px)
{
  int i;
  double v[3]; /*axis definition vector*/
  double theta; /*rotation angle*/
  double q1[3];
  double q2[3];
  double norm;
  
  theta = 2.0 * M_PI * (double) k / (double) n;

  /*define axis vector and translate in origin p -> q*/
  for(i = 0; i < 3; i++)
  {
    v[i] = p1[i] - p0[i];
    q1[i] = px[i] - p0[i];
  }

  norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  for (i = 0; i < 3; i++) v[i] /= norm;

  norm = sqrt(v[1]*v[1] + v[2]*v[2]);

 /*rotate about x axis to bring rotation axis in xz plain*/
  if (norm != 0)
  {
    q2[0] = q1[0];
    q2[1] = q1[1] * v[2] / norm - q1[2] * v[1] / norm;
    q2[2] = q1[1] * v[1] / norm + q1[2] * v[2] / norm;
  }
  else
  {
    for(i = 0; i < 3; i++) q2[i] = q1[i];
  }

  /*rotate about y axis to bring rotation axis along z*/

  q1[0] = q2[0] * norm - q2[2] * v[0];
  q1[1] = q2[1];
  q1[2] = q2[0] * v[0] + q2[2] * norm;

  /*rotate around z axis by theta angle*/
  q2[0] = q1[0] * cos(theta) - q1[1] * sin(theta);
  q2[1] = q1[0] * sin(theta) + q1[1] * cos(theta);
  q2[2] = q1[2];

  /*rotate back about y axis*/
  q1[0] = q2[0] * norm + q2[2] * v[0];
  q1[1] = q2[1];
  q1[2] = -q2[0] * v[0] + q2[2] * norm;

  /*rotate back around x axis*/
  if (norm != 0)
  {
    q2[0] = q1[0];
    q2[1] = q1[1] * v[2] / norm + q1[2] * v[1] / norm;
    q2[2] =-q1[1] * v[1] / norm + q1[2] * v[2] / norm;
  }
  else
  {
    for(i = 0; i < 3; i++) q2[i] = q1[i];
  }

  /*translate back from origin*/
  for(i = 0; i < 3; i++)
    px[i] = q2[i] + p0[i];
}

void copy_atom_data(int iat1, int iat2)
{
  int i;
  if (iat1 >= m->natom) return;
  if (iat2 >= m->natom) return;

  if (m->elem[iat2].name)
    free(m->elem[iat1].name);
  if (m->elem[iat1].name) m->elem[iat2].name = strdup(m->elem[iat1].name);
  for (i = 0; i < 4; i++) m->elem[iat2].color[i] = m->elem[iat1].color[i];
  m->elem[iat2].vdw_rad = m->elem[iat1].vdw_rad;
  m->elem[iat2].bond_rad = m->elem[iat1].bond_rad;
  m->elem[iat2].valency = m->elem[iat1].valency;
  for(i = 0; i < 3; i++)
    m->xyz[iat2][i] = m->xyz[iat1][i];
  if (m->name[iat2]) free(m->name[iat2]);
  if (m->name[iat1]) m->name[iat2] = strdup(m->name[iat1]);
  m->charge_m[iat2] = m->charge_m[iat1];
  m->charge_l[iat2] = m->charge_l[iat1];
  m->additional_numeration[iat2] = m->additional_numeration[iat1];
  m->symmetry[iat2] = m->symmetry[iat1];
}

void add_watched_coord(void)
{
  int i;
  if (!m) return;
  if (m->n_selected < 2) return;
  n_watched++;
  watch = (WATCH*) realloc(watch, sizeof(WATCH) * n_watched);
  if (m->n_selected == 2)
    watch[n_watched-1].coord_type = C_BOND;
  else if (m->n_selected == 3)
    watch[n_watched-1].coord_type = C_ANGLE;
  else if (m->n_selected == 4)
    watch[n_watched-1].coord_type = C_DIHEDRAL;
  for(i = 0; i < watch[n_watched-1].coord_type; i++)
    watch[n_watched-1].atom[i] = m->selected[i] + 1;
  return;
}

void remove_watched_coord(int icoord)
{
  int i, j;
  for(i = icoord; i < n_watched-1; i++)
  {
    watch[i].coord_type = watch[i+1].coord_type;
    for(j = 0; j < watch[i].coord_type; j++)
      watch[i].atom[j] = watch[i+1].atom[j];
  }
  n_watched--;
  watch = (WATCH*) realloc(watch, sizeof(WATCH) * n_watched);
}

void remove_watched_coords(void)
{
  if (n_watched) n_watched--;
  watch = (WATCH*) realloc(watch, sizeof(WATCH) * n_watched);
  if (n_watched == 0) watch = NULL;
}

void remove_all_watched_coords(void)
{
  free(watch);
  n_watched = 0;
  watch = NULL;
}

void add_triangle_selected(void)
{
  int i;
  if (!m) return;
  if(m->n_selected != 3) return;
  allocate_triangles(m, m->ntriangle + 1);

  for(i = 0; i < 3; i++)
  {
    m->triangle1[m->ntriangle-1][i] = m->xyz[m->selected[0]][i];
    m->triangle2[m->ntriangle-1][i] = m->xyz[m->selected[1]][i];
    m->triangle3[m->ntriangle-1][i] = m->xyz[m->selected[2]][i];
  }

  for(i = 0; i < 4; i++)
    m->triangle_color[m->ntriangle-1][i] = Input_Data.extracolor[i];
}

void add_surface_selected(void)
{
  int i;
  if (!m) return;
  if(m->n_selected != 3) return;
  allocate_surfaces(m, m->nsurf + 1);

  for(i = 0; i < 3; i++)
  {
    m->surf1[m->nsurf-1][i] = m->xyz[m->selected[0]][i];
    m->surf2[m->nsurf-1][i] = m->xyz[m->selected[1]][i];
    m->surf3[m->nsurf-1][i] = m->xyz[m->selected[2]][i];
  }

  for(i = 0; i < 4; i++)
    m->surf_color[m->nsurf-1][i] = Input_Data.extracolor[i];
}

void add_cell_selected(void)
{
  int i;
  if (!m) return;
  if(m->n_selected != 4) return;
  allocate_cells(m, m->ncells + 1);

  for(i = 0; i < 3; i++)
  {
    m->cell1[m->ncells-1][i] = m->xyz[m->selected[0]][i];
    m->cell2[m->ncells-1][i] = m->xyz[m->selected[1]][i];
    m->cell3[m->ncells-1][i] = m->xyz[m->selected[2]][i];
    m->cell4[m->ncells-1][i] = m->xyz[m->selected[3]][i];
  }

  for(i = 0; i < 4; i++)
    m->cell_color[m->ncells-1][i] = Input_Data.extracolor[i];
}

void rebond(void)
{
  /*erase all bonds*/
  deallocate_bonds(m);
  determine_bonding(m);
}

