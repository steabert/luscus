/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
/*this structure is used for storing basic orbital data*/
/*these data is stored in tha array and accessed by index*/

const gchar orb_types_1_c[] = 
{
  'U', 'F', 'I', '1', '2', '3', 'S', 'D'
};

typedef struct orb_disp
{
  gint sym_num;   /*symmetry number*/
  gint index_num; /*index number*/
  gdouble energ;
  gdouble occ;
  gint type;

  gint x;    /*x coordinate of the left upper corner on the layout*/ 
  gint y;    /*y coordinate of the left upper corner on the layout*/

  gdouble mx; /*mouse x coordinate*/
  gdouble my; /*mouse y coordinate*/

  GtkWidget *frame; /*frame that defines the area allocated to the one "orbital display"*/
  GtkWidget *draw; /*drawing area for orbitals*/
  GLint olist;

  gdouble size;

  double mvm[16]; /*modelview matrix for the orbital display*/

  GtkWidget *label;
} ORB_DISP;

gint magnification_level;

GtkWidget *multi_view_win; /*multiview window*/ 
gint number_of_orbitals;
