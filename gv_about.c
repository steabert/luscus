/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<gtk/gtk.h>
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

void luscus_gtk_callback_about(GtkWidget *widget, gpointer data)
{
  GtkWidget *dialog;
  GtkWidget *label;
  GtkWidget *vbox;
  gchar *ftext;

  ftext = g_markup_printf_escaped("<span size=\"15000\">%s %s</span>\n\ngrid and geometry viewer\n\nWritten by:\n\n<span size=\"10000\">Valera Veryazov\nGoran Kova\304\215evi\304\207</span>",PROGRAM_NAME, VERSION);

  dialog = gtk_dialog_new_with_buttons("About Luscus", GTK_WINDOW(gtk_widget_get_toplevel(widget)),
                                       GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                       GTK_STOCK_CLOSE, GTK_RESPONSE_REJECT,
                                       NULL);
  vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

  label = gtk_label_new(NULL);
  gtk_label_set_markup(GTK_LABEL(label), ftext);
  gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
  gtk_widget_show(label);

  gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);

  gtk_dialog_run(GTK_DIALOG (dialog));
  gtk_widget_destroy (dialog);
}

