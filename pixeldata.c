/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pango/pangocairo.h>
/*#include<cairo.h>*/
#include"luscus.h"
#include"gv.h"
#include"gv_functions.h"

INPUT_DATA Input_Data;
MOL *m;

void deallocate_stxt_t(STXT_T m_stxt)
{
  if (m_stxt.pixels) free(m_stxt.pixels);
  m_stxt.pixels = NULL;
  m_stxt.width = 0;
  m_stxt.height = 0;
}

void allocate_stxt_t(STXT_T *m_stxt, int width, int height)
{
  m_stxt->width = width;
  m_stxt->height = height;
  m_stxt->pixels = (unsigned char*) malloc(m_stxt->width * m_stxt->height * sizeof(unsigned char));
}

void draw_all_pixdata(void)
{
  if (!m) return;
  if (!m->natom) return;
  if (!m->pixdata) return;

  draw_pixdata_indices();
  draw_pixdata_symbols();

  if (m->ishow & HAS_ATOMNUMS) draw_pixdata_atomnums();
  if (m->ishow & HAS_ATOMNAMES) draw_pixdata_atomnames();
  if (m->ishow & HAS_MULLIKEN) draw_pixdata_mulliken();
  if (m->ishow & HAS_LOPROP) draw_pixdata_loprop();
}

void draw_pixdata_indices(void)
{
  int i, j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;
  char ctmp[10];

  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  desc = pango_font_description_from_string(Input_Data.font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  for(i = 0; i < m->natom; i++)
  {
    sprintf(ctmp, "%d", i + 1);
    deallocate_stxt_t(m->pixdata[i].pix_index);
    init_width = strlen(ctmp) * 2 * approx_glyph_size;

    for(init_height = 1, j = 0; ctmp[j] != 0; j++)
      if (ctmp[j] == 10) init_height++;
    init_height *= (4 * approx_glyph_size);

    init_width+=init_width%32;

    stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
    tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
    memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

    surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
    cr = cairo_create(surface);
    layout = pango_cairo_create_layout(cr);

    pango_layout_set_text(layout, ctmp, -1);
    pango_layout_set_font_description(layout, desc);

    pango_layout_get_size(layout, &width, &height);

    width /= PANGO_SCALE;
    height /= PANGO_SCALE;
    width += width%8;

    allocate_stxt_t(&(m->pixdata[i].pix_index), width, height);

    cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
    pango_cairo_update_layout(cr, layout);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.0);
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_surface_flush(surface);

    for(l = 0, k = m->pixdata[i].pix_index.height-1; l < m->pixdata[i].pix_index.height; l++, k--)
      for(j = 0; j < m->pixdata[i].pix_index.width; j++)
      {
        if (l * m->pixdata[i].pix_index.width + j > m->pixdata[i].pix_index.width * m->pixdata[i].pix_index.height-1) printf("ERROR0 pixdata: index out of bounds!\n");
        if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
        m->pixdata[i].pix_index.pixels[l*m->pixdata[i].pix_index.width + j] = (tmpdat[k*stride+4*j+2]);
      }

    g_object_unref(layout);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(tmpdat);
  }
  pango_font_description_free(desc);
}

void draw_pixdata_symbols(void)
{
  int i, j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;

  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  desc = pango_font_description_from_string(Input_Data.font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  for(i = 0; i < m->natom; i++)
  {
    deallocate_stxt_t(m->pixdata[i].pix_symm);
    init_width = strlen(m->elem[i].name) * 2 * approx_glyph_size;

    for(init_height = 1, j = 0; m->elem[i].name[j] != 0; j++)
      if (m->elem[i].name[j] == 10) init_height++;
    init_height *= (4 * approx_glyph_size);

    init_width+=init_width%32;

    stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
    tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
    memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

    surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
    cr = cairo_create(surface);
    layout = pango_cairo_create_layout(cr);

    pango_layout_set_text(layout, m->elem[i].name, -1);
    pango_layout_set_font_description(layout, desc);

    pango_layout_get_size(layout, &width, &height);

    width /= PANGO_SCALE;
    height /= PANGO_SCALE;
    width += width%8;

    allocate_stxt_t(&(m->pixdata[i].pix_symm), width, height);

    cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
    pango_cairo_update_layout(cr, layout);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.0);
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_surface_flush(surface);

    for(l = 0, k = m->pixdata[i].pix_symm.height-1; l < m->pixdata[i].pix_symm.height; l++, k--)
      for(j = 0; j < m->pixdata[i].pix_symm.width; j++)
      {
        if (l * m->pixdata[i].pix_symm.width + j > m->pixdata[i].pix_symm.width * m->pixdata[i].pix_symm.height-1) printf("ERROR0 pixdata: index out of bounds!\n");
        if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
        m->pixdata[i].pix_symm.pixels[l*m->pixdata[i].pix_symm.width + j] = (tmpdat[k*stride+4*j+2]);
      }

    g_object_unref(layout);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(tmpdat);
  }
  pango_font_description_free(desc);
}

void draw_pixdata_atomnums(void)
{
  int i, j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;
  char ctmp[10];

  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  if (!m->additional_numeration) return;

  desc = pango_font_description_from_string(Input_Data.font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  for(i = 0; i < m->natom; i++)
  {
    sprintf(ctmp, "%d", m->additional_numeration[i]);
    deallocate_stxt_t(m->pixdata[i].pix_addnum);
    init_width = strlen(m->elem[i].name) * 2 * approx_glyph_size;

    for(init_height = 1, j = 0; ctmp[j] != 0; j++)
      if (ctmp[j] == 10) init_height++;
    init_height *= (4 * approx_glyph_size);

    init_width+=init_width%32;

    stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
    tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
    memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

    surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
    cr = cairo_create(surface);
    layout = pango_cairo_create_layout(cr);

    pango_layout_set_text(layout, ctmp, -1);
    pango_layout_set_font_description(layout, desc);

    pango_layout_get_size(layout, &width, &height);

    width /= PANGO_SCALE;
    height /= PANGO_SCALE;
    width += width%8;

    allocate_stxt_t(&(m->pixdata[i].pix_addnum), width, height);

    cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
    pango_cairo_update_layout(cr, layout);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.00);
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_surface_flush(surface);

    for(l = 0, k = m->pixdata[i].pix_addnum.height-1; l < m->pixdata[i].pix_addnum.height; l++, k--)
      for(j = 0; j < m->pixdata[i].pix_addnum.width; j++)
      {
        if (l * m->pixdata[i].pix_addnum.width + j > m->pixdata[i].pix_addnum.width * m->pixdata[i].pix_addnum.height-1) printf("ERROR0 pixdata: index out of bounds!\n");
        if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
        m->pixdata[i].pix_addnum.pixels[l*m->pixdata[i].pix_addnum.width + j] = (tmpdat[k*stride+4*j+2]);
      }

    g_object_unref(layout);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(tmpdat);
  }
  pango_font_description_free(desc);
}

void draw_pixdata_atomnames(void)
{
  int i, j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;

  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  if (!m->name) return;

  desc = pango_font_description_from_string(Input_Data.font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  for(i = 0; i < m->natom; i++)
  {
    deallocate_stxt_t(m->pixdata[i].pix_name);
    init_width = strlen(m->name[i]) * 2 * approx_glyph_size;

    for(init_height = 1, j = 0; m->name[i][j] != 0; j++)
      if (m->name[i][j] == 10) init_height++;
    init_height *= (4 * approx_glyph_size);

    init_width+=init_width%32;

    stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
    tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
    memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

    surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
    cr = cairo_create(surface);
    layout = pango_cairo_create_layout(cr);

    pango_layout_set_text(layout, m->name[i], -1);
    pango_layout_set_font_description(layout, desc);

    pango_layout_get_size(layout, &width, &height);

    width /= PANGO_SCALE;
    height /= PANGO_SCALE;
    width += width%8;

    allocate_stxt_t(&(m->pixdata[i].pix_name), width, height);

    cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
    pango_cairo_update_layout(cr, layout);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.0);
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_surface_flush(surface);

    for(l = 0, k = m->pixdata[i].pix_name.height-1; l < m->pixdata[i].pix_name.height; l++, k--)
      for(j = 0; j < m->pixdata[i].pix_name.width; j++)
      {
        if (l * m->pixdata[i].pix_name.width + j > m->pixdata[i].pix_name.width * m->pixdata[i].pix_name.height-1) printf("ERROR0 pixdata: index out of bounds!\n");
        if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
        m->pixdata[i].pix_name.pixels[l*m->pixdata[i].pix_name.width + j] = (tmpdat[k*stride+4*j+2]);
      }

    g_object_unref(layout);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(tmpdat);
  }
  pango_font_description_free(desc);
}

void draw_pixdata_mulliken(void)
{
  int i, j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;
  char ctmp[10];

  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  if (!m->charge_m) return;

  desc = pango_font_description_from_string(Input_Data.font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  for(i = 0; i < m->natom; i++)
  {
    sprintf(ctmp, "%5.2f", m->charge_m[i]);
    deallocate_stxt_t(m->pixdata[i].pix_charge_m);
    init_width = strlen(ctmp) * 2 * approx_glyph_size;

    for(init_height = 1, j = 0; ctmp[j] != 0; j++)
      if (ctmp[j] == 10) init_height++;
    init_height *= (4 * approx_glyph_size);

    init_width+=init_width%32;

    stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
    tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
    memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

    surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
    cr = cairo_create(surface);
    layout = pango_cairo_create_layout(cr);

    pango_layout_set_text(layout, ctmp, -1);
    pango_layout_set_font_description(layout, desc);

    pango_layout_get_size(layout, &width, &height);

    width /= PANGO_SCALE;
    height /= PANGO_SCALE;
    width += width%8;

    allocate_stxt_t(&(m->pixdata[i].pix_charge_m), width, height);

    cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
    pango_cairo_update_layout(cr, layout);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.0);
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_surface_flush(surface);

    for(l = 0, k = m->pixdata[i].pix_charge_m.height-1; l < m->pixdata[i].pix_charge_m.height; l++, k--)
      for(j = 0; j < m->pixdata[i].pix_charge_m.width; j++)
      {
        if (l * m->pixdata[i].pix_charge_m.width + j > m->pixdata[i].pix_charge_m.width * m->pixdata[i].pix_charge_m.height-1) printf("ERROR0 pixdata: index out of bounds!\n");
        if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
        m->pixdata[i].pix_charge_m.pixels[l*m->pixdata[i].pix_charge_m.width + j] = (tmpdat[k*stride+4*j+2]);
      }

    g_object_unref(layout);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(tmpdat);
  }
  pango_font_description_free(desc);
}

void draw_pixdata_loprop(void)
{
  int i, j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;
  char ctmp[10];

  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  if (!m->charge_l) return;

  desc = pango_font_description_from_string(Input_Data.font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  for(i = 0; i < m->natom; i++)
  {
    sprintf(ctmp, "%5.2f", m->charge_l[i]);
    deallocate_stxt_t(m->pixdata[i].pix_charge_l);
    init_width = strlen(ctmp) * 2 * approx_glyph_size;

    for(init_height = 1, j = 0; ctmp[j] != 0; j++)
      if (ctmp[j] == 10) init_height++;
    init_height *= (4 * approx_glyph_size);

    init_width+=init_width%32;

    stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
    tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
    memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

    surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
    cr = cairo_create(surface);
    layout = pango_cairo_create_layout(cr);

    pango_layout_set_text(layout, ctmp, -1);
    pango_layout_set_font_description(layout, desc);

    pango_layout_get_size(layout, &width, &height);

    width /= PANGO_SCALE;
    height /= PANGO_SCALE;
    width += width%8;

    allocate_stxt_t(&(m->pixdata[i].pix_charge_l), width, height);

    cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
    pango_cairo_update_layout(cr, layout);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.0);
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_surface_flush(surface);

    for(l = 0, k = m->pixdata[i].pix_charge_l.height-1; l < m->pixdata[i].pix_charge_l.height; l++, k--)
      for(j = 0; j < m->pixdata[i].pix_charge_l.width; j++)
      {
        if (l * m->pixdata[i].pix_charge_l.width + j > m->pixdata[i].pix_charge_l.width * m->pixdata[i].pix_charge_l.height-1) printf("ERROR0 pixdata: index out of bounds!\n");
        if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
        m->pixdata[i].pix_charge_l.pixels[l*m->pixdata[i].pix_charge_l.width + j] = (tmpdat[k*stride+4*j+2]);
      }

    g_object_unref(layout);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    free(tmpdat);
  }
  pango_font_description_free(desc);
}

void draw_pixdata_textbox(int itb)
{
  int j, k, l;
  unsigned char *tmpdat;
  int stride;
  int init_width, init_height, width, height;
  gint approx_glyph_size;
  PangoLayout *layout;
  PangoFontDescription *desc;
  cairo_surface_t *surface;
  cairo_t *cr;
  cairo_status_t cstatus;

  if (!m) return;
  if (itb >= m->ntextboxes) return;
  if (!m->textboxes) return;
  if (itb < 0) return;
  if (!m->textboxes[itb].message) return;

#ifdef EBUG
  printf("drawing pixeldata\n"); fflush(stdout);
#endif

  desc = pango_font_description_from_string( m->textboxes[itb].font);
  approx_glyph_size = pango_font_description_get_size(desc);
  if (!pango_font_description_get_size_is_absolute(desc))
    approx_glyph_size /= PANGO_SCALE;

  if (m->textboxes[itb].pixtext.pixels) free(m->textboxes[itb].pixtext.pixels);
  init_width = strlen(m->textboxes[itb].message) * 2 * approx_glyph_size;

  for(init_height = 1, j = 0; m->textboxes[itb].message[j] != 0; j++)
    if (m->textboxes[itb].message[j] == 10) init_height++;
  init_height *= (4 * approx_glyph_size);

  init_width+=init_width%32;

  stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, init_width);
    
  tmpdat = malloc(stride * init_height*4 * sizeof(unsigned char));
  memset(tmpdat, 0, stride * init_height*4 * sizeof(unsigned char));

  surface = cairo_image_surface_create_for_data(tmpdat, CAIRO_FORMAT_ARGB32, init_width, init_height, stride);
  cr = cairo_create(surface);
  layout = pango_cairo_create_layout(cr);

  pango_layout_set_text(layout, m->textboxes[itb].message, -1);
  pango_layout_set_font_description(layout, desc);

  pango_layout_get_size(layout, &width, &height);

  width /= PANGO_SCALE;
  height /= PANGO_SCALE;
  width += width%8;

  allocate_stxt_t(&(m->textboxes[itb].pixtext), width, height);

  cairo_set_source_rgba(cr, 0.00, 0.00, 0.00, 0.00);
  pango_cairo_update_layout(cr, layout);
  cairo_paint(cr);

  cairo_set_source_rgba(cr, 1.00, 1.00, 1.00, 1.0);
  pango_cairo_update_layout(cr, layout);
  pango_cairo_show_layout(cr, layout);
  cairo_surface_flush(surface);

  for(l = 0, k = m->textboxes[itb].pixtext.height-1; l < m->textboxes[itb].pixtext.height; l++, k--)
  {
    for(j = 0; j < m->textboxes[itb].pixtext.width; j++)
    {
      if (l*m->textboxes[itb].pixtext.width + j > m->textboxes[itb].pixtext.width * m->textboxes[itb].pixtext.height - 1) printf("ERROR0 pixdata: index out of bounds!\n");
      if (k*stride+j > 4*stride * init_height-1) printf("ERROR1 pixdata: index out of bounds!\n");
      m->textboxes[itb].pixtext.pixels[l*m->textboxes[itb].pixtext.width+j] = (tmpdat[k*stride+4*j+2]);
#ifdef EBUG
      printf("%d", (unsigned char) ((float) (m->textboxes[itb].pixtext.pixels[l * m->textboxes[itb].pixtext.width + j] / 26.0)));
#endif
    }
# ifdef EBUG
    putchar(10);
#endif
  }

  g_object_unref(layout);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
  free(tmpdat);
  pango_font_description_free(desc);
}

