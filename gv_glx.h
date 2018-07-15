/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#ifdef GTK_GLEXT
#define BEGIN_OGL_STUFF \
        GdkGLContext *glcontext = gtk_widget_get_gl_context(da); \
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(da); \
        if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) return FALSE

#define END_OGL_STUFF gdk_gl_drawable_gl_end(gldrawable);
#define SWAP_OGL_BUFFERS gdk_gl_drawable_swap_buffers(gldrawable);
#else
#ifdef LINUX
typedef struct glx_conf
{
  Display *display;
  GLXContext context;
  Colormap xcolormap;
  GdkVisual* visual;
} GLXVIS;

#ifdef GTK2
#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(da); \
        int id = gdk_x11_drawable_get_xid(window); \
        if (!glXMakeCurrent(glxvis.display, id, glxvis.context)) return FALSE

/*        Display *display = gdk_x11_display_get_xdisplay(gdk_window_get_display(window)); \*/
#endif
#ifdef GTK3
#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(da); \
        int id = gdk_x11_window_get_xid(window); \
        if (!glXMakeCurrent(glxvis.display, id, glxvis.context)) return FALSE

/*        Display *display = gdk_x11_display_get_xdisplay(gdk_window_get_display(window)); \*/
#endif
#define END_OGL_STUFF

#define SWAP_OGL_BUFFERS glXSwapBuffers(glxvis.display, id);

#else
typedef struct glx_conf
{
  Display *display;
  int id;
  this_can_not_be_compiled context;
} GLXVIS;
#endif

#endif



