/*This file is a part of luscus project*/
/*Licensed under the Academic Free License version 3.0*/
#include<gtk/gtk.h>
#ifdef GTK_GLEXT
#include<gtk/gtkgl.h>
#else
#include<gdk/gdkx.h>
#include<GL/glx.h>
#endif
#include<gdk/gdkkeysyms.h>
#include<math.h>
#include<GL/gl.h>
#include"luscus.h"
#include"gv.h"
#include"gv_gtk.h"
#include"gv_functions.h"
#include"gv_gtk_multiview.h"

/*minimal size of the subwindow in pixels*/
#define MINW 250
#define MINH 250

#ifdef GTK_GLEXT

#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(widget); \
        GdkGLContext *glcontext = gtk_widget_get_gl_context(widget); \
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget); \
        if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) return FALSE

#define END_OGL_STUFF gdk_gl_drawable_gl_end(gldrawable);
#define SWAP_OGL_BUFFERS gdk_gl_drawable_swap_buffers(gldrawable);

#else

GdkVisual* gvisual;
Colormap gxcolormap;
Display *gdisplay;

#ifdef GTK2
#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(widget); \
        Display *display = gdk_x11_display_get_xdisplay(gdk_window_get_display(window)); \
        int id = gdk_x11_drawable_get_xid(window); \
        if (!glXMakeCurrent(display, id, context)) return FALSE

#endif
#ifdef GTK3
#define BEGIN_OGL_STUFF \
        GdkWindow *window = gtk_widget_get_window(widget); \
        Display *display = gdk_x11_display_get_xdisplay(gdk_window_get_display(window)); \
        int id = gdk_x11_window_get_xid(window); \
        if (!glXMakeCurrent(display, id, context)) return FALSE
#endif

#define END_OGL_STUFF
#define SWAP_OGL_BUFFERS glXSwapBuffers(display, id);

#endif

static gboolean luscus_init_orb_sub_callback(GtkWidget*, gpointer);
static gboolean luscus_init_sizes_orb_callback(GtkWidget*, GdkEventConfigure*, gpointer);
#ifdef GTK2
static gboolean draw_callback(GtkWidget*, GdkEventExpose*, gpointer);
#endif
#ifdef GTK3
static gboolean draw_callback(GtkWidget*, cairo_t*, gpointer);
#endif
static gboolean luscus_orb_sub_button_press_event(GtkWidget*, GdkEventButton*, gpointer);
static gboolean luscus_orb_sub_button_release_event(GtkWidget*, GdkEventButton*, gpointer);
static gboolean luscus_orb_sub_motion_notify_event(GtkWidget*, GdkEventMotion*, gpointer);

/*void draw_Surface(void);*/
static guint idle_id = 0;
void luscus_gtk_get_cur_color(int, GLfloat*, GLfloat*, GLfloat*);
double cum_ang = 0.;
#ifndef GTK_GLEXT
GLXContext context;
#endif

const gchar orb_typ[] = 
{
  'U', 'F', 'I', '1', '2', '3', 'S', 'D'
};

const char *molcas_icon[] = 
{
  "64 64 938 2",
"  	c None",
". 	c #020500",
"+ 	c #000400",
"@ 	c #010301",
"# 	c #010400",
"$ 	c #030600",
"% 	c #040400",
"& 	c #030200",
"* 	c #050100",
"= 	c #060100",
"- 	c #020101",
"; 	c #020201",
"> 	c #010302",
", 	c #000300",
"' 	c #050200",
") 	c #0C0200",
"! 	c #120100",
"~ 	c #170201",
"{ 	c #180300",
"] 	c #180200",
"^ 	c #140100",
"/ 	c #100100",
"( 	c #0B0000",
"_ 	c #090000",
": 	c #020202",
"< 	c #010200",
"[ 	c #010300",
"} 	c #130100",
"| 	c #1D0200",
"1 	c #240000",
"2 	c #270100",
"3 	c #2A0300",
"4 	c #280200",
"5 	c #250000",
"6 	c #220100",
"7 	c #200200",
"8 	c #1F0200",
"9 	c #070100",
"0 	c #030102",
"a 	c #000301",
"b 	c #150200",
"c 	c #2D0300",
"d 	c #330200",
"e 	c #360300",
"f 	c #380100",
"g 	c #3A0100",
"h 	c #370200",
"i 	c #350300",
"j 	c #330400",
"k 	c #2F0300",
"l 	c #290100",
"m 	c #260100",
"n 	c #210200",
"o 	c #1A0200",
"p 	c #0C0000",
"q 	c #010202",
"r 	c #260300",
"s 	c #2F0200",
"t 	c #3C0200",
"u 	c #420300",
"v 	c #480600",
"w 	c #4D0400",
"x 	c #4C0400",
"y 	c #4B0400",
"z 	c #490500",
"A 	c #480500",
"B 	c #450400",
"C 	c #410200",
"D 	c #3C0100",
"E 	c #380300",
"F 	c #330300",
"G 	c #2C0200",
"H 	c #1A0100",
"I 	c #110100",
"J 	c #040202",
"K 	c #2D0200",
"L 	c #3C0300",
"M 	c #4B0300",
"N 	c #510400",
"O 	c #550500",
"P 	c #570400",
"Q 	c #5A0800",
"R 	c #580500",
"S 	c #580300",
"T 	c #580400",
"U 	c #570600",
"V 	c #540400",
"W 	c #4F0600",
"X 	c #490300",
"Y 	c #450300",
"Z 	c #3F0200",
"` 	c #390200",
" .	c #260200",
"..	c #150100",
"+.	c #560300",
"@.	c #5F0500",
"#.	c #640600",
"$.	c #660700",
"%.	c #670500",
"&.	c #690700",
"*.	c #680600",
"=.	c #670600",
"-.	c #650600",
";.	c #620500",
">.	c #600600",
",.	c #5D0600",
"'.	c #550300",
").	c #500500",
"!.	c #490400",
"~.	c #440400",
"{.	c #3E0300",
"].	c #370300",
"^.	c #290200",
"/.	c #0D0000",
"(.	c #620400",
"_.	c #6B0600",
":.	c #730800",
"<.	c #740800",
"[.	c #790900",
"}.	c #7A0900",
"|.	c #730900",
"1.	c #720800",
"2.	c #6C0600",
"3.	c #680800",
"4.	c #640700",
"5.	c #600400",
"6.	c #5B0400",
"7.	c #4E0400",
"8.	c #470400",
"9.	c #3B0100",
"0.	c #2B0200",
"a.	c #1B0200",
"b.	c #530400",
"c.	c #7F0700",
"d.	c #850800",
"e.	c #870900",
"f.	c #8A0600",
"g.	c #860A00",
"h.	c #860900",
"i.	c #830700",
"j.	c #800700",
"k.	c #7E0800",
"l.	c #7D0900",
"m.	c #760900",
"n.	c #700500",
"o.	c #6D0700",
"p.	c #6A0700",
"q.	c #650700",
"r.	c #5E0400",
"s.	c #5A0400",
"t.	c #550200",
"u.	c #3E0100",
"v.	c #340300",
"w.	c #2A0100",
"x.	c #230100",
"y.	c #040001",
"z.	c #020600",
"A.	c #7B0700",
"B.	c #920B00",
"C.	c #970600",
"D.	c #990800",
"E.	c #970800",
"F.	c #940A00",
"G.	c #930900",
"H.	c #930800",
"I.	c #920800",
"J.	c #8F0A00",
"K.	c #8C0B00",
"L.	c #870B00",
"M.	c #840600",
"N.	c #810700",
"O.	c #7E0900",
"P.	c #770700",
"Q.	c #720700",
"R.	c #6E0600",
"S.	c #690800",
"T.	c #620600",
"U.	c #5D0500",
"V.	c #400200",
"W.	c #2B0100",
"X.	c #240100",
"Y.	c #190200",
"Z.	c #040700",
"`.	c #870800",
" +	c #A00A00",
".+	c #9F0B00",
"++	c #A10A00",
"@+	c #A40B00",
"#+	c #A30B00",
"$+	c #A00800",
"%+	c #9B0A00",
"&+	c #980800",
"*+	c #960800",
"=+	c #8F0900",
"-+	c #8B0800",
";+	c #820800",
">+	c #710700",
",+	c #6B0800",
"'+	c #5E0500",
")+	c #580700",
"!+	c #480400",
"~+	c #080000",
"{+	c #570500",
"]+	c #9E0B00",
"^+	c #A30700",
"/+	c #AA0900",
"(+	c #AD0A00",
"_+	c #AA0A00",
":+	c #A60A00",
"<+	c #A30A00",
"[+	c #A00900",
"}+	c #9C0900",
"|+	c #950800",
"1+	c #910A00",
"2+	c #8C0900",
"3+	c #820900",
"4+	c #7B0900",
"5+	c #740700",
"6+	c #6D0600",
"7+	c #680700",
"8+	c #590400",
"9+	c #500400",
"0+	c #3E0200",
"a+	c #060201",
"b+	c #140500",
"c+	c #210300",
"d+	c #220200",
"e+	c #AC0B00",
"f+	c #B40B00",
"g+	c #B00A00",
"h+	c #AB0900",
"i+	c #A80A00",
"j+	c #A70B00",
"k+	c #940800",
"l+	c #8E0900",
"m+	c #890900",
"n+	c #830800",
"o+	c #750800",
"p+	c #6C0700",
"q+	c #670700",
"r+	c #5F0400",
"s+	c #580600",
"t+	c #310300",
"u+	c #030201",
"v+	c #270300",
"w+	c #100200",
"x+	c #AB0A00",
"y+	c #B90B00",
"z+	c #B50A00",
"A+	c #B10B00",
"B+	c #AF0B00",
"C+	c #AA0700",
"D+	c #A60C00",
"E+	c #A30D00",
"F+	c #9D0A00",
"G+	c #970A00",
"H+	c #910900",
"I+	c #8D0A00",
"J+	c #820600",
"K+	c #670800",
"L+	c #560500",
"M+	c #4B0200",
"N+	c #430400",
"O+	c #380200",
"P+	c #080100",
"Q+	c #010600",
"R+	c #250800",
"S+	c #220300",
"T+	c #A50A00",
"U+	c #B80A00",
"V+	c #B80B00",
"W+	c #B50B00",
"X+	c #B30D00",
"Y+	c #B40A00",
"Z+	c #B10A00",
"`+	c #A90D00",
" @	c #A20A00",
".@	c #870700",
"+@	c #530600",
"@@	c #4A0300",
"#@	c #3F0400",
"$@	c #280100",
"%@	c #1E0200",
"&@	c #0A0000",
"*@	c #030101",
"=@	c #0A0100",
"-@	c #AC0A00",
";@	c #BB0C00",
">@	c #B70C00",
",@	c #BF0C00",
"'@	c #BA0E00",
")@	c #AD0E00",
"!@	c #970900",
"~@	c #940900",
"{@	c #850B00",
"]@	c #7C0800",
"^@	c #730700",
"/@	c #4F0500",
"(@	c #3A0200",
"_@	c #140200",
":@	c #540500",
"<@	c #590500",
"[@	c #320300",
"}@	c #040200",
"|@	c #AC0E00",
"1@	c #BF0D00",
"2@	c #C10D00",
"3@	c #BB0B00",
"4@	c #C00C00",
"5@	c #BB0E00",
"6@	c #B20D00",
"7@	c #AA0C00",
"8@	c #9E0900",
"9@	c #980A00",
"0@	c #920900",
"a@	c #8A0900",
"b@	c #820700",
"c@	c #700700",
"d@	c #4A0400",
"e@	c #340400",
"f@	c #2C0500",
"g@	c #050201",
"h@	c #AB0B00",
"i@	c #C10E00",
"j@	c #CB0C00",
"k@	c #DD1000",
"l@	c #BE0E00",
"m@	c #A90C00",
"n@	c #A10B00",
"o@	c #810900",
"p@	c #780800",
"q@	c #650800",
"r@	c #460400",
"s@	c #1F0100",
"t@	c #0E0100",
"u@	c #850900",
"v@	c #6F0600",
"w@	c #2C0400",
"x@	c #000A02",
"y@	c #000A03",
"z@	c #000803",
"A@	c #AC0C00",
"B@	c #B90D00",
"C@	c #CE0B00",
"D@	c #F10F00",
"E@	c #B90E00",
"F@	c #A20C00",
"G@	c #600500",
"H@	c #550600",
"I@	c #4A0500",
"J@	c #3F0300",
"K@	c #300400",
"L@	c #230000",
"M@	c #6E0800",
"N@	c #880A00",
"O@	c #950700",
"P@	c #8A0B00",
"Q@	c #660500",
"R@	c #470300",
"S@	c #001A07",
"T@	c #001B08",
"U@	c #001A08",
"V@	c #001507",
"W@	c #000E03",
"X@	c #010601",
"Y@	c #020200",
"Z@	c #4A0700",
"`@	c #9D0900",
" #	c #CC0D00",
".#	c #EC0F00",
"+#	c #EC1100",
"@#	c #B10D00",
"##	c #B80F00",
"$#	c #A80D00",
"%#	c #9C0B00",
"&#	c #910800",
"*#	c #8A0A00",
"=#	c #790800",
"-#	c #630600",
";#	c #5A0300",
">#	c #4C0300",
",#	c #430300",
"'#	c #170200",
")#	c #060101",
"!#	c #9B0900",
"~#	c #9D0800",
"{#	c #610600",
"]#	c #470500",
"^#	c #3C0400",
"/#	c #190300",
"(#	c #011A04",
"_#	c #003211",
":#	c #002D0E",
"<#	c #00290D",
"[#	c #00250B",
"}#	c #001E0A",
"|#	c #001205",
"1#	c #000A04",
"2#	c #060200",
"3#	c #230101",
"4#	c #C80E00",
"5#	c #E10B00",
"6#	c #CF0F00",
"7#	c #B90F00",
"8#	c #B40C00",
"9#	c #5D0400",
"0#	c #4F0300",
"a#	c #350100",
"b#	c #290300",
"c#	c #080101",
"d#	c #760B00",
"e#	c #930A00",
"f#	c #A80E00",
"g#	c #7D0700",
"h#	c #5E0700",
"i#	c #5C0300",
"j#	c #2C0300",
"k#	c #002E0B",
"l#	c #003810",
"m#	c #00370F",
"n#	c #00330F",
"o#	c #002E0D",
"p#	c #002B0D",
"q#	c #001F0A",
"r#	c #021102",
"s#	c #110501",
"t#	c #2A0001",
"u#	c #340005",
"v#	c #2D0002",
"w#	c #900900",
"x#	c #BA0B00",
"y#	c #C50D00",
"z#	c #B90C00",
"A#	c #A50D00",
"B#	c #9C0C00",
"C#	c #890A00",
"D#	c #6B0700",
"E#	c #530300",
"F#	c #3A0001",
"G#	c #190100",
"H#	c #080001",
"I#	c #B60C00",
"J#	c #A40C00",
"K#	c #640400",
"L#	c #800A00",
"M#	c #720600",
"N#	c #4B0700",
"O#	c #0B2307",
"P#	c #00390F",
"Q#	c #00390C",
"R#	c #023202",
"S#	c #013508",
"T#	c #00370E",
"U#	c #00320E",
"V#	c #032509",
"W#	c #161401",
"X#	c #270800",
"Y#	c #2C0201",
"Z#	c #1F0001",
"`#	c #0D0100",
" $	c #010500",
".$	c #BC0D00",
"+$	c #BE0D00",
"@$	c #AA0B00",
"#$	c #9D0B00",
"$$	c #950900",
"%$	c #8D0C00",
"&$	c #830500",
"*$	c #770800",
"=$	c #630700",
"-$	c #2B0300",
";$	c #020300",
">$	c #8D0900",
",$	c #B70E00",
"'$	c #900B00",
")$	c #700600",
"!$	c #511701",
"~$	c #023308",
"{$	c #003D08",
"]$	c #093A01",
"^$	c #0B3700",
"/$	c #093700",
"($	c #013805",
"_$	c #003A0E",
":$	c #062E07",
"<$	c #241C00",
"[$	c #1D1401",
"}$	c #0E1103",
"|$	c #020701",
"1$	c #C30E00",
"2$	c #BE0C00",
"3$	c #BC0B00",
"4$	c #AC0800",
"5$	c #980900",
"6$	c #940700",
"7$	c #B70A00",
"8$	c #B40F00",
"9$	c #710600",
"0$	c #9C0A00",
"a$	c #2C2701",
"b$	c #003E07",
"c$	c #054206",
"d$	c #133F00",
"e$	c #134100",
"f$	c #0E3A00",
"g$	c #053D01",
"h$	c #00390B",
"i$	c #062E0A",
"j$	c #032D0B",
"k$	c #002209",
"l$	c #001306",
"m$	c #000602",
"n$	c #B00B00",
"o$	c #A40A00",
"p$	c #920700",
"q$	c #870A00",
"r$	c #710800",
"s$	c #640500",
"t$	c #5A0200",
"u$	c #BB0D00",
"v$	c #9F0A00",
"w$	c #690900",
"x$	c #460300",
"y$	c #781400",
"z$	c #063705",
"A$	c #00470C",
"B$	c #144701",
"C$	c #184900",
"D$	c #194A00",
"E$	c #174600",
"F$	c #074001",
"G$	c #003E02",
"H$	c #003C0B",
"I$	c #003B0E",
"J$	c #003912",
"K$	c #002D0D",
"L$	c #001E09",
"M$	c #000E04",
"N$	c #C00E00",
"O$	c #B70B00",
"P$	c #B00C00",
"Q$	c #A70800",
"R$	c #720900",
"S$	c #2A0200",
"T$	c #B40D00",
"U$	c #8C0A00",
"V$	c #004912",
"W$	c #0A4D08",
"X$	c #194C00",
"Y$	c #1B4F00",
"Z$	c #1C5400",
"`$	c #174900",
" %	c #074301",
".%	c #003F04",
"+%	c #00410C",
"@%	c #004314",
"#%	c #003711",
"$%	c #00280B",
"%%	c #001807",
"&%	c #000D04",
"*%	c #C30F00",
"=%	c #BD0C00",
"-%	c #B30B00",
";%	c #A80B00",
">%	c #710900",
",%	c #560400",
"'%	c #3A0000",
")%	c #270200",
"!%	c #6C0800",
"~%	c #280300",
"{%	c #004716",
"]%	c #005817",
"^%	c #104F03",
"/%	c #1B5000",
"(%	c #1C5600",
"_%	c #1D5F00",
":%	c #1D6000",
"<%	c #1A5400",
"[%	c #114800",
"}%	c #004509",
"|%	c #004E14",
"1%	c #004D18",
"2%	c #003F12",
"3%	c #00320F",
"4%	c #00240A",
"5%	c #010502",
"6%	c #AA0D00",
"7%	c #890700",
"8%	c #7D0800",
"9%	c #0B0100",
"0%	c #A30C00",
"a%	c #004E16",
"b%	c #00631B",
"c%	c #095D0C",
"d%	c #195600",
"e%	c #1B5C00",
"f%	c #0F5F00",
"g%	c #044E00",
"h%	c #005716",
"i%	c #005A19",
"j%	c #00521A",
"k%	c #004414",
"l%	c #003911",
"m%	c #002C0D",
"n%	c #001C08",
"o%	c #000903",
"p%	c #C40E00",
"q%	c #500300",
"r%	c #420200",
"s%	c #300300",
"t%	c #005D1B",
"u%	c #006F22",
"v%	c #007321",
"w%	c #066910",
"x%	c #0E5F00",
"y%	c #065E03",
"z%	c #005B0B",
"A%	c #00641A",
"B%	c #005E1C",
"C%	c #00551B",
"D%	c #004C18",
"E%	c #003F13",
"F%	c #00310F",
"G%	c #002109",
"H%	c #001005",
"I%	c #010501",
"J%	c #AD0D00",
"K%	c #690600",
"L%	c #5C0500",
"M%	c #4C0500",
"N%	c #3B0300",
"O%	c #750900",
"P%	c #006620",
"Q%	c #007122",
"R%	c #007726",
"S%	c #007C26",
"T%	c #017922",
"U%	c #056D16",
"V%	c #145B00",
"W%	c #065B08",
"X%	c #00691B",
"Y%	c #00681E",
"Z%	c #00621F",
"`%	c #005819",
" &	c #004E18",
".&	c #004114",
"+&	c #00350F",
"@&	c #00260B",
"#&	c #001407",
"$&	c #B80D00",
"%&	c #AD0B00",
"&&	c #840900",
"*&	c #660600",
"=&	c #1C0100",
"-&	c #4F0400",
";&	c #006419",
">&	c #007424",
",&	c #007323",
"'&	c #007826",
")&	c #007C28",
"!&	c #017320",
"~&	c #125C01",
"{&	c #0B5304",
"]&	c #006215",
"^&	c #00671B",
"/&	c #005F1E",
"(&	c #004D19",
"_&	c #00360E",
":&	c #001809",
"<&	c #010803",
"[&	c #A90A00",
"}&	c #AF0D00",
"|&	c #8F0B00",
"1&	c #5F0600",
"2&	c #3A0300",
"3&	c #0A0300",
"4&	c #090200",
"5&	c #00671E",
"6&	c #007021",
"7&	c #0A5A0B",
"8&	c #114C00",
"9&	c #015A11",
"0&	c #006422",
"a&	c #005419",
"b&	c #004A17",
"c&	c #003F14",
"d&	c #00340F",
"e&	c #001305",
"f&	c #BB0A00",
"g&	c #AC0900",
"h&	c #960900",
"i&	c #860700",
"j&	c #550400",
"k&	c #420400",
"l&	c #2E0200",
"m&	c #0F0100",
"n&	c #00661A",
"o&	c #006E20",
"p&	c #007023",
"q&	c #006B1E",
"r&	c #006A20",
"s&	c #09560B",
"t&	c #194400",
"u&	c #034E0E",
"v&	c #00581A",
"w&	c #005319",
"x&	c #004B18",
"y&	c #004112",
"z&	c #00350D",
"A&	c #032D0D",
"B&	c #042408",
"C&	c #4E0B00",
"D&	c #480300",
"E&	c #350400",
"F&	c #1B0100",
"G&	c #040800",
"H&	c #005F1C",
"I&	c #00601C",
"J&	c #005F20",
"K&	c #0C4C12",
"L&	c #283000",
"M&	c #093805",
"N&	c #004B19",
"O&	c #004516",
"P&	c #003E12",
"Q&	c #003811",
"R&	c #4D0300",
"S&	c #1E0100",
"T&	c #400400",
"U&	c #1C0400",
"V&	c #0A1701",
"W&	c #170300",
"X&	c #310200",
"Y&	c #510500",
"Z&	c #4A0200",
"`&	c #250100",
" *	c #120400",
".*	c #570200",
"+*	c #270400",
"@*	c #000500",
"#*	c #030300",
"$*	c #040600",
"%*	c #1D0100",
"&*	c #4E0200",
"**	c #590900",
"=*	c #590800",
"-*	c #520500",
";*	c #430500",
">*	c #300200",
",*	c #030400",
"'*	c #030302",
")*	c #020203",
"!*	c #060000",
"~*	c #180100",
"{*	c #080300",
"]*	c #6F0700",
"^*	c #4B0600",
"/*	c #410500",
"(*	c #3C0600",
"_*	c #320400",
":*	c #2E0400",
"<*	c #290400",
"[*	c #2A0400",
"}*	c #2A0201",
"|*	c #2F0100",
"1*	c #320200",
"2*	c #3A0500",
"3*	c #3F0500",
"4*	c #440500",
"5*	c #4C0200",
"6*	c #450500",
"7*	c #410300",
"8*	c #360200",
"9*	c #7F0800",
"0*	c #840800",
"a*	c #830900",
"b*	c #800800",
"c*	c #7C0900",
"d*	c #780A00",
"e*	c #5A0500",
"f*	c #590600",
"g*	c #5C0600",
"h*	c #630800",
"i*	c #620700",
"j*	c #5B0800",
"k*	c #530500",
"l*	c #540600",
"m*	c #160200",
"n*	c #900A00",
"o*	c #8D0800",
"p*	c #8B0900",
"q*	c #8B0B00",
"r*	c #880900",
"s*	c #7F0900",
"t*	c #810A00",
"u*	c #650500",
"v*	c #520600",
"w*	c #980700",
"x*	c #9A0900",
"y*	c #9A0B00",
"z*	c #960B00",
"A*	c #960A00",
"B*	c #8C0800",
"C*	c #980B00",
"D*	c #6E0700",
"E*	c #660900",
"F*	c #5B0500",
"G*	c #340200",
"H*	c #1C0200",
"I*	c #8E0A00",
"J*	c #9C0800",
"K*	c #A00B00",
"L*	c #A20900",
"M*	c #A00C00",
"N*	c #9F0D00",
"O*	c #9A0A00",
"P*	c #880700",
"Q*	c #850A00",
"R*	c #730600",
"S*	c #5B0300",
"T*	c #5A0600",
"U*	c #320100",
"V*	c #9F0C00",
"W*	c #A20B00",
"X*	c #A90800",
"Y*	c #A90700",
"Z*	c #AD0800",
"`*	c #AB0D00",
" =	c #A70C00",
".=	c #A50C00",
"+=	c #A00700",
"@=	c #8F0800",
"#=	c #880B00",
"$=	c #7D0A00",
"%=	c #170100",
"&=	c #040101",
"*=	c #8E0800",
"==	c #A10C00",
"-=	c #AE0E00",
";=	c #B01000",
">=	c #B30E00",
",=	c #BC0C00",
"'=	c #B40E00",
")=	c #AE0C00",
"!=	c #A60900",
"~=	c #9F0600",
"{=	c #7A0800",
"]=	c #5C0400",
"^=	c #410400",
"/=	c #250300",
"(=	c #0F0000",
"_=	c #9E0A00",
":=	c #B20C00",
"<=	c #BA0A00",
"[=	c #BD0E00",
"}=	c #B71000",
"|=	c #B31000",
"1=	c #B50F00",
"2=	c #B60D00",
"3=	c #A10900",
"4=	c #300100",
"5=	c #230200",
"6=	c #960700",
"7=	c #9E0C00",
"8=	c #AB0800",
"9=	c #B60A00",
"0=	c #B41100",
"a=	c #AB0E00",
"b=	c #AD0F00",
"c=	c #B51000",
"d=	c #A40E00",
"e=	c #8C0D00",
"f=	c #960600",
"g=	c #A70A00",
"h=	c #BD0B00",
"i=	c #BA0F00",
"j=	c #B70F00",
"k=	c #BA0C00",
"l=	c #AE0B00",
"m=	c #A50E00",
"n=	c #990900",
"o=	c #8E0B00",
"p=	c #7A0A00",
"q=	c #560600",
"r=	c #2D0100",
"s=	c #070200",
"t=	c #AA0800",
"u=	c #B60B00",
"v=	c #B91000",
"w=	c #B80E00",
"x=	c #AE0900",
"y=	c #A40D00",
"z=	c #9A0700",
"A=	c #3B0200",
"B=	c #130000",
"C=	c #990B00",
"D=	c #B20B00",
"E=	c #B90A00",
"F=	c #C10C00",
"G=	c #BD0D00",
"H=	c #B00900",
"I=	c #7B0A00",
"J=	c #510600",
"K=	c #370100",
"L=	c #B00F00",
"M=	c #AD0C00",
"N=	c #4D0100",
"O=	c #0C0300",
"P=	c #9A0D00",
"Q=	c #AF1000",
"R=	c #860800",
"S=	c #760800",
"T=	c #120200",
"U=	c #A80800",
"V=	c #AF0E00",
"W=	c #A80900",
"X=	c #970700",
"Y=	c #8B0C00",
"Z=	c #5E0600",
"`=	c #340100",
" -	c #B30A00",
".-	c #9A0800",
"+-	c #3D0300",
"@-	c #9B0600",
"#-	c #830600",
"$-	c #610500",
"%-	c #440300",
"&-	c #970C00",
"*-	c #9B0800",
"=-	c #7A0700",
"--	c #910700",
";-	c #900800",
">-	c #800900",
",-	c #550A00",
"'-	c #7A0600",
")-	c #770601",
"!-	c #770701",
"~-	c #6A0600",
"                                                                                                                                ",
"                                                                                                                                ",
"                                                                    . + @ @ # # + $                                             ",
"                                                            % & * * * = * - ; > , + , # #                                       ",
"                                                        ' ) ! ! ~ { ] ^ ! / ( _ * : < , + [ # .                                 ",
"                                                      } | 1 2 3 3 3 4 5 6 7 8 { ! ( 9 0 a + , [ #                               ",
"                                                  b 7 c d e f f f g h i j k l m n o ^ p = q [ + , ,                             ",
"                                                r s t u v w x x y z A B C D E F G 2 6 H I 9 J [ , + [                           ",
"                                              K L M N O P Q R S T U V W x X Y Z ` F G  .7 ..p * [ + + ,                         ",
"                                            k Y +.@.#.$.%.&.*.=.-.;.>.,.T '.).!.~.{.].s ^.6 o /.* [ + # ,                       ",
"                                            ).(._.:.<.[.}.}.:.|.1.2._.3.4.5.6.U 7.8.~.9.F 0.n a./.* [ + , #                     ",
"                                          b.1.c.d.e.f.e.g.h.i.j.k.l.m.n.o.p.q.r.s.t.X B u.v.w.x.a./.y.> , [ z.                  ",
"                                          A.B.C.D.E.F.G.H.I.J.K.L.M.N.O.P.Q.R.S.T.U.U w A V.v.W.X.Y.p & , , Z.                  ",
"                                        `. +.+++@+@+#+++$+%+&+*+H.=+-+`.;+k.P.>+,+q.'+)+).!+C F 0.X.^ ~+@ , ,                   ",
"                                      {+        ]+^+/+(+_+:+<+[+}+&+|+1+2+g.3+4+5+6+7+'+8+9+8.0+v.0.7 I a+, + [                 ",
"                                b+c+d+                e+f+g+h+i+j+<+]+D.k+l+m+n+k.o+p+q+r+s+y ~.9.t+ .| ( u++ [                 ",
"                              v+].s w+                  x+y+z+A+B+C+D+E+F+G+H+I+J+c.o+6+K+r.L+M+N+O+K c+} P+> , Q+              ",
"                          R+L X N+S+                      T+U+V+W+X+Y+Z+`+ @&+G..@h.c.<.6+q.U.+@@@#@F $@%@&@*@, #               ",
"                        k !.s.R #@=@                        -@;@;@U+>@,@'@)@<+!@~@2+{@]@^@,+;.s./@Y (@c 7 _@= [ ,               ",
"                      d :@q+2.<@[@}@                          |@1@2@3@4@5@6@7@8@9@0@a@b@A.c@3.5.L+d@C e@m o ~+[ ,               ",
"                      U Q.]@^@9+f@g@              # [ [ #       h@i@2@j@k@  l@m@n@D.H+a@o@p@p+q@6./@r@E G s@t@*@,               ",
"                    :@<.g.u@v@X w@~+            x@y@z@[ < , #     A@B@C@D@    E@m@F@*+l+e.O.:.S.G@H@I@J@K@L@/ * [               ",
"                    M@N@O@P@Q@R@v.| #         S@T@U@V@W@X@Y@[ +   Z@`@ #.#+#  @###$#%#&#*#N.=#p+-#;#>#,#F m '#)#[               ",
"                  #.c.!#~#n+{#>#]#^#/#      (#_#:#<#[#}#|#1#< 2#3#    j+4#5#6#7#y+8#++E.J.M.l.>+&.9#0#B a#b#Y.c#<               ",
"                  d#e#/+f#g#'+s+h#i#7.j#'#  k#l#m#n#o#p#q#r#s#t#u#v#  w#x#y#l@E@z#V+A#B#H.C#j.<.D#>.E#!+F#w@G#H#[               ",
"                  ;+ +I#J#:.@.K#P.L#o@M#N#O#P#Q#R#S#T#U#V#W#X#Y#Z#`# $  [+,@y#.$+$;@@$#$$$%$&$*$p+=$O d@D -$G#H#;$              ",
"                  >$A@,$'$&.5.)$a@B#h@:+!$~${$]$^$/$($_$:$<$[$}$|$[ ,     |@1$2$3$I#4$F@5$J.d.]@p+#.S d@{.-$G#=                 ",
"                  6$7$8$k.<@,.9$`.0$8#`+a$b$c$d$e$f$g${$h$i$j$k$l$m$, #     >@1$3$V+n$o$D.p$q$]@r$s$t$@@0+-$] }@                ",
"                  !#u$v$w$@@x$        y$z$A$B$C$D$E$F$G$H$I$J$K$L$M$, [  $    N$2@O$P$Q$!#H.g.k.R$s$S @@^#S$'#'                 ",
"                  T+T$U$@@w@            V$W$X$Y$Z$Y$`$ %.%+%@%#%$%%%&%[ #     @#*%=%-%;%8@G.f.k.>%s$,%d@'%)%`#                  ",
"                  :+n$!%~%            {%]%^%/%(%_%:%<%[%}%|%1%2%3%4%l$5%[       i@2@-%6%$+H.7%8%c@-#:@X F x.9%                  ",
"                    0%Y               a%b%c%d%e%    f%g%h%i%j%k%l%m%n%o%[       n$p%>@A@$+k+g.4+o.>.q%r%s%H                     ",
"                    1+                t%u%v%w%x%    y%z%A%B%C%D%E%F%G%H%I%        p%;@J%v$H.q$[.K%L%M%N% .t@                    ",
"                    O%                P%Q%R%S%T%U%V%W%X%Y%Z%`% &.&+&@&#&m$        $&=%%&]+H+&&o+*&L+r@t+=&                      ",
"                    -&                ;&>&,&'&)&!&~&{&]&^&/&`%(&.&_&@&:&<&        [&i@}&F+|&O.p+1&7.2&^.3&                      ",
"                    N+4&              5&6&>&,&,&6&7&8&9&0&i%a&b&c&d&@&e&          `@f&g&h&i&o+*&j&k&l&I                         ",
"                    k&m&Q+              n&o&p&q&r&s&t&u&v&w&x&y&z&A&B&          C&l+e+F+*#=#&.P D&E&F&                          ",
"                    >#c+;$G&                H&I&J&K&L&M&N&O&P&Q&                R&]@g.;+*$7+s.d@` S&                            ",
"                  <@8+T&U&Y@+                         V&                    W&X&Y&=.2._.;.{+Z&2&`& *                            ",
"                  {#T..*~.+*w+Y@+ @*                  #*[ # # # @*@*$*}@`#%*s N+&***=*,%0#r@].l ..                              ",
"                  7+o.p+;.-*;*>*8 ! =@& #*& @ + + [ ,*,*& '*)*!*=@! F&~%s 2&~.!.M+>#M !+,#i ^.~*{*                              ",
"                  6+o+*$*$]*&.L%L+^*/*(*_*:*v+r r r 3 <*[*}*|*1*2*3*u 4*!.d@M 5*>#!.6*7*8*W.7 ( #                               ",
"                  *$9*M.0*a*b*c*d*9$6+S.-#@.'+e*f*g*h*i*j*H@k*)+s+L+P P O l*N >#w !+N+f t+`&m*}@                                ",
"                  ]@i.m+n*n*w#o*p*q*r*&&o@s*l.P.o+k.J.>$t*)$7+p.K+-.u*;.'+6.f*v*5*!+,#h S$d+/ @                                 ",
"                  c.C#H+k+w*x*y*x*D.5$z*A*e#=+B*-+w#o$o$C*n+*$=#<.D*p+p.E*G@F*U M+!+u G*4 H*_ a                                 ",
"                  n+I*G.J*]+K* @<+L*<+ @M*N*`@J*5$$$O*j+!#P*Q*0*c.c*R*6+,+#.S*T*0#!.u U*2 ] !*#                                 ",
"                  a@&#E.V*W*X*C+Y*J%C+h+Z*%&`*6% =.=+=E.~@&#@=#=d.k.$=R*o.S.(.F*9+!+u U*2 %=&=                                  ",
"                  *=G.~#==X*C+e+-=6@;=>=7$,=3$;@'@'=)=!=~=*+G.@=C#M.k.{=>+D#-.]=k*!+^=X&/=(=[                                   ",
"                  a@O@_=F@C+A@:=X+W+<=3$[=E@}=|=1=2=5@@#D+3=|+G.w#N@&$l.o+p+$.'+L+!+k&4=5=t@[                                   ",
"                    6=7=#+8=)=-%9=7$x#2$.$E@0=a=b=c=E@,$6%d=`@$$0@e=d.N.[.6+q+'+.*!+T&l&S&p ,*                                  ",
"                    f=.+g=/+B+f+W+U+x#h=4@i=j=8$1=,$'@k=l=m=M*n=H.o=q$J+p=>+q+'+q=R@{.r==&s=                                    ",
"                      V*g=t=%&f+u=u=<=3@+$+$;@v=w=i=.$9=x=y=M*z=H.J.q$J+}.n.q+'+k*R@A=b#B=,*                                    ",
"                      C=:+C+e+D=z+7$V+E=3@,@F=2$,@G=u=H=8=X*]+z=H.&#q$J+I=1.q+L%J=R@K=v+/.                                      ",
"                        #+C+A@L=W+W+7$y+E=E=3@3$3$z#L=M=t=g=]+z=H.I*q$J+[.p+$.L%N=N+X&%@O=                                      ",
"                        P=!=h+J%A+z+9=u=9=O$7$9=9=T$Q=h@C+#+]+D.H.o=R=j.S=D#(.R !.t ^.T=                                        ",
"                          _=U=h+V=P$8#z+u=u=u=9=-%}&e+/+W===]+X=H+Y=&$l.>+&.Z=-*r@`=8                                           ",
"                            ++X*t=-=L=P$ -Y+Y+A+}&e+C+X*<+V*.-H.I*d.9*}.2.-.T d@+- .                                            ",
"                              n@X*t=t=t=-@M=M=M=t=C+W=n@V*@-G.H+r*J+O.D*p.'+9+u X&                                              ",
"                                B#@+C+X*X*X*X*X*X*0%V*V*.-e#0@a@#-9*:._.>.E#Y s%                                                ",
"                                  C=7= @@+J#F@V*]+#$}+C.e#&#p*J+k.:.p.$-b.%-d                                                   ",
"                                      &-O*%+0$*-&+6$0@n*a@q$b*=-r$&.9#q%,#                                                      ",
"                                            --;-*=*=7%h.>-]@<.p+{#8+,-                                                          ",
"                                                    '-)-!-~-                                                                    ",
"                                                                                                                                "
};

ORB_DISP *od = NULL;
GtkWidget *layout;

void kill_multi_view_win(GtkWidget *widget, gpointer data)
{
  gtk_widget_destroy(multi_view_win);
  printf("KILL MULTI VIEW WIN\n");

#ifndef GTK_GLEXT
  glXDestroyContext(gdisplay, context);
  XFreeColormap(gdisplay, gxcolormap);
  g_object_unref(G_OBJECT(gvisual));
#endif

  g_free(od);
  od = NULL;
}

static gboolean luscus_init_orb_sub_callback(GtkWidget *widget, gpointer data)
{
  GLfloat lpos1[4]={-.5,-.5,1,0};
  GLfloat lpos2[4]={-1,0,0,0};
  GLfloat lwhite[4]={1.,1.,1.,1.};
  GLfloat lgray[4]={0.3,0.3,0.3,1.};
  float tmp[]={1.0, 1.0, 1.0};
  gint iiorb = GPOINTER_TO_INT(data);
  gint i, j;
  gint currentview;
  GLfloat bg0, bg1, bg2, bg3;
  GLdouble *mvm_main;
BEGIN_OGL_STUFF;

  iorb = iiorb;
  bg0=Input_Data.background_color[0];
  bg1=Input_Data.background_color[1];
  bg2=Input_Data.background_color[2];
  bg3=Input_Data.background_color[3];

  glLightfv(GL_LIGHT0,GL_POSITION,lpos1);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,lwhite);
  glLightfv(GL_LIGHT0,GL_SPECULAR,lwhite);

/*  glLightfv(GL_LIGHT1,GL_POSITION,lpos2);
  glLightfv(GL_LIGHT1,GL_DIFFUSE,lwhite);
  glLightfv(GL_LIGHT1,GL_SPECULAR,lwhite);
  glLightfv(GL_LIGHT1,GL_AMBIENT,lgray);*/

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_COLOR_MATERIAL);

  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, tmp);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);

  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  glDisable(GL_CULL_FACE);
  glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

  /*draw orbital and store it in list*/

  od[iiorb].olist = glGenLists(1);
  glNewList(od[iiorb].olist, GL_COMPILE);

/*  j = 0;
  for(i = 0; i < m->ngrids; i++)
    if(currentGridData->FltTitle[i]>0)
    {
      if(j==iorb) {currentview=i; break;}
      else j++;
    }*/

   glMatrixMode(GL_MODELVIEW);
   mvm_main = get_main_mvm();

   glLoadMatrixd(mvm_main);

  if(iiorb < m->ngrids-1)
  {
    do_remove_grid();
    do_load_grid();
/*    Do_LoadGrid(currentview,0);*/

    luscus_gtk_get_cur_color(od[iiorb].type, &bg0, &bg1, &bg2);
/*    if(Input_Data.RainBowBG) getCurIndexColor(&bg0, &bg1, &bg2);*/

    glClearColor(bg0, bg1, bg2, bg3);
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

/*    glDepthMask(GL_FALSE);*/
    glDepthMask(GL_TRUE);

    draw_molecule_for_multiview();
 
    draw_grid();
  }
  else
  {
    glClearColor(Input_Data.background_color[0], Input_Data.background_color[1], Input_Data.background_color[2], Input_Data.background_color[3]);
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glDepthMask(GL_FALSE);
  }

  glDepthMask(GL_TRUE);

  glEndList();

END_OGL_STUFF
  return TRUE;
}

void luscus_gtk_get_cur_color(int itype, GLfloat* bg0, GLfloat* bg1, GLfloat*bg2) /*this function works what getCurIndexColor does, but since there is a bug in writing "ThisTitle", this function works correctly*/
{
   switch(itype)
   {
      case 0: *bg0=0.0f;  *bg1=0.0f;  *bg2=0.0f;     break;
      case 1: *bg0=1.0f;  *bg1=0.5f;  *bg2=0.5f;     break;
      case 2: *bg0=1.0f;  *bg1=0.78f; *bg2=0.5f;     break;  
      case 3: *bg0=1.0f;  *bg1=1.0f;  *bg2=0.5f;     break;
      case 4: *bg0=0.5f;  *bg1=0.78f; *bg2=0.5f;     break;
      case 5: *bg0=0.5f;  *bg1=1.0f;  *bg2=1.0f;     break;
      case 6: *bg0=0.5f;  *bg1=0.5f;  *bg2=1.0f;     break;
      case 7: *bg0=0.78f; *bg1=0.5f;  *bg2=0.78f;     break;
      default: *bg0=0.0f;  *bg1=0.0f;  *bg2=0.0f;     break;
   }
}


gboolean luscus_gtk_multiview_redraw_list(int iorb, int type)
{
  GdkRectangle rect;
  gchar *tmpchar;
  GtkWidget *widget;
  gdouble w, h;
#ifndef GTK_OLD
  GtkAllocation allocation;
#endif

  if (od == NULL) return FALSE;
  widget = od[iorb].draw;

BEGIN_OGL_STUFF;

  od[iorb].type = type;
#ifdef GTK_OLD
  w = od[iorb].draw->allocation.width;
  h = od[iorb].draw->allocation.height;
#else
  gtk_widget_get_allocation(od[iorb].draw, &allocation);

  w = allocation.width;
  h = allocation.height;
#endif

/*  od[iorb].type = luscus_gtk_change_orbital_type(od[iorb].sym_num, od[iorb].index_num);*/

  tmpchar = g_strdup_printf("t: %c", orb_types_1_c[type/*od[iorb].type*/]);
  gtk_label_set_text(GTK_LABEL(od[iorb].label), tmpchar);
  g_free(tmpchar);

  redraw_list(iorb);

END_OGL_STUFF;

  rect.x = 2;
  rect.y = 2;
  rect.width=w;
  rect.height=h;

  gdk_window_invalidate_rect(window, &rect, FALSE);
  return TRUE;
}

void redraw_list(int iiorb)
{
  gint i, j;
  gint currentview;
  GLfloat bg0, bg1, bg2, bg3;

  bg0=Input_Data.background_color[0];
  bg1=Input_Data.background_color[1];
  bg2=Input_Data.background_color[2];
  bg3=Input_Data.background_color[3];

  glDeleteLists(od[iiorb].olist, 1);
  od[iiorb].olist = glGenLists(1);

  glNewList(od[iiorb].olist, GL_COMPILE);
  
  j = 0;
  for(i=0; i< m->ngrids; i++)
    if(m->grid_type[i] == ORBITAL)
    {
      if(j==iiorb) {currentview=i; break;}
      else j++;
    }

  if(currentview< m->ngrids-1)
  {
/*    do_remove_grid();*/
    iorb = iiorb;
    do_load_grid(); /*remove_grid; load_grid*/

   /* if(Input_Data.RainBowBG)*/
    luscus_gtk_get_cur_color(od[iiorb].type, &bg0, &bg1, &bg2);
/*    getCurIndexColor(&bg0, &bg1, &bg2);*/

    glClearColor(bg0, bg1, bg2, bg3);
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

/*    glDepthMask(GL_FALSE);*/

    draw_molecule_for_multiview();
 
    draw_grid();
  }
  else
  {
    glClearColor(Input_Data.background_color[0], Input_Data.background_color[1], Input_Data.background_color[2], Input_Data.background_color[3]);
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

/*    glDepthMask(GL_FALSE);*/
  }

/*  glDepthMask(GL_TRUE);  */

  glEndList();
}

static gboolean luscus_init_sizes_orb_callback(GtkWidget *widget, GdkEventConfigure *event, gpointer data)
{ /*this function replaces Init_Sizes and setSizes i Set_Scale*/
#ifdef GTK_OLD
  GLfloat h = (GLfloat) (widget->allocation.height);
  GLfloat w = (GLfloat) (widget->allocation.width);
#else
  GLfloat w, h;

  GtkAllocation allocation;
#endif

BEGIN_OGL_STUFF;
#ifndef GTK_OLD
  gtk_widget_get_allocation(od[iorb].draw, &allocation);

  w = allocation.width;
  h = allocation.height;
#endif

/*  GLdouble d=Calc_Diameter();

  if(d<5) d=5;*/

  glViewport (0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (w > h) glOrtho(-2.50*w/h, 2.50*w/h, -2.50 * h/w, 2.50 * h/w, -50.0, 3.0 /*GuiState.win_state.front_plane, GuiState.win_state.back_plane*/);
  else glOrtho(-2.50, 2.50, -2.50 * h/w, 2.50 * h/w, -50.0, 3.0 /*GuiState.win_state.front_plane, GuiState.win_state.back_plane*/);

  glMatrixMode(GL_MODELVIEW);

  glLoadIdentity ();

END_OGL_STUFF
  return TRUE;
}

#ifdef GTK2
static gboolean draw_callback(GtkWidget *widget, GdkEventExpose *event, gpointer data)
#endif
#ifdef GTK3
static gboolean draw_callback(GtkWidget *widget, cairo_t *cr, gpointer data)
#endif
{
  gint iorb = GPOINTER_TO_INT(data);
BEGIN_OGL_STUFF;
/*  GLfloat bg0, bg1, bg2, bg3;
  gint i, j;
  gint currentview;*/

/*  bg0=Input_Data.bgcolor[0];
  bg1=Input_Data.bgcolor[1];
  bg2=Input_Data.bgcolor[2];
  bg3=Input_Data.bgcolor[3];

  glDepthMask(GL_FALSE);*/
  glScalef(od[iorb].size, od[iorb].size, od[iorb].size);
  od[iorb].size = 1.0;
  glCallList(od[iorb].olist);

SWAP_OGL_BUFFERS
END_OGL_STUFF

  return TRUE;
}

static gboolean luscus_orb_sub_button_press_event(GtkWidget *draw, GdkEventButton *event, gpointer data)
{
  gint iorb = GPOINTER_TO_INT(data);
  gint new_orb_type;
/*  GdkGLContext *glcontext;
  GdkGLDrawable *gldrawable;
  GdkRectangle rect;
  gchar *tmpchar;*/
#ifdef GTK_OLD
  gdouble w = draw->allocation.width;
  gdouble h = draw->allocation.height;
#else
  GLfloat w, h;

  GtkAllocation allocation;

  gtk_widget_get_allocation(od[iorb].draw, &allocation);

  w = allocation.width;
  h = allocation.height;
#endif

/*luscus_gtk_search_orbital_in_list*/
  luscus_gtk_search_orbital_in_list(od[iorb].sym_num, od[iorb].index_num);
  if(event->button == 1 || event->button == 3)
  {
    od[iorb].mx = event->x * 2. / w - 1.;
    od[iorb].my = (h - event->y) * 2. / h - 1.;
  }
  else if (event->button == 2)
  {

    od[iorb].type = luscus_gtk_change_orbital_type(od[iorb].sym_num, od[iorb].index_num);
/*    luscus_gtk_multiview_redraw_list(iorb);*/

/*    od[iorb].type = luscus_gtk_change_orbital_type(od[iorb].sym_num, od[iorb].index_num);

    tmpchar = g_strdup_printf("t: %c", orb_types_1_c[od[iorb].type]);
    gtk_label_set_text(GTK_LABEL(od[iorb].label), tmpchar);
    g_free(tmpchar);

    glcontext = gtk_widget_get_gl_context(draw);
    gldrawable = gtk_widget_get_gl_drawable(draw);

    if (!gdk_gl_drawable_gl_begin(gldrawable, glcontext)) return;

    redraw_list(iorb);

    gdk_gl_drawable_gl_end(gldrawable);

    rect.x = 2;
    rect.y = 2;
    rect.width=w;
    rect.height=h;

    gdk_window_invalidate_rect(draw->window, &rect, FALSE);*/
  }

  return TRUE;
}

static gboolean luscus_orb_sub_button_release_event(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
/*  gint iorb = GPOINTER_TO_INT(data);
  printf("button released iorb = %d ccords: %f %f\n", iorb, od[iorb].mx, od[iorb].my);*/

  return FALSE;
}

static gboolean luscus_orb_sub_motion_notify_event(GtkWidget *widget, GdkEventMotion *event, gpointer data)
{
  gdouble mx;
  gdouble my;
  gdouble x_ax, y_ax, z_ax, angle;
  gdouble dx, dy;
  gint iorb = GPOINTER_TO_INT(data);
  GdkRectangle rect;

  gint currentview;

#ifdef GTK_OLD
  gdouble w = widget->allocation.width;
  gdouble h = widget->allocation.height;
#else
  gdouble w, h;

  GtkAllocation allocation;
#endif

BEGIN_OGL_STUFF;

  mx = event->x * 2. / w - 1.;
  my = (h - event->y) * 2. / h - 1.;

#ifndef GTK_OLD
  gtk_widget_get_allocation(od[iorb].draw, &allocation);

  w = allocation.width;
  h = allocation.height;
#endif

  dx = mx - od[iorb].mx;
  dy = my - od[iorb].my;

  if (event->state & GDK_BUTTON1_MASK)
  {
    x_ax = (od[iorb].my - my);
    y_ax = (mx - od[iorb].mx);
    z_ax = (od[iorb].mx * my - mx * od[iorb].my);


    angle = sqrt((x_ax*x_ax+y_ax*y_ax+z_ax*z_ax) / (mx*mx+my*my+1) * (od[iorb].mx*od[iorb].mx+od[iorb].my*od[iorb].my+1));
    angle*=(180.0e0/M_PI);
    cum_ang += fabs(angle);
    if (cum_ang > 45.0)
    {
      redraw_list(iorb); /*this automatically resorts surfaces also!*/
      cum_ang = 0.0;
    }
/*    if (angle > 1.) angle = 1.;
    if (angle < -1.) angle = -1.;*/

    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX,od[iorb].mvm);
    glLoadIdentity();
    glRotated(angle,x_ax,y_ax,z_ax);
    glMultMatrixd(od[iorb].mvm);
  }
  else if (event->state & GDK_BUTTON3_MASK)
  {
    if (2.0 * fabs(dy) < fabs(dx))
    {
      if (dx > 0.0) od[iorb].size *= 1.01;
      else od[iorb].size /= 1.01;
    }
    else
    {
      if (dy > 0.0) Input_Data.lev*=1.01;
      else Input_Data.lev*=0.99;
      make_surfaces();
      redraw_list(iorb);
    }
/*    printf("od - size = %f lev = %f\n", od[iorb].size, Input_Data.lev);*/

  }

  rect.x = 2;
  rect.y = 2;
  rect.width=w;
  rect.height=h;

  gdk_window_invalidate_rect(window, &rect, FALSE);

END_OGL_STUFF

  od[iorb].mx = mx;
  od[iorb].my = my;
  return TRUE;
}

gboolean luscus_multi_view_key_press(GtkWidget* window, GdkEventKey* event, gpointer data)
{
  printf("KEY PRESS EVENT!\n");
  if (!m->ngrids) return FALSE;
  switch(event->keyval)
  {
    case GDK_KEY_1:
#ifdef GTK_OLD
    case GDK_KP_1:
#else
    case GDK_KEY_KP_1:
#endif
    case 1:
      change_orbital_type(ORBITAL_TYPE_1);
      break;
    case GDK_KEY_2:
#ifdef GTK_OLD
    case GDK_KP_2:
#else
    case GDK_KEY_KP_2:
#endif
    case 2:
    case 'a':
      change_orbital_type(ORBITAL_TYPE_2);
      break;
    case GDK_KEY_3:
#ifdef GTK_OLD
    case GDK_KP_3:
#else
    case GDK_KEY_KP_3:
#endif
    case 3:
      change_orbital_type(ORBITAL_TYPE_3);
      break;
    case 's':
      change_orbital_type(ORBITAL_TYPE_S);
      break;
    case 'd':
      change_orbital_type(ORBITAL_TYPE_D);
      break;
    case 'f':
      change_orbital_type(ORBITAL_TYPE_F);
      break;
    case 'i':
      change_orbital_type(ORBITAL_TYPE_I);
      break;
  }
  return FALSE;
}

/*gboolean luscus_multi_view_win_resize(GtkWidget *window, GdkEvent *event, gpointer data)
{
  arrange_frames_in_layout();
}*/

void arrange_frames_in_layout(void)
{
  gint i;
  gint frames_in_row;
#ifdef GTK3
  GtkAllocation allocation;
#endif

  for(i = 0; i < number_of_orbitals; i++) gtk_widget_set_size_request(od[i].frame, MINW + magnification_level*10, MINH + magnification_level*10);

#ifdef GTK2
  frames_in_row = layout->allocation.width / (MINW + magnification_level * 10);
#endif
#ifdef GTK3
  gtk_widget_get_allocation(layout, &allocation);
  frames_in_row = allocation.width / (MINW + magnification_level * 10);
#endif
  gtk_layout_set_size(GTK_LAYOUT(layout), frames_in_row * (MINW + magnification_level * 10), (number_of_orbitals/frames_in_row + 1) * (MINH + magnification_level * 10));

  for(i = 0; i < number_of_orbitals; i++)
    gtk_layout_move(GTK_LAYOUT(layout), od[i].frame, (MINW + magnification_level * 10)*(i%frames_in_row), (MINH + magnification_level * 10)*(i/frames_in_row));

}

void magnify_frames(GtkWidget button, gpointer data)
{
  magnification_level++;
  arrange_frames_in_layout();
}

void extenuate_frames(GtkWidget button, gpointer data)
{
  magnification_level--;
  arrange_frames_in_layout();
}

void luscus_gtk_wsub_create(gint iod)
{
  GtkWidget *hbox, *vbox;
  /*GtkWidget *draw;*/
  GtkWidget *label;
#ifdef GTK_GLEXT
  GdkGLConfig *glconfig;
#else
  int attributes[] ={GLX_RGBA, 
                     GLX_ALPHA_SIZE, 1,
                     GLX_RED_SIZE, 1, 
                     GLX_GREEN_SIZE, 1, 
                     GLX_BLUE_SIZE, 1, 
                     GLX_DOUBLEBUFFER, True, 
                     GLX_DEPTH_SIZE, 12, 
                     None};

  int xscreen;
  GdkScreen *screen;
  XVisualInfo *xvisual;
  Window root;
#endif
  gchar *tmpchar;

  od[iod].frame = gtk_frame_new(NULL);
  gtk_widget_set_size_request(od[iod].frame, MINW + magnification_level*10, MINH + magnification_level*10);

#ifdef GTK2
  vbox = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox), FALSE);
#endif
  gtk_container_add(GTK_CONTAINER(od[iod].frame), vbox);

#ifdef GTK_GLEXT
  glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH | GDK_GL_MODE_DOUBLE);
#else



  gdisplay = gdk_x11_get_default_xdisplay();
  xscreen = DefaultScreen(gdisplay);
  screen = gdk_screen_get_default();
  xvisual = glXChooseVisual(gdisplay, xscreen, attributes); /*xscreen => glxvis.id */
  gvisual = gdk_x11_screen_lookup_visual(screen, xvisual->visualid); /*xscreen => glxvis.id */
  root = RootWindow(gdisplay, xscreen);
  gxcolormap = XCreateColormap(gdisplay, root, xvisual->visual, AllocNone);
  gtk_widget_set_visual(multi_view_win, gvisual);
  context = glXCreateContext(gdisplay, xvisual, NULL, TRUE);
  XFree(xvisual);
#endif

  od[iod].draw = gtk_drawing_area_new();
  gtk_widget_add_events(od[iod].draw, GDK_BUTTON1_MOTION_MASK |
                                      GDK_BUTTON2_MOTION_MASK |
                                      GDK_BUTTON3_MOTION_MASK |
                                      GDK_BUTTON_PRESS_MASK |
                                      GDK_BUTTON_RELEASE_MASK |
                                      GDK_VISIBILITY_NOTIFY_MASK);

  gtk_box_pack_start(GTK_BOX(vbox), od[iod].draw, TRUE, TRUE, 0);
#ifdef GTK_GLEXT
  gtk_widget_set_gl_capability(od[iod].draw, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE);
#endif
  g_signal_connect_after(G_OBJECT(od[iod].draw), "realize", G_CALLBACK(luscus_init_orb_sub_callback), GINT_TO_POINTER(iod));
  g_signal_connect(G_OBJECT(od[iod].draw), "configure_event", G_CALLBACK(luscus_init_sizes_orb_callback), NULL);
#ifdef GTK2
  g_signal_connect(G_OBJECT(od[iod].draw), "expose_event", G_CALLBACK(draw_callback), GINT_TO_POINTER(iod));
#endif
#ifdef GTK3
  g_signal_connect(G_OBJECT(od[iod].draw), "draw", G_CALLBACK(draw_callback), GINT_TO_POINTER(iod));
#endif

  g_signal_connect(G_OBJECT(od[iod].draw), "button_press_event", G_CALLBACK(luscus_orb_sub_button_press_event), GINT_TO_POINTER(iod));
  g_signal_connect(G_OBJECT(od[iod].draw), "button_release_event", G_CALLBACK(luscus_orb_sub_button_release_event), GINT_TO_POINTER(iod));
  g_signal_connect(G_OBJECT(od[iod].draw), "motion_notify_event", G_CALLBACK(luscus_orb_sub_motion_notify_event), GINT_TO_POINTER(iod));

/*  g_signal_connect(G_OBJECT(draw), "map_event", G_CALLBACK(luscus_orb_sub_map_event), NULL);
  g_signal_connect(G_OBJECT(draw), "unmap_event", G_CALLBACK(luscus_orb_sub_unmap_event), NULL);*/

  gtk_widget_show(od[iod].draw);

#ifdef GTK2
  hbox = gtk_hbox_new(TRUE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);
#endif
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

  tmpchar = g_strdup_printf("s: %d", od[iod].sym_num);
  label = gtk_label_new(tmpchar);
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);
  g_free(tmpchar);

  tmpchar = g_strdup_printf("i: %d", od[iod].index_num);
  label = gtk_label_new(tmpchar);
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  gtk_widget_show(label);
  g_free(tmpchar);

  tmpchar = g_strdup_printf("t: %c", orb_types_1_c[od[iod].type]);
  od[iod].label = gtk_label_new(tmpchar);
  gtk_box_pack_start(GTK_BOX(hbox), od[iod].label, FALSE, FALSE, 0);
  gtk_widget_show(od[iod].label);
  g_free(tmpchar);

  gtk_widget_show(hbox);

  gtk_widget_show(vbox);

  gtk_widget_show(od[iod].frame);
}

void luscus_gtk_show_multiorb_window(void)
{
  GtkWidget *vbox, *hbox;
  GtkWidget *scrolled_window;
  GtkWidget *button;
  GdkPixbuf *pix;
  gint lw_w, lw_h; /*layout window width and height*/
  gint i, j;
  gint norb = 0;
  gint frames_in_row;
  magnification_level = 0;

  pix = gdk_pixbuf_new_from_xpm_data(molcas_icon);
  if (GTK_IS_WIDGET(multi_view_win))
  {
    gtk_widget_destroy(multi_view_win);
    g_free(od);
    return;
  }

  for(i = 0; i < m->ngrids; i++)
    if (m->grid_type[i] == ORBITAL)
      norb++;

  number_of_orbitals = norb;
  od = (ORB_DISP*) g_malloc(sizeof(ORB_DISP) * norb);

  j = 0;
  for(i = 0; i < m->ngrids; i++)
  {
    if (m->grid_type[i] == ORBITAL)
    {
      od[j].sym_num = m->grid_symmetry[i];
      od[j].index_num = m->grid_index[i];
      od[j].energ = m->grid_energy[i];
      od[j].occ = m->grid_occ[i];
      od[j].type = m->orbital_type[i];
      j++;
    }
  }

  multi_view_win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  g_signal_connect(G_OBJECT(multi_view_win), "delete-event", G_CALLBACK(kill_multi_view_win), NULL);
  g_signal_connect(G_OBJECT(multi_view_win), "key_press_event", G_CALLBACK(luscus_multi_view_key_press), NULL);
/*  g_signal_connect(G_OBJECT(multi_view_win), "configure-event", G_CALLBACK(luscus_multi_view_win_resize), NULL);*/
  gtk_window_set_title(GTK_WINDOW(multi_view_win), "Luscus multiview");
  if (GDK_IS_PIXBUF(pix)) gtk_window_set_icon(GTK_WINDOW(multi_view_win), pix);

#ifdef GTK2
  vbox = gtk_vbox_new(FALSE, 0);
#endif
#ifdef GTK3
  vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(vbox), FALSE);  
#endif
  gtk_container_add(GTK_CONTAINER(multi_view_win), vbox);

  scrolled_window = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  gtk_widget_set_size_request(GTK_WIDGET(scrolled_window), Input_Data.init_screen_size, Input_Data.init_screen_size);
  gtk_box_pack_start(GTK_BOX(vbox), scrolled_window, TRUE, TRUE, 0);
/*  gtk_container_add(GTK_CONTAINER(multi_view_win), scrolled_window);*/

  layout = gtk_layout_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(scrolled_window), layout);

  frames_in_row = Input_Data.init_screen_size / (MINW + magnification_level * 10);
  gtk_layout_set_size(GTK_LAYOUT(layout), frames_in_row * (MINW + magnification_level * 10), (norb/frames_in_row + 1) * (MINH + magnification_level * 10));

  j = 0;
  for(i = 0; i < m->ngrids; i++)
  {
    if (m->grid_type[i] == ORBITAL)
    {
      od[j].size = 1.;
      luscus_gtk_wsub_create(j);
      gtk_layout_put(GTK_LAYOUT(layout), od[j].frame, (MINW + magnification_level * 10)*(j%frames_in_row), (MINH + magnification_level * 10)*(j/frames_in_row));
      j++;
    }
  }

  gtk_widget_show(layout);

  gtk_widget_show(scrolled_window);

#ifdef GTK2
  hbox = gtk_hbox_new(FALSE, 0);
#endif
#ifdef GTK3
  hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_box_set_homogeneous(GTK_BOX(hbox), FALSE);  
#endif
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

  button = gtk_button_new_from_stock(GTK_STOCK_ZOOM_IN);
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(magnify_frames), NULL);
  gtk_widget_show(button);

  button = gtk_button_new_from_stock(GTK_STOCK_ZOOM_OUT);
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(extenuate_frames), NULL);
  gtk_widget_show(button);

  gtk_widget_show(hbox);
  gtk_widget_show(vbox);
  gtk_widget_show(multi_view_win);
}


