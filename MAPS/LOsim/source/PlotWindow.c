/* brightnessPlot.c (derived from basicwin.c, Copyright 1989 O'Reilly and Associates, Inc.

   Read a source brightess array and display on screen
 */

#include "PlotWindow.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>

///#include "./bitmaps/icon_bitmap"

#define BITMAPDEPTH 1
#define TOO_SMALL 0
#define BIG_ENOUGH 1

/* These are used as arguments to nearly every Xlib routine, so it saves 
 * routine arguments to declare them global.  If there were 
 * additional source files, they would be declared extern there. */
Display *display; int screen_num; Window win; GC gc; unsigned int width, height; unsigned long textColor;

// Plot dimensions
#define D_X 600		// size of the plot (assumes size of display is this large)
#define D_Y D_X
#define	X0DATA 0//40	// default origin for data in the window
#define Y0DATA 0//40
#define	X0BOX 0//30 	//                    box
#define Y0BOX 0//30
// Plot colors
#define N_MAIN_COLORS 8 // see draw_graphics for the colors
#define N_COLOR_STEPS 10
	
XColor colorWanted[(N_MAIN_COLORS+1)*N_COLOR_STEPS];

// Plot parameters
double xMax,yMax; // the maximum data values
int x0Box,y0Box,x0Data,y0Data;
int xOrg,yOrg; // the coordinates of the center of the screen within the image
XFontStruct *font_normal,*font_superScript;

// Function prototypes
int getBrightnessData(void);
void Redraw(Window win, GC gc, unsigned int width, unsigned int height, unsigned long textColor, int n, double x[], double y[], int colorIx);
void Load_font(char fontName[], XFontStruct **font_info);
void getGC(Window win, GC *gc);
void draw_text(Window win, GC gc, unsigned int win_width, unsigned int win_height, unsigned long color);
void draw_graphics(Window win, GC gc, unsigned int window_width, unsigned int window_height, int n, double x[], double y[], int colorIx);
void TooSmall(Window win, GC gc);
void ColorSetup(void);

int PlotWindow(int whatToDo, int nP, double x[], double y[], int colorIx)
{
	int winPosX, winPosY; 	/* window position */
	int button;
	unsigned int border_width = 4;	/* four pixels */
	unsigned int display_width, display_height;
#define MAX_WNAMELENGTH 300 // maximum length of window title
#define MAX_INAMELENGTH 80
	char window_name[MAX_WNAMELENGTH]; char *pWindow_name=window_name;
	char icon_name[MAX_INAMELENGTH]; char *pIcon_name=icon_name;
	Pixmap icon_pixmap;
	XSizeHints *size_hints;
	XWMHints *wm_hints;
	XClassHint *class_hints;
	XTextProperty windowName, iconName;
	int count,notYetDrawn;
	XEvent report;
	char *display_name = NULL;
	int window_size = BIG_ENOUGH;	/* or TOO_SMALL to display contents */
 	int l;
	// Stuff to deal with keyboard input
	char buffer[10]; int bufsize=10; KeySym keysym; XComposeStatus compose;

	switch(whatToDo)
	{
		case(PW_ADD_POINTS):
			Redraw(win, gc, width, height, textColor,nP,x,y,colorIx); notYetDrawn=0;
			return 0;
			
		case(PW_INITIALIZE):
		xMax=x[0]; yMax=y[0]; notYetDrawn=1;
		if (!(size_hints = XAllocSizeHints())) {
		fprintf(stderr, "%s: failure allocating memory\n", "PlotWindow");
        return 0;
		}
		if (!(wm_hints = XAllocWMHints())) {
			fprintf(stderr, "%s: failure allocating memory\n","PlotWindow" );
	        return 0;
	    }
	    //INDENT FOLLOWING
	if (!(class_hints = XAllocClassHint())) {
		fprintf(stderr, "%s: failure allocating memory\n", "PlotWindow");
        return 0;
    }

	strcpy(window_name,"u, v coverage"); l=strlen(window_name);/// strncat(window_name,brightnessDataFilename,MAX_WNAMELENGTH-l-1);
	strcpy(icon_name,"u, v coverage"); l=strlen(icon_name);/// strncat(icon_name,brightnessDataFilename,MAX_INAMELENGTH-l-1);
	
	/* connect to X server */
	if ( (display=XOpenDisplay(display_name)) == NULL )
	{
		(void) fprintf( stderr, "%s: cannot connect to X server %s\n", 
				"PlotWindow", XDisplayName(display_name));
		return 0;
	}

	/* get screen size from display structure macro */
	screen_num = DefaultScreen(display);
	display_width = DisplayWidth(display, screen_num);
	display_height = DisplayHeight(display, screen_num);

	/* Note that in a real application, winPosX and winPosY would default to 0
	 * but would be settable from the command line or resource database.  
	 */
	winPosX = winPosY = 0;

	/* size window with enough room for text */
#define TOP_THICKNESS 80  /*a guess at the top thickness*/
	width = D_X;/*display_width-(border_width*2);*/ height = D_Y;/*display_height-border_width-TOP_THICKNESS;*/

	/* create opaque window */
	win = XCreateSimpleWindow(display, RootWindow(display,screen_num), 
			winPosX, winPosY, width, height, border_width, BlackPixel(display,
			screen_num), WhitePixel(display,screen_num));


///	/* Create pixmap of depth 1 (bitmap) for icon */
///	icon_pixmap = XCreateBitmapFromData(display, win, icon_bitmap_bits,
///			icon_bitmap_width, icon_bitmap_height);

	/* Set size hints for window manager.  The window manager may
	 * override these settings.  Note that in a real
	 * application if size or position were set by the user
	 * the flags would be UPosition and USize, and these would
	 * override the window manager's preferences for this window. */

	/* winPosX, winPosY, width, and height hints are now taken from
	 * the actual settings of the window when mapped. Note
	 * that PPosition and PSize must be specified anyway. */

	size_hints->flags = PPosition | PSize | PMinSize;
	size_hints->min_width = 300;
	size_hints->min_height = 200;

	/* These calls store window_name and icon_name into
	 * XTextProperty structures and set their other 
	 * fields properly. */
	if (XStringListToTextProperty(&pWindow_name, 1, &windowName) == 0) {
		(void) fprintf( stderr, "%s: structure allocation for windowName failed.\n", 
				"PlotWindow");
		return -1;
	}
		
	if (XStringListToTextProperty(&pIcon_name, 1, &iconName) == 0) {
		(void) fprintf( stderr, "%s: structure allocation for iconName failed.\n", 
				"PlotWindow");
		return -1;
	}

	wm_hints->initial_state = NormalState;
	wm_hints->input = True;
	wm_hints->icon_pixmap = icon_pixmap;
	wm_hints->flags = StateHint | IconPixmapHint | InputHint;

	class_hints->res_name = "PlotWindow";
	class_hints->res_class = "PlotWindow";

	XSetWMProperties(display, win, &windowName, &iconName, 
			NULL, 0, size_hints, wm_hints,
			class_hints);

	/* Select event types wanted */
	XSelectInput(display, win, ExposureMask | KeyPressMask | 
			ButtonPressMask | StructureNotifyMask);

	Load_font("-misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-1",&font_normal);
	Load_font("-misc-fixed-medium-r-normal--13-100-100-100-c-70-iso8859-1",&font_superScript);
	textColor=WhitePixel(display,screen_num);
	
	/* create GC for text and drawing */
	getGC(win, &gc);

	/* Display window */
	XMapWindow(display, win);
 	
	x0Data=X0DATA; y0Data=Y0DATA; // set data and box origins to defaults
	x0Box=X0BOX; y0Box=Y0BOX;	
	xOrg=yOrg=0; // start with no zoom at image center
 	
	/* get events, use first to display text and graphics */
	while (notYetDrawn)  {
		XNextEvent(display, &report);
		switch  (report.type) {
		case Expose:
			/* unless this is the last contiguous expose,
			 * don't draw the window */
			if (report.xexpose.count != 0)
				break;

			/* if window too small to use */
			if (window_size == TOO_SMALL)
				TooSmall(win, gc);
			else Redraw(win, gc, width, height, textColor,0,x,y,0); notYetDrawn=0; // erases only
			break;
		case ConfigureNotify:
			/* window has been resized, change width and
			 * height to send to draw_text and draw_graphics
			 * in next Expose */
			width = report.xconfigure.width;
			height = report.xconfigure.height;
			if ((width < size_hints->min_width) || 
					(height < size_hints->min_height))
				window_size = TOO_SMALL;
			else
				window_size = BIG_ENOUGH;
			break;
		
		// Mouse click means to zoom and position
		case ButtonPress:
			button=report.xbutton.button;
			if(button!=Button2)
			{
			}
			switch(button)
			{
				case(Button1): // left button zooms out
					break;
				case(Button2): // center button resets origin
					break;
				case(Button3): // right button zooms in
					break;
				default:
			}
			Redraw(win, gc, width, height, textColor,nP,x,y,colorIx);
			break;
		
		// Press Esc to exit; other keys do nothing
		case KeyPress:
			count=XLookupString((XKeyEvent *)(&report),buffer,bufsize,&keysym,&compose);
			if(keysym!=XK_Escape) break;
			XUnloadFont(display, font_normal->fid);
			XUnloadFont(display, font_superScript->fid);
			XFreeGC(display, gc);
			XCloseDisplay(display);
			return 0;
		default:
			/* all events selected by StructureNotifyMask
			 * except ConfigureNotify are thrown away here,
			 * since nothing is done with them */
			break;
		} /* end switch */
	} /* end while */
	
	} /* end whatToDo switch */
	
	return 0;
}

void Redraw(Window win, GC gc, unsigned int width, unsigned int height, unsigned long textColor, int n, double x[], double y[], int colorIx)
{
	if(!n)
	{
		XSetWindowBackground(display,win,BlackPixel(display,screen_num)); // make the window black
		XClearWindow(display,win);
		ColorSetup();
	}
	else
	{
		draw_text(win, gc, width, height, textColor);
		draw_graphics(win, gc, width, height,n,x,y,colorIx);
	}
	return;
}

void getGC(Window win, GC *gc)
{
	unsigned long valuemask = 0; /* ignore XGCvalues and use defaults */
	XGCValues values;
	unsigned int line_width = 6;
	int line_style = LineOnOffDash;
	int cap_style = CapRound;
	int join_style = JoinRound;

	/* Create default Graphics Context */
	*gc = XCreateGC(display, win, valuemask, &values);

	/* specify font */
	XSetFont(display, *gc, font_normal->fid);

///	XSetForeground(display, *gc, BlackPixel(display,screen_num));
///	/* specify white foreground since window background is black */
///	XSetForeground(display, *gc, WhitePixel(display,screen_num));

	/* set line attributes */
	XSetLineAttributes(display, *gc, line_width, line_style, 
			cap_style, join_style);

///	/* set dashes */
///	XSetDashes(display, *gc, dash_offset, dash_list, list_length);
}

void Load_font(char fontName[], XFontStruct **font_info)
{
    unsigned long superScriptY; Bool ok;
	/* Load font and get font information structure. */
	if ((*font_info = XLoadQueryFont(display,fontName)) == NULL)
	{
		(void) fprintf(stderr,"%s: Cannot open %s font\n",
				"PlotWindow",fontName);
		return;
	}
	ok=XGetFontProperty(*font_info,XA_SUPERSCRIPT_Y,&superScriptY);
}

void draw_text(Window win, GC gc, unsigned int win_width, unsigned int win_height, unsigned long color)
{
	char s[20];
	sprintf(s,"Hello; how are you?");
	XSetForeground(display, gc, color);
	XSetFont(display, gc, font_normal->fid);
	XDrawString(display, win, gc, D_X+80, 130,s,strlen(s));
	sprintf(s,"Fine, how are you?");
	XDrawString(display, win, gc, D_X+80, 200,s,strlen(s));
}


/*
   Display some points from the image on the screen

   ---------------------------------------------------------------------
   |                                                                    |
   |                                                                    |
   |                                                                    | <--- size of image is
   |    ------------------                                              |      nx x ny points
   |    |                |                                              |
   |    |                | <--- size of screen is                       |
   |    |                |       D_X x D_Y pixels                       |      To convert points
   |    |                |                                              |      to pixels, take a
   |    |       *        |                                              |      square ptsPerPixel
   |    |    screen      |                                              |      on a side and aver-
   |    |    origin      |                                              |      age them
   |    |  (xOrg,yOrg)   |                                              |
   |    |                |                                              |
   |    ------------------                                              |
   |                                   * image origin                   |
   |                                       (0,0)                        |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   |                                                                    |
   ---------------------------------------------------------------------

*/

void draw_graphics(Window win, GC gc, unsigned int window_width, unsigned int window_height, int nP, double x[], double y[], int iColor)
{
	int i,jx,jy,cIx;
 	
 	for (i=0; i<nP; i++)
	{
		jx=x[i]/xMax*(D_X/2)+(D_X/2); // FAKE
		jy=y[i]/yMax*(D_Y/2)+(D_Y/2); // FAKE
		cIx=iColor;
		XSetForeground(display,gc,colorWanted[cIx].pixel);
		// y-coordinate inverted because in X window y increases from top to bottom
		XDrawPoint(display,win,gc,jx+x0Data,(D_Y-jy)+y0Data);
	}
}
	
void ColorSetup(void)
{
/*
	Enclosing box and coordinate axes are drawn first.
	|---------------------- window drawing area (0,0)
	|  |------------------- enclosing box (x0Box,y0Box)
	|  |  |---------------- origin for data plot (x0Data,y0Data)
	|  |  |
	
	Scale data as follows:
	The brightest point is <10^(kMax+1)=dMax
	A point b*10^k (1.0<=b<10) will have a color index
	cIx=(kMax-k)*N_COLOR_STEPS + rounded_int(b/10*N_COLOR_STEPS)
		
		cIx	COLOR		 R   G   B
		0	white		255,255,255
		10	yellow		255,255,  0
///		20	orange		255,165,  0
		20	red			255,  0,  0
		30  magenta		255,  0,255
		40	blue		  0,  0,255
		50	blue-green	  0,255,255
		60	green		  0,255,  0
		70	background color (i e, not plotted)
	Color indices not in the table above are proportional blends of adjacent colors; e g, 13 is 70%
	cIx 10 and 30% cIx 20.
*/	
	XColor nextColor;
	Colormap default_cmap;
	char *colorName[N_MAIN_COLORS]={"White","Yellow","Red","Magenta","Blue","Cyan","Green","Black"};
	int i,cIx,ok;
	int boxWidth,boxHeight;
	int redStep,blueStep,greenStep;
	// Get the pixel values for the desired colors (at full intensity)
	default_cmap=DefaultColormap(display,screen_num);
	// First, the main colors
	for (i=0; i<N_MAIN_COLORS; i++)
	{
		cIx=i*N_COLOR_STEPS;
		ok=XParseColor(display,default_cmap,colorName[i],&colorWanted[cIx]);
		ok=XAllocColor(display,default_cmap,&colorWanted[cIx]);
	}
	// Now the blended colors
	for (i=0; i<N_MAIN_COLORS; i++)
	{
		if((i+1)<N_MAIN_COLORS) nextColor=colorWanted[(i+1)*N_COLOR_STEPS];
		else nextColor.red=nextColor.blue=nextColor.green=0; // FAKE-SHOULD BE BACKGROUND
		redStep=(nextColor.red-colorWanted[i*N_COLOR_STEPS].red)/N_COLOR_STEPS;
		blueStep=(nextColor.blue-colorWanted[i*N_COLOR_STEPS].blue)/N_COLOR_STEPS;
		greenStep=(nextColor.green-colorWanted[i*N_COLOR_STEPS].green)/N_COLOR_STEPS;
		for (cIx=i*N_COLOR_STEPS+1,ok=1; cIx<(i+1)*N_COLOR_STEPS; cIx++)
		{
			colorWanted[cIx]=colorWanted[cIx-1];
			colorWanted[cIx].red+=redStep;
			colorWanted[cIx].blue+=blueStep;
			colorWanted[cIx].green+=greenStep;
			ok&=XAllocColor(display,default_cmap,&colorWanted[cIx]);
		}
	}
	
	
	XSetForeground(display,gc,WhitePixel(display,screen_num));
	XSetLineAttributes(display,gc,0/*"thin line"*/,LineSolid,CapButt,JoinBevel);
	boxWidth=D_X+(2*(x0Data-x0Box)); boxHeight=D_Y+(2*(y0Data-y0Box));
#if 0
	// Draw the enclosing box
	XDrawRectangle(display, win, gc, x0Box, y0Box, boxWidth, boxHeight);
#endif
	
	// Draw the axes (inexact by 1/2 pixel)
	XDrawLine(display,win,gc, x0Box,y0Data+D_Y/2, x0Box+boxWidth,y0Data+D_Y/2);
	XDrawLine(display,win,gc, x0Data+D_X/2,y0Box, x0Data+D_X/2,y0Box+boxHeight);
		
	return;
}
void TooSmall(Window win, GC gc)
{
	char *string1 = "Too Small";
	int y_offset, x_offset;

	y_offset = font_normal->ascent + 2;
	x_offset = 2;

	/* output text, centered on each line */
	XDrawString(display, win, gc, x_offset, y_offset, string1, 
			strlen(string1));
}
