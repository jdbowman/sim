/**********************************************************************
 *                                                                    *
 * brigen.c                                                           *
 *                                                                    *
 * Creation of a simulated sky brightness image and saving it in fits *
 * format.                                                            *
 *                                                                    *
 * Created Jan 10th 2011 by L. Benkevitch                             *
 *                                                                    * 
 *********************************************************************/
#define FALSE 0
#define TRUE  1
#define INPUT_LINE_LENGTH 256
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <paths.h>
#include <errno.h>
#include <fitsio.h>

#include "brigen.h"
#include "visgen.h"
#include "sky_source.h"
#include "read_prj.h"

double const pi = 3.1415926535897931;
int yydebug = 0;

static void brigen_help(void);
/*void read_brigen_data(void);*/

struct bg_param *bgp = NULL;
struct vg_param *vgp = NULL;
int isrc = 0; /* Index into the array of sky sources for the parser */
int iimgf = 0; /* Index into the list of fits files with external images */

extern struct bufstack *curbs; /* Linked list of the saved inpit files */


int main(int argc, char *argv[]) 
{
  char optstring[] = "m:n:x:y:a:b:s:i:dh"; /* All 1-char options */
  int optc;                          /* Command line option character */
  struct bg_param bgparam =  /* Input parameters of brigen */
    {0, 0, 0.0, 0.0, 0.0, 0.0, NULL, 0, 0, NULL, NULL, NULL, 
     {0, 0, 0, 0, 0, 0, 0, 0}};
  struct Source source[MAX_SOURCES];
  int naxis_col = 0, naxis_row = 0;
  double xsize = 0., ysize = 0.;
  double xcenter = 0., ycenter = 0.;
  int i, nfilearg, iparse;
  int ncol, nrow; /* Image grid dimensions */
  char deb = FALSE; 
  char *inpfilename = NULL, *fitsfilename = NULL;
  char *inpfbase = NULL;
  double *image; /* Array to contain the brightness image */
  /* int j; */
  /* FILE *fh = NULL; */
  
  /*
   * Initialization 
   */
  bgparam.source = source;
  bgp = &bgparam; /* Make bgparam visible for lexer and parser */
  //vgp = &vgparam; /* Make vgparam visible for lexer and parser */
  

  if (argc == 1) {
    brigen_help(); /* Command line empty: print out help page */
    exit(0);
  }

  /* 
   * Parse command line arguments 
   */
  while (1) {
    static struct option long_options[] =
      {
	/* These options set a flag:
	 * {"verbose", no_argument,       &verbose_flag, 1},
	 * {"brief",   no_argument,       &verbose_flag, 0},
	 */
	/* These options don't set a flag.
	   We distinguish them by their indices. */
	{"mgrid",          required_argument, 0, 'm'},
	{"ngrid",          required_argument, 0, 'n'},
	{"racenter",       required_argument, 0, 'x'},
	{"deccenter",      required_argument, 0, 'y'},
	{"rasize",         required_argument, 0, 'a'},
	{"decsize",        required_argument, 0, 'b'},
	{"srcfile",        required_argument, 0, 's'},
	{"addimage",       required_argument, 0, 'i'},
	{"debug",          no_argument,       0, 'd'},
	{"help",           no_argument,       0, 'h'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int err = 0; /* Error flag for the angle conversion routines */
    char *tail = NULL; /* String tail for the angle conversion routines */
     
    optc = getopt_long (argc, argv, optstring,	long_options, &option_index);
     
    /* Detect the end of the options. */
    if (optc == -1)
      break;

    switch (optc)
      {
      case 'm': /* Image grid column number in pixels */
	if (bgparam.is_set.ncol) {
	  printf("Extra '--mgrid' or '-m' cmdline option\n");
	  exit(1); 
	}
	naxis_col = strtol(optarg, &tail, 10); 
	if (naxis_col <= 0 || *tail != 0) {
	  printf("Wrong '--mgrid' or '-m' cmdline option format\n");
	  exit(1);
	}
	bgparam.ncol = naxis_col;
	bgparam.is_set.ncol = 1;
	break;

      case 'n': /* Image grid row number in pixels */
	if (bgparam.is_set.nrow) {
	  printf("Extra '--ngrid' or '-n' cmdline option\n");
	  exit(1); 
	}
	naxis_row = strtol(optarg, &tail, 10);
	if (naxis_row <= 0 || *tail != 0) {
	  printf("Wrong '--ngrid' or '-n' cmdline option format\n");
	  exit(1);
	}
	bgparam.nrow = naxis_row; 
	bgparam.is_set.nrow = 1;
	break;

      case 'a': /* Image RA size in angular units */
	if (bgparam.is_set.xsize) {
	  printf("Extra '--rasize' or '-a' cmdline option\n");
	  exit(1); 
	}
	printf("-a = '%s'\n", optarg);
	xsize = angle2deg(optarg, &err)*3600.0;
	if (err) {
	  xsize =  printf("Wrong '--rasize' or '-a' cmdline option format\n");
	  exit(1);
	}
	bgparam.xsize = xsize;
	bgparam.is_set.xsize = 1;
	break;

      case 'b': /* Image Dec size in angular units */
	if (bgparam.is_set.ysize) {
	  printf("Extra '--decsize' or '-b' cmdline option\n");
	  exit(1); 
	}
	printf("-b = '%s'\n", optarg);
	ysize = angle2deg(optarg, &err)*3600.0;
	if (err) {
	  printf("Wrong '--decsize' or '-b' cmdline option format\n");
	  exit(1);
	}
	bgparam.ysize = ysize;
	bgparam.is_set.ysize = 1;
	break;

      case 'x': /* Image RA center in angular units */
	if (bgparam.is_set.xcenter) {
	  printf("Extra '--racenter' or '-x' cmdline option\n");
	  exit(1); 
	}
	printf("-x = '%s'\n", optarg);
	xcenter = angle2deg(optarg, &err);
	printf("xcenter = %g\n", xcenter);
	if (err) {
	  printf("Wrong '--racenter' or '-x' cmdline option format\n");
	  exit(1);
	}
	bgparam.xcenter = xcenter;
	bgparam.is_set.xcenter = 1;
	break;

      case 'y': /* Image Dec center in angular units */
	if (bgparam.is_set.ycenter) {
	  printf("Extra '--deccenter' or '-y' cmdline option\n");
	  exit(1); 
	}
	printf("-y = '%s'\n", optarg);
	ycenter = angle2deg(optarg, &err);
	printf("ycenter = %g\n", ycenter);
	if (err) {
	  printf("Wrong '--deccenter' or '-y' cmdline option format\n");
	  exit(1);
	}
	bgparam.ycenter = ycenter;
	bgparam.is_set.ycenter = 1;
	break;

      case 's': /* Source list file name: overrides one given in project file */
	if (bgparam.srcfilename) {
	  printf("Extra '--srclist' or '-s' cmdline option\n");
	  exit(1); 
	}
	printf("-s = '%s'\n", optarg);
	bgparam.srcfilename = strdup(optarg);
	break;

      case 'i': /* Additional image fits file name (may be several) */
	if (iimgf < MAXEXTIMAGES)
	  bgparam.imgfilename[iimgf++] = strdup(optarg);
	else { /* maximum number of additioanl image files exceeded */
	  printf("Number of additional image files exceeds %d\n", iimgf);
	  printf("(command line option -i or --addimage)\n");
	  exit(1);
	}
	break;

      case 'd': deb = TRUE;
	break;

      case 'h': brigen_help(); exit(0);
	break;

      default:  ;
      }
  }
  
  /*
   * Figure out how many non-option arguments are there on the command line. 
   * Only 1 or 2 file mames are allowed. If the name is only one,
   * <base>.<ext>,  its base is used to fotm the ourput fits file name 
   * <base>.fits 
   */
  nfilearg = argc - optind; /* Number of non-option, filename arguments */
  if (nfilearg > 2) {
    printf("More than two file names on command line.\n"); 
    exit(1);
  }
  if (nfilearg == 0) {
    printf("No file name(s) on command line.\n"); 
    exit(1);
  }

  inpfilename = (char *) strdupa(argv[optind]);
  /* 
   * Assume the file on cmdline is simulation project file.
   * What it really is will be clear arter parsing it.
   */
  //bgparam.srcfilename = (char *) strdupa(argv[optind]);
  //bgparam.is_set.srclist = 1;
  
  /*
   * If on cmdline only one file name, its base is used to form
   * the name of output fits file. Otherwise the second file name
   * is assumed the output fits file name.
   */
  if (nfilearg == 1) { /* Only one file in cmdline */
    //ifile_set = TRUE;
    long i, j, slen;
    /* Name with extension, but no dirs */
    inpfbase = (char *) strdupa((char *) basename(inpfilename)); 
    slen = strlen(inpfbase);
    /* Find the dot separating extension from the file name */
    j = slen-1;
    for (i = 0; i < slen; i++) { if (inpfbase[j] == '.') break; j--; }
    if (inpfbase[j] == '.') { /* Extension exists */
      inpfbase[j] = 0; /* Cut the base name at the dot */
    }
    fitsfilename = (char *) malloc(sizeof(char)*(strlen(inpfbase)+5));
    strcpy(fitsfilename, inpfbase);
    strcat(fitsfilename, ".fits");
    bgparam.fitsfilename = (char *) strdupa(fitsfilename);
    /* We intentionally leave bgparam.is_set.fitsfilename = 0,
     * so while parsing the project file the fits filename in 
     * the project file, if present, would supersede this 
     * "assumed" filename. */
    bgparam.is_set.fitsfilename = 0;
  /* printf("Output FITS file will have the same name as "		\
     *	   "input file, '%s', with '.fits' extension.\n", inpfbase); */
  }
  else { /* Second file in cmdline is treated as output fits file name */
    fitsfilename = (char *) strdupa(argv[optind+1]);
    bgparam.fitsfilename = fitsfilename;   
    bgparam.is_set.fitsfilename = 1; /* fitsfile is in cmdline explicitly */

    printf("argv[optind] = '%s', argv[optind+1] = '%s'\n",
	   argv[optind], argv[optind+1]);

  }

  /* If the source file name is given in the cmdline with -s or -srclist 
   * option, it should be scanned first to override the sourcefile  given
   * in the simulation project file 
   */ 
  if (bgparam.srcfilename) {
    newfile(FSrclist, bgparam.srcfilename);
    iparse = yyparse();
    /* yy_flush_buffer(curbs); // discard buffer contents */
    if (iparse) exit(1); //==== Quit in case of a syntax error ============>>>
    bgparam.is_set.srclist = 1; /* Do not read the srclist any more */
  }


  /*
   * Open the first file from command line and put it on stack
   * To read the external files, like srclist, obsspec, or array)
   * the read_prj.l lexer makes more calls to newfile(), putting
   * the file info on stack, organized as linked list 
   * struct bufstack *curbs.  
   */
  newfile(FProject, inpfilename);
  
  /*
   * Parse the first file given on the command line. The parser will 
   * automatically determine what kind of input it is: the simulation 
   * project description file or just the list of celestial sources. 
   * In any case, all the data will be read and inserted into 
   * struct bg_param bgparam and struct Source source[]. 
   */
  iparse = yyparse();

  if (iparse) exit(1); //==== Quit in case of a syntax error ==============>>>

  printf("Input file: %s\n", inpfilename);
  if (bgparam.srcfilename) 
    printf("Source list file: %s\n", bgparam.srcfilename);
  printf("Output file: %s\n", bgparam.fitsfilename);

  if (bgparam.is_set.fitsfilename == 0) 
    bgparam.is_set.fitsfilename = 1; /* No imfitsout in the project */

  bgparam.nsrc = isrc;   /* Save calculated number of sources */
  bgparam.nimgf = iimgf; /* Save calculated number of external images */

  /*
   * Check if all the data have been entered correctly
   */
  printf("ncol, nrow = %d, %d\n", bgparam.ncol,  bgparam.nrow);
  printf("xsize, ysize = %g\", %g\"\n", bgparam.xsize,  bgparam.ysize);
  printf("xcenter, ycenter = %g deg, %g deg\n\n", 
	 bgparam.xcenter,  bgparam.ycenter);
  printf("List of the celestial sources:\n");
  printf("Type          Intens.   Q       U       V    SpecIx     " \
	 "x       y      Maj.ax  Min.ax  PA(deg)\n");
  for (i = 0; i < isrc; i++) {
    struct Source *s = source;
    printf("%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
	   s[i].id, s[i].fluxDensity, s[i].Q, s[i].U, s[i].V, 
	   s[i].specIndex, s[i].xOffset, s[i].yOffset, 
	   s[i].majorAxis, s[i].minorAxis, s[i].positionAngle);
  }
  printf("\n");

  /*
   * Check if all necessary data are entered
   */
  int ierr = 0;
  if (!bgparam.is_set.ncol) { printf("Error: need grid ncol size.\n"); ierr++;}
  if (!bgparam.is_set.nrow) { printf("Error: need grid nrow size.\n"); ierr++;}
  if (!bgparam.is_set.xsize) { printf("Error: need FOV RA size.\n"); ierr++; }
  if (!bgparam.is_set.ysize) { printf("Error: need FOV DEC size.\n"); ierr++; }
  if (!bgparam.is_set.xcenter) 
    { printf("Error: need RA of the FOV center.\n"); ierr++; }
  if (!bgparam.is_set.ycenter) 
    { printf("Error: need DEC of the FOV center.\n"); ierr++; }
  if (!bgparam.is_set.srclist) { printf("Error: need source list.\n"); ierr++;} 
  if (!bgparam.is_set.fitsfilename) 
    { printf("Error: need output FITs file name.\n"); ierr++;} 
  if (ierr) exit(1); 

  /*
   * Allocate the ncol x nrow  array for the output image
   */
  ncol = bgparam.ncol;
  nrow = bgparam.nrow;
  image = (double *) malloc(ncol*nrow*sizeof(double));

  /*
   * Sum up all the additional external images (if any) into one image
   */
  if (bgparam.nimgf)
      add_extimages(&bgparam, image);

  /*
   * Combine all the sources in one brightness image
   */
  sources_to_image(&bgparam, image);
  
  /* fh = fopen("imout.txt", "w"); */
  /* for (i = 0; i < nrow; i++) { */
  /*   for (j = 0; j < ncol; j++) { */
  /*     fprintf(fh, "%g ", image[i*ncol+j]); */
  /*   } */
  /*   fprintf(fh, "\n"); */
  /* } */
  /* fclose(fh); */

  /*
   * Save the brightness image as a FITs file
   */
  image_to_fits(&bgparam, image);

  /*
   * Once again inform the user
   */
  printf("\n");
  printf("Success. Output file: %s\n", bgparam.fitsfilename);

  return 0;
}


static void brigen_help(void) { /* 'static"- visible within this file only */
  printf("=== MAPS Sky Brightess Image Simulator ===\n"			\
	 "Creates simulated brightness map with sources, whose "	\
	 "coordinates are specified in infile.\n\n" \
	 "Usage: brigen [OPTIONS] infile [outfile]\n\n"	\
	 "Example: brigen skysources.txt skyimage.fits\n\n" \
	 "OPTIONS:\n"	\
	 "-m <M>,\t  --mgrid <M>\t  image RA  grid size in pixels\n" \
	 "-n <N>,\t  --ngrid <N>\t  image Dec grid size in pixels\n" \
	 "-a <A>,\t --rasize  <A> \t image RA  size in angular units\n" \
	 "-b <B>,\t --decsize <B> \t image Dec size in angular units\n" \
	 "-x <X>,\t --racenter  <X> \t image RA  center in angular units\n" \
	 "-y <Y>,\t --deccenter <X> \t image Dec center in angular units\n" \
	 "-s <file>,\t --srcfile <X> \t name of the text file with " \
	        "specification of sources\n" \	
	 "-i <file>,\t --addimage <file> \t name of a FITs file with image" \
	 " to be added to the output\n" \
	 "-d,\t --debug\t verbose output for debugging\n" \
	 "-h,\t --help\t This help page\n\n" \
	 "Angular units examples: 0.1367rad, 12deg; 460'; 20000\"; " \
	 "17:23:34.29hms; -21:52:43.765dms\n" \
	 "Values in arcminutes must be in double quotes, like \"258'\" \n" \
	 "Values in arcseconds must be in single quotes, like '85400\"' \n\n");
  return;
}
