#define MAXEXTIMAGES 16
#define MAXNAXES 16


#ifndef BGPARAM
#define BGPARAM
/*
 * brigen input parameters
 */ 

struct bg_param {
  int ncol;   /* Number of columns of brightness image in pixels */ 
  int nrow;   /* Number of rows of brightness image in pixels */ 
  double xsize;    /* RA size of brightness image in arcseconds */  
  double ysize;    /* DEC size of brightness image in arcseconds */ 
  double xcenter;    /* RA center of brightness image in degrees */  
  double ycenter;    /* DEC center of brightness image in degrees */ 
  struct Source *source; /* Link to the array of sources */
  int nsrc;              /* Number of the sources */
  int nimgf;             /* Number of fits files with external images */
  char *simfilename; /* Simulation project file name */
  char *srcfilename; /* Source list file name */
  char *fitsfilename; /* Output brightness image fits file name */
  char *imgfilename[MAXEXTIMAGES]; /* Additional image fits file names */
  /* Bit flags. 1: parameter entered in command line, 2: in project file,
   *            0: not set at all. */ 
  struct { 
    unsigned int ncol: 2;
    unsigned int nrow : 2;
    unsigned int xsize : 2;
    unsigned int ysize : 2;
    unsigned int xcenter : 2;
    unsigned int ycenter : 2;
    unsigned int srclist : 2;
    unsigned int fitsfilename : 2;
  } is_set; 
};
#endif

void sources_to_image(struct bg_param *bgparam, double *image);
void add_extimages(struct bg_param *bgparam, double *image);
void image_to_fits(struct bg_param *bgparam, double *image);
