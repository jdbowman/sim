#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <fitsio.h>
#include "brigen.h"

/*
 * add_extimage.c
 *
 * Add external image(s) from FITs file(s) to brightness array 'image'
 */

void add_extimages(struct bg_param *bgparam, double *image) {
  int result;      /* Error code */
  int naxis;       /* Number of image dimensions */
  int anynul = 0;  /* Nobody knows what is this for :) */
  long naxes[MAXNAXES];   /* RA and Dec image dimensions in pixels */
  fitsfile *infile = NULL;    /* FITs file handle */
  double *data = NULL; /* The image read from a fitsfile */
  double *dptr = NULL, *imptr = NULL; /* Pointers into data and image */
  long ncol = bgparam->ncol, nrow = bgparam->nrow;
  long imsize = ncol*nrow;
  long i, j;
  char *fname = NULL; /* holds current bgparam->imgfilename[i] */

  /*
   * Allocate the data array to acquire the images 
   * from external FITs files 
   */
  data = (double *) malloc(ncol*nrow*sizeof(double));
  if(!data)  { perror("malloc"); exit(1); }

  /*
   * In a loop read all external images 
   * and add them to the image array.
   */

  for (i = 0; i < bgparam->nimgf; i++) {
    fname = bgparam->imgfilename[i];

    fits_open_image(&infile, fname, READONLY, &result);
    if (result !=0) {
      printf("Failed to open FITs file '%s'; Error code: %d\n", 
	     fname, result);
      exit(1);
    }

    fits_get_img_dim(infile, &naxis, &result);
    if (result != 0) {
      printf("Failed to get image naxis for FITs file '%s'; "		\
	     "error code %d\n", fname, result);
      exit(1);
    }
    if (naxis != 2) {
      printf("The FITs image in '%s' is not 2-dimensional; naxis = %d\n", 
	     fname, naxis);
      /* exit(1); */
    }

    fits_get_img_size(infile, naxis, naxes, &result);
    if (result != 0) {
      printf("Failed to get axis sizes for file '%s'; error code %d\n",
	     fname, result);
      exit(1);
    }
    
    /*
     * Check if the dimensions above 2 are 
     * degenerate: naxes[3:naxis] == 1.
     * If not, exit.
     */
    if (naxis > 2) {
      int degenerate = 1; /* Assume naxes above 2 are degenerate */ 
      for (j = 2; j < naxis; j++) 
	if (naxes[j] != 1) degenerate = 0;
      if (!degenerate) exit(1); /* Message has been issued already */
    }
    if (naxes[0] != ncol || naxes[1] != nrow) {
      printf("The dimensions of external FITs image '%s', %ld x %ld, \n" \
	     "are not equal to those of output image, %ld x %ld\n", 
	     fname, naxes[0], naxes[1], ncol, nrow);
      exit(1);
    }
    
    fits_read_img(infile, TDOUBLE, 1, imsize, NULL, data, 
		  &anynul, &result );
    if (result !=0) {
      printf("Reading FITS data from '%s' failed with error code %d\n", 
	     fname, result);
      exit(1);
    }

    /* Add the image from FITs file to the image array */
    imptr = image;
    dptr = data;
    for (j = 0; j < imsize; j++) *(imptr++) += *(dptr++);

    fits_close_file(infile, &result);
    if (result !=0) 
      printf("WARNING: fits_close_file() failed with code %d.\n", result);
  }
  free(data);
}
