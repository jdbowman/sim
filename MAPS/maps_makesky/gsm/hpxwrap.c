/* ------------------------------------------------------------------------------
 *
 * Copyright (C) 2007 Stephen Ord
 *
 * ------------------------------------------------------------------------------ */


/*! hpxwrap :   Simple routine to transfer a ring format Healpix sky into a FITs file
  		as a binary table. Written simply to parse the files created by the GSM code. 
	       
		Strongly advise against general use of this routine as it makes assumptions
		specific to the GSM file format - namely NSIDE and the COORDS.

		More general (and accurate) transformation routines - utilising the WCS library 
		will be employed later.

		In fact this is designed to provide an interface between M. Calabretta's HPXcvt
		and the GSM

		Stephen Ord (October 2007) <sord@cfa.harvard.edu	
*/

/* The input arguments:
 *
 *    signal   : signal in each pixel
 *    nside    : parameter to define simulation resolution 
 *    filename : name of the fits file written
 *    nest     : structure type of the input HealPix signal (RING/NESTED) (0/1)
 *
*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fitsio.h"
#include "utils.h"


int hpxwrap (float *signal, long nside, const char *filename, int nest, float frequency) {

  int tfields,status;
  long nrows,firstrow,firstelem,nelements,colnum;

  char extname[] = "GSM_SKY";
  char *ttype[] = {"Flux"};
  char *tform[] = {"1E"};
  char *tunit[] = {" "};

   fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */

  /* create the fits file */

  status = 0;
  if (remove(filename)) { /* clobber existing file, no questions asked */
    fprintf(stderr,"%s (%d): cannot remove file %s\n",__FILE__, __LINE__,filename);
  }
  if (fits_create_file(&fptr, filename, &status)) {
    fprintf(stderr, "%s (%d): Could not create new fits file called %s. status: %d\n",__FILE__, __LINE__,filename,status);
    return EXIT_FAILURE;
  }
  /* create a new table */
 /* write the header keywords */

  status = 0;
  nrows = 1;
  tfields = 1;

  if (fits_create_tbl(fptr,BINARY_TBL,nrows,tfields,ttype,tform,tunit,extname, &status)) {
    fprintf(stderr, "%s (%d): Could not create new BINARY TABLE.\n",__FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  /* we have to write some more header information */


  /* write the floats to the first column */

  colnum = 1;
  firstrow = 1;
  firstelem = 1;
  nelements = 12 * nside * nside;

  fits_update_key(fptr,TSTRING,"COORDSYS","G","Coordinate SYSTEM",&status);
  fits_update_key(fptr,TSTRING,"ORDERING","RING","HEALPix ORDERING Scheme",&status);
  fits_update_key(fptr,TULONG,"NSIDE",&nside,"Number of pix to a side",&status);


  if (fits_write_col(fptr,TFLOAT,colnum,firstrow,firstelem,nelements,signal,&status)) {
    fprintf(stderr, "%s (%d): Could not write column.\n",__FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  if ( fits_close_file(fptr, &status) ) {       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n",__FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
} 


