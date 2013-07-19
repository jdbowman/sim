/* ------------------------------------------------------------------------------
 *
 * Copyright (C) 2007 Stephen Ord
 *
 * ------------------------------------------------------------------------------ */


/*! dat2carre : 		Stephen Ord (October 2007) <sord@cfa.harvard.edu	
  				simple routine to take a column-major format, square
				sky and add the WCS to it for a plan-caree projection
*/

/* The input arguments:
 *
 *    signal   : signal in each pixel
 *    nside    : parameter to define simulation resolution 
 *    filename : name of the fits file written
 *
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fitsio.h"
#include "utils.h"

void int_and_norm(float ** , float ** , long , long nrows);

int dat2carre(float *signal, long nside, const char *filename, float freq) {

  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int ii,jj; 
  int status;
  double area = 1.0; // this needs fixing  
  long naxis   =   2;
  long naxes[] = {nside,nside}; // npix across, npix high
  long pixtoaside = nside; 
  long first[2] = {1,1};

  long npix = nside*nside;
  long in_npix = nside*nside;
  
  float **dat_map;
  float *memory;
  float *hits;
  float **hits_arr;
  
  
  double delta11,delta22,delta21,delta12,ref1,ref2;
  double pix1,pix2;
  double wcs_temp=0.0;

  naxes[0] = pixtoaside; /* RA */
  naxes[1] = pixtoaside; /* Dec */

  /* lets pretend we are doing the whole sky */
  /* these pix aren't square but hey .... */
  
  delta22= 180.0/pixtoaside;

  /* need our reference pixel and its value in the real world */

  /* pix */
  
  pix1 = pixtoaside/2.0;
  pix2 = pixtoaside/2.0;

  delta11= 360.0/(pixtoaside);
  delta22 = -180/(-1.0*pixtoaside);
  delta12 = 0.0;
  delta21 = 0.0;
 
  ref1 = 180-delta11/2.0;
  ref2 = 0+delta22/2.0;


  /* just giving a dynamic 2 dimensional image array */ 
  memory = (float *) calloc(npix,sizeof(float));
  dat_map = (float **) calloc(pixtoaside,sizeof(float *));

  if (memory == NULL) {
    fprintf(stderr,"Cannot allocate memory for array\n");
    exit(EXIT_FAILURE);
  }

  hits = (float *) calloc(npix,sizeof(float));
  hits_arr = (float **) calloc(pixtoaside,sizeof(float *));

  if (hits == NULL) {
    fprintf(stderr,"Cannot allocate memory for array\n");
    exit(EXIT_FAILURE);
  }


  jj=0;
  for (ii=0; ii<npix; ii=ii+pixtoaside) {
    dat_map[jj] = &memory[ii];
    hits_arr[jj] = &hits[ii];
    jj++;
  }	
  /* Each entry in hpx_map can be thought of a column of the image
     remeber FITs images are stored like FORTRAN ones are. */

  /* initialize status before calling fitsio routines */
  status = 0;

  /* create new FITS file */

  if (fits_create_file(&fptr, filename, &status)) 
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
	__FILE__, __LINE__);

  area=1.0; // for the moment
  for (jj = 0; jj < in_npix; jj++) {
      memory[jj] = convert_to_Jy((double)signal[jj],(double) freq,area );
      hits[jj]++;

  }
  /* Serious problem here as the HPX does not map simply on to a plane */
  /* Essentially there are holes in the pixelation which have to be interpolated across */
  /* You can see why in "Mapping on the HealPix grid by M. Calabretta" */

  int_and_norm(dat_map, hits_arr, naxes[1], naxes[0]); 


  if ( fits_create_img(fptr,  FLOAT_IMG, naxis, naxes, &status) )
    fprintf(stderr, "%s (%d): Could not create new image file.\n", 
	__FILE__, __LINE__);


  /* put a WCS in the image. Need CTYPE, CRVAL, CRPIX, CDELT */
  /* From maps_makesky */

  fits_update_key(fptr, TSTRING, "CTYPE1","RA---CAR", NULL,&status);
  fits_update_key(fptr, TSTRING, "CTYPE2","DEC--CAR", NULL,&status);
  fits_update_key(fptr, TSTRING, "CUNIT1","deg", NULL,&status);
  fits_update_key(fptr, TSTRING, "CUNIT2","deg", NULL,&status);
  fits_update_key(fptr, TSTRING, "RADECSYS","FK5", NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CRPIX1",&pix1, NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CRPIX2",&pix2, NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CRVAL1",&ref1, NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2",&ref2, NULL,&status);
  wcs_temp = 1; 
  fits_update_key(fptr, TDOUBLE, "BSCALE",&wcs_temp, NULL,&status);
  wcs_temp = 0; 
  fits_update_key(fptr, TDOUBLE, "BZERO",&wcs_temp, NULL,&status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1",&delta11,NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CDELT2",&delta22,NULL,&status);
  
  fits_update_key(fptr, TDOUBLE, "CD1_1",&delta11,NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CD2_2",&delta22,NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CD1_2",&delta12,NULL,&status);
  fits_update_key(fptr, TDOUBLE, "CD2_1",&delta21,NULL,&status);
  fits_update_key(fptr, TFLOAT, "FREQ",&freq,"Frequency in MHz",&status);


  fits_write_pix(fptr, TFLOAT, first, npix, dat_map[0], &status); 

  /*   fits_write_pix(fptr, TFLOAT, first, npix, hits_arr[0], &status); */ 

  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n", 
	__FILE__, __LINE__);

  fits_report_error(stderr,status);
  
  /*
  
     remapimage(filename);
  
   */

  return 0;

}

void int_and_norm(float ** map, float ** hits, long ncols, long nrows) {

  int col,row,scol,srow,offset,cshift,rshift;
  float newval,newhits;
  int maxoffset = 20;
  float minhits = 2;


  printf("Nearest neighbour matching."); 
  for (col = 0; col < ncols ; col++) {
    for (row = 0; row < nrows ; row++) {
      if (hits[col][row] == 0) {
	newhits = 0;
	newval = 0.0;
	offset = 1;
	/* printf("COL %d- ROW %d ...\n",col,row); */
	while (newhits < minhits && (offset < maxoffset)) {

	  srow = 0;
	  scol = 0;
	  newhits = 0;
	  newval = 0.0;
	  /* printf("OFFSET %d\n",offset); */
	  for (cshift = -1*offset; cshift < offset+1; cshift++) {
	    scol = col + cshift;
	    /* printf("SCOL: %d ---- \n",scol);*/

	    if ((scol >=0) && (scol < ncols)) {

	      for (rshift = -1*offset; rshift < offset+1; rshift++) {
		srow = row + rshift;	
		/*  printf("SROW: %d ---- ",srow); */
		if ((srow >= 0) && (srow < nrows)) {

		  if (hits[scol][srow] > 0) {

		    /*	    printf("(%d,%d) OLD -- HITS %f: VAL %f \n",scol,srow,hits[scol][srow],map[scol][srow]); */
		    /*	    newhits = newhits + hits[scol][srow]; */
		    newhits = newhits+1; 
		    newval = newval + map[scol][srow];
		    /*	    printf("(%d,%d) NEW -- HITS %f: VAL %f \n",col,row,newhits,newval); */
		  }
		  else {
		    /* printf("EMPTY\n");*/
		  }

		}
		/* printf("OUTOFRANGE\n"); */
	      }

	    }

	  }
	  offset++;
	}


	if (newhits > 0) {
	  hits[col][row] = 1.0;
	  map[col][row] = newval/newhits;
	}
      }
    }
  }
  printf(".done.\n"); 

  /* normalising does not have to be done */
  /*
  printf("Normalising by hits array.\n");

  for (col = 0; col < ncols ; col++) {
    for (row = 0; row < nrows ; row++) {

      if (hits[col][row] > 0) {
		map[col][row] = map[col][row]/hits[col][row];
      }
    }
  }
  */
  

}


