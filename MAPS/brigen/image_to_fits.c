#include <stdio.h>
#include <string.h>
#include <fitsio.h>
#include "brigen.h"

/*
 * Write resulting brightness array 'image' to file in FITS format
 */
 
void image_to_fits(struct bg_param *bgparam, double *image) {	
 
  double bscale = 1.0, bzero = 0.0, zero = 0.0;
  double crpix1 = bgparam->ncol/2 + 1; /* RA reference pixel */
  double crpix2 = bgparam->nrow/2 + 1; /* DEC reference pixel */
  double crval1 = bgparam->xcenter;  /* FOV RA  center in degrees */
  double crval2 = bgparam->ycenter;  /* FOV DEC center in degrees*/
  /* RA and DEC  increments in degrees */
  double cdelt1 = -bgparam->xsize/(bgparam->ncol - 1)/3600.0; /* degrees */
  double cdelt2 =  bgparam->ysize/(bgparam->nrow - 1)/3600.0; /* degrees */
  fitsfile *fptr = NULL;
  int status = 0;
  long startpixel = 1, naxis = 2, npixels;
  long naxes[2]; 
  char fitsfilename[256] = "!"; /* CFITSIO: '!' before filename means rm old */

  strcat(fitsfilename, bgparam->fitsfilename);
  naxes[0] = bgparam->ncol; /* Number of pixels in RA dimension */ 
  naxes[1] = bgparam->nrow; /* Number of pixels in DEC dimension */ 
  npixels = naxes[0]*naxes[1];
  
  /* create new file */
  fits_create_file(&fptr, fitsfilename, &status);   
  fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
  fits_update_key(fptr, TUINT, "NAXIS1", &(naxes[0]),
		  "X, pixels along RA dimensions", &status);
  fits_update_key(fptr, TUINT, "NAXIS2", &(naxes[1]),
		  "Y, pixels along DEC dimensions", &status);
  fits_update_key(fptr, TDOUBLE, "BSCALE", &bscale, 
		  "Brightness scale", &status);
  fits_update_key(fptr, TDOUBLE, "BZERO", &bzero, 
		  "Brightness zero level?", &status);
  fits_update_key(fptr, TSTRING, "BUNIT", "Jy/steradian", 
		  "Brightness (flux) units", &status);
  fits_update_key(fptr, TSTRING, "CTYPE1", "RA---SIN", "Projection type", 
		  &status);
  fits_update_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "X reference pixel", 
		  &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1, 
		  "FOV RA center in degrees", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1, "RA increment in degrees", 
		  &status);
  fits_update_key(fptr, TDOUBLE, "CROTA1", &zero, "Rotation", 
		  &status);
  fits_update_key(fptr, TSTRING, "CTYPE2", "DEC---SIN", "Projection type?", 
		  &status);
  fits_update_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "X reference pixel", 
		  &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2, 
		  "FOV DEC center in degrees", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2, "DEC increment in degrees", 
		  &status);
  fits_update_key(fptr, TDOUBLE, "CROTA2", &zero, "Rotation", 
		  &status);


  fits_write_img(fptr, TDOUBLE, startpixel, npixels, image, &status);
  fits_close_file(fptr, &status);     /* close the file */
  fits_report_error(stderr, status);  /* print out any error messages */

/* FITS header */
/*   iTemp = 2; */
/*   dTemp = 0.0; */
/*   dTemp = -brightnessMinX + 1;	 */
/*   // FITS numbers pixels from 1, not 0: */
/*   FITSwrite_headerLine("CRPIX1", Double, &dTemp, "");  */
/*   dTemp = fieldRA/oneDegree; */
/*   // CRVAL1 in degrees, not radians */
/*   FITSwrite_headerLine("CRVAL1", Double, &dTemp, "");  */
/*   dTemp = dx/oneDegree; */
/*   // No sign change, it seems to be corrected in visgen */
/*   FITSwrite_headerLine("CDELT1", Double, &dTemp,"");  */
/*   dTemp = 0.0; */
/*   FITSwrite_headerLine("CROTA1", Double, &dTemp, ""); */
/*   FITSwrite_headerLine("CTYPE2", CharString, "DEC--SIN", ""); */
/*   dTemp = -brightnessMinY + 1; */
/*   FITSwrite_headerLine("CRPIX2", Double, &dTemp, ""); */
/*   dTemp = fieldDec/oneDegree;	 */
/*   FITSwrite_headerLine("CRVAL2", Double, &dTemp, ""); */
/*   dTemp = dy/oneDegree; */
/*   FITSwrite_headerLine("CDELT2", Double, &dTemp, ""); */
/*   dTemp = 0.0; */
/*   FITSwrite_headerLine("CROTA2", Double, &dTemp, ""); */
/*   FITSwrite_endHeader(); */
/*   n=FITSwrite_array2D(image, brightnessXN, brightnessYN, datasize); */
/*   FITSwrite_endFile(); */
} 
