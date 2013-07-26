/* ------------------------------------------------------------------------------
 *
 * Copyright (C) 2007 Stephen Ord
 *
 * ------------------------------------------------------------------------------ */


/*! hpx2carre : Simple routine to transfer a ring format Healpix sky into a plan carre
  		sky. Written simply to parse the files created by the GSM code. 
	       
		Strongly advise against general use of this routine as it makes assumptions
		specific to the GSM file format.

		More general (and accurate) transformation routines - utilising the WCS library 
		will be employed later.

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

#undef _PLOT
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fitsio.h"
#include <wcslib/wcslib.h>
#include <wcslib/fitshdr.h>
#include <wcslib/getwcstab.h>

#include "utils.h"
#include "wcs_utils.h"

#include "novas.h"

#ifdef _PLOT
#include "cpgplot.h"
#endif 

#define DEGtoRAD (M_PI/180.0) 
#define RADtoDEG (180.0/M_PI)
#define RADtoHOURS (24.0/(2*M_PI))
#define HOURStoRAD ((2*M_PI)/24.0)
#define NOVAS 1
#define GET_1950 0
 
void int_and_norm(float ** , float ** , long , long nrows);
void nearest(float *, float *, long,long,long,long,long,long);

void gal2eq(double ll,double bb,double *ra,double *dec) {

  double sin_dec=0.0;
  double ra_1950 = 0.0;
  double dec_1950 = 0.0;
  double modified_ra = 0.0;
  double x,y;

  
  sin_dec = cos(bb) * sin(ll - (33.0*DEGtoRAD)) * cos((27.4*DEGtoRAD)) + sin(bb) * sin((27.4*DEGtoRAD));

  dec_1950 = asin(sin_dec);

  y = (cos(bb)*cos(ll-(33.0*DEGtoRAD)));
  x = (sin(bb)*cos(27.4*DEGtoRAD) - cos(bb)*sin(27.4*(DEGtoRAD))*sin(ll-(33.0*DEGtoRAD)));
  
  modified_ra = atan2(y,x);

  ra_1950 = modified_ra + (192.25*DEGtoRAD);

  ra_1950 = fmod(ra_1950,(2*M_PI));

  if (GET_1950) {
    	*ra = ra_1950;
  	*dec = dec_1950;
	return;
  }

  /* now we transform for FK4 to FK5 using the NOVAS routines */

  if (NOVAS) {
      
	char *catalog ="FK4";
	char *star_name = "Snoopy";
        long int num = 1; 
        double pm_dec=0;
	double pm_ra=0;
	double para=0.0;
        double rad_vel=0.0;
      
        cat_entry star_1950;
        cat_entry star_2000;

	ra_1950 = ra_1950 * RADtoHOURS;
        dec_1950 = dec_1950 * RADtoDEG;
	
	make_cat_entry(catalog,star_name,num,ra_1950,dec_1950,pm_ra,pm_dec,para,rad_vel,&star_1950);

	transform_cat(2,1950.0,&star_1950,2000.0,"FK5",&star_2000);

	*ra=star_2000.ra * HOURStoRAD;
	*dec=star_2000.dec * DEGtoRAD;

	return;

  }

}
	     	

  

int hpx2carre (const char *hpx, const char *filename, float frequency) {

  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  fitsfile *hpxfile;
  fitsfile *outfile;

  int i,j,k,dd,rval; 
  int status;
  int debug = 0; 
  double pixcrd[4][2];
  long naxis   =   2;
  long naxes[] = {0,0}; // npix across, npix high
  long in_naxes[] = {0,0};

  long pixtoaside=0;
  long first[] = {1,1};
  long in_npix=0;
  long npix = 0;
  
  double delta11,delta22,delta21,delta12,ref1,ref2;
  double pix1,pix2;
  double wcs_temp=0.0;

  /* need to declare some WCS info */
  struct wcsprm *hpx_wcs,*car_wcs;
  int  ncards, nreject, nwcs, stat[NWCSFIX];

  char *hpxheader,*carheader;

  double *in_corners[4];
  double *out_corners[4];
  char tempname[L_tmpnam];

  float inpix[2];

  double polyvert_x[4];
  double polyvert_y[4];
  double pixvert_x[4];
  double pixvert_y[4];

  float outpix1,outpix2;
  
  float pixrange[4];
  double area_poly = 0;
  double area_inpix = 0;
  double pix_bright_contr = 0.0;
  double contributed_area = 0.0;
  int out_offset=0;
  int in_offset=0;

  double lon,lat,ra,dec;

  double clip_x[12];
  double clip_y[12];
  int nvert;

  float nullval=0;
  int n_nulls=0;

  float *input_data=NULL;
  float *data = NULL;
  float *included = NULL;

  hpx_wcs = malloc(sizeof(struct wcsprm));
  car_wcs = malloc(sizeof(struct wcsprm));

#ifdef _PLOT
  cpgopen("?");
  cpgsubp(1,2);
  cpgpanl(1,1);
#endif
  in_corners[0] = (double *) malloc(8*sizeof(double));
  out_corners[0] = (double *) malloc(8*sizeof(double));

  for (i=1; i<4 ;i++) {

    in_corners[i] = in_corners[0] + (2*i);
    out_corners[i] = out_corners[0] + (2*i);
  } 

  rval = 0; 
  fprintf(stdout,"Opening %s\n",hpx);
  if (fits_open_file(&hpxfile,hpx,0,&rval)) {
    fprintf(stderr,"fits_open_file has failed\n");
    fits_report_error(stderr,rval);
  }

  if (fits_get_img_size(hpxfile,2,in_naxes,&rval)) {
    fprintf(stderr,"fits_get_img_dim has failed\n");
    fits_report_error(stderr,rval);
  }

  in_npix = in_naxes[0]*in_naxes[1];
  pixtoaside = in_naxes[0]/2.0;
  npix = pixtoaside*pixtoaside;

  naxes[0] = pixtoaside; /* RA */
  naxes[1] = pixtoaside; /* Dec */

  /* setting up the OUTPUT grid */

  /* lets pretend we are doing the whole sky */
  /* these pix aren't square but hey .... */
  
  delta22= 180.0/pixtoaside;

  /* need our reference pixel and its value in the real world */

  /* pix */
  
  pix1 = pixtoaside/2.0;
  pix2 = pixtoaside/2.0;

  delta11= -360.0/(pixtoaside);
  delta22 = -180/(-1.0*pixtoaside);
  delta12 = 0.0;
  delta21 = 0.0;
 
  ref1 = 180-delta11/2.0;
  ref2 = 0+delta22/2.0;

  /* create new FITS file */
  status = 0;
  tmpnam(tempname);
  remove(tempname); /* clobber existing file, no questions asked */
  if (fits_create_file(&fptr, tempname, &status)) {
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
	__FILE__, __LINE__);

    fits_report_error(stderr,status);
    exit(EXIT_FAILURE);
  }

  if ( fits_create_img(fptr,  FLOAT_IMG, naxis, naxes, &status) ) {
    fprintf(stderr, "%s (%d): Could not create new image file.\n", 
	__FILE__, __LINE__);

    fits_report_error(stderr,status);
    exit(EXIT_FAILURE);
  }

  /* put a WCS in the image. Need CTYPE, CRVAL, CRPIX, CDELT */
  /* From maps_makesky */

  fits_update_key(fptr, TSTRING, "CTYPE1","RA---CAR", NULL,&status);
  fits_update_key(fptr, TSTRING, "CTYPE2","DEC--CAR", NULL,&status);
  fits_update_key(fptr, TSTRING, "RADECSYS","FK5", NULL,&status);
  fits_update_key(fptr, TSTRING, "CUNIT1","deg", NULL,&status);
  fits_update_key(fptr, TSTRING, "CUNIT2","deg", NULL,&status);
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
  
//  fits_update_key(fptr, TDOUBLE, "CD1_1",&delta11,NULL,&status);
//  fits_update_key(fptr, TDOUBLE, "CD2_2",&delta22,NULL,&status);
//  fits_update_key(fptr, TDOUBLE, "CD1_2",&delta12,NULL,&status);
//  fits_update_key(fptr, TDOUBLE, "CD2_1",&delta21,NULL,&status);
  fits_update_key(fptr,TFLOAT,"FREQ",&frequency,"Frequency in MHz",&status);


  if ( fits_close_file(fptr, &status) )       /* close the FITS file */
    fprintf(stderr, "%s (%d): Could not close file.\n", 
	__FILE__, __LINE__);

  // This has created the WCS for the output
  // Next step is to open the two files

  fprintf(stdout,"Opening %s\n",hpx);
  if (fits_open_file(&hpxfile,hpx,0,&rval)) {
    fprintf(stderr,"fits_open_file has failed\n");
    fits_report_error(stderr,rval);
  }
  if (fits_hdr2str(hpxfile, 0, NULL, 0, &hpxheader, &ncards, &rval)) {
    fprintf(stderr,"fits_hdr2str has failed\n");
    fits_report_error(stderr,rval);
  }
  else {
    printf("readin %d cards\n",ncards);
  }  

  rval = wcspih(hpxheader, ncards, WCSHDR_all, 2, &nreject, &nwcs, &hpx_wcs);
  if (rval) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", rval,wcshdr_errmsg[rval]);
    return 1;
  }

/*
  rval = fits_read_wcstab(hpxfile, hpx_wcs->nwtb, (wtbarr *)hpx_wcs->wtb,&rval);
  if (rval) {
    fits_report_error(stderr, rval);
    return 1;
  }

  rval = wcsfix(7, 0, hpx_wcs, stat) ;
  if (rval){
    for (i = 0; i < NWCSFIX; i++) {
      if (stat[i] > 0) {
	fprintf(stderr, "wcsfix ERROR %d: %s.\n", rval,
	    wcsfix_errmsg[stat[i]]);
      }
    }

    return 1;
  }

*/

  if ((rval= wcsset(hpx_wcs))) {
    fprintf(stderr, "wcsset ERROR %d: %s.\n", rval, wcs_errmsg[rval]);
    return 1;
  }


  fprintf(stdout,"Opening %s\n",tempname);
  if (fits_open_file(&fptr,tempname,0,&rval)) {
    fprintf(stderr,"fits_open_file has failed\n");
    fits_report_error(stderr,rval);
  }
  if (fits_hdr2str(fptr, 0, NULL, 0, &carheader, &ncards, &rval)) {
    fprintf(stderr,"fits_hdr2str has failed\n");
    fits_report_error(stderr,rval);
  }
  else {
    printf("readin %d cards\n",ncards);
  }  

  rval = wcspih(carheader, ncards, WCSHDR_all, 2, &nreject, &nwcs, &car_wcs);
  if (rval) {
    fprintf(stderr, "wcspih ERROR %d: %s.\n", rval,wcshdr_errmsg[rval]);
    return 1;
  }

/*
  rval = fits_read_wcstab(fptr, car_wcs->nwtb, (wtbarr *)car_wcs->wtb,&rval);
  if (rval) {
    fits_report_error(stderr, rval);
    return 1;
  }

  rval = wcsfix(7, 0, car_wcs, stat) ;
  if (rval){
    for (i = 0; i < NWCSFIX; i++) {
      if (stat[i] > 0) {
	fprintf(stderr, "wcsfix ERROR %d: %s.\n", rval,
	    wcsfix_errmsg[stat[i]]);
      }
    }

    return 1;
  }

*/

  if ((rval= wcsset(car_wcs))) {
    fprintf(stderr, "wcsset ERROR %d: %s.\n", rval, wcs_errmsg[rval]);
    return 1;
  }
  
  // Now we need the input_data 
  // the output data
  // and included
  input_data = (float *) calloc (in_npix,sizeof(float));
  data = (float *) calloc (npix,sizeof(float));
  included = (float *) calloc (npix,sizeof(float)); 

  nullval = nanf("0");

  // read input data from input file

  if (fits_read_pix(hpxfile,TFLOAT,first,in_npix,&nullval,input_data,&n_nulls,&rval)){
    fprintf(stderr,"fits_read_pix failed\n");
  fits_report_error(stderr,rval);
    exit(EXIT_FAILURE);
  }

  // POLYSAMP

  for (i=1; i <= in_naxes[0]; i++) { /* ~RA axis */

     fprintf(stderr,"                              \r");
     fprintf(stderr,"Projecting column %d/%ld\r",i,in_naxes[0]);
    /* You may have noticesd that the i,j actualy seems to correspond to column,row 
       instead of row,column.
       This is simply because the FITS format is column major and axis1 == RA and axis2 == DEC
     */


    for (j=1; j <= in_naxes[1] ; j++) { /* ~Dec Axis */
      /* which pixel in the output image does (i,j) correspond to? */ 

      pixcrd[0][0] = (float) i;	
      pixcrd[0][1] = (float) j;	


      inpix[0] = pixcrd[0][0];
      inpix[1] = pixcrd[0][1];

      in_offset = (j-1)*in_naxes[0] + (i-1);
      if (isnan(input_data[in_offset])) {
        continue;
      }
      else {
	get_corners(hpx_wcs,inpix,in_corners,0);

	/* so we now have the sky coordinates of the input pixel vertices */
	/* in GLON and GLAT - we will have to slalib them to RA and DEC */
        
	for (k=0;k<4;k++) {

          lon = *(in_corners[0]+(2*k)) * M_PI / 180.0;
	  lat = *(in_corners[0]+(2*k+1)) * M_PI / 180.0;
	  // NO LONGER USED slaGaleq (lon,lat, &ra, &dec);
	  // fprintf(stderr,"RA(SLA) == %f, DEC(SLA) == %f\n",ra,dec);
          gal2eq(lon,lat,&ra,&dec);
	  // fprintf(stderr,"RA(NOVAS) == %f, DEC(NOVAS) == %f\n",ra,dec);
	  *(in_corners[0] + (2*k)) = ra * 180.0/M_PI;
	  *(in_corners[0] + (2*k+1)) = dec * 180.0/M_PI;
	}
	/* unfortunately we need the coordinates in the frame of the output */


	if (debug) {
	  fprintf(stderr,"IN CORNERS\n");
	  for (dd = 0; dd < 4; dd++) {
	    fprintf(stderr,"  (%7.3f,%8.3f)",
		in_corners[dd][0],  in_corners[dd][1]);
	  }
	  fprintf(stderr,"\n");
	}


	get_pixels(car_wcs,pixrange,in_corners,out_corners,0);

	/* the above gets the range over which to test the input pixel against
	   output pixels
	 */

	if (debug) {
#ifdef _PLOT

	  cpgpanl(1,2);
	  cpgeras();
	  cpgswin(pixrange[0]-0.5,pixrange[1]+0.5,pixrange[2]-0.5,pixrange[3]+0.5);
	  cpgbox("ABCN",0.0,0.0,"ABCN",0.0,0.0);
	  cpglab("HPX Pixel Number","HPX Pixel Number","Input Pixel in Output Pixel Frame");
	  cpgslw(2);

#endif

	  fprintf(stderr,"in pixel %d,%d needs touches = (%f,%f,%f,%f) \n",i,j,pixrange[0], pixrange[1],pixrange[2],pixrange[3]);

	  fprintf(stderr,"OUT CORNERS\n");

	  for (dd = 0; dd < 4; dd++) {
	    fprintf(stderr,"  (%7.3f,%8.3f)",
		out_corners[dd][0],  out_corners[dd][1]);

	    pixvert_x[dd] = out_corners[dd][0]; 
	    pixvert_y[dd] = out_corners[dd][1]; 

#ifdef _PLOT
	    inpix_x[dd] = (float)out_corners[dd][0]; 
	    inpix_y[dd] = (float)out_corners[dd][1]; 
#endif
	  }
	  fprintf(stderr,"\n");

#ifdef _PLOT
	  cpgpanl(1,2);
	  cpgsfs(0);
	  cpgpoly(4,inpix_x,inpix_y);
#endif
	}
	for (outpix1=pixrange[0];outpix1<=pixrange[1];outpix1++) {

	  for (outpix2=pixrange[2];outpix2<=pixrange[3];outpix2++) {

	    if (fabs(pixrange[1] - pixrange[0]) > 20) {
	      continue;
	    }
	    else if (fabs(pixrange[3] - pixrange[2]) > 20) {
	      continue;
	    }


	    if (outpix1 > pixtoaside || outpix2 > pixtoaside || outpix1 < 1 || outpix2 < 1) {
	      continue;
	    }

	    /* get input pixel in output pixel space */


	    if (debug) {
	      fprintf(stderr,"Output clipping xmin=%f,xmax=%f,ymin=%f,ymax=%f\n",outpix1-0.5,outpix1+0.5,outpix2-0.5,outpix2+0.5);
	      for (dd = 0; dd < 4; dd++) {
		fprintf(stderr,"  (%7.3f,%8.3f)",
		    out_corners[dd][0],  out_corners[dd][1]);
	      }
	      fprintf(stderr,"\n");
	    }
	    for (dd = 0; dd < 4; dd++) {
	      polyvert_x[dd] = out_corners[dd][0];
	      polyvert_y[dd] = out_corners[dd][1];


	    }

	    nvert = rectClip(4,polyvert_x,polyvert_y,clip_x,clip_y,outpix1-0.5,outpix2-0.5,outpix1+0.5,outpix2+0.5);
#ifdef _PLOT
	    cpgbbuf();
	    cpgsfs(0);
	    cpgpanl(1,2);	  
	    cpgswin(pixrange[0]-0.5,pixrange[1]+0.5,pixrange[2]-0.5,pixrange[3]+0.5);
	    cpgrect(outpix1-0.5,outpix1+0.5,outpix2-0.5,outpix2+0.5); 

	    cpgsfs(1);
	    if (debug) {
	      for (dd=0;dd<nvert;dd++) {
		fprintf(stderr,"clip_x %d = %f, clip_y %d = %f\n",dd,clip_x[dd], dd, clip_y[dd]);
		poly_x[dd] = (float) clip_x[dd];
		poly_y[dd] = (float) clip_y[dd];
	      }
	    }
	    cpgpoly(nvert,poly_x,poly_y);
	    usleep(125000);
#endif
	    area_poly=0;
	    area_poly = Bourke(nvert,clip_x,clip_y);
	    // area_poly = convexArea(nvert,clip_x,clip_y);
	    // area_inpix = convexArea(4,pixvert_x,pixvert_y);
	    // area_poly = area_poly/area_inpix;

	    if (debug) {
	      fprintf(stderr,"INPUT PIXEL AREA: %f\n",area_inpix);
	    }

	    if (debug && nvert > 2) {
	      fprintf(stderr,"Overlap polygon has %d vertices and area %f\n",nvert,area_poly);
	    }

	    out_offset = (rint(outpix2)-1)*pixtoaside + (rint(outpix1)-1); // same ordering as image

	    in_offset = (j-1)*in_naxes[0] + (i-1);

	    if (included[out_offset] == 0) {
	      data[out_offset] = 0.0;
	    }

	    if (included[out_offset] < 100.0) { // superfluous now

	      if (area_poly > 0.0) {

		contributed_area = included[out_offset] + area_poly;

		pix_bright_contr = input_data[in_offset]*area_poly;

		/* here we write it out */
		data[out_offset] = data[out_offset] + (pix_bright_contr);
		included[out_offset] = contributed_area;
	      }
#ifdef _PLOT
	      cpgpanl(1,1);
	      cpgswin(1,regrid->output_image->sizex,1,regrid->output_image->sizey);
	      cpgbox("ABCN",0.0,0.0,"ABCN",0.0,0.0);
	      cpglab("HPX Pixel Number","HPX Pixel Number","Output Pixel Frame");


	      cpgsci(rint(floorf(included[out_offset])+1));
	      cpgrect(outpix1-0.5,outpix1+0.5,outpix2-0.5,outpix2+0.5); 
	      cpgsci(1);
	      cpgebuf();
	      cpgupdt();
#endif
	      if (debug) {
		fprintf(stdout," IN (%d,%d) OUT (%f,%f) area %e, and total contibuted area so far %e \n",
		    i,j,outpix1,outpix2,area_poly,contributed_area);

	      }
	    }
	    if (included[out_offset] > 1.0) {
	      if (debug) {
		fprintf(stdout," ADVISORY IN (%d,%d) OUT (%f,%f) contibuted area so far %e (GREATER THAN 1.0) \n",
		    i,j,outpix1,outpix2,contributed_area);

	      }
	    }
	  }
	}
      }
    }
  }
  // unload the new one
  // quick normalisation

  for (i=0;i<npix;i++) {
    if (isnan(data[i]) || isnan(included[i])) {
      fprintf(stderr,"%d is NAN\n",i);
      data[i] = 0;
      included[i] = 1.0;
    }
    if (included[i] != 0.0) {
      data[i] = data[i]/included[i];
      included[i] = 1.0;
    }
  }

  for (i=0;i<npix;i++) {
    if (included[i] == 0) {
      /* this is going to be a bit strange */
      /* 1 d kludge */
      float pixies = (float) i;
      float pixie2 = (float) pixtoaside;

      long row =  floor (pixies/pixie2);
      long col = rint (i%pixtoaside);

      fprintf(stderr,"OFFSET = %d, ROW = %ld, COL = %ld\n",i,row,col);
      nearest(data, included, pixtoaside, pixtoaside, row, row, col, col);
    }
  } 
  


 //STEVE

  remove(filename); /* clobber existing file, no questions asked */
  if (fits_create_file(&outfile, filename, &status)) {
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
	__FILE__, __LINE__);
    fits_report_error(stderr,status);
  }

  if (fits_copy_header(fptr,outfile,&status)) {
    fprintf(stderr,"Cannot copy header\n");
    fits_report_error(stderr,status);
  }

  if (fits_write_img(outfile,TFLOAT,1,npix,data,&status)) {
    fprintf(stderr,"Cannot write img\n");
    fits_report_error(stderr,status);
  }

  if (fits_close_file(hpxfile,&status)) {
  
    fprintf(stderr,"Cannot close hpxfile\n");
    fits_report_error(stderr,status);
  }
  
  if (fits_close_file(fptr,&status)) {
    fprintf(stderr,"Cannot close tempfile\n");
    fits_report_error(stderr,status);
  }

  if (fits_close_file(outfile,&status)) {
    fprintf(stderr,"Cannot close outfile\n");
    fits_report_error(stderr,status);
  }

  remove(tempname);
  return 0;

}


int healpix2carre (float *signal, long nside, const char *filename, int nest, float frequency) {

  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int ii,jj,ni,nj; 
  int status;
  
  long naxis   =   2;
  long naxes[] = {0,0}; // npix across, npix high
  
  long pixtoaside = 10*nside;
  long first[2] = {1,1};

  long npix = 100*nside*nside;
  long in_npix = 12*nside*nside;
  
  float **hpx_map;
  float *memory;
  float *hits;
  float **hits_arr;
  
  double lat,lon,ra,dec;
  double area=1.0;

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

  delta11= -360.0/(pixtoaside);
  delta22 = -180/(-1.0*pixtoaside);
  delta12 = 0.0;
  delta21 = 0.0;
 
  ref1 = 180-delta11/2.0;
  ref2 = 0+delta22/2.0;


  /* just giving a dynamic 2 dimensional image array */ 
  memory = (float *) calloc(npix,sizeof(float));
  hpx_map = (float **) calloc(pixtoaside,sizeof(float *));

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
    hpx_map[jj] = &memory[ii];
    hits_arr[jj] = &hits[ii];
    jj++;
  }	
  /* Each entry in hpx_map can be thought of a column of the image
     remeber FITs images are stored like FORTRAN ones are. */

  /* initialize status before calling fitsio routines */
  status = 0;

  /* create new FITS file */

  remove(filename); /* clobber existing file, no questions asked */
  if (fits_create_file(&fptr, filename, &status)) 
    fprintf(stderr, "%s (%d): Could not create new fits file.\n", 
	__FILE__, __LINE__);
  /*
     for (jj = 0; jj < naxes[1]; jj++) {
     for (ii = 0; ii < naxes[0]; ii++) {
     hpx_map[jj][ii] = ii + jj;
     } 
     }
   */
  
  area = 4*M_PI/in_npix; // because we are equal area
  
  for (jj = 0; jj < in_npix; jj++) {
    if (!nest) {

      pix2ang_ring(nside,jj,&lat,&lon);

      /* this has NOT come out as Galactic Lat and Lon -- I hope
	 going to change this to FK5 via a SLALIB call */

      lat = M_PI/2.0 - lat;

      // slaGaleq (lon,lat, &ra, &dec);

      // fprintf(stderr,"RA(SLA) == %f, DEC(SLA) == %f\n",ra,dec);
      gal2eq(lon,lat,&ra,&dec);
      // fprintf(stderr,"RA(NOVAS) == %f, DEC(NOVAS) == %f\n",ra,dec);
      /* comes back in rad - need degree measure for FITS */


      ra = ra * 180.0/M_PI ;
      dec = dec * 180.0/M_PI ;

      wcs_temp = dec + 90.0;

      /* Going to drop it into plan carree as its easiest */
      /* RA increases from right to left */

      ni = pixtoaside + (ra/delta11);
      nj = (wcs_temp/delta22);
      /* 
	 printf("pixel %ld -> %ld,%ld: x=%f,y=%f : SIG: %f \n",jj,ni,nj,ra,dec,signal[jj]); 
       */

/*
      if (ni == (int) (pixtoaside/2.0)) {
	ref1 = ra;
      }
      if (nj == (int) (pixtoaside/2.0)) {
	ref2 = dec;
      }
*/
      
//      hpx_map[nj][ni] = convert_to_Iv(signal[jj],frequency,area);
      hpx_map[nj][ni] = signal[jj];
      
	 //  printf("VALUE: %g\n",hpx_map[nj][ni]);
      hits_arr[nj][ni]++;

      if (hits_arr[nj][ni] > 100) {

	printf("pixel %d -->%d,%d -- l=%f,b=%f : ra=%f, dec=%f\n",jj,ni,nj,lon,lat,ra,dec); 
      }
    }
  }
  /* Serious problem here as the HPX does not map simply on to a plane */
  /* Essentially there are holes in the pixelation which have to be interpolated across */
  /* You can see why in "Mapping on the HealPix grid by M. Calabretta" */

  /* there is a better way of doing this now but it involves a pixel re-ordering */

  int_and_norm(hpx_map, hits_arr, naxes[1], naxes[0]); 


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
  
 // fits_update_key(fptr, TDOUBLE, "CD1_1",&delta11,NULL,&status);
 // fits_update_key(fptr, TDOUBLE, "CD2_2",&delta22,NULL,&status);
 // fits_update_key(fptr, TDOUBLE, "CD1_2",&delta12,NULL,&status);
 // fits_update_key(fptr, TDOUBLE, "CD2_1",&delta21,NULL,&status);
  fits_update_key(fptr,TFLOAT,"FREQ",&frequency,"Frequency in MHz",&status);


  fits_write_pix(fptr, TFLOAT, first, npix, hpx_map[0], &status); 

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

void nearest(float * map, float *hits, long ncols, long nrows, long startrow, long stoprow, long startcol, long stopcol) {

  int col,row,scol,srow,offset,cshift,rshift;
  float newval,newhits;
  int maxoffset = 20;
  float minhits = 2;


  // printf("Nearest neighbour matching."); 
  for (row = startrow; row <= stoprow ; row++) {
    for (col = startcol; col <= stopcol ; col++) {
      int p_offset = (row * nrows) + col;
     // fprintf(stderr,"NN offset %d\n",p_offset);
      if (hits[p_offset] == 0) {
	newhits = 0;
	newval = 0.0;
	offset = 1;
//	printf("COL %d- ROW %d ...\n",col,row); 
	while (newhits < minhits && (offset < maxoffset)) {

	  srow = 0;
	  scol = 0;
	  newhits = 0;
	  newval = 0.0;
	 // printf("OFFSET %d\n",offset); 
	  for (cshift = -1*offset; cshift < offset+1; cshift++) {
	    scol = col + cshift;
	    // printf("SCOL: %d ---- \n",scol);

	    if ((scol >=0) && (scol < ncols)) {

	      for (rshift = -1*offset; rshift < offset+1; rshift++) {
		srow = row + rshift;	
	 //	printf("SROW: %d ---- ",srow); 
		if ((srow >= 0) && (srow < nrows)) {
                  int n_offset = (srow * nrows) + scol; 
		  if (hits[n_offset] > 0) {
                    
//		    printf("(%d,%d) SHIFTED -- HITS %f: VAL %f \n",scol,srow,hits[n_offset],map[n_offset]); 
//		    printf("(%d,%d) ORIG -- HITS %f: VAL %f \n",col,row,hits[p_offset],map[p_offset]); 
		    newhits = newhits + hits[n_offset]; 
		    newval = newval + map[n_offset];
//		    printf("(%d,%d) NEW -- HITS %f: VAL %f \n",col,row,newhits,newval); 
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
	  hits[p_offset] = 1.0;
	  map[p_offset] = newval/newhits;
	}
      }
    }
  }
//  printf(".done.\n"); 

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


