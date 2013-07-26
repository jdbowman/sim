#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "fitsio.h"

#define NDEBUG 0

#include <assert.h>

#include <wcslib/wcslib.h>
#include <wcslib/fitshdr.h>

/* WARNING: outpix must be a preallcated double[4][2]
            pixrange must be a float[4]

	    */

int get_pixels (struct wcsprm *wcs, float *pixrange, double *corners[4], double *outpix[4], int range) {
  int debug = 0;  
  double phi[4],imgcrd[4][2], theta[4];
  int stat[4];
  double **pixcrd; 
  float max_long,max_lat,min_long,min_lat;
  
  int status = 0, k=0;

  pixcrd =  outpix;
  
  if ((status = wcss2p(wcs, 4, 2, corners[0],phi,theta,imgcrd[0],pixcrd[0],stat))) {
    fprintf(stderr, "\n\nwcss2p ERROR %d: %s.\n", status,
	wcs_errmsg[status]);
  
    /* Invalid pixel coordinates. */

  }

  /* now find the range of pixels */
  max_long = 0;
  min_long = 0;
  max_lat = 0;
  min_lat = 0;
	
  for (k=0 ;  k<4 ; k++) {
    if (debug) {
      printf("%f:%f\n",pixcrd[k][0],pixcrd[k][1]);
    }
        
    if (max_long < pixcrd[k][0]) {
      max_long = pixcrd[k][0];
    }
    if (max_lat < pixcrd[k][1]) {
      max_lat = pixcrd[k][1];
    }
  
  }
  
  min_long = max_long;
  min_lat = max_lat;
	
  for (k=0 ;  k<4 ; k++) {
  
    if (min_long > pixcrd[k][0]) {
      min_long = pixcrd[k][0];
    }
    if (min_lat > pixcrd[k][1]) {
      min_lat = pixcrd[k][1];
    }
  
  }
  /* added the extra padding to cope with distortions */
  /* may have to reverse the logic to get this right */

  pixrange[0] = (float) rint ((min_long)-range);
  pixrange[1] = (float) rint ((max_long)+range);
  pixrange[2] = (float) rint ((min_lat)-range);
  pixrange[3] = (float) rint ((max_lat)+range);


  return status;;

}

int get_corners (struct wcsprm *wcs, float pix[2], double *corners[4], int image) {

  double phi[4], pixcrd[4][2], imgcrd[4][2], theta[4], world[4][2];
  int stat[4];

  int k=0,status = 0;
  int verbose=0;

  pixcrd[0][0] = pix[0] - 0.5; // bottom left lon
  pixcrd[0][1] = pix[1] - 0.5; // bottom left lat

  pixcrd[1][0] = pix[0] - 0.5; // top left lon
  pixcrd[1][1] = pix[1] + 0.5; // top left lat

  pixcrd[2][0] = pix[0] + 0.5; // top right lon
  pixcrd[2][1] = pix[1] + 0.5; // top right lat

  pixcrd[3][0] = pix[0] + 0.5; // bottom right lon
  pixcrd[3][1] = pix[1] - 0.5; // bottom right lat

  /* project the pixels onto the sky */

	
  if ((status = wcsp2s(wcs, 4, 2, pixcrd[0], imgcrd[0], phi, theta, world[0],stat))) {
    fprintf(stderr, "\n\nwcsp2s ERROR %d: %s.\n", status,
	wcs_errmsg[status]);

    /* Invalid pixel coordinates. */
  }

  if (status == 0) {
    if (verbose) {
    printf("\n\nCorner world coordinates:\n"
	"            Pixel           World\n");
    }
    for (k = 0; k < 4; k++) {
      if (verbose) {
      printf("  (%5.1f,%6.1f) -> (%7.3f,%8.3f)",
	  pixcrd[k][0], pixcrd[k][1], 
	  world[k][0],  world[k][1]);
      }
      if (!image) {
      corners[k][0] = world[k][0];
      corners[k][1] = world[k][1];
      }
      else {
	corners[k][0] = imgcrd[k][0];
	corners[k][1] = imgcrd[k][1];
      }

      if (stat[k] && verbose) printf("  (BAD)");

      if (verbose) {			
	printf("\n");
      }
    }
  }
  return status;
}		
