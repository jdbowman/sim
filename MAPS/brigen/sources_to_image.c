//#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "sky_source.h"
#include "brigen.h"

void sources_to_image(struct bg_param *bgparam, /* image parameters */
		      double *image             /* Image plane */) {

  double const PI = 3.1415926535897931;
  double const dtor = PI/180.0;      /* degrees-to-radians multiplier */
  double const astor = dtor/3600.0;  /* arcseconds-to-radians multiplier */
  struct Source *source = bgparam->source; /* Descriptor of a sky source */ 
  int nsrc = bgparam->nsrc;        /* Number of sky sources */ 
  int ncol = bgparam->ncol;        /* Number of image grid columns */ 
  int nrow = bgparam->nrow;        /* Number of image grid rows */ 
  double xsize = bgparam->xsize*astor;   /* Image RA extent in radians */ 
  double ysize = bgparam->ysize*astor;   /* Image DEC extent in radians */ 
  double *imptr = image; /* Points at a location within the image plane */
  int isrc, i, j;
  int ix, iy; 
  int iy0, iy1;
  int np = 11;  /* Pixel sub grid size? */
  int nph = np/2;
  double fnp = (float)np, fnp2 = pow(fnp,2);
  long ixMax, iyMax, ixNeg = -1, iyNeg = -1;
  long ncmin, ncmax; /* Column numbers wrt the center of image */ 
  long nrmin, nrmax; /* Row numbers wrt the center of image */ 
  double t1, t2, coeffA, coeffB, coeffC, k = 74.860; 
  double bx, b0, cx2, cx, c0, discr, y0, y1;
  double x, y, dx, dy, I0, xof, yof, pAng, sigX, sigY;
  double r11, r12, r21, r22, xt, yt, contrib;
  double sigXR, sigYR, xp, yp, xpp, xpp0, xpp1, ypp;
  double contribMax, contribNeg; 
  double contrib_total; 

  /* If ncol or nrow is even, the extra point is negative */
  ncmax =  (ncol - 1)/2; 
  ncmin = -(ncol - 1 - ncmax); 
  nrmax = (nrow - 1)/2; 
  nrmin = -(nrow - 1 - nrmax);
  /* memset(image, 0, ncol*nrow*sizeof(double)); // requires about 30 sec */

  dx = xsize/((double)(ncmax - ncmin));  /* Cell size along RA */
  dy = ysize/((double)(nrmax - nrmin));  /* Cell size along DEC */
  printf("Field size is %g x %g radians;\n", xsize, ysize);
  printf("Grid x: %ld to %ld  y: %ld to %ld;\n", ncmin, ncmax, nrmin, nrmax);
  printf("dx=%g rad,  dy=%g rad;\n", dx, dy);
  printf("xcenter = %g deg, ycenter = %g deg.\n\n", 
	 bgparam->xcenter, bgparam->ycenter);

  /* 
   * For each point (x,y) in the ncol x nrow grid, 
   * calculate its position (x",y") relative to the source center by
   * 1) Translate the coordinate system to make the origin and the 
   * source center coincide, using
   *
   *		x' = x - xof
   *		y' = y - yof
   *
   * 2) Rotate the coordinate system by angle a counter-clockwise, 
   * where a = p + pi/2 (p=position angle of source major
   * axis, measured counter-clockwise from +y-axis) using
   *
   *  x" = x' * cos a + y' * -sin a  = x' * -sin p + y' * -cos p
   *  y" = x' * sin a + y' *  cos a  = x' *  cos p + y' * -sin p
   *
   *	        [r11 r12]
   *  I.e., R =           , with r11=-sin p, r12=-cos p, 
   *            [r21 r22]   r21=cos p, r22=-sin p (R=rotation matrix, 
   *                        converting unrotated coordinates to 
   *                        rotated coordinates.
   */


  // Loop on sources
  for (isrc=0; isrc < nsrc; isrc++) {
    contrib_total = 0.0;
    /* Convert all angles in radians */
    xof = source[isrc].xOffset*astor;
    yof = source[isrc].yOffset*astor;
    sigX = source[isrc].majorAxis*astor; 
    sigY = source[isrc].minorAxis*astor;
    sigXR = 1.0/sigX; 
    sigYR = 1.0/sigY;
    pAng = source[isrc].positionAngle*dtor;
    /* Rotation matrix: */
    r11 = -sin(pAng); r12 = cos(pAng); 
    r21 = -cos(pAng); r22 = -sin(pAng);
    /* 
     * Normalize so \integral{I0*exp(-0.5*[x^2/sigX^2+y^2/sigY^2])dxdy} = I 
     * Compute the x-independent parts of the coefficients for the truncation
     * The quadratic coeffA*y^2 + coeffB*y + coeffC = 0 is solved for y 
     * to give the trucation ellipse at each x
     */
    I0 = source[isrc].fluxDensity/(2*PI*sigX*sigY); 
    t1 = r12*sigXR; 
    t2 = r22*sigYR; 
    coeffA = t1*t1 + t2*t2;
    bx = 2.0*(r11*r12*sigXR*sigXR + r21*r22*sigYR*sigYR); 
    // coeffB=bx*(x-xof)+b0;
    b0 = -(2.0*yof*coeffA); 
    // coeffC=cx2*(x-xof)^2 + cx*(x-xof) + c0
    t1 = r11*sigXR; 
    t2 = r21*sigYR; 
    cx2 = t1*t1 + t2*t2; 
    cx = -(bx*yof); 
    c0 = coeffA*(yof*yof) - k;
    contribMax = -1.23;
    // ixNeg, iyNeg are used to detect illegal values for brightness:
    ixNeg = ncmin - 1; 
    iyNeg = nrmin - 1; 

    for (ix = ncmin; ix <= ncmax; ix++) {
      x = ix*dx;
      xp = x - xof;  /* x' = x - x_offset   */
      xpp0 = xp*r11; /* x"_1 = -x'sin(pAng) */
      xpp1 = xp*r21; /* x"_2 = -x'cos(pAng) */

      // Calculate the intersection of x=x with the truncation ellipse
      coeffB = bx*xp + b0;
      coeffC = cx2*xp*xp + cx*xp + c0;
      discr = coeffB*coeffB - 4*coeffA*coeffC;
      /* if (discr >= 0) { --  -- use continue; instead of if-clause */
      if (discr < 0) continue; //== No contribution outside the ellipse ==>>>
      if(discr > 0) { /* Two real solutions to ay^2+by+c=0: y0 and y1 */
	y0 = (-coeffB - sqrt(discr))/(2.0*coeffA); 
	y1 = (-coeffB + sqrt(discr))/(2.0*coeffA);
      }
      else { /* i.e. if (discr == 0.0) exactly: multiple root */
	y0 = y1 = (-coeffB)/(2.0*coeffA); // discr=0 -- unlikely
      }

      // Find next further grid points
      /* We only calculate the Gaussian between the points iy0 and iy1 */
      iy0 = ((int)(y0/dy)) - 1; 
      if (iy0 < nrmin) iy0 = nrmin;
      iy1 = ((int)(y1/dy)) + 1; 
      if (iy1 > nrmax) iy1 = nrmax;

      /* imptr points to brightness[ix][iy0]; */
      imptr = image + (ix-ncmin)*nrow + (iy0-nrmin); 

      for (iy = iy0; iy <= iy1; iy++) { /* Only inside the trunc. ellipse */
	y = iy*dy;
	yp = y - yof;
	// Grid point indices (ix,iy) -> grid point (x,y)
	// Rotation, optimized by moving r*x out of iy-loop):
	xpp = xpp0 + yp*r12;  /* x" = -x'sin(pAng) + y'cos(pAng) */
	ypp = xpp1 + yp*r22;  /* y" = -x'cos(pAng) - y'sin(pAng) */
	xt = xpp*sigXR;       /* xt = x"/sigX */
	yt = ypp*sigYR;       /* yt = y"/sigY */
	contrib = 0;
	/* sub pixelise the pixel */
	for (i = 0; i < np; i++) {
	  for (j = 0; j < np; j++) {
	    double tx, ty; /* Temporary x and y */
	    tx = xt + dx*sigXR*(j - nph)/fnp;
	    ty = yt + dy*sigYR*(i - nph)/fnp;
	    // FAKE-always Gaussian: 
	    /* contrib += I0*exp(-0.5*(tx^2 + ty^2)) / np^2; */
	    contrib += I0*exp(-0.5*(tx*tx + ty*ty))/fnp2; 
	  }
	}
	/* printf("Contribution at x=%.10le, y=%10lg is %10lg.  "	\
	 *	 "xt=%10lg yt=%10lg xp=%10lg yp=%10lg\n",
	 *	 x, y, contrib, xt, yt, xp, yp); */

	/*
	 * image["x"]["y"] += contrib: 
	 */
	//*(imptr++) += contrib; 
	/* image[(ix-ncmin)*nrow + (iy-nrmin)] += contrib; */
	/* image[(iy-nrmin)*ncol + (ix-ncmin)] += contrib; */
	image[(iy-nrmin)*ncol + (ncmax-ix+1)] += contrib;

	contrib_total += contrib;

	if(contrib < 0) {
	  ixNeg = ix; iyNeg = iy; contribNeg = contrib;
	}
	if(contrib > contribMax) {
	  contribMax = contrib; ixMax = ix; iyMax = iy;
	}
      } /* for (iy = iy0; ... */ 
	/* } if (discr >= 0) { ... -- use continue; instead of if-clause */
    } /* for (ix = ncmin ... */
 
    if(ixNeg != (ncmin-1)) {
      printf("*** ERROR: Impossible brightness %15g at grid (%5ld,%5ld)",
	     contribNeg, ixNeg, iyNeg);
    }
    printf("Source %3d: I0 = %g; sigX =  %g; sigY = %g; total contrib = %g;\n", 
	   isrc+1, I0, sigX, sigY, contrib_total);
  } /* for (isrc = ... */ 
}

