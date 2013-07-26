/* polysamp
   A routine that implements the fast flux conserving resampling
   implementing the Sutherland-Hodgman clipping algorithm

 */

/* reference frame notes
   I'm going to attempt to increase speed over the "perfect" resampling
   algorithm applied in EqualArea by

   (1) Assume a 2d sky
   (2) Assume straight pixel vertices

   To this end the easiest frame to work in is the pixel frame

   By that I mean instead of working in angles we:

   (1) take the input pixel vertices (sky) and project into the 
   reference (pixel space).

   problem here is I do not want to perform too many WCS calls. 
   I have to perform 4 - per pixel. Just the same for everything 
else: but this can be stored from the ionospherc dedistortion step.
And we can perform this at the same time as the 

Determine a range of output pixels against which to 
check the input pixel. Need a fast check as to whether a point is
within a polygon.

(2) So we now have 2 sets of pixel vertices in the pixel 
space of the output frame

(3) We then clip the input polygon by the output pixel boundaries
assuming them to be straight,

(4) We must keep track of the fractional area of each input pixel
attributed to each otput pixel

 */   
#undef _PLOT

#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<math.h>
#include "wcs_utils.h"

#ifdef _PLOT
#include "cpgplot.h"
#endif 

#define true 1
#define false 0

// From the polySamp class
// The maximum number of vertices we will get if we start with
// a convex quadrilateral is 12, but we use larger
// arrays in case this routine is used is some other context.
// If we were good we'd be checking this when we add points in
// the clipping process.

// "local externals" for temporary vertices storage

// i'll remove these in the production version


double  rcX0[100];
double  rcX1[100];
double  rcY0[100];
double  rcY1[100];
     
/** Calculate the area of an arbitrary triangle.
 *  Use the vector formula
 *     A = 1/2 sqrt(X^2 Y^2 - (X-Y)^2)
 *  where X and Y are vectors describing two sides
 *  of the triangle.
 * 
 *  @param x0	x-coordinate of first vertex
 *  @param y0       y-coordinate of first vertex
 *  @param x1       x-coordinate of second vertex
 *  @param y1       y-coordinate of second vertex
 *  @param x2       x-coordinate of third vertex
 *  @param y2       y-coordinate of third vertex
 * 
 *  @return         Area of the triangle.
 */

double triangleArea(double x0, double y0, 
    double x1, double y1, 
    double x2, double y2) {

  // Convert vertices to vectors.
  double a = x0-x1;
  double b = y0-y1;
  double e = x0-x2;
  double f = y0-y2;

  double area=  (a*a+b*b)*(e*e+f*f) - (a*e+b*f)*(a*e+b*f);
  if (area <= 0) {
    return 0; // Roundoff presumably!
  } else {
    return sqrt(area)/2;
  }
}



/** Calculate the area of a convex polygon.
 * This function calculates the area of a convex polygon
 * by deconvolving the polygon into triangles and summing
 * the areas of the consituents.  The user provides the
 * coordinates of the vertices of the polygon in sequence
 * along the circumference (in either direction and starting
 * at any point).
 * 
 * Only distinct vertices should be given, i.e., the
 * first vertex should not be repeated at the end of the list.
 *
 * @param	n	The number of vertices in the polygon.
 * @param   x	The x coordinates of the vertices
 * @param   y	The y coordinates of teh vertices
 * @return		The area of the polygon.
 */
double convexArea(int n, double* x, double* y) {
	
	
  double area = 0;
  int i = 0;
  for(i=1; i<n-1; i += 1) {

    area += triangleArea(x[0],y[0], x[i], y[i], x[i+1], y[i+1]);
  }
	
  return area;
}
double Bourke(int nv, double *xval, double *yval) { 
  int i,j;
  double area = 0;
  if (nv < 3) 
    return(0.);


    for (i=0;i<nv;i++) {
      j = (i + 1) % nv;
      area += xval[i] * yval[j];
      area -= yval[i] * xval[j];
    }

    area /= 2;

    if(isnan(area))
      area = 0.;
    
    if (area < 0) {
      area = -1 * area;
    }

    return(area);

}

static int inPlane(double test, double divider, int direction) {

  // Note that since we always include
  // points on the dividing line as 'in'.  Not sure
  // if this is important though...

  if (direction) {
    if (test >= divider) {
      return 1;
    }
    else {
      return 0;
    }
  } 
  else {
    if (test <= divider) {
      return 1;
    }
    else {
      return 0;
    }
  }
}

static int lineClip(int n, double *x, double *y, double *nx, double *ny, double val, int dir) {

  int nout=0;
  int i=0;	
  // Need to handle first segment specially
  // since we don't want to duplicate vertices.

  int last = inPlane(x[n-1], val, dir);

  for (i=0; i < n; i += 1) {

    if (last) {

      if (inPlane(x[i], val, dir)) {
	// Both endpoints in, just add the new point
	nx[nout] = x[i];
	ny[nout] = y[i];
	nout    += 1;
      } else {
	double ycross;
	// Moved out of the clip region, add the point we moved out
	if (i == 0) {
	  ycross = y[n-1] + (y[0]-y[n-1])*(val-x[n-1])/(x[0]-x[n-1]);
	} else {
	  ycross = y[i-1] + (y[i]-y[i-1])*(val-x[i-1])/(x[i]-x[i-1]);
	}
	nx[nout] = val;
	ny[nout] = ycross;
	nout    += 1;
	last     = false;
      }

    } else {

      if (inPlane(x[i], val, dir)) {
	// Moved into the clip region.  Add the point
	// we moved in, and the end point.
	double ycross;
	if (i == 0) {
	  ycross = y[n-1] + (y[0]-y[n-1])*(val-x[n-1])/(x[i]-x[n-1]);
	} else {
	  ycross = y[i-1] + (y[i]-y[i-1])*(val-x[i-1])/(x[i]-x[i-1]);
	}
	nx[nout]  = val;
	ny[nout] = ycross;
	nout += 1;

	nx[nout] = x[i];
	ny[nout] = y[i];
	nout += 1;
	last     = true;

      } else {
	// Segment entirely clipped.
      }
    }
  }
  return nout;
}


/** Clip a polygon by a non-rotated rectangle.
 * 
 *  This uses a simplified version of the Sutherland-Hodgeman polygon
 *  clipping method.  We assume that the region to be clipped is
 *  convex.  This implies that we will not need to worry about
 *  the clipping breaking the input region into multiple
 *  disconnected areas.
 *    [Proof: Suppose the resulting region is not convex.  Then
 *     there is a line between two points in the region that
 *     crosses the boundary of the clipped region.  However the
 *     clipped boundaries are all lines from one of the two
 *     figures we are intersecting.  This would imply that
 *     this line crosses one of the boundaries in the original
 *     image.  Hence either the original polygon or the clipping
 *     region would need to be non-convex.]
 * 
 *  Private arrays are used for intermediate results to minimize
 *  allocation costs.
 * 
 *  @param n	Number of vertices in the polygon.
 *  @param x	X values of vertices
 *  @param y        Y values of vertices
 *  @param nx	X values of clipped polygon
 *  @param ny       Y values of clipped polygon
 * 
 *  @param          minX Minimum X-value
 *  @param		minY Minimum Y-value
 *  @param          maxX MAximum X-value
 *  @param          maxY Maximum Y-value
 * 
 *  @return		Number of vertices in clipped polygon.
 */

int rectClip(int n, double *x, double *y, double *nx, double *ny,
    double minX, double minY, double maxX, double maxY) {

  int nCurr;

  // lineClip is called four times, once for each constraint.
  // Note the inversion of order of the arguments when
  // clipping vertically.

  nCurr = lineClip(n, x, y, rcX0, rcY0, minX, true);

  if (nCurr > 0) {
    nCurr = lineClip(nCurr, rcX0, rcY0, rcX1, rcY1, maxX, false);

    if (nCurr > 0) {
      nCurr = lineClip(nCurr, rcY1, rcX1, rcY0, rcX0, minY, true);

      if (nCurr > 0) {
	nCurr = lineClip(nCurr, rcY0, rcX0, ny, nx, maxY, false);
      }
    }
  }

  // We don't need to worry that we might not have set the output arrays.
  // If nCurr == 0, then it doesn't matter that
  // we haven't set nx and ny.  And if it is then we've gone
  // all the way and they are already set.
  return nCurr;
}


