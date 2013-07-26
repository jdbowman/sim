// SumSources.c: Convert a list of input sources to a 
// brightnessXN x brightnessYN brightness file
// Currently all sources are Gaussian

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 5/8/01	Disable non-functional truncation code
// 5/24/01	Add source truncation
// 5/30/01  Correct rotation matrix sign
// 7/11/01  Array sizes variable (v1.0)
// 8/22/01	Read sources from file, not array
// 9/27/01	Correct brightness array index calculation
// 10/23/01 v1.09: output brightness file is in FITS format
// 12/2/01  v1.10: position checking
// 1/11/02  v1.13: warn on 0 sources
// 2/22/02	v1.15: correct FITS brightness file sign error

#include "SumSources.h"
#include "SourceListRead.h"
#include "Parameters.h"
#include "ScreenIO.h"
#include "SimFilename.h"
#include "GlobalReferences.h"
//#include "FITSread.h" // for FITS value types
#include "FITSwrite.h"
#include "ImagesRead.h"
#include "ImageAdd.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// Debugging options used in this function:
// debugOption_noSourceTruncate set to 1 to force calculation without 
// truncating Gaussian at small values

void Square(int x0, int y0, int s, double pArray[], int ny, double v);

void SumSources(void)
{
  int debug, noWrite; // 'debug' non-0 makes all sources the test case
  int whatToDo, eof, iI, result = 0;
  struct Source source;
  struct ImageFileParams image[MAX_IMAGES];
  char filename[101];
  size_t brightnessSize = (size_t)brightnessXN*(size_t)brightnessYN*
    sizeof(double);
  double *pBrightness0 = NULL;
  int datasize = SINGLE; /* single precision float written to FITS*/
  // C++ allocation: 
  //double *pBrightness0 = new double [brightnessXN*brightnessYN];
  long iS; 
  int ix, iy; 
  char screenLine[200];
  double x, y, dx, dy, I0, xSource, ySource, pAng, sigX, sigY;
  double sigXR, sigYR, xp, yp, xpp, xpp0, xpp1, ypp;
  double r11, r12, r21, r22, xt, yt, contrib;
  double sumContrib; 
  long timeReq, timeStartSource;
  // These declarations are here because otherwise gdb doesn't 
  // see the symbols during debugging
  double t1, t2, coeffA, coeffB, coeffC, kk=74.860; 
  // FAKE=2ln2*(m+1), m=53
  double bx, b0, cx2, cx, c0, disc, y0, y1;
  int iy0, iy1, nImages;
  double *pBrightness;
  double contribMax, contribNeg; 
  long ixMax, iyMax, ixNeg = -1, iyNeg = -1;
  long n;
  int debugShowDisc = 0;
  // 1 to display discriminant calc at every 10th x:
  int i, j, pixel_sub_grid_size = 11; 
  FILE *fh = NULL;
	
  pBrightness0 = (double *) malloc(brightnessSize);
  if (pBrightness0 == NULL) {
    fprintf(stderr, "FATAL: no malloc in SumSources\n");
    exit(1);
  }
  memset(pBrightness0, 0, brightnessSize); // requires about 30 sec

  printf("fieldXSize = %g, fieldYSize = %g\n", fieldXSize, fieldYSize); 

  dx = fieldXSize/((double)(brightnessMaxX - brightnessMinX));
  dy = fieldYSize/((double)(brightnessMaxY - brightnessMinY));
  SkipLine();
  sprintf(screenLine,"Field size is %g x %g radians.   " \
	  "Grid x: %ld to %ld  y: %ld to %ld   dx=%g  dy=%g\n",
	  fieldXSize, fieldYSize, brightnessMinX, brightnessMaxX,
	  brightnessMinY, brightnessMaxY, dx, dy);
  OutputToScreen(screenLine, 0);

  // If there are any external images, get them now
  if(ImagesRead(&nImages,image))
    {
      fprintf(stderr,"Errors reading image files; exiting\n");
      return;
    }
  if (nImages>0) fprintf(stdout,"There are %d images to read\n",nImages);
	    
  for (iI = 0; iI < nImages; iI++) 
    ImageAdd(image[iI],pBrightness0,dx,dy);
	
  /* 
   * For each point (x,y) in the brightnessXN x brightnessYN grid, 
   * calculate its position (x",y") relative to the source center by
   * 1) Translate the coordinate system to make the origin and the 
   * source center coincide, using
   *
   *		x' = x - xSource
   *		y' = y - ySource
   *
   * 2) Rotate the coordinate system by angle a counter-clockwise, 
   * where a = p + pi/2 (p=position angle of source major
   * axis, measured counter-clockwise from +y-axis) using
   *
   *  x" = x' * cos a + y' * -sin a  = x' * -sin p + y' * -cos p
   *  y" = x' * sin a + y' *  cos a  = x' *  cos p + y' * -sin p
   *
   *	     [r11 r12]
   *  I e, R =           , with r11=-sin p, r12=-cos p, 
   *           [r21 r22]  r21=cos p, r22=-sin p (R=rotation matrix, 
   *                      converting unrotated coordinates to 
   *                      rotated coordinates.
   */
  sprintf(screenLine," Count  Source ID              I0          " \
	  "sigmaX          sigmaY    Tot contrib       Time       " \
	  "Max brightness at\n");
  OutputToScreen(screenLine,0);
  // Where the source list is found
  CategoryFilename(filename, "SourceList", sourceListName, "txt"); 
	
  // Loop on sources, terminating when EOF is encountered
  for (iS=0; ; iS++) {
      whatToDo = iS == 0 ? SLR_FIRST : SLR_NEXT;
      result = SourceListRead(whatToDo, filename, &source, &eof);
      if(eof || result != 0) break; // Get a source, if any
      timeStartSource = clock();
      sumContrib = 0.0;
      xSource = source.xOffset;
      ySource = source.yOffset;
      pAng = source.positionAngle;
      debug = 0;
      // if(debug) {xSource=2.0; ySource=1.0; pAng=30.*(2.*PI/360.);}
      r11 = -sin(pAng); r12 = cos(pAng); r21 = -cos(pAng); r22 = -sin(pAng);
      sigX = source.majorAxis; sigY = source.minorAxis;
      // if(debug) {sigX=3.0; sigY=5.0;}
      sigXR = 1.0/sigX; sigYR = 1.0/sigY;
      I0 = source.fluxDensity/(2*PI*sigX*sigY); // Normalize so
      //   integral{I0*exp(-0.5*[x^2/sigX^2+y^2/sigY^2])dxdy} = I  
      // Compute the x-independent parts of the coefficients for the truncation
      // The quadratic coeffA*y^2 + coeffB*y + coeffC = 0 is solved for y 
      // to give the trucation ellipse at each x
      t1 = r12*sigXR; t2 = r22*sigYR; coeffA = t1*t1+t2*t2;
      bx = 2.0*(r11*r12*sigXR*sigXR + r21*r22*sigYR*sigYR); 
      // coeffB=bx*(x-xSource)+b0;
      b0 = -(2.0*ySource*coeffA); 
      // coeffC=cx2*(x-xSource)^2 + cx*(x-xSource) + c0
      t1 = r11*sigXR; t2 = r21*sigYR; cx2 = t1*t1 + t2*t2; 
      cx = -(bx*ySource); c0 = coeffA*(ySource*ySource) - kk;
		
      sprintf(screenLine,"%6ld %10s %15g %15g %15g\n", iS+1, source.id,
	      I0, sigX, sigY);
      OutputToScreen(screenLine, 0);
      pBrightness = pBrightness0; // We process brightness in storage order
      contribMax = -1.23;
      // ixNeg, iyNeg are used to detect illegal values for brightness:
      ixNeg = brightnessMinX - 1; 
      iyNeg = brightnessMinY - 1; 

      //printf("brightnessXN = %ld\n",  brightnessXN);
      //printf("brightnessYN = %ld\n",  brightnessYN);
      //printf("brightnessMinX = %ld\n",  brightnessMinX);
      //printf("brightnessMaxX = %ld\n",  brightnessMaxX);

      if (debugShowDisc) 
	printf("\n\n   ix         x            disc         " \
	       "y0             y1      iy0  iy1");
      for (ix = brightnessMinX; ix <= brightnessMaxX; ix++) {
	  x = ix*dx;
	  xp = x - xSource; xpp0 = xp*r11; xpp1 = xp*r21;
	  // Calculate the intersection of x=x with the truncation ellipse
	  coeffB = bx*xp + b0;
	  coeffC = cx2*xp*xp + cx*xp + c0;
	  disc = coeffB*coeffB - 4*coeffA*coeffC;
	  if (debugShowDisc && (ix%10 == 0)) 
	    printf("\n%5d %13.4e %13.4e", ix, x, disc);
	  if ((disc>=0) || debugOption_noSourceTruncate) {
	      if(disc > 0) { // real solutions to ay^2+by+c=0
		y0 = (-coeffB - sqrt(disc))/(2.0*coeffA); 
		y1 = (-coeffB + sqrt(disc))/(2.0*coeffA);
	      }
	      else 
		y0 = y1 = (-coeffB)/(2.0*coeffA); // disc=0 -- unlikely

	      // Find next further grid points
	      iy0 = ((int)(y0/dy)) - 1; 
	      if (iy0 < brightnessMinY) iy0=brightnessMinY;
	      iy1 = ((int)(y1/dy))+1; 
	      if (iy1 > brightnessMaxY) iy1=brightnessMaxY;
	      if(debugOption_noSourceTruncate) {
		iy0 = brightnessMinY; iy1 = brightnessMaxY;
	      }
	      if (debugShowDisc && (ix%10 == 0)) 
		printf("%13.4e %13.4e %5d %5d", y0, y1, iy0, iy1);
	      // points to brightness[ix][iy0];
	      pBrightness = pBrightness0 + 
		(ix - brightnessMinX)*brightnessYN + (iy0 - brightnessMinY); 
	      for (iy = iy0; iy <= iy1; iy++) {
		  y = iy*dy;
		  yp = y - ySource;
		  // Grid point indices (ix,iy) -> grid point (x,y)
		  // Rotation, optimized by moving r*x out of iy-loop):
		  xpp = xpp0 + yp*r12;
		  ypp = xpp1 + yp*r22;	
		  xt = xpp*sigXR;
		  yt = ypp*sigYR;
		  contrib = 0;
		  /* sub pixelise the pixel */
		  for (i = 0; i < pixel_sub_grid_size; i++) {
		    for (j = 0; j < pixel_sub_grid_size; j++) {
		      double temp_x, temp_y;
		      temp_x = xt + dx*sigXR*(j - pixel_sub_grid_size/2)/
			(float)pixel_sub_grid_size;
		      temp_y = yt + dy*sigYR*(i - pixel_sub_grid_size/2)/
			(float)pixel_sub_grid_size;
		      // FAKE-always Gaussian:
		      contrib += I0*exp(-0.5*(temp_x*temp_x + temp_y*temp_y))/
			(pixel_sub_grid_size*pixel_sub_grid_size); 
		    }
		  }
		  /* printf("Contribution at x=%.10le, y=%10lg is %10lg.  " \
		   *	 "xt=%10lg yt=%10lg xp=%10lg yp=%10lg\n",
		   *	 x, y, contrib, xt, yt, xp, yp); */
		  // brightness["x"]["y"] += contrib:
		  *(pBrightness++) += contrib; 
		  //debugBrightness[iy-brightnessMinY]=contrib;
		  sumContrib += contrib;
		  if(contrib < 0) {
		      ixNeg = ix; iyNeg = iy; contribNeg = contrib;
		  }
		  if(contrib > contribMax) {
		      contribMax = contrib; ixMax = ix; iyMax = iy;
		  }
	      } // iy
	  } // disc>=0
      } // ix

      timeReq = clock() - timeStartSource;
      if (sumContrib > 0.0) 
	sprintf(screenLine, "Total: %15g, Time: %10f, Max: %15g, " \
		"Range:(%5ld,%5ld)", sumContrib,
		((float)timeReq)/CLOCKS_PER_SEC, contribMax, ixMax, iyMax);
      else 
	sprintf(screenLine, "%15g ** intensity below noise at all " \
		"x-grid lines", sumContrib);
      OutputToScreen(screenLine, 1);
      if(ixNeg != (brightnessMinX-1)) {
	sprintf(screenLine, 
		"*** ERROR: Impossible brightness %15g at grid (%5ld,%5ld)",
		contribNeg, ixNeg, iyNeg);
	  OutputToScreen(screenLine, 1);
      }
      printf("\nSource %d done.\n", iS+1);
  } /* for (iS=0; ; iS++) { */
 	
  // Write resulting brightness array to file in FITS format
  SimFilename(filename, "Brightness", "fts");
  noWrite = 0; // 1 for debugging
	
  if(!iS) {
      fprintf(stderr, "\n*** WARNING: No sources in source list input file");
      // don't bother to write an all-0 brightness file, 
      // but do write if external images were added
      noWrite = !nImages; 
  }

  /* ASCII output of the brightness image */
  int ncol = brightnessXN, nrow = brightnessYN;
  fh = fopen("brimout.txt", "w");
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      fprintf(fh, "%g ", pBrightness0[i*ncol+j]);
    }
    fprintf(fh, "\n");
  }
  fclose(fh);


  if(!noWrite) {
      // FITS header
      int iTemp; 
      double dTemp;
      FITSwrite_startHeader(filename); // includes SIMPLE
      iTemp=datasize;
      FITSwrite_headerLine("BITPIX", Integer, &iTemp,
			   " IEEE 32 or 64 bit floating point");
      iTemp = 2;
      FITSwrite_headerLine("NAXIS", Integer, &iTemp, " 2-D array");
      FITSwrite_headerLine("NAXIS1", Integer, &brightnessXN, " x");
      FITSwrite_headerLine("NAXIS2", Integer, &brightnessYN, " y");
      dTemp = 1.0;
      FITSwrite_headerLine("BSCALE", Double, &dTemp, "");
      dTemp = 0.0;
      FITSwrite_headerLine("BZERO", Double, &dTemp, "");
      //FITSwrite_headerLine("BUNIT",CharString," ",""); // Changed by benkev:
      FITSwrite_headerLine("BUNIT", CharString, "Jy/steradian", "");
      FITSwrite_headerLine("CTYPE1", CharString, "RA---SIN", "");
      dTemp = -brightnessMinX + 1;	
      // FITS numbers pixels from 1, not 0:
      FITSwrite_headerLine("CRPIX1", Double, &dTemp, ""); 
      dTemp = fieldRA/oneDegree;
      // CRVAL1 in degrees, not radians
      FITSwrite_headerLine("CRVAL1", Double, &dTemp, ""); 
      dTemp = dx/oneDegree;
      // No sign change, it seems to be corrected in visgen
      FITSwrite_headerLine("CDELT1", Double, &dTemp,""); 
      dTemp = 0.0;
      FITSwrite_headerLine("CROTA1", Double, &dTemp, "");
      FITSwrite_headerLine("CTYPE2", CharString, "DEC--SIN", "");
      dTemp = -brightnessMinY + 1;
      FITSwrite_headerLine("CRPIX2", Double, &dTemp, "");
      dTemp = fieldDec/oneDegree;	
      FITSwrite_headerLine("CRVAL2", Double, &dTemp, "");
      dTemp = dy/oneDegree;
      FITSwrite_headerLine("CDELT2", Double, &dTemp, "");
      dTemp = 0.0;
      FITSwrite_headerLine("CROTA2", Double, &dTemp, "");
      FITSwrite_endHeader();
      n=FITSwrite_array2D(pBrightness0, brightnessXN, brightnessYN, datasize);
      FITSwrite_endFile();
      sprintf(screenLine,"Brightness file %s written: %ld doubles", 
	      filename, n);
      OutputToScreen(screenLine, 0);
    }
  else printf("\nNo brightness file written");
	
  free(pBrightness0);
  return;
}
		
// Functions which generate test images
void Square(int x0, int y0, int s, double pArray[], int ny, double v) {
  int s2 = s/2, ix, iy;
  for (ix = x0-s2; ix < (x0+s2); ix++) 
    for (iy = y0-s2; iy < (y0+s2); iy++) 
      pArray[ix*ny+iy]=v;
  return;
}
