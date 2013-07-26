
// ImageAdd.c: Add a user-supplied image to the 2-D brightness array

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 1/23/02	Create

#include "ImageAdd.h"
#include "FITSread.h"
#include "GlobalReferences.h"
#include "Parameters.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

int ImageAdd(struct ImageFileParams imageParams, double *pBrightness0, 
	     double dx, double dy) {
  int i,err,dataSize,nDim,dim[10],ix,iy,scaling,stretching;
  int n,nTimes, indegrees;
  char *typeName,bUnit[68+1],line[80-9+1+1];
  struct Axis axis[10];
  double bScale,bZero,*pImageD=NULL,r11,r12,r21,r22,v;
  double BMajS,BMinS,beamSize,BMaj,BMin;


  printf("ImageAdd: brightnessXN = %ld\tbrightnessYN = %ld, " \
	 "pointer: %p\n", brightnessXN, brightnessYN,pBrightness0);

  // Get the info from the FITS file header
  if((err = FITSread_primaryHeaderInfo(imageParams.filename,
				       &dataSize, &nDim, axis, 
				       &bScale,&bZero,bUnit)))
    {
      fprintf(stderr,"Error code %d in opening or " \
	      "interpreting FITS file %s",err,imageParams.filename);
      return 1;
    }
  // See if it's a file we can deal with
  for (i = nDim-1; i >= 0; --i) 
    if(axis[i].nPixels != 1) break; // find the real number of dimensions
  if((i+1) != 2) {
    fprintf(stderr,"Array dimensions=%d in FITS file %s; " \
	    "expected 2 dimensions",i+1,imageParams.filename);
    return 2;
  }
  if(!(dataSize==SINGLE || dataSize==DOUBLE))
    {
      fprintf(stderr,"In FITS file %s, array is not IEEE single " \
	      "or double precision format (BITPIX=%d)",
	      imageParams.filename, dataSize);
      return 3;
    }
	
  // Allocate a working array for the data 
  // and check to see that FITS image is same
  // size as brightness array
  for (i = 0; i < nDim; i++) 
    dim[i] = axis[i].nPixels; // copy for FITSread_array2D call
  if ((dim[0] != brightnessXN) || (dim[1] != brightnessYN)) {
    fprintf(stderr,"Brightness array size=%ldx%ld doesn't " \
	    "match FITS size=%dx%d\n",
	    brightnessXN, brightnessYN, dim[0], dim[1]);
    return 4;
  }
  if (dataSize == DOUBLE) typeName="IEEE double";
  else typeName = "IEEE float";
  pImageD = (double *) malloc(dim[0]*dim[1]*sizeof(double));
  if (pImageD == NULL) {
    fprintf(stderr,"FATAL: ImageAdd: malloc failed\n"); exit(1);
  }
  fprintf(stdout,"pImageD is: %p\n",pImageD);

  // Now read the file
  err = FITSread_primaryHeaderVerify(imageParams.filename, dataSize,
				     nDim,dim);
  err = FITSread_array2D(pImageD, dataSize, dim[0], dim[1]);
		
  // Find the beam size info as the last occurrence of BMAJS=...BMINS=...
  err = FITSfind_LastHeaderLine("HISTORY", "BMAJS", line, &nTimes);
  if (nTimes) {
    indegrees = 0;
    // these are the beam dimensions in arcsec
    sscanf(strstr(line,"BMAJS"), "BMAJS =%le BMINS =%le", &BMajS, &BMinS);
    beamSize = PI*(BMajS*oneArcsec)*(BMinS*oneArcsec); // in radians^2
  }
  else { // data in image file is assumed to be Jy/sr, not Jy/beam
    BMajS = BMinS = 0.0; beamSize=1.0;
  } 

  if (!nTimes) {
    printf ("\nLooking for BMIN and BMAJ in degrees\n");
    // Find the beam size info as the last occurrence of BMAJ=...BMIN=...
    err = FITSfind_LastHeaderLine("HISTORY", "BMAJ", line, &nTimes);
    if(nTimes) {
      indegrees = 1;
      // these are the beam dimensions in degrees
      sscanf(strstr(line,"BMAJ"),"BMAJ =%le BMIN =%le",&BMaj,&BMin); 
      beamSize=PI*(BMaj)*(BMin)*(PI/180.0)*(PI/180.0); // in radians^2
    }
    else { // data in image file is assumed to be Jy/sr, not Jy/beam
      BMaj=BMin=0.0; beamSize=1.0;} 
  }
	
  FITSread_endFile();
	
  // Debugging only: get the max, min of the data
  {
    int i,j,nx=dim[0],ny=dim[1],nZero=0,nNeg=0;
    double tD,dataMinD=FLT_MAX,dataMaxD=0;
		
    for (i=0; i<nx; i++)
      {
	for (j=0; j<ny; j++)
	  {
	    tD = pImageD[i+j*ny];
						
	    if (tD > dataMaxD) dataMaxD=tD;
	    if(tD<dataMinD) dataMinD=tD;
	    if(tD<0.0) nNeg++;
	    if(tD==0.0) nZero++;
	  }
      }
		
    printf("Data type is %s  Array is %d x %d  Min,max=%g,%g  "
	   "Zero=%d  Negative=%d\n", typeName, nx, ny, 
	   dataMinD, dataMaxD, nZero, nNeg);
    if (indegrees == 0) 
      printf("Beam is %g x %g arcsec  Beam area=%g steradians\n",
	     BMajS, BMinS, beamSize);
    else 
      printf("Beam is %g x %g degrees  Beam area=%g steradians\n", 
	     BMaj, BMin, beamSize);
    printf("Axis #  Name              Pixels    Ref pixel          " \
	   "Ref value       Pix spacing         Rotate\n");
    for(i=0; i<nDim; i++) 
      printf("%4d %-20s %5d %16.9le  %16.9le  %16.9le  %16.9le  \n",
	     i+1, axis[i].name, axis[i].nPixels, axis[i].refPix,
	     axis[i].refVal, axis[i].pixSpacing, axis[i].rotate);
    fflush(stdout);
  }

  // Convert the angles to radians
  for(i = 0; i < nDim; i++) {
    axis[i].refVal *= oneDegree; 
    axis[i].pixSpacing *= oneDegree; 
    axis[i].rotate *= oneDegree;
  }
	
  /* Relation between image pixel index and coordinate 
   * value in radians is
   *
   *  coord = refVal + (index - (refPix-1))*pixSpacing
   *  index = (coord - refVal)/pixSpacing + refPix-1
   *	
   *  taking account of refPix being 1-based, and C array being 0-based
   */
  // Decode the pixel values if necessary
  n = dim[0]*dim[1];
  if(bScale != 1.0 || bZero != 0.0) {
    for (i = 0; i < n; i++) pImageD[i]=pImageD[i]*bScale+bZero;
  }
	
  // Pixel threshholding, limiting, scaling and stretching (power law)
  // to speed up loop below:
  scaling = (imageParams.scale != 1.0); 
  stretching = (imageParams.stretch != 1.0); 
  for (i = 0; i < n; i++) {
    v = pImageD[i];
    if (v < imageParams.blank) {
      pImageD[i] = 0.0; // don't wast time scaling or stretching
      continue;
    }
    else if (v > imageParams.limit) v = imageParams.limit;
    if (scaling) v *= imageParams.scale;
    if (stretching) v = pow(v,imageParams.stretch);
    pImageD[i] = v;
  }
	
  // Convert from Jy/beam to Jy/steradian, the units 
  // the Gaussian source brightnesses are in
  for (i = 0; i < n; i++) pImageD[i] /= beamSize;
		
  r11 = r22 = 1.0; r12 = r21 = 0.0; // CHANGE if rotation implemented
	
  // Now add the result to the main brightness array

  for (ix = 0; ix < brightnessXN; ix++) {
    for (iy = 0; iy < brightnessYN; iy++) {
      // flip x axis 17Jul03 -SSD 
      pBrightness0[ix*brightnessYN+iy] += 
	pImageD[(brightnessXN-ix-1)*brightnessYN+iy];
    }
  }
	
  // Free the storage for the external image
  free(pImageD);
	
  printf("pointer: %p\n",pBrightness0);
	
  return 0;
}
