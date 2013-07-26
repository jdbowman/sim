// RowTransformBrightness.c: FFT of the zero-padded brightness rows

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// ?/?/01	create
// 7/11/01	v1.0
// 9/14/01	decrease amount of progress output
// 9/18/01	v1.04: compansate for column spacing in padded brightness grid
// 9/25/01	v1.06: recorrect the above for retro change in v1.05
// 10/9/01	a row=const x, varying y, not the reverse
// 10/19/01	v1.08; shift negative x to xMax/2
// 10/23/01 v1.09: input brightness file is in FITS format
// 12/2/01	v1.10: position checking

// Read the brightness file as input, pad with 0s, and write the row-transformed
// result as BrFFT_<simulation name>.dat

#include "RowTransformBrightness.h"
#include "LargeFiles.h"
#include "DataStructures.h"
#include "GlobalReferences.h"
#include "Parameters.h" // also includes ArrayDimensions.h
#include "SimFilename.h"
#include "FITSread.h"
#include "FFTsubs.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
// #include "../fftw-2.1.3/fftw/fftw.h"
#include "DFT.h"


void RowTransformBrightness(char inputBrightnessFilename[], char *title)
{
	double *bReal,*bImag,*sine,*tangent,*pBrightness; struct Complex *rowTr; // these will be used as arrays
	clock_t start,finish,begin,end;
	double cps=(double)CLOCKS_PER_SEC,elapsedSec;
	double z,t,colSpacing; int i,j,k,kx,iy,nw,err,progresCounter,dim[10];
	int testingTranspose=0; // set testing=1 to do FFT with package FFTW and a DFT of one row
									  // set testingTranspose=1 to use CTB as the input
        // int testing=0;
	char filename[101];
	FILE *debugOutfile; FileHandle outputRowTransformedFile;
//	fftw_complex in[1<<15], out[1<<15], outDFT[1<<15]; // PARAMETERIZE max size
//	fftw_plan plan;
	int debugAnalytic=0; // 1 to calculate the analytic soln for a single 1 Jy circular Gaussian source at the origin
	double sigmaX,sigmaY,x,dx,dy,v,dv,f1,f2,f3,diff; struct Complex analSoln[1<<15]={{0,0}};
	long halfNY;

	// Variables used in reading brightness FITS file
	struct Axis axis[10];
	double bScale,bZero,BMajS,BMinS,beamSize;
	int dataSize, nDim, nTimes, n;
	char bUnit[68+1],line[80-9+1+1];
	
#if 0
	char inputBrightnessFilenameCheck[100]; FILE *inputBrightnessFileCheck; double *bCheck;
	int firstMismatch,nMismatch,lastMismatch,nc;
#endif
	
	// Allocate the arrays
	pBrightness=(double *)malloc(brightnessXN*brightnessYN*sizeof(double));
	bReal=(double *)malloc(visibilityYN*sizeof(double));
#if 0
	bCheck=(double *)malloc(visibilityYN*sizeof(double));
#endif
	bImag=(double *)malloc(visibilityYN*sizeof(double));
	sine=(double *)malloc((visibilityYN/4+1)*sizeof(double)); // N/4+1 values are computed in SinesAndTangents
	tangent=(double *)malloc((visibilityYN/4+1)*sizeof(double));
	rowTr=(struct Complex *)malloc(visibilityYN*sizeof(struct Complex));
	
	SinesAndTangents(visibilityYN,visibilityYLn2,sine,tangent,&z);
	
	// Get the input brightness array
	if(testingTranspose) {strcpy(inputBrightnessFilename,"CTBtest_brightness.fts"); dim[0]=dim[1]=1024;}
	
#if 0
	// For checking, use old .dat file
	strcpy(inputBrightnessFilenameCheck,inputBrightnessFilename); i=strlen(inputBrightnessFilenameCheck);
	strcpy(&inputBrightnessFilenameCheck[i-3],"dat"); inputBrightnessFileCheck=fopen(inputBrightnessFilenameCheck,"rb");
#endif


	// Get the info from the FITS file header
	printf ("Input Brightness FITS File: %s\n",inputBrightnessFilename);

	if((err=FITSread_primaryHeaderInfo(inputBrightnessFilename,&dataSize,&nDim,axis,&bScale,&bZero,bUnit)))
	{
		fprintf(stderr,"Error in opening or interpreting FITS file %s",inputBrightnessFilename);
		return; 
	}
	// See if it's a file we can deal with
	for (i=nDim-1; i>=0; --i) if(axis[i].nPixels!=1) break; // find the real number of dimensions
	printf ("Number of Axes in FITS file: %d   Number of Axes with >1 pixel: %d\n",nDim, i+1);

	if((i+1)!=2)
	{
		fprintf(stderr,"Array dimensions=%d in FITS file %s; expected 2 dimensions",i+1,inputBrightnessFilename);
		return;
	}
	if(!(dataSize==SINGLE || dataSize==DOUBLE))
	{
		fprintf(stderr,"In FITS file %s, array is not IEEE single or double precision format (BITPIX=%d)",inputBrightnessFilename,dataSize);
		return;
	}
	printf("Datasize in FITS file: %d\n", dataSize);

	for (i=0; i<nDim; i++) dim[i]=axis[i].nPixels; // copy for FITSread_array2D call
	printf("Size of Array in FITS file: %dx%d\n", dim[0],dim[1]);
	printf("Size we expect: %ldx%ld\n", brightnessXN, brightnessYN);

	dim[0]=brightnessXN; dim[1]=brightnessYN;  // What we expect
    
	// Now read the file and check that expected dimensions are found
	if((err=FITSread_primaryHeaderVerify(inputBrightnessFilename,dataSize,nDim,dim)))
	{
		fprintf(stderr,"Brightness file not in correct FITS format. Error %d from FITSreadHeaderVerify",err);
		return;
	}

	if ((err=FITSread_array2D(pBrightness,dataSize,dim[0],dim[1])))
	{
	  fprintf(stderr,"Brightness file not in correct FITS format. Error %d from FITSread_array",err);
	  return;
	}
	
	// Find the beam size info as the last occurrence of BMAJS=...BMINS=...
	err=FITSfind_LastHeaderLine("HISTORY", "BMAJS",line,&nTimes);
	if(nTimes)
	{
		sscanf(strstr(line,"BMAJS"),"BMAJS =%le BMINS =%le",&BMajS,&BMinS); // these are the beam dimensions in arcsec
		beamSize=PI*(BMajS*oneArcsec)*(BMinS*oneArcsec); // in radians^2
	}
	else {BMajS=BMinS=0.0; beamSize=1.0;} // data in image file is assumed to be Jy/sr, not Jy/beam

	printf("\nBeam Size in Steradians is: %e\n", beamSize);

	FITSread_endFile();

	// Convert from Jy/Beam to Jy/SR which is what LOsim expects.
	n = dim[0]*dim[1];
	//	for (i=0; i<n; i++)  {pBrightness[i]/=beamSize;
	  // if (pBrightness[i] > 1e-5) {printf("Element %d of pBrightness was %g, set to 1.0\n",i,pBrightness[i]);
	  // pBrightness[i] = 1.0;
	  // }
	//}

#if testPosition
	{
		// Verify that the position-checking input is as expected
		int i,j,iExp,jExp; double check;
		for (iExp=0; iExp<brightnessXN; iExp++)
		{
			for (jExp=0; jExp<brightnessYN; jExp++)
			{
				check=pBrightness[iExp*brightnessYN+jExp]; i=check/bBase; j=check-(i*bBase);
					if(i!=iExp || j!=jExp || check!=(i*bBase+j))
					{
						printf("\nRowTransformBrightness position error: check=%23.16lg iExp=%d jExp=%d i=%d j=%d",check,iExp,jExp,i,j);
					}
			}
		}
	}
#endif
	
	if(debugAnalytic) debugOutfile=fopen("RowTransformDebugOutput.txt","w");
	
	// Output row-transformed brightness file
	SimFilename(filename,"BrFFT","dat");
	outputRowTransformedFile=fopenL(filename,"wb");

	begin=clock(); progresCounter=0;
	for (i=0; i<brightnessXN; i++)
	{
		// Pad with 0s at beginning and end (NOTE: assumes 0.0 is represented
		// by 8 bytes of 0; true for glibc on Intel
		memset(&bReal[                             0],0,brightnessYNzero0*sizeof(double));
		memset(&bReal[visibilityYN-brightnessYNzero1],0,brightnessYNzero1*sizeof(double));
		// Real input (OPTIMIZE-do transform in half-size complex)
		memset(bImag,0,visibilityYN*sizeof(double));
		// Copy the real brightness data
		memcpy(&bReal[brightnessYNzero0],&pBrightness[i*brightnessYN],sizeof(double)*brightnessYN);
#if 0
		nc=fread(bCheck,sizeof(double),brightnessYN,inputBrightnessFileCheck);
		
		for (j=0,firstMismatch=-1,nMismatch=0; j<brightnessYN; j++)
		{
			if(bCheck[j]!=bReal[brightnessYNzero0+j])
			{
				nMismatch++;
				if(firstMismatch<0) firstMismatch=j;
				lastMismatch=j;
			}
		}
		if(nMismatch)
		{
			j=nMismatch;
		}
#endif
					
		// Compensate for the area of the grid square
		colSpacing=fieldYSize/(brightnessYN-1); // column spacing of visibility grid (padded brightness grid) in radians
		// Note: j<brightnessXNzero0 and j>=(brightnessXNzero0+brightnessXN) were zeroed above, and
		//       imaginary part=0, so these components do not need to be adjusted
	
       // printf("colSpacing in RowTransformBrightness = %g\n", colSpacing);

		for (j=brightnessYNzero0; j<(brightnessYNzero0+brightnessYN); j++) bReal[j]*=colSpacing;
		
		// Move -A to +A/2  (A=interval) so we can do FFT of 0-(N-1) rather than -N/2 - (N/2-1)
		// OPTIMIZE by combining with readin and adjustment above; also by not switching 0s
		for (j=0; j<visibilityYN/2; j++) {k=(j+(visibilityYN/2))%visibilityYN; t=bReal[j]; bReal[j]=bReal[k]; bReal[k]=t;}
		
//		if(testing)
//		  for (j=0; j<visibilityYN; j++) {in[j].re=bReal[j]; in[j].im=bImag[j];} // copy starting padded brightness to 'in' for compare test
		
		start=clock();
		BitReverse(bReal,bImag,visibilityYLn2);
		Butterfly(bReal,bImag,sine,tangent,visibilityYN,1);
		finish=clock();
		elapsedSec=((double)(finish-start))/cps;


/*		if(testing)
		{
			// Verify the FFT by doing the same row transformation using a different package, then comparing results
			plan = fftw_create_plan(visibilityYN, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_one(plan, in, out);
			fftw_destroy_plan(plan);

			// Check both packages by doing a DFT
			if(i==15 || 1)
			{
				DFT(visibilityXN,in,outDFT); // FAKE-depends on structures fftw_complex and Complex being the same
			}
		}
*/
					
		// Put real and imaginary parts together
		for (j=0; j<visibilityYN; j++) {rowTr[j].re=bReal[j]; rowTr[j].im=bImag[j];}
		
		if(debugAnalytic)
		{
			// Analytic soln for a single 1 Jy circular Gaussian source at the origin
			halfNY=visibilityYN/2;
			dx=fieldXSize/((double)(brightnessMaxX-brightnessMinX)); dy=fieldYSize/((double)(brightnessMaxY-brightnessMinY));
		 	dv=1.0/(visibilityYN*dy);
			kx=i-(visibilityXN-1-brightnessMaxX); x=kx*dx;
			sigmaX=0.21233*oneArcsec; sigmaY=0.21233*oneArcsec;
			for(iy=0; iy<visibilityYN; iy++)
			{
				// u is ordered with the zero frequency first, then increasing positive frequencies, etc.
				if(iy<=halfNY) v=iy*dv; else v=(iy-visibilityYN)*dv; // make the ambiguous center freq positive
				f1=((1.0L)/(sqrt(2.0L*PI)*sigmaX)); f2=exp(-2.0L*PI*PI*2*sigmaY*sigmaY*v*v); f3=exp(-0.5L*(x/sigmaX)*(x/sigmaX));
				analSoln[iy].re=((1.0L)/(sqrt(2.0L*PI)*sigmaX))*exp(-2.0L*PI*PI*sigmaY*sigmaY*v*v)*exp(-0.5L*(x/sigmaX)*(x/sigmaX));
				analSoln[iy].im=0.0;
			}
			if(abs(i-(visibilityXN/2))<20) // look at the center 40 rows
			{
				fprintf(debugOutfile,"\nRow %d  kx=%d  x=%le",i,kx,x);
				for (j=0; j<40; j++) // look at the first 20 and last 20 frequencies
				{
					iy=j<20?j:visibilityYN-(40-j);
					if(analSoln[iy].re!=0.0) diff=(rowTr[iy].re-analSoln[iy].re)/analSoln[iy].re; else diff=0.0;
					fprintf(debugOutfile,"\n%5d %24.17le %24.17le %24.17le %12.9f",iy,rowTr[iy].re,rowTr[iy].im,analSoln[iy].re,diff);
				}
			}
		}
		
#if testPosition
		// Replace actual data with dummy test data
		for (j=0; j<visibilityYN; j++) {rowTr[j].re=i*bBase+j; rowTr[j].im=-(i*bBase+j);}
#endif
		nw=fwriteL(rowTr,sizeof(rowTr[0]),visibilityYN,outputRowTransformedFile);
		if((--progresCounter)<0) // give feedback every brightnessYN/32 rows, starting with row 0
		{
			printf("\nDone with row%6d  Elements written=%d  Time for this row=%6.2f",i,nw,elapsedSec);
			progresCounter=brightnessXN/32;
		}
	}
	
	end=clock(); elapsedSec=((double)(end-begin))/cps;
	printf("\n* Total time=%.2f sec",elapsedSec);
	
	// Deallocate the arrays
	free(pBrightness); free(bReal); free(bImag); free(sine); free(tangent); free(rowTr);
	
	// Close the input and output files
//	fclose(inputBrightnessFile);
	if(debugAnalytic) fclose(debugOutfile);
	fcloseL(outputRowTransformedFile);
	
	return;
}
