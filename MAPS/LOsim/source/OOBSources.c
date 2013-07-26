// OOBSources.c: Calculate the contribution to the visibility of the out-of-beam sources

// Author:	Peter Sherwood	iti@world.std.com	(617) 244-0836

// 6/29/01	First light
// 7/11/01	Variable array sizes (v1.0)
/*
The visibility of the OOB source characterized by (xSource, ySource, sigmaX, sigmaY, p) is

V(u, v) = I0·4*pi^2*sigmaX^2*sigmaY^2*
          exp{-pi(sigmaX^2[-sin p·u + cos p·v]^2+ sigmaY^2[-cos p·u - sin p·v]^2)}*
          exp(-2*pi*i*xSource*[-sin p·u + cos p·v])*exp(-2*pi*iySource[-cos p·u - sin p·v])

          i=sqrt(-1)

The OOB sources are in source[nSources]
*/

#include "OOBSources.h" // also includes LargeFiles.h
#include "RefreshRow.h"
//#include "SourceArrays.h"
#include "Parameters.h"
#include "GlobalReferences.h"
#include "ScreenIO.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define DEBUGDISPLAY 0	// non-zero to enable saving debug display

void OOBSources(char visibilityFilename[])
{
#if 0 // code below needs to be fixed to read in OOB source file, not use array
	int iS,ki,kj,kDisp;
	long i,j,halfNX=visibilityXN/2,halfNY=visibilityYN/2,iCached;
	double dx,dy,du,dv, xSource,ySource,pAng,sigX,sigY, r11,r12,r21,r22, I0,sigX2,sigY2;
	double u,v,u1,v1,upp,vpp, V0,V1,Vreal,Vimag, au,av,cu,cv, argExp;
	struct Complex *visRow; // used as an array to hold a single row from the visibility file
#if DEBUGDISPLAY // disable saving debug display
	int zeroCount;
	double z,Vmag,debugDisplayv[600],debugDisplayReal[600],debugDisplayImag[600],debugDisplayMag[600]; // store for display
	FILE *realFile,*imagFile,*magFile;
 	realFile=fopen("OOBSourceReal.txt","w"); fprintf(realFile,"# OOB source visibilities, real part\n");
 	imagFile=fopen("OOBSourceImag.txt","w"); fprintf(imagFile,"# OOB source visibilities, imag part\n");
 	magFile=fopen("OOBSourceMag.txt","w"); fprintf(magFile,"# OOB source visibilities, magnitude\n");
#endif
	FileHandle visibilityFH;
	
	// Open the visibility file to which the OOB sources are to be cumulated
	if((visibilityFH=fopenL(visibilityFilename,"rb+"))<0)
	{
		printf("Cannot open the visibility file '%s'",visibilityFilename);
		return;
	}
	
	// Allocate arrays
	visRow=(struct Complex *)malloc(visibilityXN*sizeof(struct Complex)); // holds a single row from the visibility file

	dx=fieldXSize/((double)(brightnessMaxX-brightnessMinX)); dy=fieldYSize/((double)(brightnessMaxY-brightnessMinY));
	du=1.0/(visibilityXN*dx); dv=1.0/(visibilityYN*dy);

	for (iS=0; iS<nSources; iS++)
	{
		// Source parameters
		xSource=source[iS].xOffset; ySource=source[iS].yOffset;
		pAng=source[iS].positionAngle;
		r11=-sin(pAng); r12=cos(pAng); r21=-cos(pAng); r22=-sin(pAng); // the rotation matrix
		sigX=source[iS].majorAxis; sigY=source[iS].minorAxis;
		I0=source[iS].fluxDensity/(2.0L*PI*sigX*sigY);	// normalize so
														// integral{I0*exp(-0.5*[x^2/sigX^2+y^2/sigY^2])dxdy}=I
		// Combinations of source parameters
		sigX2=sigX*sigX; sigY2=sigY*sigY;
		V0=I0*4.0L*PI*PI*sigX2*sigY2; au=PI*sigX2; av=PI*sigY2; cu=2.0L*PI*xSource; cv=2.0L*PI*ySource;
		printf("\nStart source %d:",iS); fflush(stdout);
		iCached=-1; // nothing in rowCache
		for (i=0,ki=0; i<visibilityXN; i++,ki++)
		{
			// Calculate u from the position in the array
			if(i<halfNX) u=i*du;
			else u=(i-visibilityXN)*du;
			u1=r11*u; v1=r21*u; // one part of rotated u and v, respectively
			
			for (j=0,kj=0,kDisp=0; j<visibilityXN; j++,kj++)
			{
				// Calculate v from the position in the array
				if(j<halfNY) v=j*dv;
				else v=(j-visibilityYN)*dv;
				
				upp=u1+r12*v; vpp=v1+r22*v; // rotated u,v
				V1=V0*exp(-(au*upp*upp + av*vpp*vpp)); // the real factor
				if(V1!=0.0)
				{
					argExp=cu*upp + cv*vpp; // arg for the complex exponential
					Vreal=V1*cos(argExp); Vimag=-V1*sin(argExp); // V=V1*exp(-argExp*sqrt(-1)); use exp(-iz)=cos z -i*sin z
					
					// Add these to the big visulaization array
					// The row cache holds a single row (N complex values); at the moment it contains row iCached, or -1 if nothing
					if(i!=iCached) RefreshRow(visibilityFH,visRow,&iCached,i); // write any old row and get row i
					visRow[j].re+=Vreal; visRow[j].im+=Vimag;
#if DEBUGDISPLAY // if we're saving debug display
					if(ki==54 && kj==54)
	       			{	// Save point for debugging display if the row is one we're saving
	       				if(kDisp<600)
	       				{
	       					kj=0;
		       				Vmag=sqrt(Vreal*Vreal+Vimag*Vimag);
    	    				debugDisplayv[kDisp]=v; debugDisplayReal[kDisp]=Vreal; debugDisplayImag[kDisp]=Vimag; debugDisplayMag[kDisp]=Vmag;
        					kDisp++;
	       				}
	       			}
#endif
				}
			}
			if(ki==54)
			{
	   			printf(" %ld%%",100*i/visibilityXN); fflush(stdout);
	   			ki=0;
#if DEBUGDISPLAY
	   			// Write a 600 points in ASCII format for gnuplot, suppressing more than 10 0s in a row
				for (kDisp=0,zeroCount=0; kDisp<600; kDisp++)
				{
					if((z=debugDisplayReal[kDisp])==0.0) zeroCount++;
					else zeroCount=0;
					if(zeroCount<10) fprintf(realFile,"%e %e %e\n",u,debugDisplayv[kDisp],z);
				}
				for (kDisp=0,zeroCount=0; kDisp<600; kDisp++)
				{
					if((z=debugDisplayImag[kDisp])==0.0) zeroCount++;
					else zeroCount=0;
					if(zeroCount<10) fprintf(imagFile,"%e %e %e\n",u,debugDisplayv[kDisp],z);
				}
				for (kDisp=0,zeroCount=0; kDisp<600; kDisp++)
				{
					if((z=debugDisplayMag[kDisp])==0.0) zeroCount++;
					else zeroCount=0;
					if(zeroCount<10) fprintf(magFile,"%e %e %e\n",u,debugDisplayv[kDisp],z);
				}
#if 0
    			// Write a column of 600 points in ASCII format for gnuplot
//				fprintf(realFile,"# Column %ld\n",i); // don't use this, since comment lines are NOT ignored in plot files
				for (kDisp=0; kDisp<600; kDisp++) fprintf(realFile,"%e ",debugDisplayReal[kDisp]);
    			fprintf(realFile,"\n");
				for (kDisp=0; kDisp<600; kDisp++) fprintf(imagFile,"%e ",debugDisplayImag[kDisp]);
    			fprintf(imagFile,"\n");
				for (kDisp=0; kDisp<600; kDisp++) fprintf(magFile,"%e ",debugDisplayMag[kDisp]);
    			fprintf(magFile,"\n");
#endif
#endif // DEBUGDISPLAY
    		}
    	}
	}
	// Deallocate arrays
	free(visRow);
#endif // code to be revised
	
	return;
}


