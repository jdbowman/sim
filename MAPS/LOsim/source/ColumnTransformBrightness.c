// ColumnTransformBrightness.c: Read the brightness file as input, pad with 0s, and write the row-transformed result
// Input filename: is in inputBrFFTfilename, normally <simulation name>_BrFFT.dat
// Output filename: <simulation name>_Visibility.dat (a "large file", actually stored as <simulation name>_Visibility.dat_nn, nn=00, 01, ...)

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 7/01		Create
// 9/18/01	v1.04: compansate for column spacing in padded brightness grid
// 9/25/01	v1.06: recorrect the above for retro change in v1.05
// 10/20/01	v1.08: put negative x above +x/2
// 12/2/01	v1.10: position testing; rotate imaginary component

#include "ColumnTransformBrightness.h"
#include "LargeFiles.h"
#include "Parameters.h" // also includes ArrayDimensions.h
#include "DataStructures.h"
#include "GlobalReferences.h"
#include "GetColumnBlock.h"
#include "SimFilename.h"
#include "FFTsubs.h"
#include "uvgrid.h" // Added for header support, CJL
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void ColumnTransformBrightness(char inputBrFFTfilename[])
{
	double *bReal,*bImag,*sine,*tangent; // these will be used as arrays
	clock_t start,finish,begin,end;
	double cps=(double)CLOCKS_PER_SEC,elapsedSec,elapsedHr;
	double zsc,rowSpacing; int nw; char filename[101];
	FileHandle inputRowTransformedFile,outputVisibilityFile;
	long i,j,j0,k,iPrint; struct Complex z;
	long long filePosition;
	struct Complex *s; // s->buffer of visibilityYN blocks
	int blocksPerRow,iP;
	double t;	
                                        // Added for header support, CJL
    struct uvgrid_header uvhdr;
    long hdrsize;
	
	// Allocate the arrays
	bReal=(double *)malloc(visibilityYN*sizeof(double));
	bImag=(double *)malloc(visibilityYN*sizeof(double));
	sine=(double *)malloc((visibilityYN/4+1)*sizeof(double));
	tangent=(double *)malloc((visibilityYN/4+1)*sizeof(double));
	s=(struct Complex *)malloc(visibilityYN*BLOCKSIZE); // the big memory allocation, to hold 1 block from each of visibilityYN rows
	
	SinesAndTangents(visibilityYN,visibilityYLn2,sine,tangent,&zsc);
	
	// Input row-transformed brightness file
	inputRowTransformedFile=fopenL(inputBrFFTfilename,"rb");
	printf("\nInput file: %s (%ld x %ld)",inputBrFFTfilename,brightnessXN,visibilityYN);
	
	// Output visibility file
	SimFilename(filename,"Visibility","dat");
	outputVisibilityFile=fopenL(filename,"wb");
	printf("   Output file: %s (%ld x %ld)  Process in %ld pieces\n\n",filename,visibilityXN,visibilityYN,visibilityXN/ELEMENTSPERBLOCK);
	
	begin=time(NULL);
	
	blocksPerRow=visibilityXN/ELEMENTSPERBLOCK;
	printf("PIECE    START COL-END COL                                            TIME REQD, SEC\n");
					
	// Read one block's worth of elements from each row, putting together to make columns
	for (j0=0,iP=0; j0<visibilityXN; j0+=ELEMENTSPERBLOCK,iP++)
	{
		/* Read blocks containing b[i][0], b[i][1], ... , b[i][visibilityYN-1] (b=row-transformed brightness data, complex)
		   Since each block contains ELEMENTSPERBLOCK b's, there are blocksPerRow=visibilityXN/ELEMENTSPERBLOCK to read
		   j0=first column in each block read  iP=the part (used for information only)
		*/
		printf("%5d%7ld=%04lx%7ld=%04lx",
		         iP,j0,j0,j0+ELEMENTSPERBLOCK-1,j0+ELEMENTSPERBLOCK-1);
		start=clock(); printf("   Read"); iPrint=0; fflush(stdout);
		// For this group of columns, j0-(j0+ELEMENTSPERBLOCK-1), read the corresponding block from each row, i
		for (i=0; i<visibilityYN; i++)
		{
			GetColumnBlock(inputRowTransformedFile,i,j0,&s[i*ELEMENTSPERBLOCK]);
			if(i==iPrint) {printf("."); iPrint+=(visibilityYN/8); fflush(stdout);}
		}
	
		// Copy each column into the vector (bReal, bImag) and transform, then copy back
		printf("   FFT"); iPrint=0; fflush(stdout);
 		for (j=0; j<ELEMENTSPERBLOCK; j++)
		{
		
#if testPosition
			{
				// Verify that we have column j0+j of the input data
				double expectedR,expectedI; int iRow;
				// Copy column j0+j, separating the real and imaginary parts (done again below, so can be disabled)
				for (i=0; i<visibilityXN; i++) {z=s[i*ELEMENTSPERBLOCK+j]; bReal[i]=z.re; bImag[i]=z.im;}
			
				for (iRow=0; iRow<visibilityXN; iRow++)
				{
					if(iRow<brightnessYNzero0 || iRow>=(visibilityXN-brightnessYNzero1)) {expectedR=0.0L; expectedI=0.0L;}
					else
					{
						expectedR= ((iRow-brightnessXNzero0)*bBase+(j0+j));
						expectedI=-((iRow-brightnessXNzero0)*bBase+(j0+j));
					}
					if(bReal[iRow]!=expectedR || bImag[iRow]!=expectedI)
					{
							printf("\nColumnTransformBrightness position error: iRow=%d j0=%ld j=%ld actual=(%23.16lg,%23.16lg) expected=(%23.16lg,%23.16lg)",
							  iRow,j0,j,bReal[iRow],bImag[iRow],expectedR,expectedI);
					}
				}
			}
#endif
			
			// Copy column j0+j, separating the real and imaginary parts, and compensate for the area of the grid square
			rowSpacing=fieldXSize/(brightnessXN-1); // row spacing of visibility grid in radians
            // printf("rowSpacing = %g\n",rowSpacing);
			for (i=0; i<visibilityXN; i++) {z=s[i*ELEMENTSPERBLOCK+j]; bReal[i]=z.re*rowSpacing; bImag[i]=z.im*rowSpacing;}
			
			// Move -A to +A/2  (A=interval) so we can do FFT of 0-(N-1) rather than -N/2 - (N/2-1)
			// OPTIMIZE by combining with readin and adjustment above; also by not switching 0s
			for (i=0; i<visibilityXN/2; i++)
			{
				k=(i+(visibilityXN/2))%visibilityXN;
				t=bReal[i]; bReal[i]=bReal[k]; bReal[k]=t;
				t=bImag[i]; bImag[i]=bImag[k]; bImag[k]=t;
			}
		
			// FFT the column
			BitReverse(bReal,bImag,visibilityYLn2);
			Butterfly(bReal,bImag,sine,tangent,visibilityYN,1);

#if testPosition
			// Replace the FFT result with test data again
			for (iRow=0; iRow<visibilityXN; iRow++)
			{
				expectedR= (iRow*bBase+(j0+j));
				expectedI=-(iRow*bBase+(j0+j));
				bReal[iRow]=expectedR; bImag[iRow]=expectedI;
			}
#endif
			
			// Copy back to blocks, combining the real and imaginary parts again
			for (i=0; i<visibilityYN; i++) {z.re=bReal[i]; z.im=bImag[i]; s[i*ELEMENTSPERBLOCK+j]=z;}
			if(j==iPrint) {printf("."); iPrint+=(ELEMENTSPERBLOCK/8); fflush(stdout);}
		}
 		
                                        // Write out the header information
                                        // (for now, just two integers)
                                        // Added Feb 20 2002 by CJL
         filePosition = 0;
         fseekL (outputVisibilityFile, filePosition, SEEK_SET);
         uvhdr.skygridsize = brightnessXN;
         uvhdr.uvgridsize = visibilityXN;
         hdrsize = sizeof (struct uvgrid_header);
         nw = fwriteL (&uvhdr, hdrsize, 1, outputVisibilityFile);
         if (nw != 1) printf ("Error writing header to output file");

 		// Write the blocks containing the column-transformed data
		 printf("   Write"); iPrint=0; fflush(stdout);
		 for (i=0; i<visibilityYN; i++)
		{
                                        // Adjust for header, CJL
			filePosition=((long long)i*visibilityXN+j0)*sizeof(struct Complex) + hdrsize;
			fseekL(outputVisibilityFile,filePosition,SEEK_SET);
			nw=fwriteL(&s[i*ELEMENTSPERBLOCK],sizeof(struct Complex),ELEMENTSPERBLOCK,outputVisibilityFile);
			if(nw!=ELEMENTSPERBLOCK) printf("** ERROR Wrote %d (should be %d **)",nw,ELEMENTSPERBLOCK);
			if(i==iPrint) {printf("."); iPrint+=(visibilityYN/8); fflush(stdout);}
		}
		finish=clock(); elapsedSec=((double)(finish-start))/cps;
		printf("   Done  %6.2f sec\n",elapsedSec); fflush(stdout);
	}

	// Deallocate the arrays
	free(bReal); free(bImag); free(sine); free(tangent); free(s);
	
	fcloseL(inputRowTransformedFile);
	fcloseL(outputVisibilityFile);
	
	end=time(NULL); elapsedSec=difftime(end,begin); elapsedHr=elapsedSec/3600.;
	printf("\n* %s  Total time=%.2f hr\n",ctime(&end),elapsedHr); fflush(stdout);
	return;
}
