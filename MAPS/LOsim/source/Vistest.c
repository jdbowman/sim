// VisibilityView.c: Displays values from a visibility file interactively

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 7/12/01	create
// 9/17/01	modify for centered brightness array
// 9/20/01	reverse the above modification and better bounds checking
// 11/30/01	go back to du,dv=.../N from N-1; add magnitude/phase calcs

/*
*/

#include "LargeFiles.h"
#include "SimFilename.h"
#include "Parameters.h"
#include "DataStructures.h"
#include "GlobalReferences.h"
#include "GlobalVariables.h"
#include <math.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
//#include <stdio.h>
//#include <time.h>

#define DEBUGDISPLAY 0	// non-zero to enable saving debug display

#define LIMITARG(arg,lo,hi) ((arg)<(lo)?(lo):((arg)>(hi)?(hi):(arg)))

int main(void)
{
	char visibilityDataFilename[10+1+50+1+3],line[100+1],uSign[6+1],vSign[6+1]; int b,ln,er,k,nr; long long n2,filePosition;
	long i,j,i0,j0,j1,jx,n,halfNX,halfNY;
	double u,v,dx,dy,du,dv,uGrid,vGrid,uMin,uMax,vMin,vMax,mag,phase,mageq,pheq,temp,sigma,lgmagdif,magrat,phdiff;
	struct statL fileInfo; FileHandle visibilityDataFileHandle; struct Complex visibilityData[7];	
	
	if(GlobalVariables(NULL)) exit(0);
	halfNX=visibilityXN/2; halfNY=visibilityYN/2;

	// See if the visibility file exists, and, if so, how big it is
	SimFilename(visibilityDataFilename,"Visibility","dat");
	if(statL(visibilityDataFilename,&fileInfo))
	{
		er=errno;
		if(er==ENOENT) printf("'\aFile '%s' does not exist\n",visibilityDataFilename);
		else if(er==EACCES) printf("'\aDon't have permission to access file '%s'\n",visibilityDataFilename);
		else if(er==ENAMETOOLONG) printf("'\aFilename '%s' is too long\n",visibilityDataFilename);
		else if(er==ENOTDIR) printf("'\aPart of the path '%s' is not a directory\n",visibilityDataFilename);
		else if(er==ELOOP) printf("'\aFile '%s' has too many symbolic links (possible loop)\n",visibilityDataFilename);
		else printf("'\aUnable to get information for file '%s'\n",visibilityDataFilename);
		exit(0);
	}
	// Get the size of the matrix represented by the file; it must be N x N x 16 bytes
	n2=fileInfo.st_sizeL/sizeof(struct Complex); n=(int)(sqrt((double)n2)+0.5);
	if(n2*sizeof(struct Complex)!=fileInfo.st_sizeL)
	{
		printf("File '%s' size is not a multiple of size of complex (%d)",visibilityDataFilename,sizeof(struct Complex));
		exit(0);
	}
	if(visibilityXN*visibilityYN!=n2)
	{
		printf("File '%s' size is not the correct size; did you change the description file?",visibilityDataFilename);
		exit(0);
	}
	
	// Open the visibility file, which will be read in small pieces after each interactive request
	visibilityDataFileHandle=fopenL(visibilityDataFilename,"rb");
	
	dx=fieldXSize/((double)(brightnessMaxX-brightnessMinX)); dy=fieldYSize/((double)(brightnessMaxY-brightnessMinY));
 	du=1.0/(visibilityXN*dx); dv=1.0/(visibilityYN*dy); uMax=halfNX*du; uMin=-halfNX*du; vMax=halfNY*dv; vMin=-halfNY*dv;
	printf("\nFile=%s    Visibility grid u-spacing=%le   v-spacing=%le",visibilityDataFilename,du,dv);
        printf("\nU range: %le  to  %le\n",uMin,uMax);
	printf("V range: %le  to  %le\n",vMin,vMax);
			
#if testPosition
	{
		int nChecks,iR; double expectedR,expectedI; struct Complex visData[16000];
		fprintf(stderr,"\nChecking random file elements:");
		nChecks=visibilityXN*visibilityYN/1000;
		for (iR=0; iR<nChecks; iR++)
		{
			i0=rand()%visibilityXN; j0=rand()%visibilityYN;
			filePosition=((long long)i0*visibilityYN+j0)*sizeof(struct Complex);
			fseekL(visibilityDataFileHandle,filePosition,SEEK_SET);
			nr=freadL(visData,sizeof(struct Complex),1000,visibilityDataFileHandle); // may go past EOF
			
			// Compare to the expected values
			for (k=0,i=i0,j=j0; k<1000; k++)
			{
				expectedR= (i*bBase+j);
				expectedI=-(i*bBase+j);
				if(visData[k].re!=expectedR || visData[k].im!=expectedI)
				{
					fprintf(stderr,
					  "\nVisibilityView position error: iR=%d k=%d i=%ld j=%ld actual=(%23.16lg,%23.16lg) expected=(%23.16lg,%23.16lg)",
					  iR,k,i,j,visData[k].re,visData[k].im,expectedR,expectedI);
				}
				j++; if(j>=visibilityYN) {i++; j=0;}
				if(i>=visibilityXN) break;
			}
		if(iR%10==0) fprintf(stderr," %d",iR);	
		}
	}
#endif

	sigma = 5.0/(sqrt(8*log(2)));
	printf("sigma=%18.16lg\n",sigma);
	printf("   i,   j u                    \tv                    \tmag(LOsim)           \t mag(theory)    \tlog rel rror(delmag/theory)\tphase diff \n");
	for (b=0;b< 1000; b++)
	{
		i=rand()%visibilityXN; j=rand()%visibilityYN;
		
		// Now go back and calculate the exact (u, v) for array point (i, j)
		if(i<=halfNX) uGrid=i*du; else uGrid=(i-visibilityXN)*du; // make the ambiguous center freq positive
		if(j<=halfNY) vGrid=j*dv; else vGrid=(j-visibilityYN)*dv;
		if(i!=halfNX) strcpy(uSign,""); else strcpy(uSign,"+-");
		if(j!=halfNY) strcpy(vSign,""); else strcpy(vSign,"+-");
	
	
		mageq = exp(-1.0*PI*PI*PI*PI*25.0 * 4.0*  (uGrid*uGrid + vGrid*vGrid)/(4*log(2)*360.0*360.0*3600.0*3600.0));
		pheq = (uGrid*2.0*PI/(360.0*3600.0) + vGrid*PI/(360.0*3600.0))*360.0;
		temp = fmod (pheq,360.0); /* printf("temp = %16.9lg\n",temp); */
		if (temp > 180.0) pheq=temp-360.0;
		if (temp < -180.0) pheq=temp+360.0;
		if (abs(temp)<180.0) pheq = temp;
	
			filePosition=((long long)i*visibilityYN+j)*sizeof(struct Complex);
			fseekL(visibilityDataFileHandle,filePosition,SEEK_SET);
			nr=freadL(visibilityData,sizeof(struct Complex),1,visibilityDataFileHandle);
				mag=sqrt(visibilityData[0].re*visibilityData[0].re+visibilityData[0].im*visibilityData[0].im);
				phase=atan2(visibilityData[0].im,visibilityData[0].re)*(360.0L/(2.0L*PI));
				magrat=log10(mag/mageq); if (fabs(magrat)<1.0e-30) magrat = 1.0e-17;
				temp = fabs(mag-mageq)/mageq; 
				if (fabs(temp) < 1.0e-30) temp = 1.0e-17; 
				lgmagdif = log10(temp);
				
				phdiff = phase-pheq;
				if (mageq<1.0e-17) continue; 
				printf("%4d,%4d %s%.16le\t%s%.16le\t%18.16lg\t%18.16lg\t%18.16lg\t%18.16lg\n",i,j,uSign,uGrid,vSign,vGrid,mag,mageq,lgmagdif,phdiff); 
				/* printf("%4d %4d %s%.16le\t%s%.16le\t%18.16lg\t%18.16lg\n",i,j,uSign,uGrid,vSign,vGrid,mag,mageq); */
		
	}
	
	return 0;
}


