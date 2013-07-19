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
#include "ComplexArithmetic.h"
#include "SourceListRead.h"
#include <math.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "uvgrid.h" // Added for header support, SSD 5/28/03
//#include <stdio.h>
//#include <time.h>

#define DEBUGDISPLAY 0	// non-zero to enable saving debug display

#define LIMITARG(arg,lo,hi) ((arg)<(lo)?(lo):((arg)>(hi)?(hi):(arg)))

#define ECM 20 // number of error histogram buckets
#define VCM 20 // number of value histogram buckets
struct Error
{
	long nTot;				// total # of examples (sum of n[])
	long n[VCM];			// # of examples by value
	long iMax[VCM];			// index of largest + error
	double eMax[VCM];		// the maximum error
	double vMax[VCM];		// the value at the maximum error
	long eCount[VCM][ECM];	// histogram of log(error)
};

// Global variables
double du,dv;
long halfNX,halfNY;
FILE *outFile=NULL; // dummy for ScreenIO
int nSources;
struct Source source[100];

// Local functions
void ExactUV(long i, long j, double *pu, double *pv);
void VCalc(double u, double v, struct Complex *pV);
void StatInit(struct Error *pE);
void Stat(double value, double error, struct Error *pE);
void StatPrint(struct Error *pE);

int main(void)
{
	char visibilityDataFilename[10+1+50+1+3],line[100+1],uSign[6+1],vSign[6+1]; int ln,er,k,nr; long long n2,filePosition;
	long i,j,i0,j0,j1,jx,n,hdrsize;
	double u,v,dx,dy,uGrid,vGrid,uMin,uMax,vMin,vMax,mag,phase;
	struct statL fileInfo; FileHandle visibilityDataFileHandle; struct Complex V,visibilityData[1<<15];	
	
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
	hdrsize = sizeof (struct uvgrid_header);
	n2=(fileInfo.st_sizeL - hdrsize)/sizeof(struct Complex); n=(int)(sqrt((double)n2)+0.5);
	if(n2*sizeof(struct Complex)!=(fileInfo.st_sizeL - hdrsize))
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
	{
		double fractError,magV,phaseV;
		int iS,whatToDo,eof; struct Error errorCR,errorCI,errorR,errorI;
		char filename[101];
		
		// Read in the entire source list (currently works only for 100 sources or fewer)
		CategoryFilename(filename,"SourceList",sourceListName,"txt"); // where the source list is found
		StatInit(&errorCR); StatInit(&errorCI); StatInit(&errorR); StatInit(&errorI);	
		
		// Loop on sources, terminating when EOF is encountered
		for (iS=0; ; iS++)
		{
			whatToDo=iS==0?SLR_FIRST:SLR_NEXT; SourceListRead(whatToDo,filename,&source[iS],&eof); if(eof) break; // get a source, if any
		}
		nSources=iS;
				
#if 0 // following calculates theoretical values for every point and creates error histograms
		for (i=0; i<visibilityXN; i++)
		{	
			// Get a full row of LOsim-calculated visibilities
			nr=freadL(visibilityData,sizeof(struct Complex),visibilityYN,visibilityDataFileHandle);
			for (j=0; j<visibilityYN; j++)
			{
				ExactUV(i,j,&u,&v);
				VCalc(u,v,&V);
				
				// Tabulate the LOsim approximation error
				fractError=fabs((V.re-visibilityData[j].re)/V.re); Stat(V.re,fractError,&errorR);
				fractError=fabs((V.im-visibilityData[j].im)/V.im); Stat(V.im,fractError,&errorI);
					
	#if 0 // unused code
				mag=sqrt(visibilityData[j].re*visibilityData[j].re+visibilityData[j].im*visibilityData[j].im);
				phase=atan2(visibilityData[j].im,visibilityData[j].re)*(360.0L/(2.0L*PI));
	#endif
			}
		}
		StatPrint(&errorCR); StatPrint(&errorCI); StatPrint(&errorR); StatPrint(&errorI);
#endif

	}
	
	while(1)
	{
		printf("\nu,v: "); fgets(line,sizeof(line),stdin);
		ln=strlen(line); if(ln>=1 && line[ln-1]=='\n') line[ln-1]=0; // remove the NL if we read a complete line (normal)
		n=sscanf(line,"%le,%le",&u,&v);
		if(n<2) return 0;
		u=LIMITARG(u,uMin,uMax); v=LIMITARG(v,vMin,vMax); // keep within the array
		
		// Calculate the position in the array
		if(u<0) i=visibilityXN-((int)(-u/du+0.5));
		else i=(int)(u/du+0.5);
		if(v<0) j=visibilityYN-((int)(-v/dv+0.5));
		else j=(int)(v/dv+0.5);
		// Now go back and calculate the exact (u, v) for array point (i, j)
		ExactUV(i,j,&uGrid,&vGrid);
		if(i!=halfNX) strcpy(uSign,""); else strcpy(uSign,"+-");
		if(j!=halfNY) strcpy(vSign,""); else strcpy(vSign,"+-");
		printf("  Nearest grid point is (%ld, %ld) which is at (u,v)=(%s%.8le,%s%.8le)",i,j,uSign,uGrid,vSign,vGrid);
		
		// Read (up to) 3 rows before through 3 rows after, 3+1+3 elements from each row
		// BUG: won't work if we're at halfN boundary
		j0=LIMITARG(j-1,0,visibilityYN); j1=LIMITARG(j+1,0,visibilityYN); // j0,1=starting/ending column to read
		printf("\n"); for (jx=j0; jx<=j1; jx++) printf("%34ld",jx);
		printf("\n"); for (jx=j0; jx<=j1; jx++) printf(" %16s,%16s","real","imag");
		printf("\n"); for (jx=j0; jx<=j1; jx++) printf(" %16s,%16s","magnitude","phase,degrees");
		for (k=-1; k<=1; k++)
		{
			i0=i+k; if(i0<0 || i0>=visibilityXN) continue; // i0=valid row to read, or skip this row entirely
			filePosition=((long long)i0*visibilityYN+j0)*sizeof(struct Complex) + hdrsize;
			fseekL(visibilityDataFileHandle,filePosition,SEEK_SET);
			nr=freadL(visibilityData,sizeof(struct Complex),j1-j0+1,visibilityDataFileHandle);
			printf("\n%5ld",i0); for (jx=0; jx<(j1-j0+1); jx++) printf(" %16.9lg,%16.9lg",visibilityData[jx].re,visibilityData[jx].im);
			printf("\nTheor");   for (jx=0; jx<(j1-j0+1); jx++)
			{
				ExactUV(i0,j0+jx,&uGrid,&vGrid);
				VCalc(uGrid,vGrid,&V);
				printf(" %16.9lg,%16.9lg",V.re,V.im);
			}
			printf("\nMagPh");   for (jx=0; jx<(j1-j0+1); jx++)
			{
				mag=sqrt(visibilityData[jx].re*visibilityData[jx].re+visibilityData[jx].im*visibilityData[jx].im);
				phase=atan2(visibilityData[jx].im,visibilityData[jx].re)*(360.0L/(2.0L*PI));
				printf(" %16.9lg,%16.9lg",mag,phase);
			}
            printf("\n");
		}
	}
	
	return 0;
}

void ExactUV(long i, long j, double *pu, double *pv)
{
	if(i<=halfNX) *pu=i*du; else *pu=(i-visibilityXN)*du; // make the ambiguous center freq positive
	if(j<=halfNY) *pv=j*dv; else *pv=(j-visibilityYN)*dv;
	return;
}

// Theoretical visibility for a Gaussian source
void VCalc(double u, double v, struct Complex *pV)
{
	double p,sigX,sigY,r1,r2,t1,t2,xSource,ySource,twoPi=2*PI;
	double a,sinp,cosp; // so debugger can access
	struct Complex e1,e2,e3,s,t;
	int iS;
	
//	pV->re=0.0; pV->im=0.0; // sum visibility due to each source
//	for (iS=0; iS<nSources; iS++)
//	{
//		xSource=source[iS].xOffset; ySource=source[iS].yOffset;
//		p=source[iS].positionAngle;
//		sigX=source[iS].majorAxis; sigY=source[iS].minorAxis;
//		sinp=sin(p); cosp=cos(p);
//		r1=-sin(p)*u+cos(p)*v; r2=-cos(p)*u-sin(p)*v;
//		t1=sigX*r1; t2=sigY*r2;					
//		a=2*PI*PI*(t1*t1+t2*t2);
//		s.re=source[iS].fluxDensity*exp(-2*PI*PI*(t1*t1+t2*t2)); s.im=0.0;
//		CExpi(-twoPi*xSource*r1,&e1); CExpi(-twoPi*ySource*r2,&e2);
//		CMult(e1,e2,&e3); CMult(s,e3,&t); CAdd(t,*pV,pV);
//	}
//
    p = 2.0*PI*(u*20.0/206264.806247096355 + v*20.0/206264.806247096355);
    pV->re = cos(p); pV->im = sin(p);
	return;
}

// Error statistics
void StatInit(struct Error *pE)
{
	int i,j;
	for (i=0; i<VCM; i++)
	{
		pE->n[i]=0; pE->iMax[i]=-1; pE->eMax[i]=0.0;
		for (j=0; j<ECM; j++) pE->eCount[i][j]=0;
	}
	pE->nTot=0;
	return;
}

void Stat(double value, double error, struct Error *pE)
{
	double log10Err,log10Val; int j,k;
	log10Val=log10(fabs(value)); j=(int)(log10Val+(log10Val<0?-.5:.5)); // round to nearest integer
	if(j>0) j=0; if(j<-(VCM-1)) j=-(VCM-1);
	if(error>pE->eMax[-j])
	{
		pE->iMax[-j]=pE->nTot;
		pE->eMax[-j]=error;
	}
	log10Err=log10(error); k=(int)(log10Err+(log10Err<0?-.5:.5)); // round to nearest integer
	if(k>0) k=0; if(k<-(ECM-1)) k=-(ECM-1);
	pE->eCount[-j][-k]++;
	pE->n[-j]++; pE->nTot++;
	return;
}

void StatPrint(struct Error *pE)
{
	int i,j;
	printf("\n      log10 value"); for (i=0; i<VCM; i++) printf("%9d",-i); printf(" or less");
	printf("\nlog10 fract error");
	for(j=0; j<ECM; j++)
	{
		printf("\n%9d%8s",-j,j==(ECM-1)?" or less":"");
		for (i=0; i<VCM; i++)
		{
			printf("%9ld",pE->eCount[i][j]);
		}
	}
	return;
}

