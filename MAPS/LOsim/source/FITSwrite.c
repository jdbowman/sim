// FITSwrite.c: Utilities to write Flexible Image Transport System (FITS) data files

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 10/22/01	create
// 12/2/01	array write returns # elements written

/*
	Based on the FITS standard NOST 100-2.0 (March 29, 1999)
*/

#include "FITSwrite.h"
#include "DataStructures.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static int nHeaderRecs; // # 80-byte header lines written to file so far
static FILE *outputFile;

// Create a FITS file and initialize the FITS primary header
int FITSwrite_startHeader(char filename[])
{
	outputFile=fopen(filename,"w"); nHeaderRecs=0;
	FITSwrite_headerLine("SIMPLE",Logical,"T","File conforms to FITS standard");
	return 0;
}

// Send a header line to the FITS file
int FITSwrite_headerLine(char keyword[], enum valueType type, void *pValue, char comment[])
{
	char line[80],temp[50]; int len,i0;
	memset(line,' ',80);
	
	// Keyword
	len=strlen(keyword); if(len>8) len=8;
	strncpy(line,keyword,len);
	
	// Value Indicator and value
	if(type!=None)
	{
		line[9-1]='='; i0=31; // for all the fixed-format cases
		switch(type) // put the value into 'temp'
		{
			case(Logical):
							line[30-1]=*((char *)pValue); temp[0]=0;
							break;
			case(Integer):
							sprintf(temp,"%20d",*((int *)pValue));
							break;
			case(LongInteger):
							sprintf(temp,"%20ld",*((long *)pValue));
							break;
			case(Float):
							sprintf(temp,"%20.13G",*((float *)pValue));
							break;
			case(Double):
							sprintf(temp,"%20.13lG",*((double *)pValue));
							break;
			case(ComplexInt):
#if 0 // no ComplexInt structure defined yet
							sprintf(temp,"(%d,%d)",((struct ComplexInt *)pValue)->re,((struct ComplexInt *)pValue)->im);
							i0=0;
#endif
							break;
			case(Complex): // i e, complex double
							sprintf(temp,"(%lG,%lG)",((struct Complex *)pValue)->re,((struct Complex *)pValue)->im);
							i0=0;
							break;
			case(CharString):
							sprintf(temp,"'%s'",(char *)pValue);
							break;
			case(Text):
							strcpy(&line[9],(char *)pValue);
							i0=9+strlen((char *)pValue);
							break;
			case(None): // only here to prevent compiler warning
				break;
		}
 		if(type!=Text) strncpy(&line[11-1],temp,strlen(temp));
		if(!i0) i0=strcspn(&line[11-1],")")+1+1; // skip the ) and convert to 1-based column
 	}
	else i0=9;
	
	// Comment
	if(comment!=NULL)
	{
		len=strlen(comment); if(len>(80-(i0+2)+1)) len=80-(i0+2)+1;
		line[i0+1-1]='/';
		strncpy(&line[i0+2-1],comment,len);
	}
	fwrite(line,sizeof(char),80,outputFile);
	nHeaderRecs++;
	return 0;
}

// Write the FITS primary header, including END and padding to multiple of 2880 bytes
int FITSwrite_endHeader(void)
{
	FITSwrite_headerLine("END",None,NULL,NULL);
	while(nHeaderRecs%36) FITSwrite_headerLine("",None,NULL,NULL); // all-blank lines
	return 0;
}

// Output a C 2-dimensional array of doubles to the FITS file, reversing the indexing so it is in
// FORTRAN index order, as required by FITS. Also convert little- to big-endian.
// Returns the number of array elements written (should be nx*ny)
int FITSwrite_array2D(void *pArray, size_t nx, size_t ny,size_t dsize)
{
	void *pOut; int i,j,n,nCum=0; unsigned char c;
	union U
	{
		double d;
		float f;
		unsigned char c[8];
	};
	union U t;

	// Sanity check
	if (dsize != SINGLE && dsize != DOUBLE) {
	  fprintf(stderr,"FITSwrite_array2D: data size does not seem to be float or double. exiting.\n");
	  exit(1);
	}
	
	// Allocate memory for a single output row
	pOut=malloc(nx*abs(dsize)/8);
	if (pOut==NULL) {
	  fprintf(stderr,"FITSwrite_array2D: no malloc\n");
	  exit(1);
	}
	for (j=0; j<ny; j++)
	  {
	    for (i=0; i<nx; i++)
	      {
		if(dsize==DOUBLE)
		  {
		    // Convert little-endian doubles to big-endian
		    t.d=((double *)pArray)[i*nx+j];
		    if (BYTE_ORDER == LITTLE_ENDIAN) {
		      c=t.c[7]; t.c[7]=t.c[0]; t.c[0]=c; c=t.c[6]; t.c[6]=t.c[1]; t.c[1]=c;
		      c=t.c[5]; t.c[5]=t.c[2]; t.c[2]=c; c=t.c[4]; t.c[4]=t.c[3]; t.c[3]=c;
		    }
		    ((double *)pOut)[i]=t.d;
		  }
		if(dsize==SINGLE) {
		    t.f=((double *)pArray)[i*nx+j];
		    if (BYTE_ORDER == LITTLE_ENDIAN) {
		      c=t.c[3]; t.c[3]=t.c[0]; t.c[0]=c; c=t.c[2]; t.c[2]=t.c[1]; t.c[1]=c;
		    }
		    ((float *)pOut)[i]=t.f;
		}
	      }
	    n=ftell(outputFile); // for debugging only
	    n=fwrite(pOut,abs(dsize)/8,nx,outputFile); nCum+=n;
	  }
	
	free(pOut);
	return nCum;
}

// Pad the array to a multiple of 2880 bytes and close the file
int FITSwrite_endFile(void)
{
	long n; char zeros[2880-1]; int rem;
	n=ftell(outputFile);
	if(n%2880)
	{
		rem=2880-(n%2880);
		memset(zeros,0,rem);
		fwrite(zeros,1,rem,outputFile);
	}
	fclose(outputFile);
	return 0;
}



