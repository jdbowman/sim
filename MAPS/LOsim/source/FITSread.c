// FITSread.c: Utilities to read Flexible Image Transport System (FITS) data files

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 10/28/01	create
// 1/22/02	read single-precision arrays, convert endian-ness

/*
	Based on the FITS standard NOST 100-2.0 (March 29, 1999)
*/

#include "FITSread.h"
#include "DataStructures.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int IsAllBlanks(char line[], int len);

static int nHeaderRecs; // # 80-byte header lines read from file so far
static FILE *inputFile=NULL;

// Get the info from the FITS file header
int FITSread_primaryHeaderInfo(char filename[],
  int *pDataSize, int *pNDim, struct Axis axis[], double *pBScale, double *pBZero, char bUnit[])
{
	char keyword[8+1],strValue[68+1]; int err,intValue,i; double dblValue;
	struct Axis axisDefault={"???",0,0,0,0};
	if((err=FITSread_startHeader(filename))) return 1;
	if(FITSread_headerLine("BITPIX",Integer,&intValue,NULL,1)) return 3;
	*pDataSize=intValue;
	if(FITSread_headerLine("NAXIS",Integer,&intValue,NULL,1)) return 5;
	*pNDim=intValue;
	for (i=0; i<*pNDim; i++)
	{
		axis[i]=axisDefault; // in case one or more of the optional values is not present
		sprintf(keyword,"NAXIS%d",i+1);
		if(FITSread_headerLine(keyword,Integer,&intValue,NULL,1)) return 7000+i;
		axis[i].nPixels=intValue;
	}
	
	// Optional keywords (can appear in any order)
	*pBScale=1.0; *pBZero=0.0; strcpy(bUnit,"???"); // defaults in case one or more of these values is not present
	if(!FITSfind_headerLine("BSCALE",Double,&dblValue,NULL)) *pBScale=dblValue;
	if(!FITSfind_headerLine("BZERO",Double,&dblValue,NULL)) *pBZero=dblValue;
	if(!FITSfind_headerLine("BUNIT",CharString,strValue,NULL)) strcpy(bUnit,strValue);
	
	for (i=0; i<*pNDim; i++)
	{
		sprintf(keyword,"CTYPE%d",i+1);
		if(!FITSfind_headerLine(keyword,CharString,strValue,NULL)) strcpy(axis[i].name,strValue);
		sprintf(keyword,"CRPIX%d",i+1);
		if(!FITSfind_headerLine(keyword,Double,&dblValue,NULL)) axis[i].refPix=dblValue;
		sprintf(keyword,"CRVAL%d",i+1);
		if(!FITSfind_headerLine(keyword,Double,&dblValue,NULL)) axis[i].refVal=dblValue;
		sprintf(keyword,"CDELT%d",i+1);
		if(!FITSfind_headerLine(keyword,Double,&dblValue,NULL)) axis[i].pixSpacing=dblValue;
		sprintf(keyword,"CROTA%d",i+1);
		if(!FITSfind_headerLine(keyword,Double,&dblValue,NULL)) axis[i].rotate=dblValue;
	}
	fclose(inputFile);
	inputFile=NULL;
	return 0;
}

// Open a FITS file and verify the entire primary header
// Returns 0 if no error, else file did not match specs
int FITSread_primaryHeaderVerify(char filename[], int dataSize, int nDim, int dim[])
{
	char keyword[8+1]; int err,intValue,i;
	if((err=FITSread_startHeader(filename))) return 1;
	if(FITSread_headerLine("BITPIX",Integer,&intValue,NULL,1)) return 3;
	if(intValue!=dataSize) return 4;
	if(FITSread_headerLine("NAXIS",Integer,&intValue,NULL,1)) return 5;
	if(intValue!=nDim) return 6;
	for (i=0; i<nDim; i++)
	{
		sprintf(keyword,"NAXIS%d",i+1);
		if(FITSread_headerLine(keyword,Integer,&intValue,NULL,1)) return 7000+i;
		if(intValue!=dim[i]) return 8000+i;
	}
	if(FITSread_endHeader()) return 9;
	return 0;
}

// Open a FITS file and get the FITS primary header
// Returns 0=file successfully opened and identified as FITS format
//         1=file not found, or not able to open
//         2=not a FITS file
int FITSread_startHeader(char filename[])
{
	char comment[80],logicalValue;
	if((inputFile=fopen(filename,"rb"))==NULL) return 1;
	nHeaderRecs=0;
	if(FITSread_headerLine("SIMPLE",Logical,&logicalValue,comment,1)) return 2;
	if(logicalValue!='T') return 2;
	return 0;
}

// Get the last occurrence of a keyword line containing specified text
int FITSfind_LastHeaderLine(char keyword[], char text[], char lastLine[], int *nTimes)
{
	int first; char line[80-9+1+1];
	first=!FITSfind_headerLine(keyword,Text,line,NULL); *nTimes=0;
	if(first)
	{
		while(first || !FITSread_headerLine(keyword,Text,line,NULL,0)) // read_header is not called the first time
		{
			first=0;
			if(strstr(line,text)!=NULL)
			{
				strcpy(lastLine,line);
				(*nTimes)++;
				fprintf(stderr,"\n%3d %s",*nTimes,line);
			}
		}
	}
	return (*nTimes)==0;
}

// Get the specified header line from the FITS file
// If the input keyword is null string, an all-blank line is expected
// ...find... scans the entire header, while ...read... starts at the current position
// If not found, *pValue and comment are unchanged
int FITSfind_headerLine(char keyword[], enum valueType type, void *pValue, char comment[])
{
	fseek(inputFile,0L,SEEK_SET); // rewind the file
	return FITSread_headerLine(keyword,type,pValue,comment,0);
}

int FITSread_headerLine(char keyword[], enum valueType type, void *pValue, char comment[], int mustBeNext)
{
	char line[80]; int i,j,len,i0,eof,n,noMatch;
	
	// Keyword
	len=strlen(keyword); if(len>8) len=8;
	// First, skip any all-blank lines, except in the special case when keyword is null, we're expecting an all-blank line
	do {n=fread(line,sizeof(char),80,inputFile); nHeaderRecs++;} while(!(eof=(feof(inputFile) || n!=80)) && IsAllBlanks(line,80) && len);
	if(eof) return 1; // EOF without finding this keyword [next if 'mustBeNext']

		
	// Then, either examine the next line, or read all lines looking for keyword
	while((noMatch=strncmp(line,keyword,len)) && !mustBeNext && !(eof=(feof(inputFile) || n!=80)))
	  {n=fread(line,sizeof(char),80,inputFile); nHeaderRecs++;}
	if(eof) return 1; // EOF without finding this keyword
	if(noMatch) return 4; // next keyword didn't match
	
	if(strncmp("        ",&line[len],8-len)) return 3; // characters following keyword were non-blank
			
	// Value Indicator and value
	if(type!=None)
	{
		if(type!=Text && line[9-1]!='=') return 2; // value indicator is missing
		i0=31; // for all the fixed-format cases
		switch(type) // get the value
		{
			case(Logical):
							*((char *)pValue)=line[30-1];
							break;
			case(Text):
							// Remove the leading and trailing blanks before copying
							i=strspn(&line[9-1]," ")+(9-1); // i->1st non-blank char after column 9
							j=79; while (line[j]==' ' && j>=(9-1)) --j; // j->last non-blank char, or (9-1)-1=7 if none
							strncpy((char *)pValue,&line[i],j-i+1); ((char *)pValue)[j-i+1]=0;
							i0=79; // start position for comment, which will always be null
							break;
			case(CharString):
							// A legal character string value is formatted as <blanks>'<text>'<blanks>[/comment]
							// <text> may contain '' (meaning ' in string)
							// Locate the open '
							for (i=11-1; i<80; i++) if(line[i]!=' ') break; // skip any number of blanks
							for (; i<80; i++) if(line[i]=='\'') break;
							if(i>=80) return 3; // missing string open '
							// Locate the close '
							for (j=i+1; j<80; j++)
							{
								if(line[j]=='\'')
								{
									if(j>=80 || line[j+1]!='\'') break;
									else j++; // skip the second '
								}
							}
							if(j>=80) return 4; // missing string close '
							strncpy((char *)pValue,&line[i+1],j-(i+1)); ((char *)pValue)[j-(i+1)]=0;
							i0=-(j+1); if(i0>=80) i0=79; // start position for comment
							break;
			case(Integer):
							sscanf(&line[11-1],"%d",(int *)pValue);
							break;
			case(LongInteger):
							sscanf(&line[11-1],"%ld",(long *)pValue);
							break;
			case(Float):
							sscanf(&line[11-1],"%G",(float *)pValue);
							break;
			case(Double):
							sscanf(&line[11-1],"%lG",(double *)pValue);
							break;
			case(ComplexInt):
#if 0 // no ComplexInt structure defined yet
							sscanf(&line[11-1],"(%d,%d)",&((struct ComplexInt *)pValue)->re,&((struct ComplexInt *)pValue)->im);
							i0=0;
#endif
							break;
			case(Complex): // i e, complex double
							sscanf(&line[11-1],"(%lG,%lG)",&((struct Complex *)pValue)->re,&((struct Complex *)pValue)->im);
							i0=0;
							break;
			case(None): // only here to prevent compiler warning
				break;
		}
		if(!i0) i0=strcspn(&line[11-1],")")+1+1; // skip the ) and convert to 1-based column
		if(i0<0) i0=-i0; // case of string
 	}
	else i0=9;
	
	// Comment
	if(comment!=NULL)
	{
		len=80-(i0+2)+1;
		strncpy(comment,&line[i0+2-1],len); comment[len]=0;
	}
	return 0;
}

// Read the FITS primary header, including END and padding to multiple of 2880 bytes
int FITSread_endHeader(void)
{
	FITSread_headerLine("END",None,NULL,NULL,0);
	while(nHeaderRecs%36) FITSread_headerLine("",None,NULL,NULL,1); // all-blank lines
	return 0;
}

// Get a C 2-dimensional array of doubles from the FITS file, reversing the indexing so it is in
// C index order, not the FORTRAN order required by FITS.
int FITSread_array2D(double *pArray, int dataSize, size_t nx, size_t ny)
{
	void *pIn=NULL; int i,j,elementSize; long n; unsigned char c;
	union U
	{
		double d;
		float f;
		unsigned char c[8];
	};
	union U t;
	
	// Allocate memory for a single input row
	if(dataSize==DOUBLE) {elementSize=sizeof(double); pIn=(double *)malloc(nx*elementSize);}
	else                 {elementSize=sizeof(float);  pIn=(float  *)malloc(nx*elementSize);}
	if (pIn ==NULL) {fprintf(stderr,"FATAL: FITSread_array2D: no malloc\n"); exit(1);}

	for (j=0; j<ny; j++)
	{
		n=ftell(inputFile); // for debugging only
		n=fread(pIn,elementSize,nx,inputFile);
		if(n!=nx) return 100000+j;
		if(feof(inputFile)) return 200000+j;
		for (i=0; i<nx; i++)
		{
				if(dataSize==DOUBLE)
				{
					 // Convert big-endian doubles to little-endian
					 t.d=((double *)pIn)[i];
					 if (BYTE_ORDER == LITTLE_ENDIAN) {
					   c=t.c[7]; t.c[7]=t.c[0]; t.c[0]=c; c=t.c[6]; t.c[6]=t.c[1]; t.c[1]=c;
					   c=t.c[5]; t.c[5]=t.c[2]; t.c[2]=c; c=t.c[4]; t.c[4]=t.c[3]; t.c[3]=c;
					 }
					 pArray[i*ny+j]=t.d;
				}
				else
				{
					 // Convert big-endian floats to little-endian
					 t.f=((float  *)pIn)[i];
					 if (BYTE_ORDER == LITTLE_ENDIAN) {
					   c=t.c[3]; t.c[3]=t.c[0]; t.c[0]=c; c=t.c[2]; t.c[2]=t.c[1]; t.c[1]=c;
					 }
					 pArray[i*ny+j]=t.f;
				}
		}
	}
	
	free(pIn);
	return 0;
}

// Pad the array to a multiple of 2880 bytes and close the file
int FITSread_endFile(void)
{
	long n; char zeros[2880-1]; int rem;
	n=ftell(inputFile);
	if(n%2880)
	{
		rem=2880-(n%2880);
		memset(zeros,0,rem);
		fread(zeros,1,rem,inputFile);
	}
	fclose(inputFile);
	inputFile=NULL;
	return 0;
}

// Close file
int FITSclose(void)
{
     fclose(inputFile);
     inputFile=NULL;
     return 0;
}

// See if the input is all blanks
int IsAllBlanks(char line[], int len)
{
	int i;
	for (i=0; i<len; i++) if(line[i]!=' ') return 0;
	return 1;
}
