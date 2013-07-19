// ImagesRead.c: Read the input file containing the user's external image parameters

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 1/22/02	Create

#include "ImagesRead.h" // includes DataStructures
#include "GlobalReferences.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include "SimFilename.h"
#include "Parameters.h"
#include "Utilities.h"

#define MAX_FIELDS 20 // used by following #include'd file
#include "TextFileParse.h"

// Local function
void ImageStore(struct ImageFileParams imageFrom, struct ImageFileParams *pImageTo, int nFieldNames, int rq[], int fieldSeen[]);

int ImagesRead(int *pnImages, struct ImageFileParams image[])
{
	int err;
	char temp[200];
	int i,iI,ln,nLine=0,prevFieldSeen,v;
	FILE *imFile=NULL;
	
	// Data for one image
	struct ImageFileParams im;
	im.filename[0] = '\0';
	im.limit=0;
	im.blank=0;
	im.stretch=0;
	im.scale=0;
	im.rotate=0;
	im.scaleX=0;
	im.scaleY=0;
	im.ctrX=0;
	im.ctrY=0;
	im.reflectX=0;
	im.reflectY=0;
	im.isFITS=0;

	// Create the input parsing table
	nFieldNames=0;
	// The following must be first (it signals the start of an image description)
	SET_INPUT("File containing image",1,DT_STR,"",0,0,0,"","","",(im.filename))
	
	SET_INPUT("FITS file",1,DT_BOOL,"",0,0,0,"1 if file is in FITS format, else 0","","",im.isFITS)
	SET_INPUT("Reflect through x-axis",0,DT_BOOL,"",0,0,0,"","","",im.reflectX)
	SET_INPUT("Reflect through y-axis",0,DT_BOOL,"",0,0,0,"","","",im.reflectY)
	SET_INPUT("Center x",0,DT_DBL,"",-1300.0,1300.0,0,"","","",im.ctrX)
	SET_INPUT("Center y",0,DT_DBL,"",-1300.0,1300.0,0,"","","",im.ctrY)
	SET_INPUT("Scale x",0,DT_DBL,"",0.01,100.0,0,"Scaling must be 0.01-100","","",im.scaleX)
	SET_INPUT("Scale y",0,DT_DBL,"",0.01,100.0,0,"Scaling must be 0.01-100","","",im.scaleY)
	SET_INPUT("Rotate",0,DT_DBL,"",-PI,+PI,0,"The rotation angle must be in the range -pi < angle < +pi","","",im.rotate)
	SET_INPUT("Pixel scale",0,DT_DBL,"",1e-6,1e6,0,"Pixels can be scaled by 10^+-6","","",im.scale)
	SET_INPUT("Pixel stretch power",0,DT_DBL,"",0.1L,5.0L,0,"The pixel power law exponent must be in the range 0.1-5","","",im.stretch)
	SET_INPUT("Pixel blank",0,DT_DBL,"",0,1e05,0,"Pixel blank range is 0-10^5","","",im.blank)
	SET_INPUT("Pixel limit",0,DT_DBL,"",0,1e10,0,"Pixel limit range is 0-10^10","","",im.limit)

//    printf ("filename = %s\n",filename);
	
//	if(!strlen(filename))
    if(!strlen(imageParametersName))
	{
		*pnImages=0; // # images stored
		printf("\n[no external images were specified]");	
		fflush(stdout);
		return 0;		
	}
	
	// Open the images file
	CategoryFilename(filename,"Images",imageParametersName,"txt");

    printf ("filename = %s\n",filename);

	if((imFile=fopen(filename,"r"))==NULL)
	{
		sprintf(temp,"Can't find the images parameters file '%s'\n",filename);
		ErrorPrint(temp);
		return 1;
	}
	
	// Blank all the fields to be filled in from this file
	for (i=0; i<nFieldNames; i++) fieldSeen[i]=0;
	
	nErrors=0; iI=0;
	
	// Process the lines of the images file
	while(1)
	{
		if((ln=GetLine(imFile,lineIn,line,MAX_LINE))<0) break; // get a line from the file, stopping on error or EOF
		nLine++;
		printf("%3d: %s\n",nLine,lineIn);
		
		if(!ln) continue; // skip a blank (or comment-only) line
		
		/* Lines in the images file have the format
		   <fieldname> = <value>
		   The blanks around the "=" are optional
		*/
		if((p=strchr(line,'='))==NULL)
		{
			ErrorPrint("Line format should be <fieldname> = <value>");
			continue; // ignore this line
		}
		strncpy(fieldname,line,p-line); fieldname[p-line]=0; RemoveWhitespace(fieldname);
		strcpy(value,p+1); RemoveWhitespace(value); // skip the '='
		
		// Identify the field
		for (i=0; i<nFieldNames; i++)
		{
			if(FieldMatches(fieldname,pText[i])) break;
		}
		
		if(i>0 && !fieldSeen[0])
		{
			ErrorPrint("The name of the file containing the image must come first");
		}
		
		// Get the field value and store it
		if(i<nFieldNames)
		{
			// Check for repeated field
			if(i>0 && fieldSeen[i])
			{
				sprintf(temp,"This field was already assigned in line %d",fieldSeen[i]); ErrorPrint(temp);
				continue;
			}
			prevFieldSeen=fieldSeen[i]; // non-0 if this is the second or subsequent occurrence of the field
			fieldSeen[i]=nLine; // remember where we first saw this field
			
			// Match the appropriate pattern for this kind of data
			if(*pPat[i])
			{
				if(!PatternMatch(value,pPat[i])) ErrorPrint(pErr[i][0]);
			}
			
			switch(dataType[i])
			{
				case(DT_BOOL):
					v=toupper(value[0]);
					switch(v)
					{
						case('0'): case('N'): case('F'): v=0; break;
						case('1'): case('Y'): case('T'): v=1; break;
						default: v=-1;
					}
					if(v>=0) *((int *)pResult[i])=v;
					break;
				case(DT_INT):
					IntegerValue((long *)pResult[i],(long)lim[i][0],(long)lim[i][1],value); break;
				case(DT_DBL):
					if((err=DoubleValue((double *)pResult[i],lim[i][0],lim[i][1],value)))
					{
						if(err==1) ErrorPrint(pErr[i][0]);
						else if(err==2) ErrorPrint(pErr[i][1]);
					}
					break;
				case(DT_STR):
					if(!i)
					{
						if(prevFieldSeen) // if we've already seen an image, store its data
						{
							ImageStore(im,&image[iI],nFieldNames,rq,fieldSeen);
							iI++;
						}
						for (i=1; i<nFieldNames; i++) fieldSeen[i]=0; // all but the filename
						i=0;
						if(iI>=MAX_IMAGES) {sprintf(temp,"Too many images; only %d are allowed",MAX_IMAGES); ErrorPrint(temp); --iI;}
						// Initialize all the non-required fields to their defaults
						im.reflectX=im.reflectY=0;
						im.ctrX=im.ctrY=0.0;
						im.scaleX=im.scaleY=1.0; im.rotate=0.0; im.scale=im.stretch=1.0;
						im.blank=-DBL_MAX; im.limit=DBL_MAX;
					}
					StringValue((char *)pResult[i],value,L_ONELINE);
					break;
				default: 
					break;
			}
		}
		else
		{
			sprintf(temp,"Field name '%s' is not known",fieldname);
			ErrorPrint(temp);
		}
	}
	
	if(fieldSeen[0]) {ImageStore(im,&image[iI],nFieldNames,rq,fieldSeen); iI++;} // if we've already seen an image, store its data
	*pnImages=iI; // # images stored
	fflush(stdout);
	return nErrors;		
}

void ImageStore(struct ImageFileParams imageFrom, struct ImageFileParams *pImageTo, int nFieldNames, int rq[], int fieldSeen[])
{
	int i,err; char temp[200+1];
	
	// Check that all the mandatory fields were filled in
	for (i=0,err=0; i<nFieldNames; i++)
	{
		if(rq[i] && !fieldSeen[i])
		{
			err++;
			sprintf(temp,"The field '%s' must be defined for this image",pText[i]); ErrorPrint(temp);
		}
	}
	if(!err) *pImageTo=imageFrom;
	return;
}
