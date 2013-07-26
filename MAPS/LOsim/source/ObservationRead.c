// ObservationRead.c: Read the input file containing the user's observing parameters

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/23/01	Create

#define SET_INPUT(text,r,type,pat,ver0,ver1,ver2,er0,er1,er2,variable) {pText[nFieldNames]=text; rq[nFieldNames]=r; \
  dataType[nFieldNames]=type; pResult[nFieldNames]=&variable; lim[nFieldNames][0]=ver0; lim[nFieldNames][1]=ver1; lim[nFieldNames][2]=ver2; \
  pPat[nFieldNames]=pat; pErr[nFieldNames][0]=er0; pErr[nFieldNames][1]=er1; pErr[nFieldNames][2]=er2; nFieldNames++;}

// Data types
#define DT_STR 1	// character string
#define DT_INT 2	// integer
#define DT_HMS 3	// hour:min:sec
#define DT_DEG 4	// deg, arcmin, arcsec
#define DT_DBL 8	// real number, double precision
#define DT_TIME 10	// date and time

#include "ObservationRead.h"
#include "GlobalReferences.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "SimFilename.h"
#include "Parameters.h"
#include "Utilities.h"

int ObservationRead(void)
{
#define MAX_FIELDS 10
	int nFieldNames,fieldSeen[MAX_FIELDS],rq[MAX_FIELDS],err;
	char *pText[MAX_FIELDS]; char *pPat[MAX_FIELDS]; char *pErr[MAX_FIELDS][3]; double lim[MAX_FIELDS][3];
	void *pResult[MAX_FIELDS]; int dataType[MAX_FIELDS];
	#define MAX_LINE 150
	char filename[100],lineIn[MAX_LINE],line[MAX_LINE],*p,fieldname[MAX_FIELDNAME+1],value[L_EXTENDEDLINE],temp[100];
	int i,ln,nLine=0,d,h,m,s; double sf;
	FILE *obsFile;

	// Create the input parsing table
	nFieldNames=0;
	SET_INPUT("Field of view center right ascension",1,DT_HMS,"2D\":\"2D\":\"2D",0,0,0, \
	          "RA format is hours:minutes:seconds, each 2 digits (ex: 12:34:56)","","",fieldRA)
	SET_INPUT("Field of view center declination",1,DT_DEG,"0.1\"-\"1.D\"d \"1.D\"' \"1.D0.1\".\"0.D\"\\\"\"",-(90*60*60),+(90*60*60),0, \
	          "Declination format is like 42d 7' 22\" or -42d 7' 22\"","Maximum declination is +- 90d 0' 0\"","",fieldDec)
	SET_INPUT("Field of view x-size",1,DT_DBL,"",1.0L,648000.0L,0,"","","",fieldXSize)
	SET_INPUT("Field of view y-size",1,DT_DBL,"",1.0L,648000.0L,0,"","","",fieldYSize)
	SET_INPUT("Start time",1,DT_TIME,"",0,0,0,"Date/time not valid, or out of range 1901-2099","","",startObs)
	SET_INPUT("End time",1,DT_TIME,"",0,0,0,"Date/time not valid, or out of range 1901-2099","","",endObs)
	SET_INPUT("Integration time",1,DT_DBL,"",1.0,600.0,0,"Invalid numeric format","The integration time should be 1-600 sec","",integrationTime)
	SET_INPUT("Observation frequency",1,DT_DBL,"",10.0,20000.0,0,"The observation frequency range is 10 MHz-20 GHz","","",obsFreq)
	SET_INPUT("Bandwidth",1,DT_DBL,"",0.01L,500.0L,0,"The bandwidth range is 10 KHz-500 MHz","","",obsBandwidth)
	SET_INPUT("Spectral points",1,DT_INT,"",1,MAX_SPECTRAL_POINTS,0, \
	  "There must be 1-1024 spectral points","","",nSpectralPoints) // FAKE-this should be MAX_SPECTRAL_POINTS

//    printf ("filename = %s\n", filename);
	
	// Open the observation file
	CategoryFilename(filename,"Observation",observationName,"txt");
	if((obsFile=fopen(filename,"r"))==NULL)
	{
		sprintf(temp,"Can't find the observation file '%s'",filename);
		ErrorPrint(temp);
		return 1;
	}
	
	// Blank all the fields to be filled in from this file (LATER CHANGE-set defaults)
	for (i=0; i<nFieldNames; i++) fieldSeen[i]=0;
	
	nErrors=0;
	
	// Process the lines of the observation file
	while(1)
	{
		if((ln=GetLine(obsFile,lineIn,line,MAX_LINE))<0) break; // get a line from the file, stopping on error or EOF
		nLine++;
		printf("\n%3d: %s",nLine,lineIn);
		
		if(!ln) continue; // skip a blank (or comment-only) line
		
		/* Lines in the observation file have the format
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
		
		// Get the field value and store it
		if(i<nFieldNames)
		{
			// Check for repeated field
			if(fieldSeen[i])
			{
				sprintf(temp,"This field was already assigned in line %d",fieldSeen[i]); ErrorPrint(temp);
				continue;
			}
			fieldSeen[i]=nLine; // remember where we first saw this field
			
			// Match the appropriate pattern for this kind of data
			// There is a special check for deg,min,sec, which uses " in an unbalanced fashion, and therefore
			// may still have a trailing comment
			if(dataType[i]==DT_DEG)
			{
				p=strstr(value,"//"); if(p!=NULL) {*p=0; RemoveWhitespace(value);} // terminate the string, thus removing the comment
			}
			if(*pPat[i])
			{
				if(!PatternMatch(value,pPat[i])) ErrorPrint(pErr[i][0]);
			}
			
			switch(dataType[i])
			{
			case(DT_INT):
				IntegerValue((long *)pResult[i],(long)lim[i][0],(long)lim[i][1],value); break;
			case(DT_HMS):
				sscanf(value,"%d:%d:%d",&h,&m,&s);
				if(h>23 || m>59 || s>59) ErrorPrint("Hours:minutes:seconds limited to 23:59:59");
				else
				{
					s=((h*60)+m)*60+s; *((double *)pResult[i])=s*SECONDStoRADIANS; // convert seconds to radians
				}
				break;
			case(DT_DEG):
				sscanf(value,"%dd %d' %lf\"",&d,&m,&sf);
				if(abs(d)>90 || m>59 || sf>=60.0) ErrorPrint("Degrees, minutes, seconds limited to 90, 59, 59.99...");
				else
				{
				  sf=((abs(d)*60)+m)*60+sf; if(sf<lim[i][0] || sf>lim[i][1]) ErrorPrint(pErr[i][1]);
				  else *((double *)pResult[i])=sf*oneArcsec*(d<0 ? -1.0 : 1.0); // arcseconds->radians
				}
				break;
			case(DT_DBL):
				if((err=DoubleValue((double *)pResult[i],lim[i][0],lim[i][1],value)))
				{
					if(err==1) ErrorPrint(pErr[i][0]);
					else if(err==2) ErrorPrint(pErr[i][1]);
				}
				break;
			case(DT_TIME):
				if(DateAndTime((long long *)pResult[i],value)) ErrorPrint(pErr[i][0]);
				break;
			case(DT_STR):
				StringValue((char *)pResult[i],value,L_ONELINE); break;
			}
		}
		else
		{
			sprintf(temp,"Field name '%s' is not known",fieldname);
			ErrorPrint(temp);
		}
	}
	
	// Check that all the mandatory fields were filled in
	for (i=0; i<nFieldNames; i++)
	{
		if(rq[i] && !fieldSeen[i]) {sprintf(temp,"The field '%s' must be defined in the observation file",pText[i]); ErrorPrint(temp);}
	}
	
	// Convert to radians
	fieldXSize*=oneArcsec; fieldYSize*=oneArcsec;
	// Convert to MHz
	obsFreq*=1e6; obsBandwidth*=1e6;
	
	return nErrors;		
}
