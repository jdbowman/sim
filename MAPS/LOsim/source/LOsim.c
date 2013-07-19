// LOsim.c: Low frequency array (LOFAR) simulator main program
//

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// ?/?/?	Create
// 7/6/01	v1.0
// 10/23/01	v1.09: list of functions as command-line input
// 12/2/01	v1.10: position checking
// 12/2/01	v1.11: station beam shape
// 2/22/02	v1.15: correct FITS brightness file sign error

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "SourceListRead.h"
#include "GlobalReferences.h"
#include "GlobalVariables.h"
#include "ScreenIO.h"
#include "ObservationRead.h"
#include "StationListRead.h"
#include "SumSources.h"
#include "OOBSources.h"
#include "RowTransformBrightness.h"
#include "ColumnTransformBrightness.h"
#include "Utilities.h"
#include "SimFilename.h"
#include "Parameters.h"
// #include "uvCoverage.h" ***** taken out 24sept03 
//           by ssd - problem with x11 libraries
// #include "VectorIntegrate.h" ***** not needed for LOsim now

FILE *outFile;

extern int pos;	// line number currently at the top of the screen

int main(int nArgs, char *pArgs[])
{
  char filename[100], temp[13], logFile[25], functions[20+1]; 
  time_t now; int running; char *pFilename;
  struct Source dummySource; int dummyEOF;
  double obsDuration; int nTimes;
  int year,month,day,hr,min,sec,msec; long long time0;
  int result=0;
	
  ScreenInit();
  printf("\nLOFAR simulator v1.15  2/22/02");
#if testPosition
  printf("\n*** Position testing is enabled; output files contain " \
	 "test data only");
#endif
		
  // If simulation name was not given when the pgm was invoked, 
  // ask for it now
  if(nArgs<2 || pArgs[1]==NULL) pFilename=NULL;
  else pFilename=pArgs[1];
  if(GlobalVariables(pFilename)) 
    // reads description and observation files
    {fprintf(stderr,"\nExiting\n"); exit(0);} 
	
  // Pick up list of functions to be executed
  functions[0]=0; // default is no functions specified on command line
  if(nArgs>=3 && pArgs[2]!=NULL && pArgs[2][0]=='-' && pArgs[2][1]=='f')
    {
      // don't count "-f" as part of the length
      int len=strlen(pArgs[2])-2; 
      if(len>20) 
	fprintf(stderr,"WARNING: %d functions given on command " \
		"line; only the first 20 are being used", len);
      // in case of max of 20 functions
      strncpy(functions,&pArgs[2][2],20); functions[20]=0; 
    }
	
	
  // Get a date_time suffix for the filename
  now=time(NULL); 
  strftime(temp,13,"_%y%m%d%_H%M",localtime(&now)); // _yymmdd_hhmm
  CategoryFilename(logFile,"LOsimLog",temp,"txt");
  strcpy(logFile,"/dev/null"); // TEMPORARY-eliminate log files
  // creates a new file if none exists, else overwrites existing file
  outFile=fopen(logFile,"w"); 

  OutputToScreen(ctime(&now),0);
  running=1;
  while(running){
    char funcC;
    if(strlen(functions)) {
      // remove the first char from 'functions'
      funcC=functions[0]; memmove(functions,&functions[1],strlen(functions)); 
      fprintf(stderr,"\nStarting function '%c'",funcC);
    }
    else {
      printf("\nFunction: "); 
      fgets(temp,sizeof(temp),stdin); 
      funcC=temp[0];
    }
    switch(funcC)
      {
	// Ignore blanks (useful if functions from command line)
      case(' '):
	break;
				
	// for debugging ColumnTransformBrightness
      case('C'): case('c'):
	SimFilename(filename,"BrFFT","dat");
	ColumnTransformBrightness(filename);
	fflush(outFile);
	break;
	// for debugging RowTransformBrightness
      case('W'): case('w'):
	SimFilename(filename,"Brightness","fts");
	RowTransformBrightness(filename,title);
	fflush(outFile);
	break;
			
	// Fourier transform brightness to visibilities
      case('F'): case('f'):
	SimFilename(filename,"Brightness","fts");
	RowTransformBrightness(filename,title);
	SimFilename(filename,"BrFFT","dat");
	ColumnTransformBrightness(filename);
	fflush(outFile);
	break;
			
	// Sum out-of-beam sources' contributions to brightness
      case('O'): case('o'):
	CategoryFilename(filename,"OOBSourceList",OOBSourceListName,"txt");
	SourceListRead(SLR_CHECK,filename,&dummySource,&dummyEOF);
	fflush(outFile);
	SimFilename(filename,"Visibility","dat");
	OOBSources(title);
	break;
 			
	// Check in-beam source list (not necessary for a simulation, 
	// as SumSources reads source file)
      case('R'): case('r'):
	CategoryFilename(filename,"SourceList",sourceListName,"txt");
	result = SourceListRead(SLR_CHECK,filename,&dummySource,&dummyEOF);
	fflush(outFile);
	break;

	// Read list of station and antenna locations
      case('L'): case('l'):
	CategoryFilename(filename,"StationList",stationListName,"txt");
	StationListRead(filename);
	fflush(outFile);
	break;

	// Sum in-beam sources
      case('S'): case('s'):
	SumSources();
	fflush(outFile);
	break;
			
	// Read observation parameters (only used in special cases, 
	// since GlobalVariables reads these parameters)
      case('T'): case('t'):
	ObservationRead();
	fflush(outFile);
	break;
			
	// Print and plot (u,v) coverage
      case('U'): case('u'):
	obsDuration=DiffTime(endObs,startObs); // total observation time in sec
	nTimes=ceil(obsDuration/integrationTime);
	// Since center of integration will be at 1/2 integration 
	// time + start time, make calculations at center
	UnpackTime(startObs,&year,&month,&day,&hr,&min,&sec,&msec);
	AddTime(integrationTime/2,&year,&month,&day,&hr,&min,&sec,&msec);
	PackTime(year,month,day,hr,min,sec,msec,&time0);
	/* uvCoverage taken out - problem with x11 libraries after upgrade */
	/* uvCoverage(nTimes,time0,integrationTime,-1,-1,NULL,NULL); 
	 * // print and plot */
	fflush(outFile);
	break;

	// Integrate visibility data over observation period
      case('I'): case('i'):
	obsDuration=DiffTime(endObs,startObs); // total observation time in sec
	nTimes=ceil(obsDuration/integrationTime);
	/* ***Not needed for LOsim now - this done in visgen	*/
	/* VectorIntegrate(nTimes,startObs,integrationTime); */
	fflush(outFile);
	break;

	// Exit simulator
      case('Q'): case('q'): case('X'): case('x'):
	running=0;
	fclose(outFile);
	break;
			
      default:
	fprintf(stderr,"\a   Available functions are"
		"\nR=read in-beam source list  L=read station list  " \
		"T=read observation parameters"
		"\nS=sum in-beam sources  O=add visibility from " \
		"out-of-beam sources"
		"\nF=FFT sources (same as W followed by C)  " \
		"W=row transform sources  C=column transform " \
		"row-transformed sources"
		"\nI=integrate visibility data"
		"\nX=exit");
	break;
      }
		
  }
	
  return 0;
}
