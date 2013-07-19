// GlobalVariables.c: All the global variables are stored here; 
// also invokes the description file read-in

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/4/01	create
// 8/22/01	add observation variables

#include "GlobalVariables.h"
#include "GlobalReferences.h" // includes DataStructures.h
#include "Parameters.h"
#include "DescriptionFileRead.h"
#include "SimFilename.h"
#include "ObservationRead.h"
#include <stdio.h>
#include <string.h>

//
// Global variables obtained from the description file:
//
// M short phrase describing this simulation 
// (e g, "Random 245 sources at full resolution")
char title[L_ONELINE]; 
// M single word describing this simulation (e g, "Rand245")
char simName[L_BRIEFNAME];
// mm/dd/yy hh:mm:ss when this description file was first created
long long createDateTime; 
// The name of the author of this simulation description file
char creator[L_NORMAL];
// The institutional affiliation of the author of this 
// simulation description file
char institution[L_NORMAL];
// The e-mail address of the author of this simulation description file
char email[L_NORMAL];			
char phone[L_BRIEFNAME];		// the telephone number of the author of this simulation description file
long long changeDateTime;		// similar information for the person who last modified this description file
char changer[L_NORMAL];
char description[L_EXTENDEDLINE];	// a one-line description of the simulation, providing additional information from the title
char comments[L_EXTENDEDLINE];		// any further remarks about this simulation

long brightnessXN;		// M number of points in the sky plane grid, x
long brightnessYN;		// M                                         y

double fieldXSize;		// M size of the beam, x, in radians
double fieldYSize;		// M                   y, in radians

long visibilityXLn2;	// M log2 of number of points in the visibility grid, x
long visibilityYLn2;	// M                                                  y

// Observation parameters
double fieldRA;				// right ascension of the center of the beam (field of view), in radians
double fieldDec;			// declination of the center of the beam (field of view), in radians
long long startObs,endObs;	// limits for observation time, microseconds since ???
double integrationTime;		// time over which individual visibilities are averaged, seconds
double obsFreq;				// center observation frequency, in Hz
double obsBandwidth;		// bandwidth frequency in Hz
int nSpectralPoints;		// number of sub-bandwidths the bandwidth is divided into

// Names of files containing user specifications
char sourceListName[L_BRIEFNAME];			// M in-beam source list is in file SourceList_<name>.txt
char OOBSourceListName[L_BRIEFNAME]; 		// M out-of-beam source list is in file OOBSourceList_<name>.txt
char observationName[L_BRIEFNAME];			// M observing parameters are in file Observation_<name>.txt
char stationListName[L_BRIEFNAME];			// M station list is in file StationList_<name>.txt
char ionosphereParametersName[L_BRIEFNAME];	// M parameters for model of ionosphere are in file Ionosphere_<name>.txt
char imageParametersName[L_BRIEFNAME];		// parameters for externally-supplied sky image files are in file Images_<name>.txt

// Stations and their component antennae
int nStations;
struct Station station[MAX_STATIONS];
struct Antenna antenna[MAX_ANTENNAE];

// Derived variables
long visibilityXN,visibilityYN; // number of points in the visibility grid
/* The brightness array is padded with zeros with each row and column consisting of

		brightnessXNzero0 zero elements  brightnessXN brightness elements   brightnessXNzero1 zero elements (row)
		brightnessYNzero0 zero elements  brightnessYN brightness elements   brightnessYNzero1 zero elements (column)
		
		brightnessXNzero0 + brightnessXN + brightnessXNzero1 = visibilityXN
		brightnessYNzero0 + brightnessYN + brightnessYNzero1 = visibilityYN
*/
long brightnessXNzero0, brightnessXNzero1;
long brightnessYNzero0, brightnessYNzero1;
long brightnessMinX, brightnessMaxX, brightnessMinY, brightnessMaxY;

// Miscellaneous constants
double oneDegree,oneArcsec; // in radians

// Other stuff
int nErrors; // overall number of errors in input (initialized and updated by parameter-file-read functions and their callers)

// Debugging options  All default to 0 (false) unless set explicitly by the description file
int debugOption_noSourceTruncate;

/* Inputs: mainArg==NULL means to ask for simulation name

   Output:	0=no error
   			1=simulation name null
   			2=description file not found
   			3=error reading description file
   			4=observation file not found
   			5=error reading observation file
*/

int GlobalVariables(char *mainArg)
{
	int ln,err; char filename[100];
	
	// Make the default debugging options disabled
	debugOption_noSourceTruncate=0;
	
	// Store the conversion constants (first, as they're used by DescriptionFileRead)
	oneDegree=2.0L*PI/360.0L; oneArcsec=2.0L*PI/(360.0L*60.0L*60.0L); // in radians
	
	// If simulation name was not given when the pgm was invoked, ask for it now
	if(mainArg==NULL)
	{
		printf("\n\nSimulation name: "); fgets(simName,sizeof(simName),stdin);
		ln=strlen(simName); if(ln>=1 && simName[ln-1]=='\n') simName[ln-1]=0; // remove the NL if we read a complete line (normal)
	}
	else strcpy(simName,mainArg);
	if(!strlen(simName)) return 1; // exit on ENTER (do we enter on ESCAPE?!)
	
	//  Get the data from the description file
	if((err=DescriptionFileRead(simName,filename)))
	{
		fflush(stdout);
		if(err==1) {fprintf(stderr,"\nCan't find the description file '%s'",filename); return 2;}
		else {fprintf(stderr,"\nErrors occurred while reading the description file '%s'",filename); return 3;}
	}
	
	//  Get additional data from the observation file
	CategoryFilename(filename,"Observation",observationName,"txt");
	printf("\n--- Reading observation file %s",filename);
	if((err=ObservationRead()))
	{
		if(err==1) {fprintf(stderr,"\nCan't find the observation file '%s'",filename); return 4;}
		else {fprintf(stderr,"\nError in reading the observation file '%s'",filename); return 5;}
	}
	
	// Set up the derived global variables
	visibilityXN = 1<<visibilityXLn2; 
	visibilityYN = 1<<visibilityYLn2; // number of points in the visibility grid                                               y
	brightnessXNzero0 = (visibilityXN - brightnessXN)/2; 
	brightnessXNzero1 = visibilityXN - brightnessXN - brightnessXNzero0;
	brightnessYNzero0 = (visibilityYN - brightnessYN)/2; 
	brightnessYNzero1 = visibilityYN - brightnessYN - brightnessYNzero0;

	brightnessMaxX =  (brightnessXN - 1)/2; 
	brightnessMinX = -(brightnessXN - 1 - brightnessMaxX); // if brightnessZN is even, the extra point is negative
	brightnessMaxY = (brightnessYN - 1)/2; 
	brightnessMinY = -(brightnessYN - 1 - brightnessMaxY);
	
    printf("Global Vars: visibilityXN = %ld\n", visibilityXN);
    printf("Global Vars: brightnessXN = %ld\n", brightnessXN);
    printf("Global Vars: brightnessMixX = %ld\n", brightnessMinX);
    printf("Global Vars: brightnessMaxX = %ld\n", brightnessMaxX);

    return 0;
}
