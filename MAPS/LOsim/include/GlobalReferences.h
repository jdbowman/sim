// GlobalReferences.h: Declares all the global variables (which are stored in GlobalVariables.c)

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// ?/?/01	create
// 8/22/01	add observation variables

#include "DataStructures.h"

extern char title[],simName[],creator[],institution[],email[],phone[],
            changer[],description[],comments[];
extern long long createDateTime,changeDateTime;
extern long brightnessXN,brightnessYN,visibilityXLn2,visibilityYLn2;
extern double fieldRA,fieldDec,fieldXSize,fieldYSize,obsFreq;

extern char sourceListName[],OOBSourceListName[],observationName[];
extern char stationListName[],ionosphereParametersName[],imageParametersName[];

extern long visibilityXN,visibilityYN;
extern long brightnessXNzero0,brightnessXNzero1,brightnessYNzero0,brightnessYNzero1;
extern long brightnessMinX,brightnessMaxX,brightnessMinY,brightnessMaxY;

extern long long startObs,endObs;
extern double integrationTime,obsFreq,obsBandwidth;
extern int nSpectralPoints;

extern double oneDegree,oneArcsec;

extern int nStations;
extern struct Station station[];
extern struct Antenna antenna[];

extern int nErrors;

extern int debugOption_noSourceTruncate;
