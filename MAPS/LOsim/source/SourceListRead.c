// SourceListRead.c: Get the list of input sources from a text file
// Currently all sources are Gaussian

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 6/1/01	Remove \r from input lines, just in case; 
//              more exact conversions to radians; include 1st source
// 7/11/01	Variable array sizes (v1.0)
// 8/22/01	Option to check sources, or to read in one at a time
// 1/11/02	Allow comment lines (v1.13)

#include <stdio.h>
#include <float.h>
#include <string.h>
#include "Parameters.h"
#include "GlobalReferences.h"
#include "SourceListRead.h"
//#include "SourceArrays.h"
#include "Utilities.h"
#include "ScreenIO.h"
#define MAX_LINE 300

/* Read a source file
 * If whatTodo=SLR_CHECK, open the source list file and read all entries, 
 *      checking values (but not storing)
 * If whatToDo=SLR_FIRST, open the source list file and read, check and 
 *      return the first entry
 * If whatToDo=SLR_NEXT, read, check and return the next entry
 * If the output EOF=0, the entry read is returned in 'source'; if EOF=1, 
 *        there are no more entries.
 */

/* private function prototypes */
int GetNextSource(struct Source *pSource);

// Stored when the file is opened, and preserved across 
// invocations of this function:
static FILE *inFile = NULL; 
static char screenLine[200]; // Static to prevent export to linker

int SourceListRead(int whatToDo, char sourceListFilename[], 
		   struct Source *pSource, int *pEOF) {
  int nSources;
  double fdMin, xOMin, yOMin, majAMin, minAMin, angMin, fdMax, xOMax, yOMax;
  double majAMax, minAMax, angMax;

  if (whatToDo != SLR_NEXT) {
      inFile = fopen(sourceListFilename, "r");
      if (inFile == NULL) {
	fprintf(stderr,"Failed to read SourceList file: %s\n",
		sourceListFilename);
	return 1;
      }
      sprintf(screenLine, " Source ID    Flux      Stokes    " \
	      "Spectral         xOffset         yOffset      Major" \
	      "      Minor  Position");
      sprintf(screenLine, "             density     Q U V    index" \
	      "                                          axis" \
	      "      axis   angle");
      OutputToScreen(screenLine, 0);
    }
	
  if(whatToDo == SLR_CHECK) {
      // Read and check all sources
      nSources = 0;
      fdMin = xOMin = yOMin = majAMin = minAMin = angMin = DBL_MAX; 
      fdMax = xOMax = yOMax = majAMax = minAMax = angMax = -DBL_MAX;
      while(!GetNextSource(pSource)) {
	  // Must do both tests because max<min initially
	  if (pSource->fluxDensity < fdMin) fdMin=pSource->fluxDensity; 
	  if (pSource->fluxDensity > fdMax) fdMax=pSource->fluxDensity;
	  if (pSource->xOffset < xOMin) xOMin = pSource->xOffset; 
	  if (pSource->xOffset > xOMax) xOMax = pSource->xOffset;
	  if (pSource->yOffset < yOMin) yOMin = pSource->yOffset; 
	  if (pSource->yOffset > yOMax) yOMax = pSource->yOffset;
	  if (pSource->majorAxis < majAMin) majAMin = pSource->majorAxis; 
	  if (pSource->majorAxis > majAMax) majAMax = pSource->majorAxis;
	  if (pSource->minorAxis < minAMin) minAMin = pSource->minorAxis; 
	  if (pSource->minorAxis > minAMax) minAMax = pSource->minorAxis;
	  if (pSource->positionAngle < angMin) angMin=pSource->positionAngle; 
	  if (pSource->positionAngle > angMax) angMax = pSource->positionAngle;
	  nSources++;
	}
      SkipLine();
      sprintf(screenLine, "%g <= fluxDensity <= %g", fdMin, fdMax);
      OutputToScreen(screenLine, 0);
      sprintf(screenLine,
	      "%g <= xOffset <= %g     %g <= yOffse t <= %g       " \
	      "%g <= majorAxis <= %g  %g <= minorAxis <= %g",
	      xOMin/oneArcsec, xOMax/oneArcsec, 
	      // report results in the input units, not radians
	      yOMin/oneArcsec, yOMax/oneArcsec,
	      majAMin/oneArcsec, majAMax/oneArcsec,
	      minAMin/oneArcsec, minAMax/oneArcsec);
      OutputToScreen(screenLine,0);
      sprintf(screenLine, "%g <= positionAngle <= %g",
	      angMin/oneDegree, angMax/oneDegree);
      OutputToScreen(screenLine, 0);
      *pEOF=0;
    }
  else *pEOF = GetNextSource(pSource);
	
  return 0;
}
	
// Returns 1 if no source seen before EOF, else 0
int GetNextSource(struct Source *pSource) {
  char lineIn[MAX_LINE], line[MAX_LINE], id[10+1];
  double fd, Q, U, V, six, xO, yO, majA, minA, ang; 
  int n;
	
  while(1) { // read till we get a source, or come to EOF
    // End-of-file:
      if((n = GetLine(inFile, lineIn, line, MAX_LINE)) < 0) return 1; 
      if(!n) continue; // skip a blank or comment line
      n = sscanf(line, "%10s %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
		 id, &fd, &Q, &U, &V, &six, &xO, &yO, &majA, &minA, &ang);
      if(n != 11) {
	ErrorPrint("Source lines require 11 values: ID, " \
		   "intensity, Q, U, V, spectral_index, x, y, " \
		   "major_axis, minor_axis, position_angle");
      }
      sprintf(screenLine, "%10s %15g %15g %15g %15g %15g %15g", 
	      id, fd, xO, yO, majA, minA, ang);
      OutputToScreen(screenLine, 0);
      // ADD CHECKING HERE for fd, Q, etc.
      strcpy(pSource->id, id); 
      pSource->fluxDensity = fd; 
      pSource->Q = Q; 
      pSource->U = U; 
      pSource->V = V; 
      pSource->specIx = six;
		
      // Convert units to radians
      // 1 arcsec=4.848*10^-6 radians, 
      // approx 1 degree=1.745*10^-2 radians.
      pSource->xOffset = xO*oneArcsec; 
      pSource->yOffset = yO*oneArcsec; 
      pSource->majorAxis = majA*oneArcsec; 
      pSource->minorAxis = minA*oneArcsec;
      pSource->positionAngle = ang*oneDegree;
      return 0;
    }
}
