// StationListRead.c: Get the list of stations from a text file

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/4/01	create
// 12/2/01	v1.11: station beam shape
// 1/11/02	v1.13: correct sign of expArg and args of atan2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "GlobalReferences.h"
#include "Parameters.h"
#include "SimFilename.h"
#include "StationListRead.h"
#include "Utilities.h"
#include "ScreenIO.h"
#include "julian_date.h"
#include "sidereal_time.h"
#include "ComplexArithmetic.h"

#define MAX_LINE 200

double Sx,Sy,Sz,sinLat,cosLat,sinLong,cosLong;

void StationToEarth(double *px, double *py, double *pz);

// Read a station file
int StationListRead(char stationListFilename[])
{
	FILE *inFile,*beamFile; char lineIn[MAX_LINE],line[MAX_LINE]; int i,n,nLine,nAntennae,nTotalAntennae,iS,iA,iA0;
	struct Complex wt,B,e,t;
	double x,y,z,r,rTolerance=RADIUSEARTH*0.01,w,sumWts,d0x,d0y,d0z,dx,dy,dz,jd,jd_high,jd_low,ee,gst;
	double hour,H,H0,sinH,cosH,sinL,cosL,sinDec0,cosDec0,dec,decStart,decEnd,decInc,sinDec,cosDec;
	double sinH0,cosH0,HStart,HEnd,HInc,xComp,yComp,zComp,k,expArg;
	int year,month,day,hr,min,sec,msec;
	char screenLine[200],id[20+1],lineType,beamFilename[100],coordOrigin,co;
	
	inFile=fopen(stationListFilename,"r");
	
	nStations=0; nAntennae=-1; nTotalAntennae=0; nLine=nErrors=0;
	sprintf(screenLine,"        Station ID  Coords       x        y         z  # antennae");
	OutputToScreen(screenLine,0);
	
	while((n=GetLine(inFile,lineIn,line,MAX_LINE))>=0)
	{
		nLine++; sprintf(screenLine,"%4d: %s",nLine,lineIn); OutputToScreen(screenLine,0);
		if(!n) continue; // skip a blank (or comment-only) line
		lineType=line[0]; line[0]=' '; RemoveWhitespace(line); // replace type with white space, which is then removed
		switch(lineType)
		{
			case('S'):
			case('s'):
			if(nAntennae>station[nStations-1].nAnt)
			{
				sprintf(screenLine,"The previous station was supposed to have %d antennae, but you only listed %d",
				  nAntennae,station[nStations-1].nAnt);
				ErrorPrint(screenLine);
			}
			else
			{
				n=sscanf(line,"%20s %c %lg %lg %lg %d",id,&co,&x,&y,&z,&nAntennae);
				if(n!=6)
				{
					ErrorPrint("Station lines require 6 values: station_id  coord_origin  x  y  z  number_of_antennae");
					break;
				}
				if(nAntennae<=0)
				{
					ErrorPrint("A station must have at least one antenna");
					break;
				}
				if(nAntennae>MAX_ANTENNAEPERSTATION)
				{
					sprintf(screenLine,"A station cannot have more than %d antennae",MAX_ANTENNAEPERSTATION);
					ErrorPrint(screenLine);
					break;
				}
				
				// Check for a non-null, unique station ID
				if(!strlen(id))
				{
					ErrorPrint("You must supply a station ID");
					break;
				}
				for (i=0; i<nStations; i++) if(!strcmp(station[i].id,id)) break;
				if(i<nStations) // wasn't a unique ID
				{
					sprintf(screenLine,"The ID '%s' was used for station %d already",id,i+1);
					ErrorPrint(screenLine);
					break;
				}
				
				// Make sure the coordinate system is one we recognize
				coordOrigin=toupper(co);
				if(coordOrigin!='S' && coordOrigin!='E')
				{
					sprintf(screenLine,
					  "Antenna coordinate system \"%c\" should be either S (relative to station) or E (relative to center of earth)",co);
					ErrorPrint(screenLine);
				}
				
				// Check that the station coordinates are reasonable
				// NOTE: these coordinates are overwritten by the weighted mean position of the antennae calculated below, and are
				// therefore ignored.
				r=sqrt(x*x+y*y+z*z); if(abs(r-RADIUSEARTH)>rTolerance) // r=distance of station-center of earth
				{
					sprintf(screenLine,"WARNING: The station is at an indicated elevation of %lg km",r-RADIUSEARTH);
					ErrorPrint(screenLine);
				}
				if(nStations>=MAX_STATIONS)
				{
					sprintf(screenLine,"Cannot have more than %d stations",MAX_STATIONS);
					ErrorPrint(screenLine);
					break;
				}
				
				// Store a station
				strcpy(station[nStations].id,id);
				station[nStations].x=x; station[nStations].y=y; station[nStations].z=z;
				station[nStations].longi=atan2(station[nStations].y,station[nStations].x); // -pi   < longitude <= pi
				sinLat=station[nStations].z/r; station[nStations].lat=asin(sinLat); // -pi/2 < latitude  <= pi/2
				if(coordOrigin=='S')
				{
					//Calculate the values needed to convert station-centered to earth-centered antenna coords
					Sx=x; Sy=y; Sz=z; // copy for convenience
					cosLat=cos(station[nStations].lat); sinLong=sin(station[nStations].longi); cosLong=cos(station[nStations].longi);
				}
				station[nStations].nAnt=0;				// exactly 'nAntennae' antenna lines should follow
				station[nStations].iAnt=nTotalAntennae; // where the 1st antenna will appear
				nStations++;
			}
				break;
			
			// Antenna line
			// So far, we've stored station[nStations].nAnt antennae for this station; 'nAntennae' are expected
			case('A'):
			case('a'):
				if(nAntennae<0)
				{
					ErrorPrint("A station line must come before the antenna lines");
					break;
				}
				else if(station[nStations-1].nAnt>=nAntennae)
				{
					ErrorPrint("More antenna lines seen than the number given in the preceding station line");
					break;
				}
				n=sscanf(line,"%10s %lg %lg %lg %lg,%lg",id,&x,&y,&z,&wt.re,&wt.im);
				if(n==4) {wt.re=1.0; n=5;} // weight defaults to (1,0)
				if(n==5) {wt.im=0.0; n=6;} // imaginary part of weight may be omitted
				if(n!=6)
				{
					ErrorPrint("Antenna lines require 4 values: antenna_id  x  y  z weight_real[,imag]");
					break;
				}
				// Check for a non-null, unique antenna ID
				if(!strlen(id))
				{
					ErrorPrint("You must supply an antenna ID");
					break;
				}
				// k=index in the antenna array of an antenna belonging to this station
				for (i=0; i<station[nStations-1].nAnt; i++) if(!strcmp(antenna[station[nStations-1].iAnt+i].id,id)) break;
				if(i<station[nStations-1].nAnt) // wasn't a unique ID
				{
					sprintf(screenLine,"The ID '%s' was used for antenna %d already",id,i+1);
					ErrorPrint(screenLine);
					break;
				}
				
				// Convert to earth-centered coordinates if necessary
				if(coordOrigin=='S') StationToEarth(&x,&y,&z);
				// Check that the coordinates are reasonable
				r=sqrt(x*x+y*y+z*z); if(abs(r-RADIUSEARTH)>rTolerance)
				{
					sprintf(screenLine,"WARNING: The antenna is at an indicated elevation of %lg m",r-RADIUSEARTH);
					ErrorPrint(screenLine);
				}
				// Store an antenna; we already checked (in the station line) that max # antennae was not exceeded
				strcpy(antenna[nTotalAntennae].id,id);
				antenna[nTotalAntennae].x=x; antenna[nTotalAntennae].y=y; antenna[nTotalAntennae].z=z; antenna[nTotalAntennae].wt=wt;
				station[nStations-1].nAnt++; // count an antenna seen for this station
				nTotalAntennae++;
				break;
			default:
				ErrorPrint("Lines in the station list file must start with \"S\" or \"A\"");
		}
				
	}
	
	// If antenna coordinates are relative to the station center, convert them to earth-centered
	if(coordOrigin=='S')
	{
		for (iS=0; iS<nStations; iS++)
		{
			for (iA=iA0; iA<(iA0+station[iS].nAnt); iA++)
			{
				//???d0x+=w*antenna[iA].x; d0y+=w*antenna[iA].y; d0z+=w*antenna[iS].z;
			}
		}
	}
	
	// Open the output beam shape file if wanted
	SimFilename(beamFilename,"BeamShape","txt");
	beamFile=fopen(beamFilename,"w");
	
	// Do the station beam shape calculations for each station
	for (iS=0; iS<nStations; iS++)
	{
		// Calculate the reference station (zero-pad) position as the weighted mean position of all its antennae
 		iA0=station[iS].iAnt;
 		sumWts=d0x=d0y=d0z=0.0;
		for (iA=iA0; iA<(iA0+station[iS].nAnt); iA++)
		{
			w=antenna[iA].wt.re; sumWts+=w;
			d0x+=w*antenna[iA].x; d0y+=w*antenna[iA].y; d0z+=w*antenna[iS].z;
		}
		d0x/=sumWts; d0y/=sumWts; d0z/=sumWts;
		station[iS].x=d0x; station[iS].y=d0y; station[iS].z=d0z; // save for future use
		r=sqrt(station[iS].x*station[iS].x+station[iS].y*station[iS].y+station[iS].z*station[iS].z); // distance of antenna-center of earth
		station[iS].longi=atan2(station[iS].y,station[iS].x); station[iS].lat=asin(station[iS].z/r); // -pi<longitude<=pi  -pi/2<latitude=<pi/2
		fprintf(beamFile,"Station %d:%s  Phase center=(%le,%le,%le)  Longitude=%le  Latitude=%le",
		  iS,station[iS].id,station[iS].x,station[iS].y,station[iS].z,station[iS].longi,station[iS].lat);
		
		// Calculate the hour angle of the center of the beam at the start of observation
		UnpackTime(startObs,&year,&month,&day,&hr,&min,&sec,&msec);
		hour=((msec/1000.+sec)/60.+min)/60.+hr;
		jd=julian_date((short int)year,(short int)month,(short int)day,hour); jd_low=modf(jd,&jd_high); // split Julian date into integer
																										// and fractional parts for sidereal_time
		ee=0.0; // FAKE-need to use EarthTilt to calculate
		sidereal_time(jd_high,jd_low,ee,&gst); // gst <- Greenwich apparent sidereal time, in hours
		H0=gst*HOURStoRADIANS+station[iS].longi-fieldRA; // hour angle in radians
		
		// Form the beam using Shep memo #? eq 13
		k=2.0L*PI*obsFreq/cLight;
		sinH0=sin(H0); cosH0=cos(H0); sinL=sin(station[iS].lat); cosL=cos(station[iS].lat);
		sinDec0=sin(fieldDec); cosDec0=cos(fieldDec);
		HStart=H0-0.1e0; HEnd=H0+0.1e0; HInc=1e-2; // FAKE: 0.2 radian by 0.01=21 values
		decStart=fieldDec-0.1e0; decEnd=fieldDec+0.1e0; decInc=1e-2; // FAKE:  0.2 radian by 0.01=21 values
		for (H=HStart; H<=HEnd; H+=HInc)
		{
			sinH=sin(H); cosH=cos(H);
			for (dec=decStart; dec<=decEnd; dec+=decInc)
			{
				sinDec=sin(dec); cosDec=cos(dec);
				xComp=cosL*(sinDec-sinDec0)-sinL*(cosDec*cosH-cosDec0*cosH0);
				yComp=cosDec*sinH-cosDec0*sinH0;
				zComp=sinL*(sinDec-sinDec0)-cosL*(cosDec*cosH-cosDec0*cosH0);
				B.re=0; B.im=0;
				for (iA=iA0; iA<(iA0+station[iS].nAnt); iA++)
				{
					dx=antenna[iA].x-d0x; dy=antenna[iA].y-d0y; dz=antenna[iS].z-d0z;
					expArg=k*(dx*xComp+dy*yComp+dz*zComp);
					CExpi(-expArg,&e); CMult(antenna[iA].wt,e,&t); CAdd(t,B,&B);
				}
				fprintf(beamFile,"\n%le %le %le %le",H,dec,B.re,B.im);
			}
		}
	}

	
	sprintf(screenLine,"* %d stations and %d antennae processed",nStations,nTotalAntennae); OutputToScreen(screenLine,0);
	SkipLine();
	return 0;
}
				
/*
	Convert coordinates given relative to the center of the station (x=north,
	-y=East, z=elevation above earth surface), on a plane tangent to the earth
	at the station center, to (x,y,z) relative to the center of the earth. The position of the station
	is given by sinLat,cosLat,sinLong,cosLong (east Longitude).
*/


void StationToEarth(double *px, double *py, double *pz)
{
	double x,y,z;
	x= - cosLong*sinLat*(*px) - sinLong*-(*py) + cosLong*cosLat*(*pz) + Sx;
	y= - sinLong*sinLat*(*px) + cosLong*-(*py) + sinLong*cosLat*(*pz) + Sy;
	z=   cosLat        *(*px)                 + sinLat        *(*pz) + Sz;
	*px=x; *py=y; *pz=z;
	return;
}
