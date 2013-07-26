// uvCoverage.c: Compute the values of u,v as functions of time and baseline

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/23/01	Create

/*	Inputs:		nTimes		calculate and optionally return (u,v) for a start time 'time0' and that time + k*increment, k=1,...,nTimes-1
				time0		starting time, in same format as startObs
				timeInc		time increment, in seconds
				iStation0,1	<0: al possible pairs of stations
							>=0: single baseline only
				uOut,vOut	if non-NULL: addresses of (2*nSpectralPoints)-element arrays of doubles to receive results
							if NULL, results are printed and plotted
				
	Outputs:	uOut, vOut[i][j], i=0,1, ..., nTimes-1; j=0, 1, ..., nSpectralPoints-1
	
*/
#include "float.h"
#include "uvCoverage.h"
#include "GlobalReferences.h"
#include "Parameters.h"
#include <stdio.h> // must precede Utilities.h
#include "PlotWindow.h"
#include "Utilities.h"
#include "julian_date.h"
#include "sidereal_time.h"
#include <math.h>

void uvCoverage(int nTimes, long long time0, double timeInc, int iStation0, int iStation1, double uOut[], double vOut[])
{

/*
Reference: "Synthesis Imaging in Radio Astronomy II", ed G B Taylor, C L Carilli & R A Perley, eqn 2-30

	[u]     1    [      sin H           cos H          0   ] [Lx]
	[v] = ------ [ -sin dec cos H   sin dec sin H   cos dec] [Ly]
	[w]   lamdba [  cos dec cos H  -cos dec sin H   sin dec] [Lz]

dec = declination of the center of the field of view
H = angular distance from the Greenwich meridian to the right ascension of the source (hour angle) = GMST - field center right ascension
Lx, Ly, Lz = coordinate differences for stations comprising a baseline
lambda = wavelength corresponding to center frequency of this frequency channel (spectral point)
*/

	double Lx,Ly,Lz,H,lambda[MAX_SPECTRAL_POINTS],channelWidth,obsFreqLo,obsFreqHi,freqLo,freqHi,freq0,jd,jd_high,jd_low,ee,gst,hour;
	double u,v,w, a11,a12,a13,a21,a22,a23,a31,a32,a33, sinH,cosH,sinDec,cosDec,plotX[10],plotY[10],uMax,vMax,r,s,t;
	int year,month,day,hr,min,sec,msec,iTime,iS,jS,iF,iColor,printing,plotting,iOut,s00,s01,ds,s11;
	
	printing=uOut==NULL; plotting=iStation0<0; // FAKE
	
	// FOLLOWING SHOULD ALL BE PRE-CALCULATED
	channelWidth=obsBandwidth/nSpectralPoints; obsFreqLo=obsFreq-(obsBandwidth/2.); obsFreqHi=obsFreq+(obsBandwidth/2.);

	// Terms that remain constant over the entire observation period
	sinDec=sin(fieldDec); cosDec=cos(fieldDec);
	
	if(plotting)
	{
		// Get the bounding box for all the (u,v) tracks to be able to scale the plot
		// Each ellipse has its center at (0,Lz*cos(dec)/lambda) and "u" axis sqrt(Lx^2+Ly^2)/lambda, "v" axis sqrt(Lx^2+Ly^2)/(sin(dec)*lambda)
		// Reference: TCP eq 2-31 and fig 2-12
		uMax=vMax=-DBL_MAX;
		for (iS=0; iS<nStations; iS++)
		{
			for (jS=iS+1; jS<nStations; jS++)
			{
				Lx=station[iS].x-station[jS].x; Ly=station[iS].y-station[jS].y; Lz=station[iS].z-station[jS].z;
				r=Lz*cosDec; s=sqrt(Lx*Lx+Ly*Ly); t=s/sinDec;
				if(s>uMax) uMax=s; if((fabs(r+t))>vMax) vMax=fabs(r+t);
			}
		}
		plotX[0]=uMax; plotY[0]=vMax; PlotWindow(PW_INITIALIZE,0,plotX,plotY,0);
	}	
	
	UnpackTime(time0,&year,&month,&day,&hr,&min,&sec,&msec);
	
	// Loop on calendar observation time
	for(iTime=0; iTime<nTimes; iTime++)
	{
		if(printing) printf("\n%4d-%02d-%02d %02d:%02d:%02d.%03d",year,month,day,hr,min,sec,msec);
		
		 // Compute the hour angle at this observation time
		hour=((msec/1000.+sec)/60.+min)/60.+hr;
		jd=julian_date((short int)year,(short int)month,(short int)day,hour); jd_low=modf(jd,&jd_high); // split Julian date into integer
																										// and fractional parts for sidereal_time
	
		ee=0.0; // FAKE-need to use EarthTilt to calculate
		sidereal_time(jd_high,jd_low,ee,&gst); // gst <- Greenwich apparent sidereal time, in hours
		H=gst*HOURStoRADIANS-fieldRA; // Greenwich hour angle in radians
//NEED HOUR ANGLE AT STATION!
	    if(printing) printf("  %lg",H);
	
		// Some trigonometric values we'll need
		sinH=sin(H); cosH=cos(H);
	
		// Tabulate terms of the matrix (these don't vary with the station)
		a11=sinH;         a12=cosH;         a13=0.0;
		a21=-sinDec*cosH; a22=sinDec*sinH;  a23=cosDec;
		a31=cosDec*cosH;  a32=-cosDec*sinH; a33=sinDec;
	
		// Calculate all the wavelengths we'll use
		for (iF=0; iF<nSpectralPoints; iF++)
		{
			freqLo=obsFreqLo+(iF*channelWidth); freqHi=obsFreqLo+((iF+1)*channelWidth); // frequency limits and center for this frequency channel
			freq0=(freqLo+freqHi)/2.; lambda[iF]=cLight/freq0; // center for this frequency channel, and its wavelength
		}
		
		// Finally, calculate the positions in the (u,v) plane for each pair of stations
		iColor=0; iOut=iTime*nSpectralPoints;
		if(iStation0<0) {s00=0;         s01=nStations-1; ds=1;                   s11=nStations;}   // all station pairs
		else            {s00=iStation0; s01=iStation0+1; ds=iStation1-iStation0; s11=iStation1+1;} // single station pair
		for (iS=s00; iS<s01; iS++)
		{
			for (jS=iS+ds; jS<s11; jS++)
			{
				Lx=station[iS].x-station[jS].x; Ly=station[iS].y-station[jS].y; Lz=station[iS].z-station[jS].z;
				if(printing) printf("\n%3d %20s %3d %20s %14g %14g %14g",iS,station[iS].id,jS,station[jS].id,Lx,Ly,Lz);
				if(plotting)
				{
					plotX[0]=(a11*Lx+a12*Ly+a13*Lz); plotX[1]=-plotX[0]; plotY[0]=(a21*Lx+a22*Ly+a23*Lz); plotY[1]=-plotY[0];
					PlotWindow(PW_ADD_POINTS,2,plotX,plotY,(iColor%7)*10); iColor++;
				}
				for (iF=0; iF<nSpectralPoints; iF++)
				{
					u=(a11*Lx+a12*Ly+a13*Lz)/lambda[iF]; v=(a21*Lx+a22*Ly+a23*Lz)/lambda[iF]; w=(a31*Lx+a32*Ly+a33*Lz)/lambda[iF];
					if(printing) printf("  %14g %14g %14g",u,v,w);
					if(uOut!=NULL) {uOut[iOut+iF]=u; vOut[iOut+iF]=v;}
				}
			}
		}
		AddTime(timeInc,&year,&month,&day,&hr,&min,&sec,&msec); // advance to next observation time
	}
	
	return;
}
