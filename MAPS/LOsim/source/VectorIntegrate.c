// VectorIntegrate.c: Vector integrate corrupted visibility data over (u,v,w) region representing sample data

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 9/4/01	Create

/*	This is box 15 on the master diagram

	Inputs:		
				
	Outputs:	
	
*/

#include "VectorIntegrate.h"
#include "GlobalReferences.h"
#include "uvCoverage.h"
#include "Parameters.h"
#include <stdio.h> // must precede Utilities.h
#include "Utilities.h"
#include <math.h>

void VectorIntegrate(int nTimes, long long time0, double timeInc)
{
	int iS,jS,iTime,iF,year,month,day,hr,min,sec,msec,iuMin,ivMin,iuMax,ivMax; long long time;
	double u[2*MAX_SPECTRAL_POINTS]={0.0},v[2*MAX_SPECTRAL_POINTS]={0.0}; // conceptually u,v[2][MAX_SPECTRAL_POINTS]
	double du,dv,uMin,vMin,uMax,vMax;
	
	du=brightnessXN/(visibilityXN*fieldXSize); dv=brightnessYN/(visibilityYN*fieldYSize); // grid spacings in (u,v) plane
	
	// Loop on baselines
	for (iS=0; iS<(nStations-1); iS++)
	{
		for (jS=iS+1; jS<nStations; jS++)
		{
			// Loop on calendar observation time
			for(iTime=0,time=time0; iTime<nTimes; iTime++)
			{
				// Loop on center frequencies
				for (iF=0; iF<nSpectralPoints; iF++)
				{
					uvCoverage(2,time,timeInc,iS,jS,u,v);
					// Get next grid points outside area to be integrated; these are points < smaller and > larger of two endpoints
					if(u[0*nSpectralPoints+iF]<u[1*nSpectralPoints+iF]) {uMin=u[0*nSpectralPoints+iF]; uMax=u[1*nSpectralPoints+iF];}
					else                                                {uMin=u[1*nSpectralPoints+iF]; uMax=u[0*nSpectralPoints+iF];}
					if(v[0*nSpectralPoints+iF]<v[1*nSpectralPoints+iF]) {vMin=v[0*nSpectralPoints+iF]; vMax=v[1*nSpectralPoints+iF];}
					else                                                {vMin=v[1*nSpectralPoints+iF]; vMax=v[0*nSpectralPoints+iF];}
					iuMin=(int)floor(uMin/du); iuMax=(int)ceil(uMax/du); ivMin=(int)floor(vMin/dv); ivMax=(int)ceil(vMax/dv);
					// Array positions from grid positions
				}
				
				// Advance to next observation time
				UnpackTime(time,&year,&month,&day,&hr,&min,&sec,&msec);
				AddTime(timeInc,&year,&month,&day,&hr,&min,&sec,&msec);
				PackTime(year,month,day,hr,min,sec,msec,&time);
			}
		}
	}
	return;
}
