// SinesAndTangents.c: Calculates the sines and tangents needed for an FFT

/*
Translated from James S. Walker "Fast Fourier Transforms" (2nd ed), procedure SINESTANS (p 378).

Uses Buneman's method to calculate sin(2pi*m/N) recursively from the identity

	sin theta = [sin(theta+phi) + sin(theta-phi)]/(2*cos phi)

starting with sin(0)=0 and sin(N/4)=1, phi=pi/4
*/

#include "FFTsubs.h"
#include <math.h>

void SinesAndTangents(int N, int R, double sines[], double tangents[], double *pcosPhi)
{
	int N2,H,H9,L,N4,D,j,k; double cosPhi;
	cosPhi = 0.0L; N2 = N/2; H = N2/2; sines[0]=0.0L; sines[H]=1.0L; D = 1;
	for (j = 0; j<=R - 3; j++)
	{
		cosPhi = sqrt(2.0L + cosPhi); H9 = H; H = H/2; L = H;
		for (k = 1; k<=D; k++)
		{
			sines[L] = (sines[L + H] + sines[L - H])/cosPhi;
			L = H9 + L;
		}
		D<<=1;
	}
	N4 = N/4;
	tangents[0] = 0.0; tangents[N4] = 1.0;
	for (j = 1; j<N4; j++)
	{
		tangents[j] = (1.0L - sines[N4 - j])/sines[j];
	}
	*pcosPhi=cosPhi;
	return;
}
