// DFT.c: Calculate the discrete Fourier transform of a complex sequence

#include "DFT.h"
#include "Parameters.h"
#include <math.h>

int DFT(int N, struct Complex h[], struct Complex H[])
{
	struct Complex W[1<<15],w; double x; int j,k,n;
	
	if(N>(1<<15)) return 1;
	
	// Calculate the N possible values of W^nk 
	// using x=2*pi*j/N, exp{ix}=cos x + i*sin x
	for (j=0; j<N; j++)
	{
		x=2.0L*PI*j/N; W[j].re=cos(x); W[j].im=sin(x);
	}
	
	//
	for (n=0; n<N; n++)
	{
		H[n].re=H[n].im=0.0L;
		for (k=0; k<N; k++)
		{
			j=(n*k)%N; w=W[j];
			H[n].re+=w.re*h[k].re-w.im*h[k].im;
			H[n].im+=w.re*h[k].im+w.im*h[k].re;
		}
	}
	return 0;
}
