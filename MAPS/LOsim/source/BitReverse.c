// BitReverse.c: Puts the arrays F and G into bit-reversed order

// Author:	James S. Walker, in "Fast Fourier Transforms" (2nd edition), CRC Press 1996.
// Adapted from procedure BITREV (appendix B, p 377) by	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 6/01		translated from QuickBASIC
// 9/14/01	correct scope of R mod 2 IF

#define DEBUG 0 // set to 1 to test correctness of this function

#include <math.h>
#if DEBUG
	#include <stdio.h>
#endif
#include "FFTsubs.h"

#define MAX_HALFBITLENGTH 8
#define MAX_SQRTSIZE (1<<MAX_HALFBITLENGTH)

// Returns 0 if no error, 1 if R out of range
int BitReverse(double F[], double G[], int R)
{
	long M; int i,bitRev[MAX_SQRTSIZE],N,N1,N9,C,R2,M7,L,L9,bitLength,N2,N6,N8,K,M6; double T,Y;
#define DEBUG 0 // set to 1 to test correctness of this function
#if DEBUG
	long long b1,b2,b2rev,Binary(int x, int reverse); int bad=0; char check[65536]={0}; // debug stuff
#endif
	if(R<0 || R>(2*MAX_HALFBITLENGTH)) return 1;
	M=1<<R;
	N2 = M / 2; N = M; M7 = 0; R2 = R / 2; N1 = N; Y = 0; C = 1;
	if(R%2==1) {N1 = N2; R2 = (R - 1) / 2; C = 0;}
	N9 = sqrt(N1);
	if(R%2==1) N8 = N9 + N9; else N8 = N9;
	
	
	/* Create the square-root-sized array bitRev by Buneman's method
	   If bitRev[0, ..., N-1] are the bit-reversed integers 0, ..., N-1, then bitRev[0, ..., 2N-1] can be computed by
	   putting 2*bitRev[i] in the first N places, and 2*bitRev[i]+1 in the last N.
	   Example: for N=4        i   0 1 2 3
	                    bitRev[i]  0 2 1 3            bit-reversed 4-bit integers
	   so for N=8              i   0 1 2 3    4 5 6 7
                        bitRev[i]  0 2 4 6    1 3 5 7 bit-reversed 8-bit integers
	
	*/
	bitRev[0]=0; bitRev[1]=1; L9 = 2;
	for (bitLength=2; bitLength<=R2; bitLength++)
	{
		for (i=0; i<L9; i++)
		{
			bitRev[i]<<=1; bitRev[i+L9]=bitRev[i]+1;
		}
		L9<<=1;
	}
	
	while(C < 2)
	{
		for (L = 1; L<N9; L++)
		{
			M7 = M7 + N8; N6 = bitRev[L] + Y;
			T = G[M7]; G[M7] = G[N6]; G[N6] = T;
			T = F[M7]; F[M7] = F[N6]; F[N6] = T;
#if DEBUG
			b1=Binary(M7,0); b2=Binary(N6,0); b2rev=Binary(N6,R); printf("\n%X=%lld <-> %X=%lld %d",M7,b1,N6,b2,b1==b2rev);
			bad+=(b1!=b2rev); if(check[M7] || check[N6]) bad++; check[M7]=1; check[N6]=1;
#endif
			for (K=1; K<L; K++)
			{
				M6 = M7 + bitRev[K]; N6 = N6 + N8;
				T = G[M6]; G[M6] = G[N6]; G[N6] = T;
				T = F[M6]; F[M6] = F[N6]; F[N6] = T;
#if DEBUG
				b1=Binary(M6,0); b2=Binary(N6,0); b2rev=Binary(N6,R); printf("\n%X=%lld <-> %X=%lld %d",M6,b1,N6,b2,b1==b2rev);
				bad+=(b1!=b2rev); if(check[M6] || check[N6]) bad++; check[M6]=1; check[N6]=1;
#endif
			}
		}
		Y = N9; M7 = N9; C = C + 1;
	}
#if DEBUG
	for (i=0; i<M; i++)
	{
		if(!check[i] && (Binary(i,0)!=Binary(i,R))) bad++;
	}
	if(bad) printf("\n** ERROR in BitReverse (count=%d)",bad);
#endif
	return 0;
}

// reverse=bit length of x (used only if reversed)
long long Binary(int x, int reverse)
{
	long long b=0,p=1; int t=x,d,i=0;
	while(t!=0)
	{
		d=t&1; t>>=1; // get ls bit
		if(reverse) {b*=10; b+=d; i++;}
		else {b=b+(p*d); p*=10;}
	}
	for(; i<reverse; i++) b*=10;
	return b;
}
