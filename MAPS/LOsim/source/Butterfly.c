//DEFINT A-D, H-N, P-R, Y
//DEFDBL E-G, O, S-X, Z
//DECLARE SUB FTBFLY (fReal#(), fImag#(), SA#(), TA#(), N%, stepSize%)

#include "FFTsubs.h"

void Butterfly(double fReal[], double fImag[], double sine[], double tangent[], int N, int stepSize)
{
	int K,N2, L, M2, M4,Q,D,Q2,D2,A,A1,B1,C,D4,AV,BV; double S,T,T1,T2,T3,T4,V,SL,TL;
	int nMult=0;
	N2=N/2; L=0;  M2=(stepSize*N)/2; M4=M2/2;
	for (K=0; K<N2; K++)
	{
		T1=fReal[L]+fReal[L+1];  T2=fImag[L]+fImag[L+1];
		fReal[L+1]=fReal[L]-fReal[L+1];  fImag[L+1]=fImag[L]-fImag[L+1];
		fReal[L]=T1; fImag[L]=T2; L=L+2;
	}
	Q=stepSize*N2; D=2; Q2=N2;
	while (D<N)
	{
		Q2=Q2/2; Q=Q/2; D2=D; D=D2+D2; A1=0; L=0; A=A1; D4=D2/2;
		for (C=1; C<=Q2; C++)
		{
			B1=A+D2;
			T=fReal[B1]; S=fImag[B1]; fReal[B1]=fReal[A]-T; fImag[B1]=fImag[A]-S;
			fReal[A]=fReal[A]+T; fImag[A]=fImag[A]+S;
			AV=A+D4; BV=B1+D4;
			T=fReal[BV]; S=fImag[BV]; fReal[BV]=fReal[AV]+S; fImag[BV]=fImag[AV]-T;
			fReal[AV]=fReal[AV]-S; fImag[AV]=fImag[AV]+T; A=A+D;
		}
		
		L=L+Q; A1=A1+1; A=A1;
		for (K=2; K<=D4; K++)
		{
			SL=sine[L]; TL=tangent[L];
			for (C=1; C<=Q2; C++)
			{
				B1=A+D2; V=fImag[B1]+TL*fReal[B1]; T3=fReal[B1]-V*SL; nMult+=2;
				T4=T3*TL+V; fReal[B1]=fReal[A]-T3; fImag[B1]=fImag[A]-T4; nMult+=1;
				fReal[A]=fReal[A]+T3; fImag[A]=fImag[A]+T4; AV=A+D4;
				BV=B1+D4; V=fReal[BV]-TL*fImag[BV]; T3=-fImag[BV]-V*SL; nMult+=2;
				T4=T3*TL+V; fReal[BV]=fReal[AV]-T3; fImag[BV]=fImag[AV]-T4; nMult+=1;
				fReal[AV]=fReal[AV]+T3; fImag[AV]=fImag[AV]+T4; A=A+D;
			}
			L=L+Q; A1=A1+1; A=A1;
		}
	}
	return;
}
