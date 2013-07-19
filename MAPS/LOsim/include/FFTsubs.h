// FFTSubs.h: Function prototypes for the Walker FFT routines

int BitReverse(double F[], double G[], int rB);
void SinesAndTangents(int nSt, int rSt, double sines[], double tangents[], double *pZ);
void Butterfly(double fReal[], double fImag[], double sine[], double tangent[], int nB, int stepSize);
