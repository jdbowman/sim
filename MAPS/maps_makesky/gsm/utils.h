#ifndef __UTILS_H
#define __UTILS_H

#define BOLT 1.3806503E-23
#define SOL 299792458
#define JANSKY 1E26

#ifdef __cplusplus
extern "C" {
#endif
	long npix2nside ( const long npix );
	void pix2ang_ring( long nside, long ipix, double *theta, double *phi);
	int hpx2carre (const char *hpxname, const char *outname,float freq);
	int healpix2carre (float *signal, long nside, const char *filename, int nest,float freq);
	int dat2carre (float *signal, long nside, const char *filename,float freq);
	
	int hpxwrap (float *signal, long nside, const char *filename, int nest,float freq);

	double convert_to_Jy(double temp,double frequency,double area);
	double convert_to_Iv(double temp,double frequency,double area);
	double get_pix_area(double dx,double dy, double dec);


/* from simple poly */
	int rectClip(int,double *,double *,double *,double *,double,double,double,double);
	double Bourke(int nv, double *xval, double *yval) ;

#ifdef __cplusplus
}
#endif

#endif


