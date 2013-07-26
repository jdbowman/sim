#ifndef ANGULAR_CORRELATION_TOOLS_H
#define ANGULAR_CORRELATION_TOOLS_H

#include <gsl/gsl_spline.h>

#include "healpix_map_fitsio.h"
#include "alm.h"
#include "healpix_map.h"
#include "alm_map_tools.h"

void angcor2alm(double theta[], double w[], int size, int nside, Alm<xcomplex<double> > &alm)
{
  // This routine gets lmax and mmax from the passed alm.
  // It uses a dumb method to transform w(theta) into C_l by populating
  // a map with the w(theta) values then calling map2alm to decompose
  // it into a_lm.
  // NOTE: This routine assumes w(theta) is given where theta is in
  // units of degrees. And that the theta,w pairs in the rows of the
  // data file are log10(theta), log10(w)
  
  Healpix_Map<double> map(nside,RING,SET_NSIDE);
  uint j=0;
  uint niter=0;
  double mintheta = pow(10, theta[0]);
  double maxtheta = pow(10, theta[size-1]);

  // Do a quick sanity check on the angular correlation function
  if (size < 2)
    printf("Warning Angular correlation function only contains %d entries\n", size);

  if (mintheta > (180 / 3.14 * sqrt(4.0*3.14/map.Npix())))
    printf("Warning: Angular correlation function is only provided for theta > %g\n", mintheta);

  if (maxtheta < 180)
    printf("Warning: Angular correlation function is only provided for theta < %g\n", maxtheta);
  
  // Setup the GSL interpolation
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_linear, size);
	gsl_spline_init(spline, theta, w, size);

  // Loop over the pixels and populate with w(theta)
  for (j=0; j<(uint)map.Npix(); j++)
  {
    if ((map.pix2ang(j).theta >= mintheta) & (map.pix2ang(j).theta <= maxtheta)) 
    {
      map[j] = pow(10, gsl_spline_eval(spline, log10(180 / 3.141592653 * map.pix2ang(j).theta), acc));
    }
    else 
    {
      map[j] = 0;
    }
  }

  // Free the GSL interpolation 
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  // Transform to a_lm
  arr<double> ring_weights(2*map.Nside());
  ring_weights.fill(1.0);
  map2alm_iter(map, alm, niter, ring_weights);

//  fitshandle out;
//  out.create("myangcormap.fits");
//  write_Healpix_map_to_fits(out, map, FITSUTIL<double>::DTYPE);
//	out.close();


}




void read_xydata_from_file(const string &strFilename, vector<double> &x, vector<double>& y)
{

    // This routine reads a simple text file containing two columns
    // and returns two vectors corresponding to x-y pairings of each
    // row.  No header lines or comments are allowed in the file.
		FILE* fid=NULL;
		double xi=0;
		double yi=0;

		// Read flux CDF from file
		fid = fopen(strFilename.c_str(), "r");
		while ( fscanf(fid, "%lg, %lg", &xi, &yi) == 2)
		{
			x.push_back(xi);
			y.push_back(yi);
		}
		fclose(fid);
}






#endif
