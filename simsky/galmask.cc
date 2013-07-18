#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"

using namespace std;


int main (int argc, const char **argv)
{
	int out_nside = atoi(argv[1]);
  int factor = atoi(argv[2]);
  string outfile(argv[3]);
  fitshandle out;

  printf("-------------------------------\n");
  printf("            GALMASK\n");
  printf("-------------------------------\n");
  printf("out_nside: %d\n", out_nside);
  printf("outfile: %s\n", outfile.c_str());

	// Setup the GSL random number generator
	gsl_rng_env_setup();
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, 1);

	// Allocate an empty map
	Healpix_Map<double> map(out_nside/factor,RING,SET_NSIDE);
  
	// Populate map with projected density drawn from Poisson distribution
  double Nsources = 2000.0;
	double threshold = Nsources / map.Npix();
  printf("Nsources=%g,  Npix=%d,  threshold=%g\n", Nsources, map.Npix(), threshold);

	for (uint j=0; j<(uint)map.Npix(); j++)
	{
		map[j] = gsl_rng_uniform(rng) < threshold ? 0.0 : 1.0;
	}

  gsl_rng_free(rng);

  Healpix_Map<double> out_map(out_nside, RING, SET_NSIDE);
  out_map.Import_upgrade(map);

  cout << "galmask: writing mask to fits" << endl;
  out.create(outfile.c_str());
  write_Healpix_map_to_fits(out, out_map, FITSUTIL<double>::DTYPE);
	out.close();

	return(0);
}
