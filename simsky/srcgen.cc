#include "xcomplex.h"
#include "cxxutils.h"
#include "healpix_data_io.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "alm_map_tools.h"
#include "fitshandle.h"
#include "powspec.h"
#include "powspec_fitsio.h"
#include "SourceCatalog.h"
#include "SourceConfusionMap.h"
#include "SourceModel.h"
#include "angular_correlation_tools.h"

extern "C" {
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
}


using namespace std;

 
void map2alm_iter3 (const Healpix_Map<double> &map, Alm<xcomplex<double> > &alm, double err_rel)
{
  arr<double> wgt(2*map.Nside());
  wgt.fill(1);
  Healpix_Map<double> map2(map);
  alm.SetToZero();
  while(true)
  {
    map2alm(map2,alm,wgt,true);
    alm2map(alm,map2);
    double errmeasure=0;
    for (int m=0; m<map.Npix(); ++m)
    {
      if (map[m] != 0) errmeasure = max(errmeasure,abs(map[m]-map2[m])/map[m]);
      map2[m] = map[m]-map2[m];
    }
    cout << "map error measure: " << errmeasure << endl;
    if (errmeasure<err_rel) break;
  }
}


int main (int argc, const char **argv)
{
  // Input parameters
  int nside = atoi(argv[1]);
  double nsources = atof(argv[2]);
  int seed = atoi(argv[3]);
  double fluxcut = atof(argv[4]);  // Jy
  string prefix(argv[5]);

  // Initialization
  int j = 0;
  int k = 0;
  int lmax = 4*nside;
  int mmax = 4*nside;
  int niter = 0;
  int nspec = 2;
  double center_freq = 150; // MHz
  double freqmin = 100; // MHz
  double freqmax = 200; // MHz

  // Require input data
  string strFluxCDF("/home/piyanat/local/src/simsky/data/src_flux_cdf.dat");
  string strIndexCDF("/home/piyanat/local/src/simsky/data/src_index_cdf.dat");
  string strAngCorFile("/home/piyanat/local/src/simsky/data/src_angcor.dat");

  // Output files
  string strCatalogFile(prefix+"_catalog.dat");
  string strConfusionBase(prefix+"_confusion_map_");
  string outfile_map(prefix+"_map.fits");
  string outfile_pow(prefix+"_pow.fits");
  string basemap(prefix+"_basmap.fits");
  
  printf("-------------------------------\n");
  printf("            SRCGEN\n");
  printf("-------------------------------\n");
  printf("nside: %d\nnsources: %g\nseed: %d\n", nside, nsources, seed);
  printf("outfile_map: %s\n", outfile_map.c_str());
  printf("outfile_pow: %s\n", outfile_pow.c_str());

  // Declare our fitsio writer helper
  fitshandle out;
 
	// -------------------------------------------------------------------
	//
	//	Populate the sky with sources according to Poisson statistics
	//
	// -------------------------------------------------------------------
  
	// Setup the GSL random number generator
	gsl_rng_env_setup();
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);

	// Allocate an empty map
	Healpix_Map<double> map(nside,RING,SET_NSIDE);
	
	// Populate map with projected density drawn from Poisson distribution
	double pixmean = nsources / map.Npix();
	for (j=0; j<map.Npix(); j++)
	{
		map[j] = (double) gsl_ran_poisson(rng, pixmean)/pixmean - 1.0;
	}

/*
  Alm<xcomplex<double> > almT(lmax, mmax);
  for (int l=0; l<=lmax; ++l)
  {
    almT(l,0).re=gsl_ran_gaussian(rng, 1.0); 
    almT(l,0).im=0.0;
  }
	
  for (int m=1; m<=mmax; ++m)
  {
    for (int l=m; l<=lmax; ++l)
    {
      almT(l,m).re=gsl_ran_gaussian(rng, 1.0); 
      almT(l,m).im=gsl_ran_gaussian(rng, 1.0); 
    }
  }

  alm2map(almT, map);
*/

	// Clean up the GSL random number generator
	gsl_rng_free(rng);

  // Print out some sanity check that the map was populated
  printf("map[100]: %g\n", map[100]);
  out.create(basemap);
  write_Healpix_map_to_fits(out, map, FITSUTIL<double>::DTYPE);
	out.close();

  // print out some statistics about our final source count map
  double pmin=0, pmax=0;
  map.minmax(pmin, pmax);
  printf("Poisson-only source count map...\n");
  printf("--Average pixel: %g\n", map.average());
  printf("--Min/max pixel: %g/%g\n", pmin, pmax);
  printf("--RMS: %g\n", map.rms());

	// -------------------------------------------------------------------
	//
	//	Modify this intial map to add clustering
	//
	// -------------------------------------------------------------------

	// Get the Alms and power spectrum of the initial map
  arr<double> ring_weights(2*map.Nside());
  ring_weights.fill(1.0);

  Alm<xcomplex<double> > alm(lmax, mmax);

  double avg = map.average();
  map.add(-avg);
  map2alm_iter(map, alm, 3, ring_weights);
  alm(0,0) += avg*sqrt(4.0*3.141592653);
  map.add(avg);

  PowSpec pow_poisson(1, lmax);
  extract_powspec(alm, pow_poisson);

  printf("pow_poisson.tt(10): %g\n", pow_poisson.tt(10));
  printf("pow_poisson.tt(200): %g\n", pow_poisson.tt(200));
  printf("pow_poisson.tt(300): %g\n", pow_poisson.tt(300));
  printf("pow_poisson.tt(400): %g\n", pow_poisson.tt(400));
  printf("pow_poisson.tt(500): %g\n", pow_poisson.tt(500));

	// Get the clustering power spectrum
  vector<double> w_x;
  vector<double> w_y;
  read_xydata_from_file(strAngCorFile, w_x, w_y);
  Alm<xcomplex<double> > w_alm(lmax, mmax);
  angcor2alm(&w_x[0], &w_y[0], w_x.size(), map.Nside(), w_alm);

  PowSpec pow_cluster(1, lmax);
	extract_powspec(w_alm, pow_cluster);

  arr<double> myarr(lmax, 0);
  pow_cluster.Set(myarr);

  printf("pow_cluster.tt(10): %g\n", pow_cluster.tt(10));
  printf("pow_cluster.tt(200): %g\n", pow_cluster.tt(200));
  printf("pow_cluster.tt(300): %g\n", pow_cluster.tt(300));
  printf("pow_cluster.tt(400): %g\n", pow_cluster.tt(400));
  printf("pow_cluster.tt(500): %g\n", pow_cluster.tt(500));

	// Define the combined power spectrum:
	// pow_final = sqrt(pow_cluster/pow_poisson + 1)
	PowSpec pow_final(1, lmax);
	for (j=0; j<=lmax; j++)
	{
		pow_final.tt(j) = sqrt(pow_cluster.tt(j) / pow_poisson.tt(j) + 1);
	}
  pow_final.tt(0) = 1.0;

  printf("pow_final.tt(10): %g\n", pow_final.tt(10));
  printf("pow_final.tt(200): %g\n", pow_final.tt(200));
  printf("pow_final.tt(300): %g\n", pow_final.tt(300));
  printf("pow_final.tt(400): %g\n", pow_final.tt(400));
  printf("pow_final.tt(500): %g\n", pow_final.tt(500));

	// Scale the alms
	//alm.ScaleL(pow_final.tt());

	// Convert back into a map
  double offset = alm(0,0).real()/sqrt(4.0*3.141592653);
  alm(0,0) = 0;
  alm2map(alm, map);
  map.add(offset);

  // print out some statistics about our clustered density map
  map.minmax(pmin, pmax);
  printf("Poisson+clustering source density map...\n");
  printf("--Average pixel: %g\n", map.average());
  printf("--Min/max pixel: %g/%g\n", pmin, pmax);
  printf("--RMS: %g\n", map.rms());
  

  // Undo the density normalization to get back to the number of
  // sources per pixel (and round to nearest integer)
	for (j=0; j<map.Npix(); j++)
	{
		map[j] = floor((map[j]+1.0)*pixmean + 0.5);
    if (map[j] < 0)
      map[j] = 0;
	}

	// Save the map
  out.create(outfile_map.c_str());
  write_Healpix_map_to_fits(out, map, FITSUTIL<double>::DTYPE);
	out.close();

	// Save the final power spectrum
	out.create(outfile_pow.c_str());
	write_powspec_to_fits(out, pow_final, 1);
	out.close();

	// -------------------------------------------------------------------
	//
	//	Run through each pixel and randomly populate with sources drawn
  //  from the appropriate flux distribution
	//
	// -------------------------------------------------------------------
  SourceCatalog catalog;
  SourceConfusionMap confusion(freqmin, freqmax, center_freq, nside, nspec);
  SourceModel model(seed);
  model.LoadCDFs(strFluxCDF, strIndexCDF);
  double flux = 0;
  double index = 0;
  
  printf("\nPopulating with sources...\n\n");

  // Loop over all map pixels
  for (j=0; j<map.Npix(); j++)
  {
    if (j%10000 == 0)
      printf("\rCurrently on pixel: %d of %d   ", j, map.Npix());
      cout << flush;

    pointing ang = map.pix2ang(j);

    // Loop over all sources in pixel
    for (k=0; k<map[j]; k++)
    {
      // Randomly select flux and spectral index
      model.DrawSource(flux, index);

      if (flux > fluxcut)
      {
        // If it is a bright source, add it to the catalog
        catalog.Add(ang.theta, ang.phi, flux, index);
      }
      else
      {
        // If it is a faint source, add it to the diffuse map
        confusion.Add(j, flux, index);

      }
    }
  }

	catalog.Randomize(sqrt(4.0 * 3.141592653 / map.Npix()), seed);
	catalog.WriteToFile(strCatalogFile);

	confusion.ConvertToLogPoly();
	confusion.WriteToFiles(strConfusionBase);



	return(0);
}
