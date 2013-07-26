#include <string>
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "powspec_fitsio.h"
#include "fitshandle.h"

#include "MapInpaint.h"
#include "MapAddSmallScales.h"

using namespace std;


int main (int argc, const char **argv)
{
	int out_nside = atoi(argv[1]);
	int llow = atoi(argv[2]);
	int lhigh = atoi(argv[3]);
	int seed = atoi(argv[4]);
  int inpaint_niter = atoi(argv[5]);
	string infile_map(argv[6]);
  string infile_mask(argv[7]);
	string outfile_painted(argv[8]);
  string outfile_final(argv[9]);
  fitshandle out;

  printf("-------------------------------\n");
  printf("            GALFIX\n");
  printf("-------------------------------\n");
  printf("out_nside: %d\nseed: %d\ninpaint_niter: %d\nllow, lhigh: %d, %d\n", out_nside, seed, inpaint_niter, llow, lhigh);
  printf("infile_map: %s\n", infile_map.c_str());
  printf("infile_mask: %s\n", infile_mask.c_str());
  printf("outfile_painted: %s\n", outfile_painted.c_str());
  printf("outfile_final: %s\n", outfile_final.c_str());

  // -------------------------------------------------------------------
	//
	//	Read in the large-scale structure and fill in masked out regions
	//
	// -------------------------------------------------------------------
  MapInpaint painter;
  painter.Load(infile_map, infile_mask);
  painter.Inpaint(inpaint_niter, painter.Map().Nside()*2);

  cout << "galfix: writing inpainted map to fits" << endl;
  out.create(outfile_painted.c_str());
  write_Healpix_map_to_fits(out, painter.Painted(), FITSUTIL<double>::DTYPE);
	out.close();

  // For diagnostics, let's also save the difference between the input map
  // and the output inpainted map.
  Healpix_Map<double> diff(painter.Map().Nside(), RING, SET_NSIDE);  
  for (uint j=0; j<(uint)painter.Map().Npix(); j++)
  {
    diff[j] = painter.Painted()[j] - painter.Map()[j];
  }
  out.create("diff.fits");
  write_Healpix_map_to_fits(out, diff, FITSUTIL<double>::DTYPE);
	out.close();

  // And the input map x mask
  for (uint j=0; j<(uint)painter.Map().Npix(); j++)
  {
    diff[j] = painter.Map()[j]*painter.Mask()[j];
  }
  out.create("masked.fits");
  write_Healpix_map_to_fits(out, diff, FITSUTIL<double>::DTYPE);
	out.close();

	// -------------------------------------------------------------------
	//
	//	Extrapolate power-law to lmax and replace alms with random values
	//  above lhigh
	//
	// -------------------------------------------------------------------
	
  MapAddSmallScales extender;
  extender.Load(outfile_painted);
  extender.Fit(llow, lhigh);
  extender.Replace(lhigh, out_nside, seed);
 
  cout << "galfix: writing final map to fits" << endl;
  out.create(outfile_final.c_str());
  write_Healpix_map_to_fits(out, extender.MapExtended(), FITSUTIL<double>::DTYPE);
	out.close();

  out.create(string("initial_pow.fits").c_str());
  write_powspec_to_fits(out, extender.Cls(), extender.Cls().Num_specs());
	out.close();

  out.create(string("final_pow.fits").c_str());
  write_powspec_to_fits(out, extender.ClsExtended(), extender.ClsExtended().Num_specs());
  out.close();

	return(0);
}
