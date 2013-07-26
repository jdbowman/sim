#ifndef MAPADDSMALLSCALES_H
#define MAPADDSMALLSCALES_H

#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>


#include "alm.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "powspec.h"

using namespace std;


/*
	This class


*/
class MapAddSmallScales
  {
  private:

    Healpix_Map<double> m_map_before;
    Healpix_Map<double> m_map_after;
    Alm<xcomplex<double> > m_alm_before;
    Alm<xcomplex<double> > m_alm_after;
    PowSpec m_pow_before;
    PowSpec m_pow_after;
    double m_c0;
    double m_c1;

  public:


    MapAddSmallScales() : m_c0(0), m_c1(0)
    {
    }


    ~MapAddSmallScales()
    {
    }

    void Load(const string& strMapFile)
    {
      cout << "MapAddSmallScales: Loading..." << endl;

      double dmin=0, dmax=0;

      // Read map
      read_Healpix_map_from_fits(strMapFile.c_str(), m_map_before);
      if (m_map_before.Scheme() != RING)
      {
        cout << "--swapping map scheme to RING" << endl;
        m_map_before.swap_scheme();
      }
      m_map_before.Map().minmax(dmin, dmax);
      cout << "--map [min, max]=" << dmin << ", " << dmax << endl;
    }


    Healpix_Map<double> & Map() { return m_map_before; }

    Healpix_Map<double> & MapExtended() { return m_map_after; }

    PowSpec & Cls() { return m_pow_before; }
    
    PowSpec & ClsExtended() { return m_pow_after; }


    void Fit(uint fitlow, uint fithigh)
    {
      // Fit a power-law profile to the angular structure
      cout << "MapAddSmallScales: Fitting..." << endl;

      // Convert the input map to alms and Cls
      cout << "--calculating alms" << endl;
	    arr<double> ring_weights(2*m_map_before.Nside());
      ring_weights.fill(1.0);

      uint lmax = m_map_before.Nside()*2;
      m_alm_before.Set(lmax, lmax);

      double avg = m_map_before.average();
      m_map_before.add(-avg);
	    map2alm_iter(m_map_before, m_alm_before, 3, ring_weights);
      m_alm_before(0,0) += avg*sqrt(4.0*3.141592653);
      m_map_before.add(avg);

	    extract_powspec(m_alm_before, m_pow_before);

      // Fit the power-law over the specified range
      cout << "--fitting power law between l=[" << fitlow << "-" << fithigh << "]" << endl;
      uint l=0;
      double cov00=0, cov01=0, cov11=0, sumsq=0;
      uint len = fithigh - fitlow + 1;
      double x[len];
      double y[len]; 
       
      for (l=0; l<len; l++)
      {
        x[l] = log10((double) (fitlow + l));
        y[l] = log10(m_pow_before.tt(fitlow + l));
      }
      gsl_fit_linear (x, 1, y, 1, len, &m_c0, &m_c1, &cov00, &cov01, &cov11, &sumsq);

    }


    void Replace(uint replaceabove, uint out_nside, int seed)
    {
      cout << "MapAddSmallScales: Replacing..." << endl;

      uint lmax = out_nside*2;
      uint l=0;
      uint m=0;

      // Extrapolate the fit to smaller scales
      cout << "--extrapolating fit to small scales" << endl;
      PowSpec newpow(1, lmax);
      for (l=0; l<replaceabove; l++)
      {
        newpow.tt(l) = m_pow_before.tt(l);
      }
      for (l=replaceabove; l<=lmax; l++)
      {
        newpow.tt(l) = pow(10.0, m_c0 + m_c1*log10((double) l));
      }

      // Replace the alms for the new small-scale structure
      cout << "--generating new alms" << endl;
      gsl_rng_env_setup();
      gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(rng, seed);

      double hsqrt2 = 1/sqrt(2.0);

      m_alm_after.Set(lmax, lmax);
      for (l=0; l<replaceabove; l++)
      {
        for (m=0; m<=l; m++)
        {
          m_alm_after(l,m) = m_alm_before(l,m);
        }
      }
      for (l=replaceabove; l<=lmax; l++)
      {
        double rms = sqrt(newpow.tt(l));
        m_alm_after(l,0) = gsl_ran_gaussian(rng, rms);
     
		    for (m=1; m<=l; m++)
		    {
			    m_alm_after(l,m).Set(gsl_ran_gaussian(rng, rms*hsqrt2), gsl_ran_gaussian(rng, rms*hsqrt2));
		    }
	    }
      
      gsl_rng_free(rng);

      // Get the final map and power spectrum
      cout << "--converting coefficients back to map with nside=" << out_nside << endl;
      extract_powspec(m_alm_after, m_pow_after);

      m_map_after.SetNside(out_nside, RING);
      double offset = m_alm_after(0,0).real()/sqrt(4.0*3.141592653);
      m_alm_after(0,0) = 0;
  	  alm2map(m_alm_after, m_map_after);
      m_map_after.add(offset);

    }

  };

#endif
