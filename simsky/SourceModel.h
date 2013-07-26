#ifndef SOURCEMODEL_H
#define SOURCEMODEL_H

#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include "angular_correlation_tools.h"

using namespace std;


/*
	This class uses cummulative probability functions for both
	the flux and spectral index distribution of a model source
	population in order to draw random samples upon request.

  This class assumes that the data file for flux is in xy 
  pairs of log10(flux [Jy]), CDF

  However it assumes that the data file for spectral index
  is in xy pairs of alpha, CDF

  The DrawSource returns the flux in Jy and the spectral 
  index as alpha (flux density units, not temperature units)

*/
class SourceModel
  {
  private:

    gsl_rng* m_rng;
    gsl_spline* m_spline_flux;
    gsl_spline* m_spline_index;
    gsl_interp_accel* m_acc_flux;
    gsl_interp_accel* m_acc_index;
    int m_seed;
    vector<double> m_flux_x;
    vector<double> m_flux_cdf;
    vector<double> m_index_x;
    vector<double> m_index_cdf;

  public:


    SourceModel(int seed) : m_seed(seed)
    {
      // Setup the GSL random number generator
      gsl_rng_env_setup();
      m_rng = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(m_rng, m_seed);

      // Clear the GSL interps
	    m_acc_flux = NULL;
      m_acc_index = NULL;
      m_spline_flux = NULL;
      m_spline_index = NULL;
    }


    ~SourceModel()
    {
      // Clean up the GSL random number generator
      gsl_rng_free(m_rng);

      // Free the GSL spline interpolators
      gsl_spline_free(m_spline_flux);
      gsl_spline_free(m_spline_index);
      gsl_interp_accel_free(m_acc_flux);
      gsl_interp_accel_free(m_acc_index);
    }


	void LoadCDFs(const string& strFluxFile, const string& strIndexFile)
	{
    // Clear any existing CDF state  
    m_flux_x.clear();
    m_flux_cdf.clear();
    m_index_x.clear();
    m_index_cdf.clear();

    if (m_spline_flux)
    {
      gsl_spline_free(m_spline_flux);
      m_spline_flux = NULL;
    }

    if (m_spline_index)
    {
      gsl_spline_free(m_spline_index);
      m_spline_flux = NULL;
    }

    if (m_acc_flux)
    {
      gsl_interp_accel_free(m_acc_flux);
      m_acc_flux = NULL;
    }

    if (m_acc_index)
    {
      gsl_interp_accel_free(m_acc_index);
      m_acc_index = NULL;
    }

    
		// Read flux CDF from file
    read_xydata_from_file(strFluxFile, m_flux_x, m_flux_cdf);

		// Read index CDF from file
		read_xydata_from_file(strIndexFile, m_index_x, m_index_cdf);

    // Setup the GSL routines
    m_acc_flux = gsl_interp_accel_alloc();
		m_spline_flux = gsl_spline_alloc(gsl_interp_linear, m_flux_x.size());
		gsl_spline_init(m_spline_flux, &m_flux_cdf[0], &m_flux_x[0], m_flux_x.size());

    m_acc_index = gsl_interp_accel_alloc();
		m_spline_index = gsl_spline_alloc(gsl_interp_linear, m_index_x.size());
		gsl_spline_init(m_spline_index, &m_index_cdf[0], &m_index_x[0], m_index_x.size());
	}


    void DrawSource(double &flux, double &index)
    {
      // This function returns a random flux [Jy] and spectral
      // index [spectral-density units] drawn from the CDFs
      // provided by the data files loaded in LoadCDFs
      flux = pow(10, gsl_spline_eval(m_spline_flux, gsl_rng_uniform(m_rng), m_acc_flux));
      index = gsl_spline_eval(m_spline_index, gsl_rng_uniform(m_rng), m_acc_index);
    }

  };

#endif
