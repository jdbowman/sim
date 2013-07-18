#ifndef SOURCECONFUSIONMAP_H
#define SOURCECONFUSIONMAP_H

#include <vector>
#include <string>
#include <sstream>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include "healpix_map.h"

using namespace std;

class SourceConfusionMap
  {

  private:

    vector<Healpix_Map<double> > m_maps;
    vector<double> m_freqs;
    double m_freqmin; // MHz
    double m_freqmax; // MHz
    double m_center; // MHz
    double m_pixelarea; // str
    uint m_nspec;
    gsl_multifit_linear_workspace* m_work;


  public:


    SourceConfusionMap(double freqmin, double freqmax, double center, uint nside, uint nspec) :
      m_freqmin(freqmin),
      m_freqmax(freqmax),
      m_center(center),
      m_nspec(nspec)
    {
      // Calculate our map frequencies
      double df = (freqmax - freqmin) / (double) nspec;

      // Populate the maps and freqs vectors
      for (uint j=0; j<=nspec; j++)
      {
        m_maps.push_back(Healpix_Map<double>(nside, RING, SET_NSIDE));
        m_freqs.push_back((double) j * df + freqmin);
      }

      m_pixelarea = 4.0 * 3.141592653 / m_maps[0].Npix();

      // Setup the GSL fitting
      m_work = gsl_multifit_linear_alloc(nspec+1, nspec+1);


    }


    ~SourceConfusionMap()
    {
		// Clean up the GSL fitting
		gsl_multifit_linear_free(m_work);

	}


    void Add(int p, double flux, double index)
    {
      // Convert flux to brightness temperature
      double temp = (flux/m_pixelarea) * 1e-26 * 2.99792e8*2.99792e8 / (2.0 * 1.38e-23 * m_center*m_center*1e12);

      // Calculate temp at each map frequency and add
      for (uint j=0; j<m_maps.size(); j++)
      {
        (m_maps[j])[p] += temp * pow(m_center/m_freqs[j], index-2.0);
      }
    }


    void ConvertToLogPoly()
    {
		// Fit a polynomial of order nspec in log-log space
		// for each pixel and replace temperature maps with
		// polynomial coefficient maps

		uint j=0;
		uint k=0;
		double chisq=0;

		gsl_matrix *X = gsl_matrix_alloc(m_nspec+1, m_nspec+1);
		gsl_matrix *cov = gsl_matrix_alloc(m_nspec+1, m_nspec+1);
		gsl_vector *y = gsl_vector_alloc(m_nspec+1);
		gsl_vector *c = gsl_vector_alloc(m_nspec+1);

		// We only have to build the X matrix once
		for (j=0; j<=m_nspec; j++)
		{
			for (k=0; k<=m_nspec; k++)
			{
        gsl_matrix_set(X, j, k, pow(log10(m_freqs[j]/m_center), k));
			}
		}

		// Loop over pixels
		for (j=0; j < (uint) m_maps[0].Npix(); j++)
		{
			// Populate the y vector for this pixel
			for (k=0; k<m_freqs.size(); k++)
			{
				gsl_vector_set(y, k, log10(m_maps[k][j]));
			}

			// Do the fit now
			gsl_multifit_linear(X, y, c, cov, &chisq, m_work);

			// Replace the map values with these new coefficients
			for (k=0; k<m_freqs.size(); k++)
			{
				m_maps[k][j] = gsl_vector_get(c,k);
			}
		}

    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);

	}



	void WriteToFiles(const string& strBase)
	{
		fitshandle out;
		uint j=0;

		for (j=0; j<m_maps.size(); j++)
		{
			ostringstream os;
      os << strBase << j << ".fits";

			out.create(os.str());
  			write_Healpix_map_to_fits(out, m_maps[j], FITSUTIL<double>::DTYPE);
  			out.close();
		}
	}

  };

#endif
