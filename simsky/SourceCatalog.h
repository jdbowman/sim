#ifndef SOURCECATALOG_H
#define SOURCECATALOG_H

#include <vector>
#include <string>
#include <gsl/gsl_rng.h>

using namespace std;


/*
  This class stores a catalog of point sources with their
  theta, phi, flux, and spectral index.  It assumes theta
  and phi arein units of radians.  It has no knowledge of 
  units for the flux or spectral index.

  This class can dump the catalog to a simple text file
  with one column for each property and as many rows as
  sources in the catalog.  In the output text file, the
  positions are expressed as DEC [deg], RA [deg], flux,
  and spectral index.

  The only added funtionality is the ability to randomly
  reposition the sources within a specified radius of
  their existing positions in the catalog.
*/
class SourceCatalog
  {
  private:

    vector<double> m_theta;
    vector<double> m_phi;
    vector<double> m_flux;
    vector<double> m_index;

  public:
    
    SourceCatalog() {};
    ~SourceCatalog() {};

    void Add(double theta, double phi, double flux, double index) 
    { 
      m_theta.push_back(theta);
      m_phi.push_back(phi);
      m_flux.push_back(flux);
      m_index.push_back(index);
    }

    void Randomize(double rad, int seed)
    {
      // This function will randomize the theta and phi
      // coordinates for each source, replacing the existing
      // values with new ones.  The new position will be
      // drawn randomly with uniform weighting from the
      // solid angle centered on the current position, and
      // with radius given by "rad".  This randomization is
      // done correctly for the surface of sphere (it assumes
      // no small angle approximation).  Hence, "rad" can
      // range from 0 to PI radians.
      double dr=0, dang=0; 
      double x,y,z;
      double u,v,w;
      double x1,y1,z1;
      double theta, phi;

      // Setup the GSL random number generator
      gsl_rng_env_setup();
      gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(rng, seed);

      // Loop over all sources and change positions
      for (uint n=0; n<m_index.size(); n++)
      {
        dr = acos(1.0 - gsl_rng_uniform(rng) * (1.0 - cos(rad)));
        dang = 2.0 * 3.141592653 * gsl_rng_uniform(rng);
        theta = m_theta[n];
        phi = m_phi[n];

        // A vector to the current position of the source
        u = cos(phi)*sin(theta);
        v = sin(phi)*sin(theta);
        w = cos(theta);

        // A vector that is dr distance away in DEC from 
        // the current position of the source
        x = cos(phi)*sin(theta-dr);
        y = sin(phi)*sin(theta-dr);
        z = cos(theta-dr);

        // Now rotate the new vector around the old vector
        // by dang radians. Equations from:
        // http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
        x1 = u*(u*x+v*y+w*z) + (x*(v*v+w*w)+u*(-v*y-w*z))*cos(dang) + (-w*y+v*z)*sin(dang);
        y1 = v*(u*x+v*y+w*z) + (y*(u*u+w*w)+v*(-u*x-w*z))*cos(dang) + (w*x-u*z)*sin(dang);
        z1 = w*(u*x+v*y+w*z) + (z*(u*u+v*v)+w*(-u*x-v*y))*cos(dang) + (-v*x+u*y)*sin(dang);  

        // Calculate the new theta, phi and replace the old
        m_theta[n] = acos(z1);
        m_phi[n] = atan2(y1, x1);       

      }
    }

    void WriteToFile(const string &strFilename)
    {
      // Dump the catalog to a text file with the following columns:
      // DEC [deg -90->90], RA [deg 0->360], flux, spectral index

      FILE* fid=NULL;
printf("HELLO!!\n");
printf("%s\n", strFilename.c_str());
      fid = fopen(strFilename.c_str(), "w");
      if (fid != NULL)
      {
	printf("here!\n");
        for (uint n=0; n<m_index.size(); n++)
        {
          fprintf(fid, "%g, %g, %g, %g\n", 90.0 - 180.0 / 3.141592653 * m_theta[n], fmod(180.0 / 3.141592653 * m_phi[n] + 360.0, 360.0), m_flux[n], m_index[n]);
        }
        fclose(fid);
      }
    }


  };

#endif
