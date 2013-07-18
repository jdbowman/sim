#ifndef MAPINPAINT_H
#define MAPINPAINT_H

#include <vector>
#include <string>
#include "alm.h"
#include "alm_healpix_tools.h"
#include "alm_map_tools.h"
#include "alm_powspec_tools.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "powspec.h"


using namespace std;


/*
	This class

  The algorithm used here is derived from:

    CMB data analysis and sparsity
    Abrial P., Moudden Y., Starck J.-L., Fadili J., 
    Delabrouille J., Nguyen M.K.

*/
class MapInpaint
  {
  private:

    Healpix_Map<double> m_map;
    Healpix_Map<double> m_mask;
    Healpix_Map<double> m_sol;
    uint m_lmax;
    uint m_niter;

  public:


    MapInpaint() : m_lmax(0), m_niter(0)
    {
    }


    ~MapInpaint()
    {
    }


    void Load(const string& strMapFile, const string& strMaskFile)
    {
      cout << "MapInpaint: Loading..." << endl;
      
      double dmin=0, dmax=0;

      // Read map
      read_Healpix_map_from_fits(strMapFile.c_str(), m_map);
      if (m_map.Scheme() != RING)
      {
        cout << "--swapping map scheme to RING" << endl;
        m_map.swap_scheme();
      }
      m_map.Map().minmax(dmin, dmax);
      cout << "--map [min, max]=" << dmin << ", " << dmax << endl;

      // Read mask
      read_Healpix_map_from_fits(strMaskFile.c_str(), m_mask);
      if (m_mask.Scheme() != RING)
      {
        cout << "--swapping mask scheme to RING" << endl;
        m_mask.swap_scheme();
      }
      m_mask.Map().minmax(dmin, dmax);
      cout << "--mask [min, max]=" << dmin << ", " << dmax << endl;

      // Initialize solution
      m_sol.SetNside(m_map.Nside(), RING);
    }


    Healpix_Map<double> & Map() { return m_map; }

    Healpix_Map<double> & Mask() { return m_mask; }

    Healpix_Map<double> & Painted() { return m_sol; }


    void Inpaint(const uint niter, const uint lmax)
    {
      cout << "MapInpaint: Inpainting..." << endl;

      double lambda_max = 0;  // This gets set to the max alm on the first iteration
      double lambda_min = 0;  // Setting this to zero reproduces the data exactly (within the limits of map2alm)
      double lambda = 0;
      Alm<xcomplex<double> > alm(lmax, lmax);
      
      arr<double> ring_weights(2*m_map.Nside());
      ring_weights.fill(1.0);

      m_sol.fill(0.0);
      m_lmax = lmax;
      m_niter = niter;

      for (uint n=0; n<niter; n++)
      {
        cout << "--iteration=" << n << " of " << niter << flush; 

        // Update our best guess in image domain
        for (uint j=0; j<(uint)m_map.Npix(); j++)
        {
          m_sol[j] = m_mask[j]*m_map[j] - m_mask[j]*m_sol[j] + m_sol[j];
        }

        // Get alms for new best guess
        double avg = m_sol.average();
        m_sol.add(-avg);
	      map2alm_iter(m_sol, alm, (n==(niter-1) ? 3 : 1), ring_weights);
        alm(0,0) += avg*sqrt(4.0*3.141592653);
        
        // Update the threshold condition
        if (n==0)
        {
          // Start the iteration by setting the threshold to the 
          // largest alm coefficient
          double maxalm = 0;
          double rmsalm = 0;
          for (uint l=1; l<=m_lmax; l++)
          { 
            for (uint m=0; m<=l; m++)
            {
              maxalm = ( alm(l,m).norm() > maxalm ? alm(l,m).norm() : maxalm ); 
              rmsalm += sqrt(alm(l,m).norm());
            }
          }

          rmsalm = rmsalm / alm.Num_Alms(m_lmax, m_lmax);
          //lambda_max = sqrt(maxalm);
          lambda_max = rmsalm;
          lambda = lambda_max;
        }
        else
        {
          lambda = lambda_max  - n * (lambda_max - lambda_min) / (niter - 1.0);
        }

        cout << " lambda=" << lambda << " lambda_max=" << lambda_max;

        // Apply threshold condition using "hard thresholding" 
        // to discard insignificant alm coefficients
        uint count = 0;
        double lstep = (double) m_lmax / (double) niter;
        uint lstart = 1 + floor(lstep * (n+1.0));
        for (uint l=lstart; l<=m_lmax; l++)
        { 
          for (uint m=0; m<=l; m++)
          {
            //if (sqrt(alm(l,m).norm()) < lambda) 
            {
              count++;
              alm(l,m) = 0;
            }
            
          }
        }

        cout << " lstart=" << lstart << " count=" << count << endl;

        // Get the new estimate by reconstructing from the selected coefficients
        double offset = alm(0,0).real()/sqrt(4.0*3.141592653);
        alm(0,0) = 0;
        alm2map(alm, m_sol);
        m_sol.add(offset);

      } 

      cout << "--done" << endl;
    }

  };

#endif
