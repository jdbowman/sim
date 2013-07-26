/*===========================================================+
 * sky_source.h                                              |
 * Defines data structs for a celestial source parameters    |
 * used in brigen.                                           |
 * Derived from DataStructures by                            |
 * Peter Sherwood, 8/22/2001, sherwood@computer.org,         |
 *	                      (617) 244-0836                 |
 * 8/22/01	added Stokes parameters and spectral density |
 *              to Source                                    |
 * Leonid Benkevitch, 02 Feb 2011                            |
 *===========================================================+
 */

#define MAX_SOURCES 500
#define MAX_SRCID_LEN 10

struct Complex {
  double re;
  double im;
};

// All sources are Gaussian and elliptical
struct Source {
  double fluxDensity;	// total flux
  double Q, U, V;	// Stokes parameters
  double specIndex;     // spectral index 
  double xOffset;       /* arcsec */
  double yOffset;       /* arcsec */
  double majorAxis;
  double minorAxis;
  double positionAngle;
  char id[MAX_SRCID_LEN+1];
};
