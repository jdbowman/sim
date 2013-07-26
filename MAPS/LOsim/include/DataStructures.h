// DataStructures.h: Defines special data types (structs) used in LOsim

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/22/01	add Stokes parameters and spectral density to Source

#ifndef dataStructuresHaveBeenDefined
#include <sys/stat.h>

struct Complex
{
	double re;
	double im;
};

// File status
struct statL
{
	struct stat s;
	long long st_sizeL;
};

// All sources are Gaussian and elliptical
struct Source
{
	double fluxDensity;	// total flux
	double Q,U,V;		// Stokes parameters
	double specIx;		// spectral index
	double xOffset;
	double yOffset;
	double majorAxis;
	double minorAxis;
	double positionAngle;
	char id[10+1];
};

// User-supplied sky images
struct ImageFileParams
{
	char filename[100+1];	// name of the file containing the image
	int isFITS;				// 1=file is in FITS format
	
	// Transformations to be applied to the image, in order of application:
	int reflectX,reflectY;	// reflect image in x- or y-axis, respectively
	double ctrX,ctrY;		// where in the sky plane the image center should appear
	double scaleX,scaleY;	// factor to change dimensions of image
	double rotate;			// radians by which image is to be rotated
	
	// Transformations to be applied to individual pixel values:
	double scale;			// a multiplicative factor applied to all pixels
	double stretch;			// value in each pixel is raised to this power (power law)
	double blank;			// all pixels below this amount are set to zero
	double limit;			// pixels above this amount are set to this value
};

// FITS parameters from external images
struct Axis
{
	char name[68+1];	// CTYPEn
	int nPixels;		// NAXISn
	double refPix;		// CRPIXn	0-based pixel index for the following value
	double refVal;		// CRVALn	coordinate axis value at the above pixel
	double pixSpacing;	// CDELTn	distance between pixels, in axis units (assumed degrees)
	double rotate;		// CROTAn	rotation angle in degrees (not in use)
};

// A station consists of 1-MAX_ANTENNAEPERSTATION or more dipole antennae
struct Station
{
	double x,y,z;				// Cartesian coords relative to center of earth, in m
	double longi,lat;			// -pi<longitude<=pi (+=E)  -pi/2<latitude<=pi/2 (+=N)
	int nAnt;					// # antennae for this station
	unsigned short iAnt;		// the index of the 1st of the nAnt antennae in the antenna array
	char id[20+1];
};

struct Antenna
{
	struct Complex wt;			// weight
	double x,y,z;				// Cartesian coords relative to center of earth, in m
	char id[10+1];
};

#define dataStructuresHaveBeenDefined 1
#endif
