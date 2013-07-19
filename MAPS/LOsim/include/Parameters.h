// Parameters.h: Constants and simulator parameters

#define BLOCKSIZE 4096				// FAKE-Linux filesystem block size in bytes

#define ELEMENTSPERBLOCK ((BLOCKSIZE)/sizeof(struct Complex))
					// number of matrix elements in a single block of the 1-D
					// row-transformed matrix and the output visibility matrix
					
#define cLight 299792458 // speed of light in m/sec

#define PI 3.14159265358979323846

#define RADIUSEARTH 6378388.0			// in m
#define ONEDAYHRS  24.0					// hours in a day
#define ONEDAYSECS (ONEDAYHRS*60L*60L)	// seconds in a day
#define HOURStoRADIANS (2.0*PI/ONEDAYHRS)		// 1 hour, in radians
#define SECONDStoRADIANS (2.0*PI/ONEDAYSECS) // 1 second, in radians

#define MAX_IMAGES 3					// number of externally-supplied sky images
#define MAX_STATIONS 500
#define MAX_ANTENNAEPERSTATION 200
#define MAX_ANTENNAE (MAX_STATIONS*MAX_ANTENNAEPERSTATION)

// Lengths of the global character variables
#define L_ONELINE (80+1)
#define L_BRIEFNAME (20+1)
#define L_NORMAL (50+1)
#define L_EXTENDEDLINE (132+1)

#define MAX_FIELDNAME 50

#define MAX_SPECTRAL_POINTS 1024

// Debugging parameters
#define testPosition 0 // if this is non-zero, dummy files are written and position of each element is checked
#define bBase 1000000.0

