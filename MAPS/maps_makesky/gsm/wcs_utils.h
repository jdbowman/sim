#ifndef __WCS_UTILS_H
#define __WCS_UTILS_H
/* Some simple functions utilising the WCS library */

#include <wcslib/wcslib.h>

/*Warning: pixrane must be atleast a float[4] and outpix must be at least a double[4][2]*/

int get_pixels (struct wcsprm *wcs, float *pixrange, double **corners, double **outpix,int range);
int get_corners (struct wcsprm *wcs, float pix[2], double **corners, int image);

#endif
