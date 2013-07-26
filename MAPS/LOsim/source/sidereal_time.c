// sidereal_time.c: This is the function sidereal_time from NOVAS-C Version 2.0.1 (10 Dec 99)

/*
   Naval Observatory Vector Astrometry Subroutines
   C Version

   U. S. Naval Observatory
   Astronomical Applications Dept.
   3450 Massachusetts Ave., NW
   Washington, DC  20392-5420
*/

#include "sidereal_time.h"
#define _CONSTS_ 1 // prevents novascon.c from including novascon.h
#include "novascon.c" // defines T0
#include <math.h>

void sidereal_time (double jd_high, double jd_low, double ee, double *gst)

/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the Greenwich apparent sidereal time, at Julian date
      'jd_high' + 'jd_low'.

   REFERENCES:
      Aoki, et al. (1982) Astronomy and Astrophysics 105, 359-361.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      jd_high (double)
         Julian date, integral part.
      jd_low (double)
         Julian date, fractional part.
      ee (double)
         Equation of the equinoxes (seconds of time). [Note: this
         quantity is computed by function 'earthtilt'.]

   OUTPUT
   ARGUMENTS:
      *gst (double)
         Greenwich apparent sidereal time, in hours.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0

   FUNCTIONS
   CALLED:
      fmod    math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C programing standards.
      V1.2/03-98/JAB (USNO/AA) Expand documentation.
      V1.3/08-98/JAB (USNO/AA) Match flow of the Fortran counterpart.

   NOTES:
      1. The Julian date may be split at any point, but for highest
      precision, set 'jd_high' to be the integral part of the Julian
      date, and set 'jd_low' to be the fractional part.
      2. For Greenwich mean sidereal time, set input variable 'ee'
      equal to zero.
      3. This function is based on Fortran NOVAS routine 'sidtim'.

------------------------------------------------------------------------
*/

{
   double t_hi, t_lo, t, t2, t3, st;

   t_hi = (jd_high -  T0) / 36525.0;
   t_lo = jd_low / 36525.0;
   t = t_hi + t_lo;
   t2 = t * t;
   t3 = t2 * t;

   st =  ee - 6.2e-6 * t3 + 0.093104 * t2 + 67310.54841
      + 8640184.812866 * t_lo
      + 3155760000.0 * t_lo
      + 8640184.812866 * t_hi
      + 3155760000.0 * t_hi;

   *gst = fmod ((st / 3600.0), 24.0);

   if (*gst < 0.0)
      *gst += 24.0;

   return;
}
