// EarthTilt.c: This is the function earthtilt.c from NOVAS-C Version 2.0.1 (10 Dec 99)

/*
   Naval Observatory Vector Astrometry Subroutines
   C Version

   U. S. Naval Observatory
   Astronomical Applications Dept.
   3450 Massachusetts Ave., NW
   Washington, DC  20392-5420
*/

void EarthTilt (double tjd,

                double *mobl, double *tobl, double *eq, double *dpsi,
                double *deps)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes quantities related to the orientation of the Earth's
      rotation axis at Julian date 'tjd'.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97,
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.
      Transactions of the IAU (1994). Resolution C7; Vol. XXIIB, p. 59.
      McCarthy, D. D. (ed.) (1996). IERS Technical Note 21. IERS
         Central Bureau, Observatoire de Paris), pp. 21-22.

   INPUT
   ARGUMENTS:
      tjd (double)
         TDB Julian date of the desired time

   OUTPUT
   ARGUMENTS:
      *mobl (double)
         Mean obliquity of the ecliptic in degrees at 'tjd'.
      *tobl (double)
         True obliquity of the ecliptic in degrees at 'tjd'.
      *eq (double)
         Equation of the equinoxes in seconds of time at 'tjd'.
      *dpsi (double)
         Nutation in longitude in arcseconds at 'tjd'.
      *deps (double)
         Nutation in obliquity in arcseconds at 'tjd'.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      PSI_COR, EPS_COR, DEG2RAD

   FUNCTIONS
   CALLED:
      nutation_angles  novas.c
      fund_args        novas.c
      fabs             math.h
      pow              math.h
      cos              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V1.1/06-97/JAB (USNO/AA) Incorporate IAU (1994) and IERS (1996)
                               adjustment to the "equation of the
                               equinoxes".
      V1.2/10-97/JAB (USNO/AA) Implement function that computes
                               arguments of the nutation series.
      V1.3/07-98/JAB (USNO/AA) Use global variables 'PSI_COR' and
                               'EPS_COR' to apply celestial pole offsets
                               for high-precision applications.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'etilt'.
      2. Values of the celestial pole offsets 'PSI_COR' and 'EPS_COR'
      are set using function 'cel_pole', if desired.  See the prolog
      of 'cel_pole' for details.

------------------------------------------------------------------------
*/
{
   static double tjd_last = 0.0;
   static double t, dp, de;
   double d_psi, d_eps, mean_obliq, true_obliq, eq_eq, args[5];

/*
   Compute time in Julian centuries from epoch J2000.0.
*/

  t = (tjd - T0) / 36525.0;

/*
   Compute the nutation angles (arcseconds) from the standard nutation
   model if the input Julian date is significantly different from the
   last Julian date.
*/

  if (fabs (tjd - tjd_last) > 1.0e-6)
      nutation_angles (t, &dp,&de);

/*
   Apply observed celestial pole offsets.
*/

   d_psi = dp + PSI_COR;
   d_eps = de + EPS_COR;

/*
   Compute mean obliquity of the ecliptic in arcseconds.
*/

   mean_obliq = 84381.4480 - 46.8150 * t - 0.00059 * pow (t, 2.0)
      + 0.001813 * pow (t, 3.0);

/*
   Compute true obliquity of the ecliptic in arcseconds.
*/

   true_obliq = mean_obliq + d_eps;

/*
   Convert obliquity values to degrees.
*/

   mean_obliq /= 3600.0;
   true_obliq /= 3600.0;

/*
   Compute equation of the equinoxes in seconds of time.

   'args[4]' is "omega", the longitude of the ascending node of the
   Moon's mean orbit on the ecliptic in radians.  This is also an
   argument of the nutation series.
*/

   fund_args (t, args);

   eq_eq = d_psi * cos (mean_obliq * DEG2RAD) +
      (0.00264  * sin (args[4]) + 0.000063 * sin (2.0 * args[4]));

   eq_eq /= 15.0;

/*
   Reset the value of the last Julian date and set the output values.
*/

   tjd_last = tjd;

   *dpsi = d_psi;
   *deps = d_eps;
   *eq = eq_eq;
   *mobl = mean_obliq;
   *tobl = true_obliq;

   return;
}
