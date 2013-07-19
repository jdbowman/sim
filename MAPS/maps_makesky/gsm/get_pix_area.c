#include <stdio.h>
#include <math.h>
 
double get_pix_area (double dxin_deg, double dyin_deg, double dec) {

  // sin.theta d.theta d. phi

  double dxin_rad,dyin_rad,decin_rad;
  double area;

  dxin_rad = dxin_deg * M_PI/180.0;
  dyin_rad = dyin_deg * M_PI/180.0;
  decin_rad = (dec+90) * M_PI/180.0;

  area = fabs(sin(decin_rad)) * fabs(dxin_rad) * fabs(dyin_rad);
//  printf("dx = %g, dy=%g, dec=%g AREA: %g\n",dxin_rad,dyin_rad,decin_rad,area);
  return area;
  
}

