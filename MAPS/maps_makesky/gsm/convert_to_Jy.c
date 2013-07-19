#include <stdio.h>
#include <math.h>
#include "utils.h"

// Converts a pixel Brightness Temp to a pixel Jy/Sr

double convert_to_Jy(double temperature, double frequency, double area) {
  
  double lambda,scaling,flux,temp;

  
  lambda = SOL / (frequency*1E6);
  
  temp = (2 * BOLT * area)/(pow(lambda,2) * 4 * M_PI); 

  scaling = temp*temperature*JANSKY;

  return scaling;

  
}

