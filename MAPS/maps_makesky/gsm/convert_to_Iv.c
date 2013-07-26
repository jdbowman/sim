#include <stdio.h>
#include <math.h>
#include "utils.h"

// Converts a pixel Brightness Temp to a pixel spectral brightness (I_v) 
// in Jy/sr

double convert_to_Iv(double temperature, double frequency, double area) {
  
  double lambda,scaling,flux,temp;

  lambda = SOL / (frequency*1E6);
  
  temp = (2 * BOLT)/(pow(lambda,2)); 

  scaling = temp*temperature*JANSKY;

  return scaling;

  
}

