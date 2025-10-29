#ifndef FUNCIONAL_HOMOGENEO
#define FUNCIONAL_HOMOGENEO

#include "comum.hpp"

real_1d_array homogeneo(double h, double alpha, double b, double t,const real_1d_array &y0);
double funcional_homogeneo(const real_1d_array &y,double H, double alpha, double B, double T);


#endif