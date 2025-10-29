#ifndef FUNCIONAL_HELICOIDAL
#define FUNCIONAL_HELICOIDAL

#include "comum.hpp"

real_1d_array helicoidal(double h, double alpha, double b, double t, const real_1d_array &y0);
double funcional_helicoidal(const real_1d_array &y,double H, double alpha, double B, double T);


#endif