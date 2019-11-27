#ifndef OPSP_LOSSY_LINE_DEVICE_H
#define OPSP_LOSSY_LINE_DEVICE_H

#include "../device.hpp"
#include <iostream>

#define RATIO 1

using namespace std;
class LossyLine : public Device
{
private:
  const double g;
  const double b;
  //const double C;
  
  const double M;
  const double L;
  const double p; // y-intercept
     
public:
  LossyLine(double g, double b, double C, int len) : 
    Device(2, len), g(g), b(b),
    M(C/b), L(1.0 - sqrt(1.0 - M*M)), // choose the negative part of the sqrt
    p(L - 0.5*(2+RATIO)*(L-1)/*1.5 - 0.5*L*/)
    { }
  
  virtual void solve(const double &rho, const double &rho_old);

};

#endif
