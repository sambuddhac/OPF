#ifndef OPSP_LOSSY_DC_LINE_DEVICE_H
#define OPSP_LOSSY_DC_LINE_DEVICE_H

#include "../device.hpp"

class LossyDCLine : public Device
{
private:
  const double g;
  //const double b;
  //const double C;
  
  const double M;
  const double beta_square;
  const double L;
  
  const double BOUND;
   
public:
  LossyDCLine(double g, double b, double C, int len) : 
    Device(2, len), g(g), 
    M(C/(2.0*g)), beta_square((b*b)/(g*g)), L(1.0 - sqrt(1.0 - M*M/beta_square)),
    BOUND((1.0 - 1.0/beta_square)*M*L - M)
  { }
  
  virtual void solve(const double &rho, const double &rho_old);

};

#endif


