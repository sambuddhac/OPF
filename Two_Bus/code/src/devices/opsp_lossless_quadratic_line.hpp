#ifndef OPSP_LOSSLESS_QUADRATIC_LINE_DEVICE_H
#define OPSP_LOSSLESS_QUADRATIC_LINE_DEVICE_H

#include "../device.hpp"

class LosslessQuadraticLine : public Device
{
private:
  double delta; // quadratic penalty (typically small)
  
public:
  LosslessQuadraticLine(double delta, int len) : Device(2, len), delta(delta) { }
  
  // minimize rho/2(p1 - v1)^2 + rho/2(p2 - v2)^2 + delta/2(p1^2 + p2^2)
  // s.t. p1 + p2 = 0
  
  virtual void solve(const double &rho, const double &rho_old);
};

#endif