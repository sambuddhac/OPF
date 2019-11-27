#ifndef OPSP_LOSSLESS_LINE_DEVICE_H
#define OPSP_LOSSLESS_LINE_DEVICE_H

#include "../device.hpp"

class LosslessLine : public Device
{
public:
  LosslessLine(int len) : Device(2, len) { }
  
  virtual void solve(const double &rho, const double &rho_old);
};

#endif