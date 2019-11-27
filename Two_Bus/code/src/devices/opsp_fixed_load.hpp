#ifndef OPSP_FIXEDLOAD_DEVICE_H
#define OPSP_FIXEDLOAD_DEVICE_H

#include "../device.hpp"

/*
 * FIXED LOAD DEVICE
 */
class FixedLoad : public Device
{
private:
  double *__restrict__ loss_profile;
public:
  FixedLoad(double *__restrict__ l, int len);
  virtual void solve(const double &rho, const double &rho_old);
  
  ~FixedLoad()
  {
    delete[] loss_profile;
  }
};

#endif