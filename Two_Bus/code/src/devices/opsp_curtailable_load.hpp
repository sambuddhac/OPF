#ifndef OPSP_CURTAILABLE_LOAD_DEVICE_H
#define OPSP_CURTAILABLE_LOAD_DEVICE_H

#include "../device.hpp"
//#include "../curtailable_load/solver.h"

/*
 * CURTAILABLE LOAD DEVICE
 */
// XXX: use analytical solution (because one exists)
class CurtailableLoad : public Device
{
private:
  const double pdes;
  const double gamma;
public:
  CurtailableLoad(double pdes, double gamma, int len) : 
    Device(1, len), pdes(pdes), gamma(gamma)
  { }

  virtual void solve(const double &rho, const double &rho_old);
  
  virtual const double objective() const;
};

#endif