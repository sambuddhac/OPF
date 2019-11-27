#ifndef OPSP_DEFERRABLE_LOAD_DEVICE_H
#define OPSP_DEFERRABLE_LOAD_DEVICE_H

/*
 * DEFERRABLE LOAD DEVICE
 */
#include "../device.hpp"
#include "../deferrable_load/solver.h"

class DeferrableLoad : public Device
{
private:
  CVX_DeferrableLoad *def;
public:
  DeferrableLoad(double E, double pmin, double pmax, int Tstart, int Tend, int len);
  
  virtual void solve(const double &rho, const double &rho_old);

  ~DeferrableLoad()
  {
    delete def;
  }
};

#endif