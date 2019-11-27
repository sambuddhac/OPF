#ifndef OPSP_GENERATING_LOAD_DEVICE_H
#define OPSP_GENERATING_LOAD_DEVICE_H

#include "../device.hpp"
//#include "../curtailable_load/solver.h"

/*
 * GENERATING LOAD DEVICE
 */
// XXX: use analytical solution (because one exists)
class GeneratingLoad : public Device
{
private:
  double *__restrict__ pdes;
  const double gamma;
public:
  GeneratingLoad(double *__restrict__ l, double gamma, int len);

  virtual void solve(const double &rho, const double &rho_old);
  
  virtual const double objective() const;

  ~GeneratingLoad()
  {
    delete[] pdes;
  }
};

#endif