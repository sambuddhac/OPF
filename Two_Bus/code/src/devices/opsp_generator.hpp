#ifndef OPSP_GENERATOR_DEVICE_H
#define OPSP_GENERATOR_DEVICE_H

/*
 * GENERATOR DEVICE
 */
#include "../device.hpp"
#include "../generator/solver.h"

class Generator : public Device
{
private:
  CVX_Generator *gen;
public:
  Generator(double alpha, double beta, double pmax, double pmin, double slew, int len);
  
  virtual void solve(const double &rho, const double &rho_old);
  
  virtual const double objective() const;
  
  ~Generator()
  {
    delete gen;
  }
};

#endif // OPSP_GENERATOR_DEVICE_H