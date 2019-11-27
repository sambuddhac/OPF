#ifndef OPSP_BATTERY_DEVICE_H
#define OPSP_BATTERY_DEVICE_H

/*
 * BATTERY DEVICE
 */
 
#include "../device.hpp"
#include "../battery/solver.h"

class Battery : public Device
{
private:
  CVX_Battery *bat;
  
public:
  Battery(double q0, double qf, double P, double Q, int len);
  
  virtual void solve(const double &rho, const double &rho_old);
  
  ~Battery()
  {
    delete bat;
  }
};


#endif
