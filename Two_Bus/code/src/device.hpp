#ifndef OPSP_DEVICE_H
#define OPSP_DEVICE_H
#include <iostream>

#include "terminal.hpp"
#include "globals.h"

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

class Device
{
protected:
  const int len;  // length of each terminal
  const int N;    // total number of terminals
  Terminal ** terminals; // each terminal can only be in at most one net
  // TODO: terminals need a way of knowing if they are in a net
  
public:
  // what happens when N = 0?
  Device(int N, int len);
  
  virtual ~Device();
  
  virtual void solve(const double &rho, const double &rho_old) = 0;
  
  virtual const double objective() const { return 0.0; }
  
  inline int num_terms() const;
  
  inline const Terminal & terminal(int i) const;
  inline Terminal *mutable_terminal(int i);
  
  void print();
  
private:
  // disallow copy constructor
  Device(const Device&);
  Device& operator=(const Device&);
};

int Device::num_terms() const
{
  return N;
}

const Terminal & Device::terminal(int i) const
{
#ifdef DEBUG
  assert(i < N && i >= 0);
#endif
  return *terminals[i];
}

Terminal *Device::mutable_terminal(int i)
{
#ifdef DEBUG
  assert(i < N && i >= 0);
#endif
  return terminals[i];
}


#include "devices/opsp_generator.hpp"
#include "devices/opsp_fixed_load.hpp"
#include "devices/opsp_battery.hpp"
#include "devices/opsp_deferrable_load.hpp"
#include "devices/opsp_curtailable_load.hpp"
#include "devices/opsp_lossless_line.hpp"
#include "devices/opsp_lossless_quadratic_line.hpp"
#include "devices/opsp_lossy_line.hpp"
#include "devices/opsp_lossy_dc_line.hpp"
#include "devices/opsp_gen_load.hpp"

#endif // OPSP_DEVICE_H
