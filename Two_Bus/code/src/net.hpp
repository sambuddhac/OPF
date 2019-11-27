#ifndef OPSP_NET_H
#define OPSP_NET_H

#include <cstdlib>

#include "terminal.hpp"
#include "globals.h"

// XXXX: probaby could use templates in this code... 
// XXX: probably rewrite? so reduce collects dc_pbar and ac_pbar based on terminal "type"
class Net
{
private:  
  int N;
  Terminal **terminals;

  bool address_exists(const Terminal *t, int low, int high);
  
public:
  Net() : N(0), terminals(NULL) {}
  
  // XXX: shouldn't be const...
  // adding the terminal is a const operation, but i do "non-const" stuff with it
  void add_terminal(const Terminal *t);

  void reduce(double &primal_residual, double &dual_residual);
    
  void print();
  
  ~Net()
  {
    if(N > 0 && terminals) free(terminals);
  }
};

#endif // OPSP_NET_H
