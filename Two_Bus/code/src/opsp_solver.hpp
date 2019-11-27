#ifndef OPSP_SOLVER_H
#define OPSP_SOLVER_H

#include "power_network.pb.h"

#include "device.hpp"
#include "terminal.hpp"
#include "net.hpp"

#include <fstream>


class OPSPSolver
{
public:
  enum Problem { LOSSLESS_STAR, LOSSLESS_QUADRATIC, FULL };
  
  // TODO: needs to take rho as a parameter  
  OPSPSolver(OPF::Network& network, int T, int max_iters, double rho,
             Problem p = FULL, bool verbose = true, bool warmstart = false, double tolerance = ABSTOL) : 
    _num_devices(0), _num_nets(0), _num_terminals(0),
    T(T), verbose(verbose), warmstart(warmstart), max_iters(max_iters), rho(rho),
    _devices(NULL), _nets(NULL), network(network), tolerance(tolerance)
  { 
    switch(p) {
      case LOSSLESS_STAR:
        setup_lossless_star_network();
        break;
      case LOSSLESS_QUADRATIC:
        setup_lossless_quadratic_network();
        break;
      case FULL:
        setup_network();
        break;
      default:
        setup_network();
    }
  }
  
  int solve(bool use_stopping);
  
  double rho_val() { return rho; }
  int num_devices() { return _num_devices; }
  const Device **devices() { return (const Device **) _devices; }
  
  ~OPSPSolver()
  {
    if(_devices) {
      for(int i = 0; i < _num_devices; i++)
        if(_devices[i]) delete _devices[i];
      delete[] _devices;
    }
    
    if(_nets) {
      for(int i = 0; i < _num_nets; i++)
        if(_nets[i]) delete _nets[i];
      delete[] _nets;
    }
  }
  
private:  
  int _num_devices;
  int _num_nets;
  int _num_terminals;
  
  const int T;
  const bool verbose;
  const bool warmstart;
  const int max_iters;
  const double tolerance;
  double rho;
  
  char buffer[64];
  
  
  Device **_devices; 
  Net **_nets;
  
  OPF::Network& network;
  
  void setup_lossless_star_network();
  void setup_lossless_quadratic_network();
  void setup_network();
  
  void read_device(const char *file, OPF::Device& device)
  {
    // Read the network.
    std::fstream input(file, std::ios::in | std::ios::binary);
    if (!device.ParseFromIstream(&input)) {
      std::cerr << "Failed to parse device" << std::endl;
      exit(-1);
    }
  }
  
  void warmstart_device(Device *dev, int k);
  
  double evaluate_objective();
  
  Device* add_device(const OPF::Bus &bus);
  
  int num_distribution_bus();
};

#endif
