#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#include "utils.hpp"
#include "power_network.pb.h"

// #include "terminal.hpp"
// #include "device.hpp"
// #include "net.hpp"
#include "opsp_solver.hpp"

using namespace std;

typedef struct  {
  const char *file;
  double rho;
  int iters;
  int threads;
  bool verbose;
  bool warmstart;
  bool save;
} opf_options;

// function prototypes
opf_options set_options(int n, char *inputs[]);
void read_network(const char *file, OPF::Network& network);
void write_network(const char *file, OPF::Network &network);
void write_device(const char *file, OPF::Device &device);

// main function
int main(int argc, char *argv[]) {
  // TODO: introduce a vector object so i can have access to length?
  
  // TODO: create a struct or class that contains pointers to the relevant  
  // "math" data for nodes and lines
  // 
  // we'll create a "network" class that consists of the arrays that contain
  // the needed information
  //
  // the "read" function will use YAJL to parse JSON in to the arrays
  //
  // after reading, the JSON will disappear--only the relevant data remains
  
  // optionally, it will do a "lookup" to find the associated prox functions
  // for each element in the network
  //
  // vectorization will only support quadratics (can't guarantee that other 
  // indices use the same prox function, so can't vectorize)
  
  //printf("%d\n", omp_get_max_threads());
  
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  opf_options O = set_options(argc, argv);

  omp_set_num_threads(O.threads);

  OPF::Network network;

  cout << "Reading network... ";
  read_network(O.file, network);
  cout << "done!" << endl;
  
  // set up solvers
  OPSPSolver admm_full(network, MAX_TIME, O.iters, O.rho, OPSPSolver::FULL, O.verbose, O.warmstart);
  
  opf_utils::print_stats(network);

  
  // pre-solve with stacked generators
  // min sum f_i(p_i)
  // s.t. sum p_i = 0
  
  //cout << "Running presolver" << endl;
  //admm_lossless.solve(O.iters, O.rho);
  //cout << endl;
  
  // solve with full topology
  cout << "Running solver" << endl;
  admm_full.solve(true);
  
  // XXX: File IO only occurs if we try to save for warmstarting
  // write_network(O.file, network);
    
  if(O.save) {
    const Device **devs = admm_full.devices();
    
    // write out device results
    char buffer[64];
    for(int i = 0; i < admm_full.num_devices(); i++)
    {
      // write out device data for warmstart
      sprintf(buffer, "devices/dev%06d",i);
      OPF::Device device;
      for(int j = 0; j < devs[i]->num_terms(); j++) {
        OPF::Device::Terminal* term = device.add_terminals();
        for(int k = 0; k < MAX_TIME; k++)
        {
          term->add_p(devs[i]->terminal(j).p(k));
          term->add_u(devs[i]->terminal(j).u(k)*admm_full.rho_val());
          term->add_theta(devs[i]->terminal(j).theta(k));
          term->add_v(devs[i]->terminal(j).v(k)*admm_full.rho_val());
        }
      }  
      write_device(buffer, device);
    }
  }
  
    
  // Optional:  Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary(); 
  return 0;
}

opf_options set_options(int n, char *inputs[])
{
  int opt;
  // default options
  opf_options O = {"network", 1.0, 2000, omp_get_num_procs(), true, false, false};

  while((opt = getopt(n, inputs, "swqhf:r:n:t:")) != -1) {
    switch(opt) {
      case 'h': opf_utils::usage(); exit(0); break;
      case 'f': O.file = optarg; break;
      case 'r': O.rho = atof(optarg); break;
      case 'n': O.iters = atoi(optarg); break;
      case 't': O.threads = atoi(optarg); break;
      case 'q': O.verbose = false; break;
      case 'w': O.warmstart = true; break;
      case 's': O.save = true; break;
      default: opf_utils::usage(); exit(0); break;
    }
  }
  
  return O;
}

void read_network(const char *file, OPF::Network& network)
{
  // Read the network.
  fstream input(file, ios::in | ios::binary);
  if (!network.ParseFromIstream(&input)) {
    cerr << "Failed to parse network" << endl;
    exit(-1);
  }
  input.close();
}

void write_network(const char *file, OPF::Network &network)
{
  // Write the network back to disk.
  fstream output(file, ios::out | ios::trunc | ios::binary);
  if (!network.SerializeToOstream(&output)) {
    cerr << "Failed to write network." << endl;
    exit(-1);
  }
  output.close();
}

void write_device(const char *file, OPF::Device &device)
{
  // Write the network back to disk.
  fstream output(file, ios::out | ios::trunc | ios::binary);
  if (!device.SerializeToOstream(&output)) {
    cerr << "Failed to write devices." << endl;
    exit(-1);
  }
  output.close();
}
