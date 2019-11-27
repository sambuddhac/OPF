#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>


#include "utils.hpp"
#include "power_network.pb.h"


using namespace std;

typedef struct  {
  const char *file;
  double scale;
} opf_options;

opf_options set_options(int n, char *inputs[]);

inline double rand_d(double fmin, double fmax);
double rand_norm(double mu, double sigma);
double log_norm(double mu, double sigma);

void read_network(const char *file, OPF::Network& network);
void write_network(const char *file, OPF::Network &network);


int main(int argc, char *argv[]) {

	
  srand(123);
  
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  
  opf_options O = set_options(argc, argv);

  OPF::Network network;

  read_network(O.file, network);

  // perturb load profiles with lognormal scaling
  for(int i = 0; i < network.bus_size(); i++)
  {
    OPF::Bus *bus = network.mutable_bus(i);
    if(bus->type() == OPF::Bus::FIX)
    {
      double s = log_norm(0, O.scale);
      for(int j = 0; j < bus->loadcost().l_size(); j++) {
        double val = bus->loadcost().l(j);
        bus->mutable_loadcost()->set_l(j, s*val);
      }
    }
  }
  
  write_network(O.file, network);

  // Optional:  Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();

  
  return 0;
}

opf_options set_options(int n, char *inputs[])
{
  int opt;
  // default options
  opf_options O = {"network", 0.05};

  while((opt = getopt(n, inputs, "f:s:")) != -1) {
    switch(opt) {
      case 'f': O.file = optarg; break;
      case 's': O.scale = atof(optarg); break;
      default: exit(0); break;
    }
  }
  
  return O;
}

inline double rand_d(double fmin, double fmax)
{
  double f = (double) rand() / RAND_MAX;
  return fmin + f*(fmax - fmin);
}

double rand_norm(double mu, double sigma)
{
  double standard_norm = 0;
  for(int i = 0; i < 12; i++)
    standard_norm += rand_d(0,1);
  standard_norm -= 6;
  
  return mu + sigma*standard_norm;
}

double log_norm(double mu, double sigma)
{
  return exp(mu + sigma*rand_norm(0,1));
}

void read_network(const char *file, OPF::Network& network)
{
  // Read the network.
  fstream input(file, ios::in | ios::binary);
  if (!network.ParseFromIstream(&input)) {
    cerr << "Failed to parse network" << endl;
    exit(-1);
  }
}

void write_network(const char *file, OPF::Network &network)
{
  // Write the network back to disk.
  fstream output(file, ios::out | ios::trunc | ios::binary);
  if (!network.SerializeToOstream(&output)) {
    cerr << "Failed to write network." << endl;
    exit(-1);
  }
}



