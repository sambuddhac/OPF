#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <omp.h>

#include "opsp_solver.hpp"

#include "utils.hpp"
#include "power_network.pb.h"

using namespace std;

// this file constructs a random network of size N

typedef struct  {
  const char *file;
  int threads;  // number of threads to use
  int N;    // number of nodes
  double p; // probability of link formation
  double d; // distance within which to form link
} gen_options;

inline double rand_d(double fmin, double fmax);

gen_options set_options(int n, char *inputs[]);
void write_network(const char *file, OPF::Network &network);
void write_device(const char *file, OPF::Device &device);

void make_connected(OPF::Network &network);

// XXX. time to twiddle with parameters!
void random_generator(OPF::Generator* gen) {
  //XXX: GENERATOR changed
  int gen_size = rand() % 3;
  switch(gen_size) {
    case 0:
      // small
      gen->set_alpha(0.02);
      gen->set_beta(1);
      gen->set_pmax(10);
      gen->set_pmin(0.01);//rand_d(0,1));
      gen->set_slew(10);//10
      break;
    case 1:
      // medium
      gen->set_alpha(0.005);
      gen->set_beta(0.2);
      gen->set_pmax(20);
      gen->set_pmin(0.01);//rand_d(1,8));
      gen->set_slew(5);//5
      break;
    case 2:
      // large
      gen->set_alpha(0.001);
      gen->set_beta(0.1);
      gen->set_pmax(50);
      gen->set_pmin(0.01);//rand_d(4,12));
      gen->set_slew(3);//3
      break;
    default:
      // make it small
      gen->set_alpha(0.02);
      gen->set_beta(1);
      gen->set_pmax(20);
      gen->set_pmin(0.01);//rand_d(0,1));
      gen->set_slew(2);
  } 
  // gen->set_alpha(rand_d(0.001, 0.002));
  // gen->set_beta(rand_d(0.05,0.15));
  // gen->set_pmax(rand_d(96,104));
  // gen->set_pmin(0);
  // gen->set_slew(200);
}

void random_fixload(OPF::FixedLoad* fix) {
  double a = rand_d(0.5,1); // amplitude
  double c = a + rand_d(0,0.1); // idle power
  static const double f = 2.0*M_PI/MAX_TIME;
  int t0 = 60 + rand() % 13;  // unif_int(24,36)
  
  for(int i = 0; i < MAX_TIME; i++) {
    fix->add_l(c + a*sin(f*(i - t0)));
  }
}

void random_battery(OPF::Battery* bat) {
  // XXX: battery changed!
  double qmax = rand_d(20,50);
  bat->set_q(qmax);
  bat->set_q0(0);
  bat->set_qf(0);  // this is ignored in the solver
  bat->set_p(rand_d(5,10));
}

void random_deferrable(OPF::DeferrableLoad* def) {
  // start in the first half of time
  static const int min_len = 10;
  int start = rand() % ((MAX_TIME - min_len)/2);
  int end = start + min_len + rand() % ((MAX_TIME - start - min_len)/2);
  //printf("!!!!!!!!!! %d %d\n", start, end);
  // XXX: this has changeD!!!
  def->set_tstart(start);
  def->set_tend(end);
  def->set_e(rand_d(5,10));//25,50));
  def->set_pmin(0);
  def->set_pmax(5.0*def->e()/(end - start+1));
}

void random_curtailable(OPF::CurtailableLoad* cur) {
  cur->set_gamma(rand_d(0.1,0.2));
  cur->set_pdes(rand_d(5,15));
}


// double gPmax[3] = {0.7867, 0.6087, 0.0222};
// double gPmin[3] = {-0.4662, -0.8677, -0.4536};
// 
// double qq[3] = {0.3344, -0.1199, 0.8040};
// double pp[3] = {0.5719, 0.8215, 0.2388};

void add_random_bus(OPF::Bus* bus, int i, double sqrt_N) {
  char buffer[33];
  sprintf(buffer, "%d", i);
  bus->set_name(buffer); // doesn't matter
  
  // // roll dice
  // int elim[] = {3,4,5};
  //   int roll = elim[0];
  //   
  //   // disallow generation of hubs
  //   while(roll == elim[0]) {
  //     roll = (rand() % NUMTYPES);
  //     for(int i = 0; i < 3; i++)
  //       if(roll == elim[i]) roll = elim[0];
  //   }
  
  int type1 = 0;
  double roll = rand_d(0,1);
  // XXX: RATIOS have changed too!
  // 0 - distribution, 1 - generator, 2 - fixed load, 3 - battery, 4 - deferrable load, 5 - curtailable load
  type1 = (roll <= 0.4) ? 0 : type1;  
  type1 = (0.4 <= roll && roll <= 0.8) ? 1 : type1;
  //type1 = (0.79 <= roll && roll <= 0.8) ? 2 : type1;
  type1 = (0.8 <= roll && roll <= 0.85) ? 3 : type1;
  type1 = (0.85 <= roll && roll <= 0.9) ? 4 : type1;
  type1 = (0.9 <= roll) ? 5 : type1;
    
    
  OPF::Bus::BusType type = OPF::Bus::BusType(type1);
  if(OPF::Bus::BusType_IsValid(type)) {
    bus->set_type(type);
    
    switch(type) {
      case OPF::Bus::DIS:
        break;
      case OPF::Bus::GEN:
        random_generator(bus->mutable_gencost());
        break;
      case OPF::Bus::FIX:
        random_fixload(bus->mutable_loadcost());
        break;
      case OPF::Bus::BAT:
        random_battery(bus->mutable_batcost());        
        break;
      case OPF::Bus::DEF:
        random_deferrable(bus->mutable_defcost());        
        break;
      case OPF::Bus::CUR:
        random_curtailable(bus->mutable_curcost());        
        break;
      default:
        // do nothing
        break;
    }
  } else {
    printf("ERROR!\n");
    exit(1);
  }
  
  // random location
  OPF::Bus::Point *p = bus->mutable_location();
  p->set_x(rand_d(0,sqrt_N));
  p->set_y(rand_d(0,sqrt_N));
}

void add_random_link(uint i, uint j, OPF::Network &network) {
  #pragma omp critical
  {
  OPF::Line *line = network.add_line();
  OPF::Network::Pair *map = network.add_map();
  
  map->set_bus1(i); map->set_bus2(j);
  line->set_name("abc_abc_link"); // doesn't matter 
  }
}

void set_line_capacity(uint i, OPF::Network &network, const double& C)
{
  // this C is Cmax |p1-p2|/2 <= Cmax
 
  // set line capacity
  double max_phase = rand_d(1,5)*M_PI/180.0;	// maximum phase difference between 1 and 5 degrees
  
  double loss = C*rand_d(0.01,0.03)*sqrt(2);   // max loss percentage
  
  // p_1 + p_2 (loss) = 2*g*(1 - cos(max_phase))
  double g = loss/(2.0*(1.0 - cos(max_phase)));

  // C = b*sin(max_phase)
  double b = C / sin(max_phase);
  
  //#pragma omp critical
  //{ 
  cout << g << " + j" << b << " : " << C/b << " and " << C << endl;
 
  OPF::Line *line = network.mutable_line(i);
    
  if(rand_d(0,1) <= 0.1)
    line->set_type(OPF::Line::DC);
  else
    line->set_type(OPF::Line::AC);
  line->set_g(g);
  line->set_b(b);
  line->set_m(C);

  //}
}

inline double bus_distance(const OPF::Bus &bus1, const OPF::Bus &bus2)
{
  double x = bus1.location().x() - bus2.location().x();
  double y = bus1.location().y() - bus2.location().y();
  return sqrt(x*x + y*y);
}

// main function
//static unsigned char fileData[65536];

int main(int argc, char *argv[]) {	
  //srand(3);//time(NULL));
  srand(time(NULL));
  
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  
  gen_options O = set_options(argc, argv);
  omp_set_num_threads(O.threads);
  
  OPF::Network network;

  int repeat = 1;
  int l;
  char buffer[64];

  OPSPSolver *admm_lossless = NULL;

  // create a device directory (if there isn't already one)
  printf("creating directory\n");
    
  system("rm -rf devices");
  system("mkdir devices");
   
  // XXX: write buses and lines to separate files in a folder
  // the file "network" contains all the info
   
  while(repeat) {
    if(admm_lossless) {
	delete admm_lossless;
    }

    network.Clear();
    // add random buses
    printf("adding buses\n");
    for(int i = 0; i < O.N; i++) {
      add_random_bus(network.add_bus(), i, sqrt(O.N));
    }

    // O.p is connection probability within some radius 
    
    printf("adding lines\n");
    int *neighbors = (int *) calloc(network.bus_size(), sizeof(int));
    #pragma omp parallel for
    for(int i = 0; i < O.N; i ++) {
      double closest_dist = O.N;
      int closest = 0;
      const OPF::Bus &bus1 = network.bus(i);
      
      for(int j = 0; j < O.N; j++) {
        if(j == i) continue;
      
        const OPF::Bus &bus2 = network.bus(j);
      
        double dist = bus_distance(bus1, bus2);
      
        closest = (dist <= closest_dist) ? j : closest;
        closest_dist = MIN(dist, closest_dist);
      
        // if within some radius, add with some fixed probability
        // if outside some radius, add with probability proportional to 1/r^2
        if(j < i) {
          if( (dist <= O.d && rand_d(0,1) <= O.p) ||
              (dist > O.d && rand_d(0,1) <= O.p*(O.d*O.d)/(dist*dist)) ) 
          {
              add_random_link(i,j,network);
              #pragma omp atomic
              neighbors[i]++; 
              #pragma omp atomic
	      neighbors[j]++;
          }
        }
      }
        
      // if i'm not connected to anybody, connect the closest neighbor
      if(neighbors[i] == 0) {
        add_random_link(MIN(i, closest), MAX(i, closest), network);
        #pragma omp atomic
        neighbors[i]++; 
        #pragma omp atomic
        neighbors[closest]++;
      }
    }
    free(neighbors);
  
    // make network connected by connecting connected components
    printf("connecting network\n");
    make_connected(network);
    
    // determine capacities by solving lossless problem with small quadratic cost
  
    // XXX: can use this to warm start
    // XXX: should also provide step length
    
    printf("solving for capacities\n");
    admm_lossless = new OPSPSolver(network, MAX_TIME, 2000, 1, OPSPSolver::LOSSLESS_QUADRATIC, true, false, 0.5*ABSTOL);
    repeat = admm_lossless->solve(true);
  } 


  const Device **devs = admm_lossless->devices();
 
  printf("writing files\n");
  // set line capacities
  l = 0;
  #pragma omp parallel for
  for(int i = 0; i < admm_lossless->num_devices(); i++)
  {
    // write out device data for warmstart
    OPF::Device device;
    for(int j = 0; j < devs[i]->num_terms(); j++) {
      OPF::Device::Terminal* term = device.add_terminals();
      for(int k = 0; k < MAX_TIME; k++)
      {
        term->add_p(0.0);
        term->add_u(0.0);
      }
    } 
    #pragma omp critical
    {
      sprintf(buffer, "devices/dev%06d",i);
      printf("writing device %i\n", i);
      write_device(buffer, device);
    }
    // set line capacities
    if(devs[i]->num_terms() == 2) {
      #pragma omp critical
      {
        printf("setting line capacity for device %i\n", i);
      }
      double pmax = 0;
      for(int j = 0; j < MAX_TIME; j++) {
        double p1 = devs[i]->terminal(0).p(j);
        double p2 = devs[i]->terminal(1).p(j);
        pmax = fabs(p1 - p2)/(2.0) >= pmax ? fabs(p1 - p2)/(2.0) : pmax;          
      }
      #pragma omp critical
      { 
      set_line_capacity(l, network, MAX(30.0, int(10.0*pmax+0.5)));        
      l++;
      }
    }
  }
  
   
  // show degree of graph
  cout << "Average degree: " << 2.0*network.line_size()/network.bus_size() << endl;
  //delete admm_lossless;

 
  printf("writing network.\n"); 
  write_network(O.file, network);

  // Optional:  Delete all global objects allocated by libprotobuf.
  google::protobuf::ShutdownProtobufLibrary();
  
  // write out everything to a file for use in matlab
  printf("writing nodes and edges for matlab.\n");
  FILE *f = fopen("nodes", "w+");
  for(int i = 0; i < network.bus_size(); i ++) {
    const OPF::Bus &bus = network.bus(i);
    fprintf(f, "%g %g %d\n", bus.location().x(), bus.location().y(), bus.type());
  }
  fclose(f);
  
  f = fopen("edges", "w+");
  for(int i = 0; i < network.line_size(); i ++) {
    const OPF::Network::Pair &map = network.map(i);
    const OPF::Line &line = network.line(i);
    fprintf(f, "%d %d %f\n", map.bus1(), map.bus2(), line.m());
  }
  fclose(f);

  if(admm_lossless)
    delete admm_lossless;

  return 0;
}

inline double rand_d(double fmin, double fmax)
{
  double f = (double) rand() / RAND_MAX;
  return fmin + f*(fmax - fmin);
}

void make_connected(OPF::Network &network)
{
  // IMPORTANT: assumes vertex list is sorted in first element (then in second)
  
  // STEP 1: check if our network is connected
  int **memlocs = (int **) calloc(network.line_size(), sizeof(int*));
  int **reachables = (int **) calloc(network.bus_size(), sizeof(int*));
  
  int count = 0;
  bool connected = false;
  
  int last_connected = 0;
  
  // use pointer arithmetic to traverse |E| exactly once  
  for(int i = 0; i < network.line_size(); i++) 
  {
    const OPF::Network::Pair& map = network.map(i);
    //cout << map.bus1() << " " << map.bus2() << endl;
    int j = (map.bus1() < map.bus2()) ? map.bus1() : map.bus2();
    int k = (map.bus1() > map.bus2()) ? map.bus1() : map.bus2();
    
    if(reachables[j] == NULL && reachables[k] == NULL)
    {
      reachables[j] = (int *) malloc(sizeof(int));
      memlocs[count] = reachables[j]; // keep track of memory location to free them
      count++;
      
      *reachables[j] = count;
      reachables[k] = reachables[j];
      last_connected = j;
    }
    else if(reachables[k] == NULL && reachables[j] != NULL)
      reachables[k] = reachables[j];
    else if(reachables[k] != NULL && reachables[j] == NULL)
      reachables[j] = reachables[k];
    else if (reachables[k] != NULL && reachables[j] != NULL && *reachables[j] <= *reachables[k])
      *reachables[k] = *reachables[j];
    else
      *reachables[j] = *reachables[k];
  }
  
  
  // the sum total should be |V|
  int total = 0;
  for(int i = 0; i < network.bus_size(); i++)
    if(reachables[i] && *reachables[i] == 1) total++;
  
  if(total == network.bus_size())
    connected = true;
  
  // STEP 2: if not connected, randomly connect components
  if(!connected) {
    int last_component = 0;
    for(int i = 0; i < network.bus_size(); i++) {
      if(!reachables[i]) {
        add_random_link(i, last_connected, network);
        last_connected = i;
      } else {
        if(last_component > 0 && last_component != *reachables[i]) {
          add_random_link(i, last_connected, network);
          *reachables[i] = last_component;
        }
        last_component = *reachables[i];
        last_connected = i;
      }
    }
  }

  // free memory locations storing component labels
  for(int i = 0; i < count; i++)
    free(memlocs[i]);
  
  free(memlocs);
  free(reachables);
  
  //return connected;
}

gen_options set_options(int n, char *inputs[])
{
  int opt;
  // default options
  gen_options O = {"network", 10, omp_get_num_procs(), 0.8, 0.15};

  while((opt = getopt(n, inputs, "hn:f:p:d:t:")) != -1) {
    switch(opt) {
      case 'h': opf_utils::generator_usage(); exit(0); break;
      case 'f': O.file = optarg; break;
      case 'n': O.N = atoi(optarg); break;
      case 'p': O.p = atof(optarg); break;
      case 'd': O.d = atof(optarg); break;
      case 't': O.threads = atoi(optarg); break;
      default: opf_utils::generator_usage(); exit(0); break;
    }
  }
  
  return O;
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


