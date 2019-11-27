#include "opsp_solver.hpp"
//#include "utils.hpp"
#include <assert.h>
#include <iomanip>

static clock_t tic_timestart;

void tic(void) {
  tic_timestart = clock();
}

float toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.2f.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

float tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

void OPSPSolver::warmstart_device(Device *dev, int k)
{
  sprintf(buffer, "devices/dev%06d",k);
  OPF::Device device;
  read_device(buffer, device);
  assert(device.terminals_size() == dev->num_terms());

  for(int l = 0; l < device.terminals_size(); l++) {
    OPF::Device::Terminal term = device.terminals(l);
    for(int j = 0; j < term.u_size(); j++) {
      dev->mutable_terminal(l)->oldp()[j] = term.p(j);
      dev->mutable_terminal(l)->p()[j] = term.p(j);
      dev->mutable_terminal(l)->u()[j] = term.u(j)/rho;
            
      dev->mutable_terminal(l)->oldtheta()[j] = term.theta(j);
      dev->mutable_terminal(l)->theta()[j] = term.theta(j);
      dev->mutable_terminal(l)->v()[j] = term.v(j)/rho;
    }
  }
}

void OPSPSolver::setup_lossless_star_network()
{
  _num_devices = network.bus_size() - num_distribution_bus();
  _num_terminals = _num_devices;
  _num_nets = 1;

  _devices = new Device*[_num_devices];
  _nets = new Net*[_num_nets];

  for(int i = 0; i < _num_nets; i++)
    _nets[i] = new Net();

  int k = 0;
  // single-terminal devices
  for(int i = 0; i < network.bus_size(); i++)
  {
    const OPF::Bus &bus = network.bus(i);
    Device* tmp = add_device(bus);
    if(tmp != NULL) {
      _devices[k] = tmp;
      if(warmstart) warmstart_device(_devices[k], k);
      _nets[0]->add_terminal(&_devices[k]->terminal(0));
      k++;
    }
  }
}

void OPSPSolver::setup_lossless_quadratic_network()
{
  int num_one_term_devices = network.bus_size() - num_distribution_bus();
  
  _num_devices = num_one_term_devices + network.line_size();
  _num_terminals = num_one_term_devices + 2*network.line_size();
  _num_nets = network.bus_size();

  _devices = new Device*[_num_devices];
  _nets = new Net*[_num_nets];

  for(int i = 0; i < _num_nets; i++)
    _nets[i] = new Net();


  int k = 0;
  for(int i = 0; i < network.bus_size(); i++)
  {
    const OPF::Bus &bus = network.bus(i);
    Device* tmp = add_device(bus);
    if(tmp != NULL) {
      _devices[k] = tmp;
      if(warmstart) warmstart_device(_devices[k], k);
      _nets[i]->add_terminal(&_devices[k]->terminal(0));
      k++;
    }
  }
        
  for(int i = 0; i < network.line_size(); i++)
  {
    const OPF::Network::Pair& map = network.map(i);
  
    k = num_one_term_devices + i;
    
    _devices[k] = new LosslessQuadraticLine(EPSILON,T);
    
    // warmstart
    if(warmstart) warmstart_device(_devices[k], k);
    
    _nets[map.bus1()]
      ->add_terminal(&_devices[k]->terminal(0));
    _nets[map.bus2()]
      ->add_terminal(&_devices[k]->terminal(1));
  }
}

void OPSPSolver::setup_network()
{
  int num_one_term_devices = network.bus_size() - num_distribution_bus();
  
  _num_devices = num_one_term_devices + network.line_size();
  _num_terminals = num_one_term_devices + 2*network.line_size();  // lines have 2 terminals
  _num_nets = network.bus_size();

  _devices = new Device*[_num_devices];
  _nets = new Net*[_num_nets];

  for(int i = 0; i < _num_nets; i++)
    _nets[i] = new Net();

  int k = 0;
  for(int i = 0; i < network.bus_size(); i++)
  {
    const OPF::Bus &bus = network.bus(i);
    Device* tmp = add_device(bus);
    if(tmp != NULL) {
      _devices[k] = tmp;
      if(warmstart) warmstart_device(_devices[k], k);
      _nets[i]->add_terminal(&_devices[k]->terminal(0));
      k++;
    }
  }
        
  int dc_count = 0, ac_count = 0;
  for(int i = 0; i < network.line_size(); i++)
  {
    const OPF::Network::Pair& map = network.map(i);
    const OPF::Line& line = network.line(i);
  
    k = num_one_term_devices + i;
    // 0 - DC line, 1 - AC line
    if(line.type() == OPF::Line::AC) {
      _devices[k] = new LossyLine(line.g(), line.b(), line.m(), T);
      ac_count++;
    } else {
      _devices[k] = new LossyDCLine(line.g(), line.b(), line.m(), T);
      dc_count++;
    }

    // warmstart
    if(warmstart) warmstart_device(_devices[k], k);
    
    _nets[map.bus1()]
      ->add_terminal(&_devices[k]->terminal(0));
    _nets[map.bus2()]
      ->add_terminal(&_devices[k]->terminal(1));
  }
  printf("!!!!!!!!! dc: %d, ac %d\n", dc_count, ac_count);
}

double OPSPSolver::evaluate_objective()
{
  double obj = 0;
  for(int i =0; i < _num_devices; i++) {
    obj += _devices[i]->objective();
  }
  return obj;
}

int OPSPSolver::solve(bool use_stopping)
{
  // return code
  int code = 1;
  bool once = true;
  
  //float timings[_num_devices];
  
  // for rho updates
  double v = rho - 1;
  double rho_old = rho;
  
  double primal_residual, dual_residual;

  for(int iter = 0; iter < max_iters; ++iter) {
    primal_residual = 0.0;
    dual_residual = 0.0;
    
    #pragma omp parallel for
    for(int i = 0; i < _num_devices; ++i) {
      _devices[i]->solve(rho, rho_old);
    }
    
    // float max_time = 0;
    //     for(int i = 0; i < _num_devices; i++) {
    //       max_time = timings[i] >= max_time ? timings[i] : max_time;
    //     }
    //     printf("TMAX: %0.15f\n", max_time/1000);
    
    //#pragma omp parallel for
    for(int i = 0; i < _num_nets; ++i) {
      _nets[i]->reduce(primal_residual, dual_residual);  // modifies values on terminals
    }
    
    //printf("%f %f\n", sqrt(primal_residual), sqrt(dual_residual));
    primal_residual = sqrt(primal_residual)/sqrt(2*_num_terminals*T);
    dual_residual = rho*sqrt(dual_residual)/sqrt(2*_num_terminals*T);
    
    if(verbose) {
      std::cout << iter << " ";
      std::cout << std::setprecision(15) << evaluate_objective() << " ";
      std::cout << primal_residual << " ";
      std::cout << dual_residual << " ";
      std::cout << rho << " ";
      std::cout << v << std::endl;
    }
    
    // stopping criterion
    if(primal_residual < this->tolerance && dual_residual < this->tolerance) {
      code = 0;
      if(once) {
        std::cout << "done" << std::endl;
        once = false;
      }
      if(use_stopping) break;
    }
    
    // rho update
/*
    double vold = v;
    v = (rho*primal_residual/dual_residual) - 1.0;
    
    rho_old = rho;
    rho = exp(LAMBDA*v + MU*(v - vold))*rho;
    
    if(rho > 1.0/this->tolerance) rho = 1.0/this->tolerance;
    if(rho < this->tolerance) rho = this->tolerance;
*/
  }
  
  //rho = rho_old;
  // copy results in to network
  int k = 0;
  // for(int i = 0; i < _num_nets; i++) {
  //   _nets[i]->print();  // modifies values on terminals
  // }
  
  // for(int i = 0; i < network.bus_size(); i++) {
  //   OPF::Bus *bus = network.mutable_bus(i);
  //   if(bus->type() != OPF::Bus::DIS) {
  //     for(int j = 0; j < T; j++) {
  //       double pval = _devices[k]->terminal(0).p(j);
  //       double uval = _devices[k]->terminal(0).u(j)*rho;
  //       bus->set_p(j, pval);
  //       bus->set_u(j, uval);
  //     }
  //     k++;
  //   }
  // }
  // 
  // for(int i = 0; i < network.line_size(); i++) {
  //   OPF::Line *line = network.mutable_line(i);
  //   for(int j = 0; j < T; j++) {
  //     double pval1 = _devices[k]->terminal(0).p(j);
  //     double uval1 = _devices[k]->terminal(0).u(j)*rho;
  //     double pval2 = _devices[k]->terminal(1).p(j);
  //     double uval2 = _devices[k]->terminal(1).u(j)*rho;
  //     line->set_p1(j, pval1);
  //     line->set_u1(j, uval1);
  //     line->set_p2(j, pval2);
  //     line->set_u2(j, uval2);
  //   }
  //   k++;
  // }
  
  // trying to double-check (not sure if it's right..)
  // k = network.bus_size();
  // for(int i = 0; i < network.line_size(); i++) {
  //   OPF::Line *line = network.mutable_line(i);
  //   //for(int j = 0; j < T; j++) {
  //     double theta1 = _devices[k]->terminal(2).p(0);
  //     //double uval1 = _devices[k]->terminal(0).u(j)*rho;
  //     double theta2 = _devices[k]->terminal(3).p(0);
  //     //double uval2 = _devices[k]->terminal(1).u(j)*rho;
  //     std::cout << "line " << k - network.bus_size() << std::endl;
  //     std::cout << theta1 << " " << theta2 << std::endl;
  //   //}
  //   k++;
  // }
  
  // *****other dianostics***** //
  double obj = evaluate_objective();
  double gen = 0;
  double consumed = 0;
  k = 0;
  for(int i =0; i < network.bus_size(); i++) {  
    const OPF::Bus& bus = network.bus(i);
    if(bus.type() != OPF::Bus::DIS) {
      for(int j = 0; j < T; j++) {
        double val = _devices[k]->terminal(0).p(j);
        gen += (val <= 0) ? val : 0;
        consumed += (val >= 0) ? val : 0;
      }
    }
    if(bus.type() != OPF::Bus::DIS) k++;
  }
  
  // for(int i = 0; i < _num_devices; i++) {
  //   for(int j = 0; j < _devices[i]->num_terms(); j++) {
  //     for(int k = 0; k < T; k++) {
  //       double val = _devices[i]->terminal(j).p(k);
  //       std::cout << val << ", ";
  //     }
  //   std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  if(verbose) {
    std::cout << "--> " << obj << " " << -gen << " " << -consumed << std::endl;
    std::cout << std::endl;
    std::cout << _num_terminals << std::endl;
  }
  return code;
}

Device* OPSPSolver::add_device(const OPF::Bus &bus)
{
  Device *dev = NULL;
  
  switch(bus.type()) {
    case OPF::Bus::DIS:
      break;
    case OPF::Bus::GEN:
    {
      dev = new Generator(
        bus.gencost().alpha(), bus.gencost().beta(),
        bus.gencost().pmax(), bus.gencost().pmin(),
        bus.gencost().slew(), T
      );
      break;
    }
    case OPF::Bus::FIX:
    {
      double load_profile[T];
      for(int j = 0; j < bus.loadcost().l_size(); j++)
        load_profile[j] = bus.loadcost().l(j);

      dev = new FixedLoad(load_profile, T);
      break;
    }
    case OPF::Bus::BAT:
    {
      dev = new Battery(
        bus.batcost().q0(), bus.batcost().qf(),
        bus.batcost().q(), bus.batcost().p(), T
      );
      break;
    }
    case OPF::Bus::DEF:
    {
      dev = new DeferrableLoad(bus.defcost().e(),
        bus.defcost().pmin(), bus.defcost().pmax(),
        bus.defcost().tstart(), bus.defcost().tend(), T
      );
      break;
    }
    case OPF::Bus::CUR:
    {
      dev = new CurtailableLoad(bus.curcost().pdes(), 
          bus.curcost().gamma(), T);
      break;
    }
    default:
      break;// nothing
  }

  return dev;
}

int OPSPSolver::num_distribution_bus()
{
  int count = 0;
  for(int i = 0; i < network.bus_size(); i++)
  {
    const OPF::Bus &bus = network.bus(i);
    if(bus.type() == OPF::Bus::DIS)
      count++;
  }
  return count;
}
