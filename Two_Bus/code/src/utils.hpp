#ifndef OPF_UTILS_H
#define OPF_UTILS_H

#include <time.h>

#include "globals.h"
#include "power_network.pb.h"

using namespace std;
namespace opf_utils {
  
  static char buffer[1024];
  
  typedef struct
  {
    double max;
    double ave;
    double total;
    char *to_str()
    {
      sprintf(buffer, "%#0.6g | %#0.6g | %#0.6g", max, ave, total);
      return buffer;
    }
  } Stats;
  
  void usage() {
    printf("OPF usage is as follows\n");
    printf("\n");
    printf("./opf [-h] [-f path] [-r double]\n");
    printf("\n");
    printf("OPTIONS\n");
    printf("\t-h\n");
    printf("\t\tPrints this help document.\n");
    printf("\n");
    printf("\t-f path\n");
    printf("\t\tSpecifies the path to a configuration file. Defaults\n");
    printf("\t\tto\"conf.xml\".\n");
    printf("\n");
    printf("\t-r double\n");
    printf("\t\tSpecifies value of rho. Defaults to 0.5.\n");
  }
  
  void generator_usage() {
    printf("CREATE_NETWORK usage is as follows\n");
    printf("\n");
    printf("./create_network [-h] [-n number nodes] [-f file] [-p percentage]\n");
  }
  
  void print_stats(const OPF::Network& network) {
    static const char linesep[] = "------------------------------------------------------------------------";
    
    static const int N = network.bus_size();
    static const int M = network.line_size();
    
    int stats[NUMTYPES];
    double ptot[MAX_TIME];
    double load[MAX_TIME];
    
    for(int i = 0; i < NUMTYPES; i++)
      stats[i] = 0;
    for(int i = 0; i < MAX_TIME; i++) {
      ptot[i] = 0;
      load[i] = 0;
    }
    
    Stats gen_stats = {0,0,0};
    Stats load_stats = {0,0,0};
    Stats cur_stats = {0,0,0};
    Stats def_stats = {0,0,0};
    Stats bat_stats = {0,0,0};
    
    for(int i = 0; i < N; i++) {
      const OPF::Bus& bus = network.bus(i);

      switch(bus.type()) {
        case OPF::Bus::DIS:
          stats[OPF::Bus::DIS]++;
          break;
        case OPF::Bus::GEN:
          stats[OPF::Bus::GEN]++;
          
          gen_stats.ave += (bus.gencost().pmax() + bus.gencost().pmin())/(2.0); // average over all time
          gen_stats.max += bus.gencost().pmax();  // max at some time
          gen_stats.total += bus.gencost().pmax()*MAX_TIME; // total over all time
          
          break;
        case OPF::Bus::FIX:
          stats[OPF::Bus::FIX]++;
          
          for(int j = 0; j < MAX_TIME; j++) {
            load[j] += bus.loadcost().l(j);
//            load_stats.max = MAX(load_stats.max, bus.loadcost().l(j));
            load_stats.ave += bus.loadcost().l(j)/MAX_TIME;
            load_stats.total += bus.loadcost().l(j);
          }
                    
          break;
        case OPF::Bus::BAT:
          stats[OPF::Bus::BAT]++;

          bat_stats.max += bus.batcost().q();
          bat_stats.ave += bus.batcost().q()/(2.0);
          bat_stats.total += bus.batcost().q()*MAX_TIME;
          
          break;
        case OPF::Bus::DEF:
        {
          stats[OPF::Bus::DEF]++;
          
          double p = bus.defcost().e()/(bus.defcost().tend()-bus.defcost().tstart()+1);
          //def_stats.max = MAX(def_stats.max, p);
          
          for(int j = 0; j < MAX_TIME; j++) {
            if(bus.defcost().tstart() <= j && j <= bus.defcost().tend())
              ptot[j] += p;
          }          
          
          def_stats.ave += bus.defcost().e()/MAX_TIME;
          def_stats.total += bus.defcost().e();
          
          break;
        }
        case OPF::Bus::CUR:
          stats[OPF::Bus::CUR]++;
          
          cur_stats.max += bus.curcost().pdes();
          cur_stats.ave += bus.curcost().pdes()/2.0;
          cur_stats.total += bus.curcost().pdes()*MAX_TIME;
          
          break;
        default:
          stats[OPF::Bus::DIS]++;
      }
    }
    
    for(int i = 0; i < MAX_TIME; i++) {
      load_stats.max = MAX(load_stats.max, load[i]);
      def_stats.max = MAX(def_stats.max, ptot[i]);
    }
          
    printf("\nNetwork with %d nodes and %d edges.\n\n", N, M);
    printf("%-20s\t%-4s\t%-7s | %-7s | %-7s\n","", "num", "  max", "  ave", " total");
    printf("%s\n", linesep);
    printf("%-20s\t%-4d\n", "Disribution nodes:", stats[OPF::Bus::DIS]);
    printf("%-20s\t%-4d\t%s\n", "Generator nodes:", stats[OPF::Bus::GEN], gen_stats.to_str());
    printf("%-20s\t%-4d\t%s\n", "Battery nodes:", stats[OPF::Bus::BAT], bat_stats.to_str());
    
    printf("%-20s\t%-4d\t%s\n", "Fixed loads:", stats[OPF::Bus::FIX], load_stats.to_str());
    printf("%-20s\t%-4d\t%s\n", "Deferrable loads:", stats[OPF::Bus::DEF], def_stats.to_str());
    printf("%-20s\t%-4d\t%s\n", "Curtailable loads:", stats[OPF::Bus::CUR], cur_stats.to_str());
        
    printf("%s\n", linesep);
    printf("%-20s\t%5g\n\n", "Average degree:", (double)(M*2)/N);
    
  }
  
  
  void print_generator(const OPF::Generator& gen) {
    cout << "Generator:        ";
    cout << gen.alpha() << " ";
    cout << gen.beta() << " ";
    cout << gen.pmax() << " ";
    cout << gen.pmin() << " ";
    cout << gen.slew() << endl;
  }

  void print_fixload(const OPF::FixedLoad& fix) {
    cout << "Fixed load:       ";
    for(int i = 0; i < fix.l_size(); i++)
      cout << fix.l(i) << " ";
    cout << endl;
  }

  void print_battery(const OPF::Battery& bat) {
    cout << "Battery:          ";
    cout << bat.q0() << " ";
    cout << bat.qf() << " ";
    cout << bat.q() << " ";
    cout << bat.p() << endl;
  }

  void print_deferrable(const OPF::DeferrableLoad& def) {
    cout << "Deferrable load:  ";
  
    cout << def.tstart() << " ";
    cout << def.tend() << " ";
    cout << def.e() << " ";
    cout << def.pmin() << " ";
    cout << def.pmax() << endl;
  }

  void print_curtailable(const OPF::CurtailableLoad& cur) {
    cout << "Curtailable load: ";
    cout << cur.gamma() << " ";
    cout << cur.pdes() << endl;
  }
  
  void print_network(const OPF::Network& network)
  {
    for(int i = 0; i < network.bus_size(); i++) {
      const OPF::Bus& bus = network.bus(i);
      
      std::cout << bus.name() << "\t ";
      switch(bus.type()) {
        case OPF::Bus::DIS:
          printf("Distribution node\n");
          break;
        case OPF::Bus::GEN:
          print_generator(bus.gencost());
          break;
        case OPF::Bus::FIX:
          print_fixload(bus.loadcost());
          break;
        case OPF::Bus::BAT:
          print_battery(bus.batcost());
          break;
        case OPF::Bus::DEF:
          print_deferrable(bus.defcost());
          break;
        case OPF::Bus::CUR:
          print_curtailable(bus.curcost());
          break;
        default:
        printf("Distribution node\n");
      }
    }
    
    // printf("\n\n");
    //     for(int i = 0; i < network.line_size(); i++) {
    //       const OPF::Line& line = network.line(i);
    //       const OPF::Network::Pair& map = network.map(i);
    //       std::cout << line.name() << " connects " << map.bus1() << " " << map.bus2() << std::endl;
    //     }
  }
}

#endif
