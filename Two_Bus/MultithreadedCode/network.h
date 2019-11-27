// Network class definition.
#ifndef NETWORK_H
#define NETWORK_H

// include the definitions of classes, generator, load, transmission line, node
#include "generator.h"
#include "load.h"
#include "transl.h"
#include "node.h"
#include <vector>

using namespace std;

class Network {

public:
	Network( int ); // constructor
	~Network(); // destructor
	void setNetworkVariables( int ); // sets variables of the network
	void runSimulation(); // runs the OPF simulation
	void generatorThread( int, int, int, double [], double [], double, double [], double [], double [], int ); // operator overloading definition for running Generators' distributed optimization problems

private:
	// Define the Network
	int networkID; // ID number of the network instance
	bool Verbose; // Parameter to decide whether to display intermediate results or not
	const int genNumber, genFields; // Number of Generators & Fields
	const int loadNumber, loadFields; // Number of Loads & Fields
	const int translNumber, translFields; // Number of Transmission Lines & Fields
	const int deviceTermCount; // Number of device terminals
	const int nodeNumber; // Number of Nodes	
	double Rho; // ADMM Tuning Parameter
	//struct generatorCalculator; // structure for multithreaded calculation of generators'optimization problem
	//generatorCalculator genCluster; // instantiate an object of the struct, generatorCalculator
	// Create vectors of Generators, Loads, Transmission lines and Nodes
	vector< Generator > genObject;
	vector< Load > loadObject;
	vector< transmissionLine > translObject;
	vector< Node > nodeObject;

}; // end class Network

#endif // NETWORK_H
