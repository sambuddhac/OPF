// Network class definition.
#ifndef NETWORK_H
#define NETWORK_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
// include the definitions of classes, generator, load, transmission line, node
#include "generator.h"
#include "load.h"
#include "transl.h"
#include "node.h"
#include <vector>
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
using namespace std;

class Network {

public:
	Network( int, int ); // constructor
	~Network(); // destructor
	void setNetworkVariables( int ); // sets variables of the network
	void runSimulation(GRBEnv*); // runs the OPF simulation
	void runSimulationCentral(GRBEnv*); // runs the centralized GUROBI based OPF simulation

private:
	// Define the Network
	int networkID; // ID number of the network instance
	int genNumber, genFields; // Number of Generators & Fields
	int loadNumber, loadFields; // Number of Loads & Fields
	int translNumber, translFields; // Number of Transmission Lines & Fields
	int deviceTermCount; // Number of device terminals
	int nodeNumber; // Number of Nodes	
	double Rho; // ADMM Tuning Parameter
	bool Verbose; // Parameter to decide whether to display intermediate results or not
	int solverChoice; // 1 for GUROBI, 2 for CVXGEN, 3 for centralized
	// Names of output files string variables
	vector< int > connNodeNumList; // List of identifiers of nodes for the 300 bus system
	vector< int > nodeValList; // List of assigned numerical ranking/order/serial numbers to the nodes for the 300 bus system
	int assignedNodeSer=0; // assigned node serial for nodes of 330 bus system, initialized to zero
	string matrixResultString;
	string devProdString;
	string iterationResultString; 
	string lmpResultString;
	string objectiveResultString;
	string primalResultString;
	string dualResultString;
	// Create vectors of Generators, Loads, Transmission lines and Nodes
	vector< Generator > genObject;
	vector< Load > loadObject;
	vector< transmissionLine > translObject;
	vector< Node > nodeObject;

}; // end class Network

#endif // NETWORK_H
