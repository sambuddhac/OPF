// Main function for implementing ADMM Based Proximal Message Passing Algorithm for the OPF case in serial mode
#include <iostream>

using namespace std;
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#include "network.h" // Network class definition

int main() // function main begins program execution
{
	int netID; // Network ID number to indicate the type of the system with specifying the number of buses
	cout << "\nEnter the number of nodes to initialize the network. (Allowed choices are 3, 4, 5, 14, 30, 57, 118, and 300 Bus IEEE Test Bus Systems as of now. So, please restrict yourself to one of these)\n";
	cin >> netID;
	int solverChoice; // 1 for GUROBI, 2 for CVXGEN
	cout << "\nEnter the choice of the solver, 1 for GUROBI, 2 for CVXGEN, and 3 for GUROBI centralized" << endl;
	cin >> solverChoice; 
	GRBEnv* environmentGUROBI = new GRBEnv("GUROBILogFile.log"); // GUROBI Environment object for storing the different optimization models
	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;

	Network network( netID, solverChoice ); // create the network instance

	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS ***\n" << endl;
	if ((solverChoice==1) || (solverChoice==2)) {
		cout << endl << "\n*** ADMM BASED PROXIMAL MESSAGE PASSING OPF SIMULATION (SERIAL IMPLEMENTATION) BEGINS ***\n" << endl << endl;
		cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE WAIT. PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE ... ***\n" << endl << endl;
		network.runSimulation(environmentGUROBI); // start simulation
		cout << "\n*** OPF SIMULATION ENDS ***\n" << endl;
	}
	else {
		cout << endl << "\n*** CENTRALIZED OPF SIMULATION BEGINS ***\n" << endl << endl;
		cout << endl << "\n*** SIMULATION IN PROGRESS; PLEASE WAIT. PLEASE DON'T CLOSE ANY WINDOW OR OPEN ANY OUTPUT FILE ... ***\n" << endl << endl;
		network.runSimulationCentral(environmentGUROBI); // start simulation

		cout << "\n*** OPF SIMULATION ENDS ***\n" << endl;	
	}
	delete environmentGUROBI; // Free the memory of the GUROBI environment object
	return 0; // indicates successful program termination

} // end main
