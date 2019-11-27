// Main function for implementing ADMM Based Proximal Message Passing Algorithm for the simplest two bus OPF case in serial mode
#include <iostream>

using namespace std;

#include "network.h" // Network class definition
//extern double Rho; // ADMM tuning parameter; Global Variable
int main() // function main begins program execution
{
	int netID;
	
	cout << "\nEnter the ID number to initialize the network.\n";
	cin >> netID;

	cout << endl << "\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n" << endl << endl;

	Network network( netID ); // create the network instance

	cout << "\n*** NETWORK INITIALIZATION STAGE ENDS ***\n" << endl;

	cout << endl << "\n*** ADMM BASED PROXIMAL MESSAGE PASSING OPF SIMULATION (SERIAL IMPLEMENTATION) BEGINS ***\n" << endl << endl;
	
	network.runSimulation(); // start simulation

	cout << "\n*** OPF SIMULATION ENDS ***\n" << endl;

	return 0; // indicates successful program termination

} // end main
