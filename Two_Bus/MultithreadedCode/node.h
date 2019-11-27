// Node class definition.
// Member functions defined in node.cpp
#ifndef NODE_H
#define NODE_H

//#include <vector>

// forward declarations
//class Generator;
//class Load;
//class transmissionLine; 

class Node {

public:
	Node( int ); // constructor
	~Node(); // destructor
	int getNodeID(); // returns ID of the node to the caller
	double npinitMessage( double ); //const; // calculates initial average power
	double devpinitMessage() const; // returns the initial average node power imbalance to the connected devices before the iterations begin
	void powerangleMessage( double, double, double ); //const; // gets the power, angle and angle price from the devices
	double PavMessage() const; // function to return the average power
	double ThetaavMessage() const; // function to return the average angle 
	double uMessage(); //const; // Real power balance lagrange multiplier iterate
	double vavMessage() const; // Function to return the average price of the voltage angle constraint
	void setgConn(); // Function to set number of generators connected
	void settConn(); // Function to set number of transmission lines connected
	void setlConn(); // Function to set number of loads connected
	void reset(); // resets the P_avg, Theta_avg, v_avg to zero after each iteration

private:
	int nodeID; // Node object id number
	int gConnNumber, tConnNumber, lConnNumber; // number of generators, transmission lines and loads connected to the node
	double P_avg, Theta_avg; // average values of power and voltage angle
	double u; // Lagrange multiplier corresponding to power balance
	double v_avg; // average Lagrange multiplier corresponding to voltage angle constraint
	int PDevCount; // Number of devices connected to the particular node object
	double Pinitavg; // initial average power
	int nodeFlag; // node flag to indicate whether u has been calculated in a particular iteration

	// Create vectors of Generators, Loads, Transmission lines pointers connected to a particular node object
	//vector< Generator* > genObjectPtr;
	//vector< Load* > loadObjectPtr;
	//vector< transmissionLine* > translObjectPtr;

}; //end class Node

#endif // NODE_H
	
