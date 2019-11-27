// Node class definition.
// Member functions defined in node.cpp
#ifndef NODE_H
#define NODE_H
#include <vector>
using namespace std;

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
	void setgConn(int); // Function to set number of generators connected
	void settConn(int, int, double, int); // Function to set number of transmission lines connected
	void setlConn(int, double); // Function to set number of loads connected
	void reset(); // resets the P_avg, Theta_avg, v_avg to zero after each iteration
	int getGenLength(); // returns the number of connected generators
	int getGenSer(int); // returns the serial number of the generators connected to the node
	int getConNodeLength(); // Returns the number of intra zonal nodes connected to the node
	int getConnSer(int); // Returns the serial number of the intra zonal nodes connected to this node
	double getConnReact(int); // Returns the reactance of the connected internal node
	double getToReact(); // Sum of reciprocal of reactances of all intra and shared exiisting lines for which this is to node
	double getFromReact(); // Sum of reciprocal of reactances of all intra and shared exiisting lines for which this is from node
	double getLoadVal(); // Returns the value of the connected load
private:
	int nodeID; // Node object id number
	int gConnNumber, tConnNumber, lConnNumber; // number of generators, transmission lines and loads connected to the node
	double P_avg, Theta_avg; // average values of power and voltage angle
	double connLoadVal=0; // values of connected load to this node, default value is zero
	double u; // Lagrange multiplier corresponding to power balance
	double v_avg; // average Lagrange multiplier corresponding to voltage angle constraint
	int PDevCount; // Number of devices connected to the particular node object
	double Pinitavg; // initial average power
	int nodeFlag; // node flag to indicate whether u has been calculated in a particular iteration
	vector< int > genSerialNum; // vector consisting of the serial numbers of generators connected to a particular node
	double fromReact; // Sum of reciprocals of reactances of lines for which this is the from node
	double toReact; // Sum of reciprocals of reactances of lines for which this is the to node
	vector< int > connNodeList; // List of intra-zonal nodes that are directly connected to this node via transmission lines
	vector< double > connReactRec; // List of reciprocals of reactances of the intra zone lines connected to the node
	vector< int > tranFromSerial; // vector consisting of the transmission lines for which the node is from node
	vector< int > tranToSerial; // vector consisting of the transmission lines for which the node is to node
	vector< int > loadSerialNum; // vector consisting of the serial numbers of loads connected to a particular node
}; //end class Node

#endif // NODE_H
	
