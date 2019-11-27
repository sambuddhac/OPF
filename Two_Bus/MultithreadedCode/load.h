// Load class definition.
// Member functions defined in load.cpp
#ifndef LOAD_H
#define LOAD_H
#include "loadsolver.h" // definition of Loadsolver class
// Include definition of Node class 
//#include "node.h"

class Node; // forward declaration

class Load {

public:
	Load( int, Node *, Loadsolver &, double ); // constructor
	~Load(); // destructor
	int getLoadID(); // returns the ID of the Load
	int getLoadNodeID(); // returns the ID number of the node to which the load is connected
	//double lpowerMessage() const; // Real power iterate
	double pinitMessage(); //const; // Calculated initial Pav
	void lpowerangleMessage( double, double, double, double ); //const; // Real power and Voltage angle iterate
	void setLoadData(); // Function to set load data values
	double calcPtilde(); //const; // Calculates the difference between Power iterate and average power
	double calcPavInit() const; // gets the Ptilde before iterations start from node
	double calcvtilde() const; // Calculates the difference between v and average v
	double getu() const; // Gets the value of price of real power from node
	double calcThetatilde(); //const; // calculates the difference between voltage angle and average voltage angle
	double getv(); //const; // Gets the value of price of voltage angle

private:
	int loadID; // Load object id number
	double Pl; // MW consumption
	double Thetal; // Voltage angle
	//int nl; // node id of the connection node
	Node *connNodelPtr; // connection node object
	double v; // Voltage angle constraint Lagrange Multiplier
	Loadsolver loadSolver; // Load solver object

}; //end class Load

#endif // LOAD_H
	
