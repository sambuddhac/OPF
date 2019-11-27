// Generator class definition.
// Member functions defined in generator.cpp
#ifndef GENERATOR_H
#define GENERATOR_H
#include "gensolver.h" // definition of Gensolver class
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
// Include definition of Node class 
//#include "node.h"

class Node; // forward declaration

class Generator {

public:
	Generator( int, Node *, Gensolver & ); // constructor
	~Generator(); // destructor
	int getGenID(); // returns the ID of the Generator
	int getGenNodeID(); // returns the ID of the Node to which the Generator is connected
	void gpowerangleMessage( double, double, double, double, double, double, double ); //const; // Real power and angle iterate
	void gpowerangleMessageGUROBI( double, double, double, double, double, double, double, GRBEnv* ); //const; // Real power and angle iterate
	void setGenData(); // Function to set values of generator data
	double genPower(); //const; // Function to return the value of the power iterate
	double calcPtilde(); //const; // Calculates the difference between Power iterate and average power
	double calcPavInit() const; // gets the Ptilde before iterations start from node
	double calcvtilde() const; // Calculates the difference between v and average v
	double objectiveGen(); // Calculates the objective function value after each iteration of ADMM
	double objectiveGenGUROBI(); // Objective from GUROBI ADMM 
	double getu() const; // Gets the value of price of real power from node
	double calcThetatilde(); //const; // calculates the difference between voltage angle and average voltage angle
	double getv(); // Gets the value of price of voltage angle
	double getPMax();
	double getPMin();
	double getQuadCost();
	double getLinCost();
	double getNLCost();	

private:
	int genID; // Generator object id number
	double Pg, Thetag; // Power and angle iterates
	Node *connNodegPtr; // connection node object
	double v; // Voltage angle constraint Lagrange Multiplier
	Gensolver genSolver; // Generator Solver object	
	double objOpt; // Optimal value of objective

}; //end class Generator

#endif // GENERATOR_H
	
