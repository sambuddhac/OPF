// Member functions for class Generator.
#include <iostream>
//#include <iomanip>

// include Generator class definition from generator.h
#include "generator.h"
// Include definition of Node class 
#include "node.h"
// Include Generator solver class defintion
#include "gensolver.h"

using namespace std;

Generator::Generator( int idOfGen, Node *nodeConng, Gensolver &paramOfGen ) // constructor begins
	: genID( idOfGen ),
	  connNodegPtr( nodeConng ),
	  genSolver( paramOfGen )
{
	cout << "\nInitializing the parameters of the generator with ID: " << genID << endl;
	connNodegPtr->setgConn(); // increments the generation connection variable to node
	setGenData(); // calls setGenData member function to set the parameter values

} // constructor ends

Generator::~Generator() // destructor
{
	//cout << "\nThe generator object having ID " << genID << " have been destroyed.\n";

} // end of destructor

int Generator::getGenID() // function getGenID begins
{
	return genID; // returns the ID of the generator object
} // end of getGenID function

int Generator::getGenNodeID() // function getGenNodeID begins
{
	return connNodegPtr->getNodeID(); // returns the ID number of the node to which the generator object is connected
} // end of getGenNodeID function

void Generator::setGenData() // start setGenData function
{
	Pg = 0.0; // Initialize power iterate
	Thetag = 0.0; // Initialize angle iterate
	v = 0.0; // Initialize the Lagrange multiplier corresponding voltage angle constraint to zero
	//double pavi = 0.0;
	
} // end of setGenData function

void Generator::gpowerangleMessage( double gRho, double Pprevit, double Pnetavg, double uprev, double vprevavg, double Aprevavg, double vprev ) //const // function gpowerangleMessage begins
{
	genSolver.mainsolve( gRho, Pprevit, Pnetavg, uprev, vprevavg, Aprevavg, vprev ); // calls the Generator optimization solver
	Pg = genSolver.getPSol(); // get the Generator Power iterate
	Thetag = genSolver.getThetaSol(); // get the Generator voltage angle iterate
	connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle and v
} // function gpowerangleMessage ends

double Generator::genPower() //const // function genPower begins
{
	return Pg; // returns the Pg iterate
} // function genPower ends

double Generator::objectiveGen() // function objectiveGen begins
{
	return genSolver.getObj(); //returns the evaluated objective
} // function objectiveGen ends

double Generator::calcPtilde() //const // function calcPtilde begins
{
	double P_avg = connNodegPtr->PavMessage(); // Gets average power from the corresponding node object
	double Ptilde = Pg - P_avg; // calculates the difference between power iterate and average
	return Ptilde; // returns the difference
} // function calcPtilde ends

double Generator::calcPavInit() const // function calcPavInit begins
{
	return connNodegPtr->devpinitMessage(); // seeks the initial Ptilde from the node
} // function calcPavInit ends

double Generator::getu() const // function getu begins
{
	double u = connNodegPtr->uMessage(); // gets the value of the price corresponding to power balance from node
	//cout << "u: " << u << endl;
	return u; // returns the price
} // function getu ends

double Generator::calcThetatilde() //const // function calcThetatilde begins
{
	//cout << "Thetag: " << Thetag << endl;
	double Theta_avg = connNodegPtr->ThetaavMessage(); // get the average voltage angle at the particular node
	//cout << "Theta_avg: " << Theta_avg << endl;
	double Theta_tilde = Thetag - Theta_avg; // claculate the deviation between the voltage angle of the device and the average
	return Theta_tilde; // return the deviation
} // function calcThetatilde ends

double Generator::calcvtilde() const // function calcvtilde begins
{
	double v_avg = connNodegPtr->vavMessage(); // get the average of the Lagrange multiplier corresponding to voltage angle balance
	//cout << "v_avg: " << v_avg << endl;
	double v_tilde = v - v_avg; // calculate the deviation of the node Lagrange multiplier to the average
	return v_tilde; // return the deviation
} // function calcvtilde ends

double Generator::getv() // function getv begins
{
	//cout << "v_initial: " << v << endl;
	v = v + calcThetatilde(); // Calculate the value of the Lagrange multiplier corresponding to angle constraint
	//cout << "v_final: " << v << endl;
	return v; // Calculate the value of the Lagrange multiplier corresponding to angle constraint
} // function getv ends		
	
