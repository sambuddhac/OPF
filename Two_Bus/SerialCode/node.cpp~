// Member functions for class Node.
#include <iostream>
//#include <iomanip>

// include Node class definition from node.h
#include "node.h"
//#include "generator.h"
//#include "load.h"
//#include "transl.h"

using namespace std;

Node::Node( int idOfNode ) // constructor begins
	: nodeID( idOfNode )
{
	cout << "\nInitializing the parameters of the node with ID: " << nodeID << endl;

	// initialize the connected devices to zero for node
	gConnNumber = 0; // number of generators connected to a particular node
	tConnNumber = 0; // number of transmission lines connected to a particular node
	lConnNumber = 0; // number of loads connected to a particular node
	nodeFlag = 0; // flag to indicate if a particular node has been accounted for by any one device connected to it for calculation of u 

	PDevCount = 0; // initialize number of devices connectedto a node to zero
	P_avg = 0.0; // Initialize average power to zero
	Theta_avg = 0.0; // initialize average angle to zero
	u = 0.0; // initialize power balance price to zero
	v_avg = 0.0; // initialize average value of voltage angle price to zero
	Pinitavg = 0.0; // initialize initial average power to zero

} // constructor ends

Node::~Node() // destructor
{
	//cout << "\nThe node object having ID " << nodeID << " have been destroyed.\n";

} // end of destructor

int Node::getNodeID() // function getNodeID begins
{
	return nodeID; // returns node ID to the caller
} // end of function getNodeID

void Node::setgConn()
{
	++gConnNumber; // increment the number of generators connected by one whenever a generator is connected to the node

}

void Node::settConn()
{
	++tConnNumber; // increment the number of txr lines connected by one whenever a txr line is connected to the node

}

void Node::setlConn()
{
	++lConnNumber; // increment the number of loads connected by one whenever a load is connected to the node

}

double Node::npinitMessage( double Pload ) //const// function npinitMessage begins
{
	Pinitavg = Pinitavg + Pload / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power 
	return Pinitavg; // return initial average power
} // function npinitMessage ends

double Node::devpinitMessage() const// function devpinitMessage begins
{
	return Pinitavg; // return the initial average node power imbalance to the devices
} // function devpinitMessage ends

void Node::powerangleMessage( double Power, double AngPrice, double Angle ) //const // function powerangleMessage begins
{
	P_avg = P_avg + Power / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average power 
	//**v_avg = v_avg + AngPrice / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average voltage angle price
	Theta_avg = Theta_avg + Angle / ( gConnNumber + tConnNumber + lConnNumber ); // calculate average voltage angle
	PDevCount++; // increment device count by one indicating that a particular device connected to the node has been taken into account
} // function powerangleMessage ends

double Node::PavMessage() const // function PavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) )  // if all the devices are taken care of return the average power
		return P_avg;
} // function PavMessage ends

double Node::uMessage() //const // function uMessage begins
{
	if ( nodeFlag != 0 ) {
		//cout << nodeFlag << endl;
		return u;
	} 
	else {
		if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) ) {
			u = u + P_avg; 		
			//cout << nodeFlag << endl;
			++nodeFlag; // this node has already been accounted for
			return u;  // if all the devices are taken care of calculate and return the power price
		}
	}

} // function uMessage ends

double Node::ThetaavMessage() const // function ThetaavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) ) 
		return Theta_avg;  // if all the devices are taken care of return the average angle
} // function ThetaavMessage ends

double Node::vavMessage() const // function vavMessage begins
{
	if ( PDevCount == ( gConnNumber + tConnNumber + lConnNumber ) )  
		return v_avg;  // if all the devices are taken care of return the average angle price
} // function vavMessage ends

	
void Node::reset() // function reset begins
{
	PDevCount = 0;
	P_avg = 0.0;
	v_avg = 0.0;
	Theta_avg = 0.0;
	nodeFlag = 0;
} // function reset ends

		

