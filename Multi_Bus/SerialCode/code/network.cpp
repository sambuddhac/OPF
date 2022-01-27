// Member functions for class Network
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
// include definitions for classes generator, load, transmission line, network and node
#include "generator.h"
#include "transl.h"
#include "load.h"
#include "node.h"
#include "network.h"
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file
#define MAX_ITER 80002
#define LINE_CAP 100.00
using namespace std;

Network::Network( int val, int solverChoice )
	: networkID( val ), // constructor begins; initialize networkID  and Rho through constructor initializer list
	  Rho( 1.0 ), 
	  solverChoice(solverChoice)
{
	setNetworkVariables( networkID ); // sets the variables of the network

} // end constructor

// destructor
Network::~Network()
{	
	cout << "\nNetwork instance: " << networkID << " for this simulation destroyed. You can now open the output files to view the results of the simulation\n" << endl;

} // end destructor

void Network::setNetworkVariables( int networkID ) // Function setNetworkVariables starts to initialize the parameters and variables
{
	Verbose = false; // disable intermediate result display. If you want, make it "true"
	do {

		nodeNumber = networkID; // set the number of nodes of the network
		char genFile[ 15 ];
		char tranFile[ 15 ];
		char loadFile[ 15 ];
		// Nodes
		for ( int l = 0; l < nodeNumber; ++l ) {
			//cout << "\nCreating the " << l + 1 << " -th Node:\n";
		
			Node nodeInstance( l + 1 ); // creates nodeInstance object with ID l + 1

			nodeObject.push_back( nodeInstance ); // pushes the nodeInstance object into the vector

		} // end initialization for Nodes
		
		switch ( nodeNumber ) { // switch structure determines which system to simulate for

			case 14: // 14 Bus case
				
				strcpy( genFile, "Gen14.txt" );			
				strcpy( tranFile, "Tran14.txt" );
				strcpy( loadFile, "Load14.txt" );				
				break; // exit switch

			case 30: // 30 Bus case

				strcpy( genFile, "Gen30.txt" );			
				strcpy( tranFile, "Tran30.txt" );
				strcpy( loadFile, "Load30.txt" );	
				break; // exit switch

			case 57: // 57 Bus case

				strcpy( genFile, "Gen57.txt" );			
				strcpy( tranFile, "Tran57.txt" );
				strcpy( loadFile, "Load57.txt" );	
				break; // exit switch

			case 118: // 118 Bus case

				strcpy( genFile, "Gen118.txt" );			
				strcpy( tranFile, "Tran118.txt" );
				strcpy( loadFile, "Load118.txt" );	
				break; // exit switch

			case 300: // 300 Bus case

				strcpy( genFile, "Gen300.txt" );			
				strcpy( tranFile, "Tran300.txt" );
				strcpy( loadFile, "Load300.txt" );	
				break; // exit switch

			case 3: // 3 Bus case

				strcpy( genFile, "Gen3.txt" );			
				strcpy( tranFile, "Tran3.txt" );
				strcpy( loadFile, "Load3.txt" );	
				break; // exit switch

			case 4: // 3 Bus case variant

				strcpy( genFile, "Gen3A.txt" );			
				strcpy( tranFile, "Tran3A.txt" );
				strcpy( loadFile, "Load3A.txt" );	
				break; // exit switch

			case 5: // 5 Bus case

				strcpy( genFile, "Gen5.txt" );			
				strcpy( tranFile, "Tran5.txt" );
				strcpy( loadFile, "Load5.txt" );	
				break; // exit switch

			case 48: // 48 Bus case

				strcpy( genFile, "Gen48.txt" );			
				strcpy( tranFile, "Tran48.txt" );
				strcpy( loadFile, "Load48.txt" );	
				break; // exit switch

			
			
			default: // catch all other entries

				cout << "Sorry, invalid case. Can't do simulation at this moment.\n" << endl;
				break; // exit switch

		} // end switch

		// Generators
		ifstream matrixFirstFile( genFile, ios::in ); // ifstream constructor opens the file of Generators

		// exit program if ifstream could not open file
		if ( !matrixFirstFile ) {
			cerr << "\nFile for Generators could not be opened\n" << endl;
			exit( 1 );
		} // end if

		matrixFirstFile >> genNumber >> genFields; // get the dimensions of the Generator matrix
		double matrixGen[ genNumber ][ genFields ]; // Generator matrix
		for ( int i = 0; i < genNumber; ++i ) {
			for ( int j = 0; j < genFields; ++j ) {
				matrixFirstFile >> matrixGen[ i ][ j ]; // read the Generator matrix
			}
		}

		// Transmission Lines
		ifstream matrixSecondFile( tranFile, ios::in ); // ifstream constructor opens the file of Transmission lines

		// exit program if ifstream could not open file
		if ( !matrixSecondFile ) {
			cerr << "\nFile for Transmission lines could not be opened\n" << endl;
			exit( 1 );
		} // end if
	
		matrixSecondFile >> translNumber >> translFields; // get the dimensions of the Transmission line matrix
		double matrixTran[ translNumber ][ translFields ]; // Transmission line matrix
		for ( int i = 0; i < translNumber; ++i ) {
			for ( int j = 0; j < translFields; ++j ) {
				matrixSecondFile >> matrixTran[ i ][ j ]; // read the Transmission line matrix
			}
		}

		// Create Transmission Lines
		for ( int k = 0; k < translNumber; ++k ) {
			int l = 0;
			int tNodeID1, tNodeID2; // node object IDs to which the particular transmission line object is connected
			int tNodeID1300, tNodeID2300; // node object IDs to which the particular transmission line object is connected for 300 bus case
			do {
				if (nodeNumber==300) { // Since for IEEE 300 bus system, the nodes are not serialy numbered, but name-numbered instead, below is conversion code
					tNodeID1300 = matrixTran[ k ][ 0 ]; //From end node identifier
					tNodeID2300 = matrixTran[ k ][ 1 ]; //To end node identifier
					if (std::find(connNodeNumList.begin(), connNodeNumList.end(), tNodeID1300) != connNodeNumList.end()) { // If node identifier value for this particular node is present in the list
						auto pos = std::find(connNodeNumList.begin(), connNodeNumList.end(), tNodeID1300) - connNodeNumList.begin(); // find the position of the node identifier in the chart of node identifiers
						tNodeID1 = nodeValList[pos]; // Get the serial number of the node from the nodeValList
						//cout << "For line " << k+1 << " Identifier of the From Node: " << tNodeID1300 << " From Node assigned Serial: " << tNodeID1 << " REPEATED" << endl;
					}
					else {
						connNodeNumList.push_back(tNodeID1300); // For a new node identifier
						nodeValList.push_back(++assignedNodeSer); // Assign the node serial
						tNodeID1 = assignedNodeSer; // Get the serial number of the node from the nodeValList
						//cout << "For line " << k+1 << " Identifier of the From Node: " << tNodeID1300 << " From Node assigned Serial: " << tNodeID1 << " FRESH" << endl;
					}
					if (std::find(connNodeNumList.begin(), connNodeNumList.end(), tNodeID2300) != connNodeNumList.end()) { // If node identifier value for this particular node is present in the list
						auto pos = std::find(connNodeNumList.begin(), connNodeNumList.end(), tNodeID2300) - connNodeNumList.begin(); // find the position of the node identifier in the chart of node identifiers
						tNodeID2 = nodeValList[pos]; // Get the serial number of the node from the nodeValList
						//cout << "For line " << k+1 << " Identifier of the To Node: " << tNodeID2300 << " To Node assigned Serial: " << tNodeID2 << " REPEATED" << endl;
					}
					else {
						connNodeNumList.push_back(tNodeID2300); // For a new node identifier
						nodeValList.push_back(++assignedNodeSer); // Assign the node serial
						tNodeID2 = assignedNodeSer; // Get the serial number of the node from the nodeValList
						//cout << "For line " << k+1 << " Identifier of the To Node: " << tNodeID2300 << " To Node assigned Serial: " << tNodeID2 << " FRESH" << endl;
					}
				}
				//*cout << "Stuck while creating nodes of transmission line: " << ( k + 1 ) << endl;
				//node IDs of the node objects to which this transmission line is connected.
				else {
					tNodeID1 = matrixTran[ k ][ 0 ]; //From end
					tNodeID2 = matrixTran[ k ][ 1 ]; //To end
				}
			} while ( ( tNodeID1 <= 0 ) || ( tNodeID1 > nodeNumber ) || ( tNodeID2 <= 0 ) || ( tNodeID2 > nodeNumber ) || ( tNodeID1 == tNodeID2) ); // validity check
			double resT, reacT, ptMax, ptMin; // Parameters for Transmission Line
			do {
				//Resistance:
				resT = matrixTran[ k ][ ( l + 2 ) ];
				//Reactance:
				//reacT = matrixTran[ k ][ ( l + 3 ) ];
				double reacT_int = matrixTran[ k ][ ( l + 3 ) ];
				reacT = (pow(reacT_int, 2)+pow(resT, 2))/reacT_int;
				//values of maximum allowable power flow on line in the forward and reverse direction:
				//Forward direction:
				//cout << "\nEnter the pu value of the line flow limit (Assumed to be the same for all the lines)\n";
				ptMax = LINE_CAP;
				//ptMax = matrixTran[ k ][ ( l + 4 ) ]/100;
				ptMin = -ptMax; //Reverse direction
			} while ( (resT < 0 ) || ( reacT <= 0 ) || ( ptMax <= ptMin ) ); // check the bounds and validity of the parameter values

			Transolver tlineParam( ptMin, ptMax, reacT ); // Instantiate the copy constructor for the Transmission Line solver object
			
			// creates transLineInstance object with ID k + 1
			transmissionLine transLineInstance( k + 1, &nodeObject[ tNodeID1 - 1 ], &nodeObject[ tNodeID2 - 1 ], tlineParam, ptMax, reacT, resT ); 

			translObject.push_back( transLineInstance ); // pushes the transLineInstance object into the vector

		} // end initialization for Transmission Lines
		
		// Create Generators
		for ( int i = 0; i < genNumber; ++i ) {
			int j = 0; 
			int gNodeID; // node object ID to which the particular generator object is connected
			int gNodeID300; // node object ID to which the particular generator object is connected for 300 bus system
			do {
				if (nodeNumber==300) { // Since for IEEE 300 bus system, the nodes are not serialy numbered, but name-numbered instead, below is conversion code
					gNodeID300 = matrixGen[ i ][ 0 ]; //Generator node identifier
					if (std::find(connNodeNumList.begin(), connNodeNumList.end(), gNodeID300) != connNodeNumList.end()) { // If node identifier value for this particular node is present in the list
						auto pos = std::find(connNodeNumList.begin(), connNodeNumList.end(), gNodeID300) - connNodeNumList.begin(); // find the position of the node identifier in the chart of node identifiers
						gNodeID = nodeValList[pos]; // Get the serial number of the node from the nodeValList
						//cout << "For Generator " << i+1 << " Identifier of the Conn Node: " << gNodeID300 << " Conn Node assigned Serial: " << gNodeID << " REPEATED" << endl;
					}
					else {
						connNodeNumList.push_back(gNodeID300); // For a new node identifier
						nodeValList.push_back(++assignedNodeSer); // Assign the node serial
						gNodeID = assignedNodeSer; // Get the serial number of the node from the nodeValList
						//cout << "For Generator " << i+1 << " Identifier of the Conn Node: " << gNodeID300 << " Conn Node assigned Serial: " << gNodeID << " FRESH" << endl;
					}
				}
				else {
					gNodeID = matrixGen[ i ][ 0 ];
				}
			} while ( ( gNodeID <= 0 ) || ( gNodeID > nodeNumber ) ); // validity check

			double c2, c1, c0, PgMax, PgMin; // Parameters for Generator
			do {
				//Quadratic Coefficient: 
				c2 = matrixGen[ i ][ ( j + 1 ) ];
				//Linear coefficient: 
				c1 = matrixGen[ i ][ ( j + 2 ) ];///100
				//Constant term: 
				c0 = matrixGen[ i ][ ( j + 3 ) ];
				//Maximum Limit: 
				PgMax = matrixGen[ i ][ ( j + 4 ) ];//*100;
				//Minimum Limit: 
				PgMin = matrixGen[ i ][ ( j + 5 ) ];
			} while ( (c2 < 0 ) || ( c1 < 0 ) || ( PgMax <= 0 ) || ( PgMin < 0 ) || ( PgMax <= PgMin ) ); 
			// check the bounds and validity of the parameter values

			Gensolver genParam( c2, c1, c0, PgMax, PgMin ); // Instantiate the copy constructor for the generator solver object
	
			Generator generatorInstance( i + 1, &nodeObject[ gNodeID - 1 ], genParam ); // creates generatorInstance object with ID number i + 1

			genObject.push_back( generatorInstance ); // pushes the generatorInstance object into the vector

		} // end initialization for Generators

		// Loads
		ifstream matrixThirdFile( loadFile, ios::in ); // ifstream constructor opens the file of Loads

		// exit program if ifstream could not open file
		if ( !matrixThirdFile ) {
			cerr << "\nFile for Loads could not be opened\n" << endl;
			exit( 1 );
		} // end if

		matrixThirdFile >> loadNumber >> loadFields; // get the dimensions of the Load matrix
		double matrixLoad[ loadNumber ][ loadFields ]; // Load matrix
		for ( int i = 0; i < loadNumber; ++i ) {
			for ( int j = 0; j < loadFields; ++j ) {
				matrixThirdFile >> matrixLoad[ i ][ j ]; // read the Load matrix
			}
		}
		// Create Loads
		for ( int j = 0; j < loadNumber; ++j ) {
			//cout << "\nEnter the parameters of the " << j + 1 << " -th Load:\n";
			int k = 0;
			int lNodeID, lNodeID300; // node object ID to which the particular load object is connected
			do {
				if (nodeNumber==300) { // Since for IEEE 300 bus system, the nodes are not serialy numbered, but name-numbered instead, below is conversion code
					lNodeID300 = matrixLoad[ j ][ 0 ]; //Load node identifier
					if (std::find(connNodeNumList.begin(), connNodeNumList.end(), lNodeID300) != connNodeNumList.end()) { // If node identifier value for this particular node is present in the list
						auto pos = std::find(connNodeNumList.begin(), connNodeNumList.end(), lNodeID300) - connNodeNumList.begin(); // find the position of the node identifier in the chart of node identifiers
						lNodeID = nodeValList[pos]; // Get the serial number of the node from the nodeValList
						//cout << "For Load " << j+1 << " Identifier of the Conn Node: " << lNodeID300 << " Conn Node assigned Serial: " << lNodeID << " REPEATED" << endl;
					}
					else {
						connNodeNumList.push_back(lNodeID300); // For a new node identifier
						nodeValList.push_back(++assignedNodeSer); // Assign the node serial
						lNodeID = assignedNodeSer; // Get the serial number of the node from the nodeValList
						//cout << "For Load " << j+1 << " Identifier of the Conn Node: " << lNodeID300 << " Conn Node assigned Serial: " << lNodeID << " FRESH" << endl;
					}
				}
				else {
					//node ID of the node object to which this load object is connected.
					lNodeID = matrixLoad[ j ][ 0 ]; 
				}
			} while ( ( lNodeID <= 0 ) || ( lNodeID > nodeNumber ) ); // validity check

			double P_Load; // Parameters for Load
			do {
				//value of allowable power consumption capability of load with a negative sign to indicate consumption:
				//Power Consumption:
				P_Load = matrixLoad[ j ][ ( k + 1 ) ];//*100;
			} while ( -P_Load <= 0 ); // check the bounds and validity of the parameter values

			Loadsolver loadParam( P_Load );	// Instantiate the copy constructor for the load solver object
			Load loadInstance( j + 1, &nodeObject[ lNodeID - 1 ], loadParam, P_Load ); // creates loadInstance object object with ID number j + 1

			loadObject.push_back( loadInstance ); // pushes the loadInstance object into the vector

		} // end initialization for Loads

	} while ( (genNumber <= 0 ) || ( nodeNumber <= 0 ) || ( loadNumber <= 0 ) || ( translNumber <= 0 ) || (genFields <= 0 ) || ( loadFields <= 0 ) || ( translFields <= 0 ) );
	// check the bounds and validity of the parameter values
	
	deviceTermCount = genNumber + loadNumber + 2 * translNumber; // total number of device-terminals

} // end setNetworkVariables function

// runSimulation function definition
void Network::runSimulation(GRBEnv* environmentGUROBI) //Function runSimulation begins
{
	// Declaration of intermerdiate variables and parameters for running the simulation
	int iteration_count = 1; // iteration counter

	double dualTol = 1.0; // initialize the dual tolerance
	double primalTol, PrimalTol; // primal tolerance
	double ptolsq = 0.0; // initialize the primal tolerance square
	
	vector< int > iterationGraph; // vector of iteration counts
	vector< double > primTolGraph; // vector of primal tolerance
	vector< double > PrimTolGraph;
	vector< double > dualTolGraph; // vector of dual tolerance
	vector< double > objectiveValue; // vector of objective function values

	int bufferIndex; // index of the buffer to store past values of voltage iterate, power and angle iterate

	double V_avg[ nodeNumber ]; // array of average node angle imbalance price from last to last iterate
	double vBuffer1[ nodeNumber ]; // intermediate buffer for average node angle price from last to last iterate
	double vBuffer2[ nodeNumber ]; // intermediate buffer for average node angle price from last iterate

	double angleBuffer[ nodeNumber ]; // buffer for average node voltage angles from present iterate
	double angleBuffer1[ nodeNumber ]; // buffer for average node voltage angles from last iterate
	double angtildeBuffer[ deviceTermCount ]; // Thetatilde from present iterate

	double powerBuffer[ deviceTermCount ]; // Ptilde from present iterate
	double powerBuffer1[ deviceTermCount ]; // Ptilde from last iterate
	double pavBuffer[ nodeNumber ]; // Pav from present iterate
	double ptildeinitBuffer[ deviceTermCount ]; // Ptilde before iterations begin
	//int firstIndex = ( MAX_ITER / 100 ) + 1;
	double uPrice[ deviceTermCount ]; // u parameter from previous iteration
	double vPrice[ deviceTermCount ]; // v parameter from previous iteration
	double LMP[ nodeNumber ]; // vector of LMPs

	double Rho1 = 1.0; // Previous value of Rho from previous iteration
	double W, Wprev; // Present and previous values of W for the PID controller for modifying Rho
	double lambdaAdap = 0.0001; // Parameter of the Proportional (P) controller for adjusting the ADMM tuning parameter
	double muAdap = 0.0005; // Parameter of the Derivative (D) controller for adjusting the ADMM tuning parameter
        double xiAdap = 0.0000; // Parameter of the Integral (I) controller for adjusting the ADMM tuning parameter
        double controllerSum = 0.0; // Integral term of the PID controller
	int setTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho

	// Set the type of tuning
	cout << "Enter the tuning mode; Enter 1 for maintaining Rho * primTol = dualTol; 2 for primTol = dualTol; anything else for Adaptive Rho (with mode-1 being implemented for the first 3000 iterations and then Rho is held constant).\n" << endl;
	cin >> setTuning;

	// Calculation of initial value of Primal Tolerance before the start of the iterations
	vector< Load >::iterator loadIterator;	
	for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
		ptolsq = ptolsq + pow( loadIterator->pinitMessage(), 2.0 ); // calls the node to divide by the number of devices connected
	}
	primalTol = sqrt( ptolsq ); // initial value of primal tolerance to kick-start the iterations
	PrimalTol = primalTol;
	// Calculation of initial value of Ptilde before the iterations start
	vector< Generator >::iterator generatorIterator;	
	for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
		bufferIndex = generatorIterator->getGenID() - 1;
		ptildeinitBuffer[ bufferIndex ] = -generatorIterator->calcPavInit();
	}

	for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
		bufferIndex = genNumber + ( loadIterator->getLoadID() - 1 );
		ptildeinitBuffer[ bufferIndex ] = loadIterator->calcPavInit();
	}

	int temptrans1 = 0; // counter to make sure that two values of Ptilde are accounted for each line
	vector< transmissionLine >::iterator translIterator;	
	for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
		bufferIndex = genNumber + loadNumber + ( translIterator->getTranslID() - 1 ) + temptrans1;
		ptildeinitBuffer[ bufferIndex ] = -translIterator->calcPavInit1(); // Ptilde corresponding to 'from' end
		ptildeinitBuffer[ ( bufferIndex + 1 ) ] = -translIterator->calcPavInit2(); // Ptilde corresponding to 'to' end
		temptrans1++;
	}

	if (solverChoice == 1) {
		matrixResultString = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/Summary_of_Result_Log.txt";
		devProdString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/powerResult.txt";
		iterationResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/itresult.txt"; 
		lmpResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/LMPresult.txt";
		objectiveResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/objective.txt";
		primalResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/primresult.txt";
		dualResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_GUROBI/dualresult.txt";
	}
	else {
		matrixResultString = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/Summary_of_Result_Log.txt";
		devProdString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/powerResult.txt";
		iterationResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/itresult.txt"; 
		lmpResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/LMPresult.txt";
		objectiveResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/objective.txt";
		primalResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/primresult.txt";
		dualResultString="/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/dualresult.txt";
	}
	ofstream matrixResultOut( matrixResultString, ios::out ); // create a new file result.txt to output the results	
	// exit program if unable to create file
	if ( !matrixResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	if ( Verbose ) {
		matrixResultOut << "\nThe initial value of primal tolerance to kick-start iterations is: " << primalTol << "\nThe initial value of dual tolerance to kick-start iterations is: " << dualTol << endl;
	}
	
	clock_t start_s = clock(); // begin keeping track of the time
	//int first = 0;
	// Starting of the ADMM Based Proximal Message Passing Algorithm Iterations
	for ( iteration_count = 1; ( iteration_count < MAX_ITER ); iteration_count++ ) {
	//while( ( primalTol >= 0.006 ) || ( dualTol >= 0.006 ) ) { // ( iteration_count <= 122 )
	
		if ( Verbose ) {
			matrixResultOut << "\nThe value of primal tolerance before this iteration is: " << primalTol << "\nThe value of dual tolerance before this iteration is: " << dualTol << endl;
			matrixResultOut << "\n**********Start of " << iteration_count << " -th iteration***********\n";
		}		
		// Recording data for plotting graphs
		
		iterationGraph.push_back( iteration_count ); // stores the iteration count to be graphed
		primTolGraph.push_back( primalTol ); // stores the primal tolerance value to be graphed
		PrimTolGraph.push_back( PrimalTol ); 
		dualTolGraph.push_back( dualTol ); // stores the dual tolerance value to be graphed
		//Initialize the average node angle imbalance price (v) vector from last to last interation, V_avg
		//**if ( iteration_count <= 2 ) {
			for ( int i = 0; i < nodeNumber; i++ )
				V_avg[ i ] = 0.0; // initialize to zero for the first and second iterations if the initial values are zero
		//**}
		//**else {
			//**for ( int j = 0; j < nodeNumber; j++ )
				//**V_avg[ j ] = vBuffer1[ j ]; // initialize to the average node v from last to last iteration for 3rd iteration on
		
		//**}
		// Initialize average v, average theta, ptilde, average P before the start of a particular iteration
		if ( iteration_count >= 2 ) {
			angleBuffer1[ 0 ] = 0.0; // set the first node as the slack node, the average voltage angle is always zero
			for ( int i = 0; i < nodeNumber; i++ ) {
				//**vBuffer1[ i ] = vBuffer2[ i ]; // Save to vBuffer1, the average v from last iteration for use in next iteration
				angleBuffer1[ i ] = angleBuffer[ i ]; // Save to angleBuffer1, the average node voltage angle from last iteration
			}

			for ( int j = 0; j < deviceTermCount; j++ )
				powerBuffer1[ j ] = powerBuffer[ j ]; // Save to powerBuffer1, the Ptilde for each device term. from last itern

		}
		
		else {
			Wprev = 0.0; // for the first iteration
			for ( int i = 0; i < nodeNumber; i++ ) {
			
				angleBuffer1[ i ] = 0.0; // Set average node voltage angle to zero for 1st iteration
			}

			vector< Node >::iterator nodeIterator;
			for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
				bufferIndex = nodeIterator->getNodeID() - 1;
				pavBuffer[ bufferIndex ] = nodeIterator->devpinitMessage(); // Average node power injection before 1st iteration
			}
			for ( int j = 0; j < deviceTermCount; j++ )
				powerBuffer1[ j ] = ptildeinitBuffer[ j ]; // Save to powerBuffer1, the Ptilde before the 1st iteration
		}

		//vector< Generator >::const_iterator generatorIterator; // Distributed Optimizations; Generators' Opt. Problems
		double calcObjective = 0.0;	// initialize the total generator cost for this iteration
		for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
			double Pgit, PowerPrice, APrice; // Generator Power, Power Price, & Angle Price iterates from last iterations
			bufferIndex = generatorIterator->getGenID() - 1;
			int gnid = generatorIterator->getGenNodeID() - 1; // gets the ID number of connection node
			if ( iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				Pgit = generatorIterator->genPower();
				PowerPrice = uPrice[ bufferIndex ];
				
				if ( gnid == 0 ) {
					APrice = 0.0; // Consider node-1 as the slack node, the angle price is zero always
				}

				else {
					APrice = vPrice[ bufferIndex ];
				}
			}
			else { // If 1st iteration, initialize to zero
				Pgit = 0.0;
				PowerPrice = 0.0;
				APrice = 0.0; 
			}
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Generator Optimization Iterations for Generator " << bufferIndex + 1 << "\n";
				matrixResultOut << "Previous power iterate (MW/pu)\n" << Pgit << "\nPrevious average power (MW/pu) for this node\n" << pavBuffer[ gnid ] << "\nPrevious power price (scaled LMP)\n" << PowerPrice << "\nAngle price from last to last iterate (scaled)\n" << V_avg[ gnid ] << "\nAngle value from last iterate\n" << angleBuffer1[ gnid ] << "\nPrevious angle price (scaled)\n" << APrice << endl;
			}
			if (solverChoice == 1) {
				generatorIterator->gpowerangleMessageGUROBI( Rho, Pgit, pavBuffer[ gnid ], PowerPrice, V_avg[ gnid ], angleBuffer1[ gnid ], APrice, environmentGUROBI ); // Solve the Optimization Problem
				calcObjective = calcObjective + generatorIterator->objectiveGenGUROBI(); // calculate the total objective after this iteration
			}
			else {
				generatorIterator->gpowerangleMessage( Rho, Pgit, pavBuffer[ gnid ], PowerPrice, V_avg[ gnid ], angleBuffer1[ gnid ], APrice ); // Solve the Optimization Problem
				calcObjective = calcObjective + generatorIterator->objectiveGen(); // calculate the total objective after this iteration
			}
		}
		//vector< Load >::const_iterator loadIterator;	// Distributed Optimizations; Loads' Optimization Problems
		for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
			double APrice, PPrice; // Load Power Price and Angle Price from last iterations
			bufferIndex = genNumber + ( loadIterator->getLoadID() - 1 );
			int lnid = loadIterator->getLoadNodeID() - 1; // gets ID number of connection node
			if ( iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				
				if ( lnid == 0 ) {
					APrice = 0.0; // Consider node-1 as the slack node, the angle price is zero always
				}

				else {
					APrice = vPrice[ bufferIndex ];
				}
				PPrice = uPrice[ bufferIndex ];
			}
			else 
				APrice = 0.0; // If 1st iteration, initialize to zero
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Load Optimization Iterations for Load " << loadIterator->getLoadNodeID() << "\n";
				matrixResultOut << "\nAngle price from last to last iterate (scaled)\n" << V_avg[ lnid ] << "\nAngle value from last iterate\n" << angleBuffer1[ lnid ] << "\nPrevious angle price (scaled)\n" << APrice << endl;
			}
			loadIterator->lpowerangleMessage( Rho, V_avg[ lnid ], angleBuffer1[ lnid ], APrice ); // Solve the Optimization Problem
		}
		//vector< transmissionLine >::const_iterator translIterator;// Distributed Optimizations; TLine' Optimization Problems
		int temptrans2 = 0;	
		for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
			double Ptit1, Ptit2, PowerPrice1, PowerPrice2, APrice1, APrice2; // Tline Power, Power price, Angle price at both ends
			bufferIndex = genNumber + loadNumber + ( translIterator->getTranslID() - 1 ) + temptrans2;
			int tnid1 = translIterator->getTranslNodeID1() - 1; // gets ID number of first conection node
			int tnid2 = translIterator->getTranslNodeID2() - 1; // gets ID number of second connection node
			if (iteration_count > 1 ) { // If 2nd or higher iterations, initialize to previous iterate values
				Ptit1 = translIterator->translPower1();
				Ptit2 = translIterator->translPower2();
				PowerPrice1 = uPrice[ bufferIndex ];
				PowerPrice2 = uPrice[ ( bufferIndex + 1 ) ];
				
				if ( tnid1 == 0 ) {
					APrice1 = 0.0; // Consider node-1 as the slack node, the angle price is zero always
				}

				else {
					APrice1 = vPrice[ bufferIndex ];
				}
				
				if ( tnid2 == 0 ) {
					APrice2 = 0.0; // Consider node-1 as the slack node, the angle price is zero always
				}

				else {
					APrice2 = vPrice[ ( bufferIndex + 1 ) ];
				}
			}
			else { // If 1st iteration, initialize to zero
				Ptit1 = 0.0;
				Ptit2 = 0.0;
				PowerPrice1 = 0.0;
				PowerPrice2 = 0.0;
				APrice1 = 0.0;
				APrice2 = 0.0;
			}
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Transmission Line Optimization Iterations for Transmission line " << translIterator->getTranslID() << "\n";
				matrixResultOut << "Previous power iterate (MW/pu) for end-1\n" << Ptit1 << "\nPrevious average power (MW/pu) for end-1\n" << pavBuffer[ tnid1 ] << "\nPrevious power price (scaled LMP) for end-1\n" << PowerPrice1 << "\nAngle price from last to last iterate for end-1 (scaled)\n" << V_avg[ tnid1 ] << "\nAngle value from last iterate for end-1\n" << angleBuffer1[ tnid1 ] << "\nPrevious angle price for end-1 (scaled)\n" << APrice1 << "\nPrevious power iterate (MW/pu) for end-2\n" << Ptit2 << "\nPrevious average power (MW/pu) for end-2\n" << pavBuffer[ tnid2 ] << "\nPrevious power price (scaled LMP) for end-2\n" << PowerPrice2 << "\nAngle price from last to last iterate for end-2 (scaled)\n" << V_avg[ tnid2 ] << "\nAngle value from last iterate for end-2\n" << angleBuffer1[ tnid2 ] << "\nPrevious angle price for end-2 (scaled)\n" << APrice2 << endl;	
			}			
			translIterator->tpowerangleMessage( Rho, Ptit1, pavBuffer[ tnid1 ], PowerPrice1, V_avg[ tnid1 ], angleBuffer1[ tnid1 ], APrice1, Ptit2, pavBuffer[ tnid2 ], PowerPrice2, V_avg[ tnid2 ], angleBuffer1[ tnid2 ], APrice2 ); // Solve the Opt. Problem
			temptrans2++; 
		}
		
		
		if ( setTuning == 1 ) {
			W = ( Rho1 ) * ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with Rho1 * primalTol = dualTol
		}
		else {
			if ( setTuning == 2 ) {
				W = ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with primalTol = dualTol
			}
			else {
	 			//W = 0.0; // Definition of W for fixed Rho
				if ( iteration_count <= 3000 ) {
					W = ( Rho1 ) * ( primalTol / dualTol ) - 1; // Definition of W for adaptive Rho with Rho1 * primalTol = dualTol
				}
				else {
					W = 0.0; // Definition of W for fixed Rho
				}
			}
		}
		// Calculation of Adaptive Rho
		controllerSum = controllerSum + W;
		Rho1 = Rho; // Store previous Rho
		Rho = ( Rho1 ) * ( exp( ( lambdaAdap * W ) + ( muAdap * ( W - Wprev ) ) + ( xiAdap * controllerSum  ) ) ); // Next iterate value of Rho
		Wprev = W; // Buffering
		
		if ( Verbose ) {
			matrixResultOut << "\n*********Starting of Gather Operation************\n";
		}
		vector< Node >::iterator nodeIterator; // Distributed Optimizations; Nodes' Optimization Problem; Gather Operation
		for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
			bufferIndex = nodeIterator->getNodeID() - 1;
			//**vBuffer2[ bufferIndex ] = ( Rho1 / Rho ) * ( nodeIterator->vavMessage() ); // Gather & Calculate average v after present iteration/node
			if ( bufferIndex == 0 ) {
				angleBuffer [ bufferIndex ] = 0.0; // consider node 1 as slack node; average voltage angle always zero
			}
			else {
				angleBuffer[ bufferIndex ] = nodeIterator->ThetaavMessage(); // Calculate average angle after present iteration/node
			}
			pavBuffer[ bufferIndex ] = nodeIterator->PavMessage(); // Calculate average power after present iteration/node
			if ( Verbose ) {
				matrixResultOut << "\nNode Number: " << bufferIndex + 1 /*<< "\nV_avg = " << vBuffer2[ bufferIndex ] */<< "\nTheta_avg = " << angleBuffer[ bufferIndex ] << "\nP_avg = " << pavBuffer[ bufferIndex ] << endl;
			}
		}

		if ( Verbose ) {
			matrixResultOut << "\n*******Starting of Broadcast Operation*******\n";
		}
		// vector< Generator >::const_iterator generatorIterator;	// Broadcast to Generators
		for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
			bufferIndex = generatorIterator->getGenID() - 1;
			if ( Verbose ) {
				matrixResultOut << "\n***Generator: " << bufferIndex + 1 << " results***\n" << endl;
			}
			powerBuffer[ bufferIndex ] = generatorIterator->calcPtilde();
			uPrice[ bufferIndex ] = ( Rho1 / Rho ) * ( generatorIterator->getu() );
			angtildeBuffer[ bufferIndex ] = generatorIterator->calcThetatilde();
			//generatorIterator->calcvtilde();
			vPrice[ bufferIndex ] = ( Rho1 / Rho ) * ( generatorIterator->getv() );
			if ( Verbose ) {
				matrixResultOut << "\nPower price after this iteration ($/MWh, LMP) is: " << ( Rho / 100 ) * uPrice[ bufferIndex ] << "\nAngle price after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ] << "\nPtilde after this iteration is: " << powerBuffer[ bufferIndex ] << "\nThetatilde at the end of this iteration is: " << angtildeBuffer[ bufferIndex ] << endl;
			}
		}

		// vector< Load >::const_iterator loadIterator;	// Broadcast to Loads
		for ( loadIterator = loadObject.begin(); loadIterator != loadObject.end(); loadIterator++ ) {
			bufferIndex = genNumber + ( loadIterator->getLoadID() - 1 );
			if ( Verbose ) {
				matrixResultOut << "\n***Load: " << loadIterator->getLoadID() << " results***\n" << endl;
			}
			powerBuffer[ bufferIndex ] = loadIterator->calcPtilde();
			uPrice[ bufferIndex ] = ( Rho1 / Rho ) * ( loadIterator->getu() );
			angtildeBuffer[ bufferIndex ] = loadIterator->calcThetatilde();
			//loadIterator->calcvtilde();
			vPrice[ bufferIndex ] = ( Rho1 / Rho ) * ( loadIterator->getv() );
			if ( Verbose ) {
				matrixResultOut << "\nPower price after this iteration ($/MWh, LMP) is: " << ( Rho / 100 ) * uPrice[ bufferIndex ] << "\nAngle price after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ] << "\nPtilde after this iteration is: " << powerBuffer[ bufferIndex ] << "\nThetatilde at the end of this iteration is: " << angtildeBuffer[ bufferIndex ] << endl;
			}
		}

		int temptrans = 0; // temporary count of transmission lines to account for both the ends // Broadcast to Transmission Lines
		// vector< transmissionLine >::const_iterator translIterator;	
		for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
			bufferIndex = genNumber + loadNumber + ( translIterator->getTranslID() - 1 ) + temptrans;
			if ( Verbose ) {
				matrixResultOut << "\n***Transmission Line: " << translIterator->getTranslID() << " results***\n" << endl;
			}
			powerBuffer[ bufferIndex ] = translIterator->calcPtilde1();
			uPrice[ bufferIndex ] = ( Rho1 / Rho ) * ( translIterator->getu1() );
			angtildeBuffer[ bufferIndex ] = translIterator->calcThetatilde1();
			//translIterator->calcvtilde1();
			vPrice[ bufferIndex ] = ( Rho1 / Rho ) * ( translIterator->getv1() );
			powerBuffer[ ( bufferIndex + 1 ) ] = translIterator->calcPtilde2();
			uPrice[ ( bufferIndex + 1 ) ] = ( Rho1 / Rho ) * ( translIterator->getu2() );
			angtildeBuffer[ ( bufferIndex + 1 ) ] = translIterator->calcThetatilde2();
			//translIterator->calcvtilde2();
			vPrice[ ( bufferIndex + 1 ) ] = ( Rho1 / Rho ) * ( translIterator->getv2() );
			temptrans++;
			if ( Verbose ) {
				matrixResultOut << "\nPower price ($/MWh, LMP at end-1) after this iteration is: " << ( Rho / 100 ) * uPrice[ bufferIndex ] << "\nAngle price (end-1) after this iteration is: " << ( Rho ) * vPrice[ bufferIndex ] << "\nPtilde (end-1) after this iteration is: " << powerBuffer[ bufferIndex ] << "\nThetatilde (end-1) at the end of this iteration is: " << angtildeBuffer[ bufferIndex ] << "\nPower price ($/MWh, LMP at end-2) after this iteration is: " << ( Rho / 100 ) * uPrice[ ( bufferIndex + 1 ) ] << "\nAngle price (end-2) after this iteration is: " << ( Rho ) * vPrice[ ( bufferIndex + 1 ) ] << "\nPtilde (end-2) after this iteration is: " << powerBuffer[ ( bufferIndex + 1 ) ] << "\nThetatilde (end-2)  at the end of this iteration is: " << angtildeBuffer[ ( bufferIndex + 1 ) ] <<endl;
			}
		}

		//if ( ( iteration_count >= 100 ) && ( ( ( iteration_count % 100 ) == 0 ) || ( iteration_count == MAX_ITER - 1 ) ) ) {
			int i = 0;
			for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
				LMP[ i ] = ( Rho / 100 ) * nodeIterator->uMessage(); // record the LMP values; rescaled and converted to $/MWh
				//nodeIterator->reset(); // reset the node variables that need to start from zero in the next iteration
				++i;
			}
			//++first;
		//}
	
		for ( nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); nodeIterator++ ) {
			nodeIterator->reset(); // reset the node variables that need to start from zero in the next iteration
		}

		// Calculation of Primal Tolerance, primalTol at the end of this particular iteration
		double primsum = 0.0;
		double Primsum = 0.0;
		for ( int i = 0; i < nodeNumber; i++ ) {
			primsum = primsum + pow( pavBuffer[ i ], 2.0 );
			Primsum = Primsum + pow( pavBuffer[ i ], 2.0 );
		}
		for ( int j = 0; j < deviceTermCount; j++ )
			primsum = primsum + pow( angtildeBuffer[ j ], 2.0 );
		primalTol = sqrt( primsum );
		PrimalTol = sqrt( Primsum );
		if ( Verbose ) {
			matrixResultOut << "\nPrimal Tolerance at the end of this iteration is: " << primalTol << endl;
		}
		// Calculation of Dual Tolerance, dualTol at the end of this particular iteration
		double sum = 0.0;
		if ( iteration_count > 1 ) {
			for ( int k = 0; k < deviceTermCount; k++ ) {
				sum = sum + pow( ( powerBuffer[ k ] - powerBuffer1[ k ] ), 2.0 ); 
				//matrixResultOut << "\npowerBuffer: " << powerBuffer[ k ] << "\npowerBuffer1: " << powerBuffer1[ k ] << endl;
			}
			for ( int i = 0; i < nodeNumber; i++ ) {
				sum = sum + pow( ( angleBuffer[ i ] - angleBuffer1[ i ] ), 2.0 );
				//matrixResultOut << "\nangleBuffer: " << angleBuffer[ i ] << "\nangleBuffer1: " << angleBuffer1[ i ] << endl;
			}
		}
		else {
			for ( int i = 0; i < nodeNumber; i++ )
				sum = sum + pow( ( angleBuffer[ i ] ), 2.0 ); 
			for ( int k = 0; k < deviceTermCount; k++ )
				sum = sum + pow( ( powerBuffer[ k ] - ptildeinitBuffer[ k ] ), 2.0 );
		}
		
		dualTol = ( Rho1 ) * sqrt( sum );
		//matrixResultOut << sqrt( sum ) << endl;
		if ( Verbose ) {
			matrixResultOut << "\nDual Tolerance at the end of this iteration is: " << dualTol << endl;
			matrixResultOut << "\nObjective value at the end of this iteration is ($): " << calcObjective << endl;
			matrixResultOut << "\n****************End of " << iteration_count << " -th iteration***********\n";
		}
		objectiveValue.push_back( calcObjective ); // record the objective values

		//iteration_count++;
		//cout << iteration_count << endl;

	} // end of one iteration
	
	clock_t stop_s = clock();  // end
	matrixResultOut << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;
	matrixResultOut << "\nLast value of dual residual / Rho = " << dualTol / Rho1 << endl;
	matrixResultOut << "\nLast value of primal residual = " << primalTol << endl;
	matrixResultOut << "\nLast value of Rho = " << Rho1 << endl;
	matrixResultOut << "\nLast value of dual residual = " << dualTol << endl;
	matrixResultOut << "\nTotal Number of Iterations = " << iteration_count - 1 << endl;	
	cout << "\nExecution time (s): " << static_cast<double>( stop_s - start_s ) / CLOCKS_PER_SEC << endl;

	/**PRINT MW**/
	ofstream devProdOut( devProdString, ios::out ); // create a new file powerResult.txt to output the results
	// exit program if unable to create file
	if ( !devProdOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	devProdOut << "Gen#" << "\t" << "Conn." << "\t" << "MW" << endl;
	for ( generatorIterator = genObject.begin(); generatorIterator != genObject.end(); generatorIterator++ ) {
		devProdOut << generatorIterator->getGenID() << "\t" << generatorIterator->getGenNodeID() << "\t" <<    generatorIterator->genPower() * 100 << endl;
	}
	devProdOut << "T.line#" << "\t" << "From" << "\t" << "To" << "\t" << "From MW" << "\t" << "To MW" << endl;
	for ( translIterator = translObject.begin(); translIterator != translObject.end(); translIterator++ ) {
		devProdOut << translIterator->getTranslID() << "\t" << translIterator->getTranslNodeID1() << "\t" << translIterator->getTranslNodeID2() << "\t" << translIterator->translPower1() * 100 << "\t" << translIterator->translPower2() * 100 << endl;
	}

	/**PRINT ITERATION COUNTS**/
	ofstream iterationResultOut( iterationResultString, ios::out ); // create a new file itresult.txt to output the results	
	// exit program if unable to create file
	if ( !iterationResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	iterationResultOut << "\nIteration Count: " << endl;
	vector< int >::iterator iterationCountIterator; 
	for ( iterationCountIterator = iterationGraph.begin(); iterationCountIterator != iterationGraph.end(); iterationCountIterator++ )  		{
		iterationResultOut << *iterationCountIterator << endl;
	}

	/**PRINT LMPs**/
	ofstream lmpResultOut( lmpResultString, ios::out ); // create a new file itresult.txt to output the results	
	// exit program if unable to create file
	if ( !lmpResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	lmpResultOut << "\nLocational Marginal Prices for Real Power at nodes ($/MWh): " << endl;
	
	//for ( int j = 0; j < firstIndex; ++j ) {
		//lmpResultOut << "After " << ( j + 1 ) * 100 << " iterations, LMPs are:" << endl;
		for ( int i = 0; i < nodeNumber; ++i ) {
			lmpResultOut << i + 1 << "\t" << LMP[ i ] << endl; // print the LMP values
		}
	//}
	
	/**PRINT OBJECTIVE VALUES**/
	ofstream objectiveResultOut( objectiveResultString, ios::out ); // create a new file objective.txt to output the results
	// exit program if unable to create file
	if ( !objectiveResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	objectiveResultOut << "\nObjective value: " << endl;
	vector< double >::iterator objectiveIterator; 
	for ( objectiveIterator = objectiveValue.begin(); objectiveIterator != objectiveValue.end(); objectiveIterator++ )  {
		objectiveResultOut << *objectiveIterator << endl;
	}
	matrixResultOut << "\nLast value of Objective = " << *(objectiveIterator-1) << endl;

	/**PRINT PRIMAL RESIDUAL**/
	ofstream primalResultOut( primalResultString, ios::out ); // create a new file primresult.txt to output the results
	// exit program if unable to create file
	if ( !primalResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	primalResultOut << "\nPrimal Residual: " << endl;
	vector< double >::iterator primalToleranceIterator;
	for ( primalToleranceIterator = primTolGraph.begin(); primalToleranceIterator != primTolGraph.end(); primalToleranceIterator++ )  		{
		primalResultOut << *primalToleranceIterator << endl;
	}
	
	/**PRINT DUAL RESIDUAL**/
	ofstream dualResultOut( dualResultString, ios::out ); // create a new file dualresult.txt to output the results
	// exit program if unable to create file
	if ( !dualResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	dualResultOut << "\nDual Residual: " << endl;
	vector< double >::iterator dualToleranceIterator;
	for ( dualToleranceIterator = dualTolGraph.begin(); dualToleranceIterator != dualTolGraph.end(); dualToleranceIterator++ )  		
	{
		dualResultOut << *dualToleranceIterator << endl;
	}
} // end runSimulation

void Network::runSimulationCentral(GRBEnv* environmentGUROBI)
{	// CREATION OF THE MIP SOLVER INSTANCE //
	clock_t begin = clock(); // start the timer
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Generator>::iterator genIterator; // Iterator for Powergenerator objects
	vector<transmissionLine>::iterator tranIterator; // Iterator for Transmission line objects
	vector<Load>::iterator loadIterator; // Iterator for load objects
	vector<Node>::iterator nodeIterator; // Iterator for node objects
	string outSummaryFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/Centralized_Gurobi/Summary_of_Result_Log.txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}

        int dimRow = (2 * genNumber + 2 * translNumber + nodeNumber); // Total number of rows of the A matrix (number of structural constraints of the QP) first term to account for lower and upper generating limits, second term for lower and upper line limits for transmission lines, the third term to account for nodal power balance constraints
        int dimCol = (genNumber+nodeNumber); // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for nodes
	outPutFile << "\nTotal Number of Structural Constraints (Rows) is: " << dimRow << endl;
	outPutFile << "\nTotal Number of Decision Variables (Columns) is: " << dimCol << endl;
	// Instantiate GUROBI Problem model
	GRBModel *modelCentQP = new GRBModel(*environmentGUROBI);
	cout << "\nGurobi model created" << endl;
    	modelCentQP->set(GRB_StringAttr_ModelName, "assignment");
	cout << "\nGurobi model created and name set" << endl;
	GRBVar decvar[dimCol+1];
	cout << "\nGurobi decision variables created" << endl;
	double z; // variable to store the objective value

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	cout << "\nGurobi decision variables to be assigned" << endl;
	decvar[0] = modelCentQP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	cout << "\nGurobi dummy decision variable created" << endl;
	int colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		decvar[colCount] = modelCentQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Power Generation continuous variables for different generators: " << colCount << endl;

	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//	
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		decvar[colCount] = modelCentQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for different intrazonal nodes: " << colCount << endl;
	outPutFile << "\nTotal Number of columns for generation, angles: " << colCount-1 << endl;
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount-1 << endl;
	//Setting Objective//
	GRBQuadExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		obj += (genIterator->getQuadCost())*(decvar[colCount])*(decvar[colCount])+(genIterator->getLinCost())*(decvar[colCount])+(genIterator->getNLCost());
		++colCount;
	}
	//Columns corresponding to Voltage Phase Angles continuous variables for different intrazonal nodes//	
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		obj += 0*(decvar[colCount]);	
		++colCount;
	}

	modelCentQP->setObjective(obj, GRB_MINIMIZE);
	//Row Definitions: Specification of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	// Constraints corresponding to supply-demand balance
	string outPGenFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/Centralized_Gurobi/PgenFile.txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	// Coefficients for the supply-demand balance constraints
	outPutFile << "\nNon-zero elements of A matrix" << endl;
	outPutFile << "\nRow Number\tColumn Number\tNon-zero Entry\tFrom Reactance\tToReactance" << endl;
	outPutFile << "\nCoefficients for the supply-demand balance constraints" << endl;
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelCentQP->addConstr(lhs[0], GRB_EQUAL, 0);
	int rCount = 1; // Initialize the row count
	vector<int> busCount; // vector for storing the node/bus serial
	outPutFile << "Constraints corresponding to Supply-Demand Balance right hand side" << endl;
	for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
		outPutFile << "\nGeneration\t" << rCount << "\n";
		int genListLength = (nodeIterator)->getGenLength(); // get the number
		lhs[rCount]=0;
		for (int cCount = 1; cCount <= genListLength; ++cCount){
			lhs[rCount] += 1*(decvar[(nodeIterator)->getGenSer(cCount)]);
			outPutFile << "\n" << rCount << "\t" << (nodeIterator)->getGenSer(cCount) << "\t" << 1.0 << endl;
		}
		outPutFile << "\nIntrazonal Node Angles\t" << rCount << "\n";
		lhs[rCount] += (((nodeIterator)->getToReact())-((nodeIterator)->getFromReact()))*(decvar[genNumber+rCount]);
		outPutFile << "\n" << rCount << "\t" << genNumber+rCount << "\t" << -((nodeIterator)->getToReact())-((nodeIterator)->getFromReact()) << "\t" << -((nodeIterator)->getFromReact()) << "\t" << -((nodeIterator)->getToReact()) << endl;
		outPutFile << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
		int connNodeListLength = (nodeIterator)->getConNodeLength(); // get the number of intra-zonal nodes connected to this node
		for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
			if (((nodeIterator)->getConnReact(cCount))<=0)
				lhs[rCount] -= (((nodeIterator)->getConnReact(cCount)))*(decvar[genNumber+((nodeIterator)->getConnSer(cCount))]);
			else
				lhs[rCount] += (((nodeIterator)->getConnReact(cCount)))*(decvar[genNumber+((nodeIterator)->getConnSer(cCount))]);
			outPutFile << "\n" << rCount << "\t" << genNumber+((nodeIterator)->getConnSer(cCount)) << "\t" <<  (-((nodeIterator)->getConnReact(cCount))) << "\n";

		}
		busCount.push_back(rCount);
		if (((nodeIterator)->getLoadVal())==0) {
			modelCentQP->addConstr(lhs[rCount], GRB_EQUAL, ((nodeIterator)->getLoadVal()));
		}
		else {
			modelCentQP->addConstr(lhs[rCount], GRB_EQUAL, -((nodeIterator)->getLoadVal()));
		}
		outPutFile << "Connected load to node " << rCount << " is " << (nodeIterator)->getLoadVal()*100 << " MW" << endl;
		outPutFile << rCount << "\t";
		if (((nodeIterator)->getLoadVal())==0)
			outPutFile << ((nodeIterator)->getLoadVal())*100 << " MW" << endl;
		else
			outPutFile << -((nodeIterator)->getLoadVal())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next node object
	}
	// Coefficients corresponding to lower generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - nodeNumber];
		modelCentQP->addConstr(lhs[rCount] >= ((genIterator)->getPMin()));
		outPutFile << rCount << "\t" << (rCount - nodeNumber) << "\t" << 1.0 << "\t" << (genIterator)->getPMin() << endl;
		outPutFile << rCount << "\t";
		outPutFile << ((genIterator)->getPMin())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next generator object
	}
	// Coefficients corresponding to upper generation limits
	outPutFile << "\nCoefficients corresponding to upper generation limits\n";
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		lhs[rCount] = 0;
		lhs[rCount] += decvar[rCount - (genNumber + nodeNumber)];
		modelCentQP->addConstr(lhs[rCount] <= ((genIterator)->getPMax()));
		outPutFile << rCount << "\t" << (rCount - (genNumber + nodeNumber)) << "\t" << 1.0 << "\t" << ((genIterator)->getPMax()) << endl;
		outPutFile << rCount << "\t";
		outPutFile << ((genIterator)->getPMax())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next generator object
	}
	// Coefficients corresponding to intra-zone Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Forward Flow Limit Constraints\n";
	for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
		lhs[rCount] = 0;
		lhs[rCount] += (1/((tranIterator)->getReactance()))*(decvar[genNumber + (tranIterator)->getTranslNodeID1()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (tranIterator)->getTranslNodeID1() << "\t" << 1/((tranIterator)->getReactance()) << "\t" << 1/((tranIterator)->getReactance()) << "\n";
		lhs[rCount] += (-1/((tranIterator)->getReactance()))*(decvar[genNumber + (tranIterator)->getTranslNodeID2()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (tranIterator)->getTranslNodeID2() << "\t" << -1/((tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((tranIterator)->getReactance()) << "\n";
		modelCentQP->addConstr(lhs[rCount] <= ((tranIterator)->getFlowLimit()));
		outPutFile << rCount << "\t";
		outPutFile << ((tranIterator)->getFlowLimit())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next transmission line object		
	}
	// Coefficients corresponding to intra-zone Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Reverse Flow Limit Constraints\n";
	for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
		lhs[rCount] = 0;
		lhs[rCount] += (1/((tranIterator)->getReactance()))*(decvar[genNumber + (tranIterator)->getTranslNodeID1()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (tranIterator)->getTranslNodeID1() << "\t" << 1/((tranIterator)->getReactance()) << "\t" << 1/((tranIterator)->getReactance()) << "\n";
		lhs[rCount] += (-1/((tranIterator)->getReactance()))*(decvar[genNumber + (tranIterator)->getTranslNodeID2()]);
		outPutFile << "\n" << rCount << "\t" << genNumber + (tranIterator)->getTranslNodeID2() << "\t" << -1/((tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((tranIterator)->getReactance()) << "\n";
		modelCentQP->addConstr(lhs[rCount] >= -((tranIterator)->getFlowLimit()));
		outPutFile << rCount << "\t";
		outPutFile << -((tranIterator)->getFlowLimit())*100 << " MW" << endl;
		++rCount; // Increment the row count to point to the next transmission line object
	}	
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rCount - 1 << endl;
	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	modelCentQP->optimize(); // Solves the optimization problem
	int stat = modelCentQP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelCentQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_INF_OR_UNBD) {
		outPutFile << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelCentQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_UNBOUNDED) {
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelCentQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_OPTIMAL) {
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;
		cout << "\nThe solution to the problem is OPTIMAL." << endl;

		//Get the Optimal Objective Value results//
		z = modelCentQP->get(GRB_DoubleAttr_ObjVal);

		// Open separate output files for writing results of different variables
		string outIntAngFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/Centralized_Gurobi/AngleResult.txt";
		string outTranFlowFileName = "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/Centralized_Gurobi/TranFlow.txt";
		ofstream internalAngleOut(outIntAngFileName, ios::out); //switchStateOut
		ofstream tranFlowOut(outTranFlowFileName, ios::out);
		outPutFile << "\nThe Optimal Objective value (Generation Dispatch cost) is: " << z << endl;
		powerGenOut << "\nThe Optimal Objective value (Generation Dispatch cost) is: " << z << endl;
		cout << "\nThe Optimal Objective value (Generation Dispatch cost) is: " << z << endl;
		vector<double> x; // Vector for storing decision variable output 
		x.push_back(0); // Initialize the decision Variable vector

		//Display Power Generation
		powerGenOut << "\n****************** GENERATORS' POWER GENERATION LEVELS (MW) *********************" << endl;
		powerGenOut << "GENERATOR ID" << "\t" << "GENERATOR MW" << "\n";
		int arrayInd = 1;
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			powerGenOut << (genIterator)->getGenID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
			++arrayInd;
		}
		powerGenOut << "Finished writing Power Generation" << endl;

		// Display Internal node voltage phase angle variables
		internalAngleOut << "\n****************** INTERNAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
		internalAngleOut << "NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
			internalAngleOut << (nodeIterator)->getNodeID() << "\t" << ((decvar[arrayInd]).get(GRB_DoubleAttr_X)) << endl;		
			++arrayInd;			
		}
		internalAngleOut << "Finished writing Internal Node Voltage Phase Angles" << endl;
		// Display Internal Transmission lines' Flows
		tranFlowOut << "\n****************** INTERNAL TRANSMISSION LINES FLOWS *********************" << endl;
		tranFlowOut << "TRANSMISSION LINE ID" << "\t" << "MW FLOW" << "\n";
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			tranFlowOut << (tranIterator)->getTranslID() << "\t" << (1/((tranIterator)->getReactance()))*((decvar[genNumber +(tranIterator)->getTranslNodeID1()]).get(GRB_DoubleAttr_X)-(decvar[genNumber + (tranIterator)->getTranslNodeID2()]).get(GRB_DoubleAttr_X))*100 << " MW" << endl;
		}
		tranFlowOut << "Finished writing Internal Transmission lines' MW Flows" << endl;
		delete modelCentQP; // Free the memory of the GUROBI Problem Model
		clock_t end2 = clock(); // stop the timer
		double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
		outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
		cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
		internalAngleOut.close();
		tranFlowOut.close();
	}
	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl;
}
