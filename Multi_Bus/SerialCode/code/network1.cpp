// Member functions for class Network
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
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
//#define MAX_ITER 80002
#define LINE_CAP 1.00
using namespace std;

Network::Network( int val )
	: networkID( val ), // constructor begins; initialize networkID  and Rho through constructor initializer list
	  Rho( 1.0 )
{
	setNetworkVariables( networkID ); // sets the variables of the network

} // end constructor

// destructor
Network::~Network()
{	
	cout << "\nNetwork instance: " << networkID << " for this simulation destroyed. You can now open the output files to view the results of the simulation\n" << endl;

} // end destructor

void Network::setNetworkVariables( int networkID ) // Function setNetworkVariables startsto initialize the parametersand variables
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

			case 48: // 48 Bus case

				strcpy( genFile, "Gen48.txt" );			
				strcpy( tranFile, "Tran48.txt" );
				strcpy( loadFile, "Load48.txt" );	

			case 5: // 48 Bus case

				strcpy( genFile, "Gen5.txt" );			
				strcpy( tranFile, "Tran5.txt" );
				strcpy( loadFile, "Load5.txt" );	
				break; // exit switch

			case 4: // 3 Bus case variant

				strcpy( genFile, "Gen3A.txt" );			
				strcpy( tranFile, "Tran3A.txt" );
				strcpy( loadFile, "Load3A.txt" );	
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
		
		// Create Generators
		for ( int i = 0; i < genNumber; ++i ) {
			int j = 0; 
			int gNodeID; // node object ID to which the particular generator object is connected
			do {
				gNodeID = matrixGen[ i ][ j ];
			} while ( ( gNodeID <= 0 ) || ( gNodeID > nodeNumber ) ); // validity check

			double c2, c1, c0, PgMax, PgMin; // Parameters for Generator
			do {
				//Quadratic Coefficient: 
				c2 = matrixGen[ i ][ ( j + 1 ) ];
				//Linear coefficient: 
				c1 = matrixGen[ i ][ ( j + 2 ) ];
				//Constant term: 
				c0 = matrixGen[ i ][ ( j + 3 ) ];
				//Maximum Limit: 
				PgMax = matrixGen[ i ][ ( j + 4 ) ];
				//Minimum Limit: 
				PgMin = matrixGen[ i ][ ( j + 5 ) ];
			} while ( (c2 < 0 ) || ( c1 < 0 ) || ( PgMax <= 0 ) || ( PgMin < 0 ) || ( PgMax <= PgMin ) ); 
			// check the bounds and validity of the parameter values

			Gensolver genParam( c2, c1, c0, PgMax, PgMin ); // Instantiate the copy constructor for the generator solver object
	
			Generator generatorInstance( i + 1, &nodeObject[ gNodeID - 1 ], genParam ); // creates generatorInstance object with ID number i + 1

			genObject.push_back( generatorInstance ); // pushes the generatorInstance object into the vector

		} // end initialization for Generators

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
			do {
				//node IDs of the node objects to which this transmission line is connected.
				tNodeID1 = matrixTran[ k ][ l ]; //From end
				tNodeID2 = matrixTran[ k ][ ( l + 1 ) ]; //To end
			} while ( ( tNodeID1 <= 0 ) || ( tNodeID1 > nodeNumber ) || ( tNodeID2 <= 0 ) || ( tNodeID2 > nodeNumber ) || ( tNodeID1 == tNodeID2) ); // validity check
			double resT, reacT, ptMax, ptMin; // Parameters for Transmission Line
			do {
				//Resistance:
				resT = matrixTran[ k ][ ( l + 2 ) ];
				//Reactance:
				reacT = matrixTran[ k ][ ( l + 3 ) ];
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
			int lNodeID; // node object ID to which the particular load object is connected
			do {
				//node ID of the node object to which this load object is connected.
				lNodeID = matrixLoad[ j ][ k ]; 
			} while ( ( lNodeID <= 0 ) || ( lNodeID > nodeNumber ) ); // validity check

			double P_Load; // Parameters for Load
			do {
				//value of allowable power consumption capability of load with a negative sign to indicate consumption:
				//Power Consumption:
				P_Load = matrixLoad[ j ][ ( k + 1 ) ];
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
void Network::runSimulation() //Function runSimulation begins
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
	vector< double > objectiveValue; // vector of objective function value

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
	double lambdaAdap = 0.0001; // Parameters of the PID controller for adjusting the ADMM tuning parameter
	double muAdap = 0.0005; // Parameters of the PID controller for adjusting the ADMM tuning parameter
	int setTuning; // parameter to select adaptive rho, fixed rho, and type of adaptive rho

	// Set the type of tuning
	cout << "Enter the tuning mode; Enter 1 for maintaining Rho * primTol = dualTol; 2 for primTol = dualTol; anything else for fixed Rho\n" << endl;
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

	ofstream matrixResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/Summary_of_Result_Log.txt", ios::out ); // create a new file result.txt to output the results
	ofstream tranResultOut( "Tie_Flow.txt", ios::out); // create a new file to output the results of tie line flow
	
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
	//for ( iteration_count = 1; ( iteration_count < 101 ); iteration_count++ ) {
	while( ( primalTol >= 0.006 ) || ( dualTol >= 0.006 ) ) { // ( iteration_count <= 122 )
	
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
			generatorIterator->gpowerangleMessage( Rho, Pgit, pavBuffer[ gnid ], PowerPrice, V_avg[ gnid ], angleBuffer1[ gnid ], APrice ); // Solve the Optimization Problem
			calcObjective = calcObjective + generatorIterator->objectiveGen(); // calculate the total objective after this iteration
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
		double tieFlow = 0.0;	
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
			if ( ( ( translNumber - temptrans2 ) <= 3 ) )
				tieFlow = tieFlow + Ptit2;
			if ( Verbose ) {
				matrixResultOut << "\nStarting of Transmission Line Optimization Iterations for Transmission line " << translIterator->getTranslID() << "\n";
				matrixResultOut << "Previous power iterate (MW/pu) for end-1\n" << Ptit1 << "\nPrevious average power (MW/pu) for end-1\n" << pavBuffer[ tnid1 ] << "\nPrevious power price (scaled LMP) for end-1\n" << PowerPrice1 << "\nAngle price from last to last iterate for end-1 (scaled)\n" << V_avg[ tnid1 ] << "\nAngle value from last iterate for end-1\n" << angleBuffer1[ tnid1 ] << "\nPrevious angle price for end-1 (scaled)\n" << APrice1 << "\nPrevious power iterate (MW/pu) for end-2\n" << Ptit2 << "\nPrevious average power (MW/pu) for end-2\n" << pavBuffer[ tnid2 ] << "\nPrevious power price (scaled LMP) for end-2\n" << PowerPrice2 << "\nAngle price from last to last iterate for end-2 (scaled)\n" << V_avg[ tnid2 ] << "\nAngle value from last iterate for end-2\n" << angleBuffer1[ tnid2 ] << "\nPrevious angle price for end-2 (scaled)\n" << APrice2 << endl;	
			}			
			translIterator->tpowerangleMessage( Rho, Ptit1, pavBuffer[ tnid1 ], PowerPrice1, V_avg[ tnid1 ], angleBuffer1[ tnid1 ], APrice1, Ptit2, pavBuffer[ tnid2 ], PowerPrice2, V_avg[ tnid2 ], angleBuffer1[ tnid2 ], APrice2 ); // Solve the Opt. Problem
			temptrans2++; 
		}
		tranResultOut << abs( ( 116.4 - (tieFlow * 100 ) ) / 116.4 ) << endl;
		
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
		Rho1 = Rho; // Store previous Rho
		Rho = ( Rho1 ) * ( exp( ( lambdaAdap * W ) + ( muAdap * ( W - Wprev ) ) ) ); // Next iterate value of Rho
		Wprev = W; // Buffering
		//cout << "\nThe value of Rho and W after iteration " << iteration_count << " are " << Rho << " and " << W << endl;
		//if ( ( iteration_count >= 2900 ) && ( iteration_count <= 2910 ) ) {
			//cout << "\nThe values of Primal and Dual Tolerances are " << primalTol << " and " << dualTol << endl;
		//}

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

		iteration_count++;

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
	ofstream devProdOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/powerResult.txt", ios::out ); // create a new file powerResult.txt to output the results
	
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
	ofstream iterationResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/itresult.txt", ios::out ); // create a new file itresult.txt to output the results
	
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
	ofstream lmpResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/LMPresult.txt", ios::out ); // create a new file itresult.txt to output the results
	
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
	ofstream objectiveResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/objective.txt", ios::out ); // create a new file objective.txt to output the results
	
	// exit program if unable to create file
	if ( !objectiveResultOut ) {
		cerr << "File could not be opened" << endl;
		exit( 1 );
	}
	
	objectiveResultOut << "\nObjective value: " << endl;
	vector< double >::iterator objectiveIterator; 
	for ( objectiveIterator = objectiveValue.begin(); objectiveIterator != objectiveValue.end(); objectiveIterator++ )  {
		//objectiveResultOut << abs( ( 113000 - *objectiveIterator ) / 113000 ) << endl;
		objectiveResultOut << *objectiveIterator << endl;
	}
	matrixResultOut << "\nLast value of Objective = " << *(objectiveIterator-1) << endl;

	/**PRINT PRIMAL RESIDUAL**/
	ofstream primalResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/primresult.txt", ios::out ); // create a new file primresult.txt to output the results
	ofstream PrimalResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/loadBalanceresult.txt", ios::out );	

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
	
	PrimalResultOut << "\nLoad Balance: " << endl;
	vector< double >::iterator PrimalToleranceIterator;
	for ( PrimalToleranceIterator = PrimTolGraph.begin(); PrimalToleranceIterator != PrimTolGraph.end(); PrimalToleranceIterator++ )  		{
		PrimalResultOut << *PrimalToleranceIterator << endl;
	}
	/**PRINT DUAL RESIDUAL**/
	ofstream dualResultOut( "/home/samie/code/ADMM_Based_Proximal_Message_Passing_Distributed_OPF/OPF/Multi_Bus/SerialCode/output/ADMM_PMP_CVXGEN/dualresult.txt", ios::out ); // create a new file dualresult.txt to output the results
	
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

	




	

	
