// Member functions for class Generator.
#include <iostream>
#include <iomanip>
// include Generator class definition from generator.h
#include "generator.h"
// Include definition of Node class 
#include "node.h"
// Include Generator solver class defintion
#include "gensolver.h"
#include "gurobi_c++.h" // includes definition of the GUROBI solver header file

using namespace std;

Generator::Generator( int idOfGen, Node *nodeConng, Gensolver &paramOfGen ) // constructor begins
	: genID( idOfGen ),
	  connNodegPtr( nodeConng ),
	  genSolver( paramOfGen )
{
	//cout << "\nInitializing the parameters of the generator with ID: " << genID << endl;
	connNodegPtr->setgConn(genID); // increments the generation connection variable to node
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
	
} // end of setGenData function

void Generator::gpowerangleMessage( double gRho, double Pprevit, double Pnetavg, double uprev, double vprevavg, double Aprevavg, double vprev ) //const // function gpowerangleMessage begins
{
	genSolver.mainsolve( gRho, Pprevit, Pnetavg, uprev, vprevavg, Aprevavg, vprev ); // calls the Generator optimization solver
	Pg = genSolver.getPSol(); // get the Generator Power iterate
	Thetag = genSolver.getThetaSol(); // get the Generator voltage angle iterate
	connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle and v
} // function gpowerangleMessage ends

void Generator::gpowerangleMessageGUROBI( double gRho, double Pprevit, double Pnetavg, double uprev, double vprevavg, double Aprevavg, double vprev, GRBEnv* environmentGUROBI ) //const // function gpowerangleMessage begins
{
	// CREATION OF THE MIP SOLVER INSTANCE //
        int dimRow = 2; // Total number of rows of the A matrix (number of structural constraints of the QP): first term for the upper generation limit, the next term for the lower generation limit
        int dimCol = 2; // Total number of columns of the QP (number of Decision Variables) first term to account for power generation MW outputs, second term for voltage phase angles for generation node
	// Instantiate GUROBI Problem model
	GRBModel *modelGenQP = new GRBModel(*environmentGUROBI);
    	modelGenQP->set(GRB_StringAttr_ModelName, "assignment");
	modelGenQP->set(GRB_IntParam_OutputFlag, 0);
	GRBVar decvar[dimCol+1];
	double z; // variable to store the objective value

	// SPECIFICATION OF PROBLEM PARAMETERS //
	// Dummy Decision Variable //
	decvar[0] = modelGenQP->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
	//Decision Variable Definitions, Bounds, and Objective Function Co-efficients//
	int colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	decvar[colCount] = modelGenQP->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	++colCount;
	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
	decvar[colCount] = modelGenQP->addVar((0), (44/7), 0.0, GRB_CONTINUOUS);
	//Setting Objective//
	GRBQuadExpr obj = 0.0;
	// Objective Contribution from Dummy Decision Variable //
	obj += 0*(decvar[0]);
	colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	obj += (genSolver.getQuadCost())*(decvar[colCount])*(decvar[colCount])+(genSolver.getLinCost())*(decvar[colCount])+(genSolver.getNLCost())+(gRho/2)*(decvar[colCount]-Pprevit+Pnetavg+uprev)*(decvar[colCount]-Pprevit+Pnetavg+uprev);
	++colCount;
	//Columns corresponding to Voltage Phase Angles continuous variables for different nodes//
	obj += (gRho/2)*(decvar[colCount]-vprevavg-Aprevavg+vprev )*(decvar[colCount]-vprevavg-Aprevavg+vprev );

	modelGenQP->setObjective(obj, GRB_MINIMIZE);
	//Row Definitions: Specification of b<=Ax<=b//
	GRBLinExpr lhs[dimRow+1];
	//Row Definitions and Bounds Corresponding to Constraints/
	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	// Dummy Constraint //
	lhs[0] = 0*(decvar[0]);
	modelGenQP->addConstr(lhs[0], GRB_EQUAL, 0);
	int rCount = 1; // Initialize the row count
	// Coefficients corresponding to lower generation limits
	lhs[rCount] = 0;
	lhs[rCount] += decvar[rCount];
	modelGenQP->addConstr(lhs[rCount] >= (genSolver.getPMin()));
	++rCount; // Increment the row count to point to the next generator object
	// Coefficients corresponding to upper generation limits
	lhs[rCount] = 0;
	lhs[rCount] += decvar[rCount - 1];
	modelGenQP->addConstr(lhs[rCount] <= (genSolver.getPMax()));
	// RUN THE OPTIMIZATION SIMULATION ALGORITHM //
	modelGenQP->optimize(); // Solves the optimization problem
	int stat = modelGenQP->get(GRB_IntAttr_Status); // Outputs the solution status of the problem 

	// DISPLAY THE SOLUTION DETAILS //
	if (stat == GRB_INFEASIBLE){
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_INF_OR_UNBD) {
		cout << "\nNO FEASIBLE or BOUNDED solution to the problem exists." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_UNBOUNDED) {
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		delete modelGenQP; // Free the memory of the GUROBI Problem Model
	} else if (stat == GRB_OPTIMAL) {
		//Get the Optimal Objective Value results//
		z = modelGenQP->get(GRB_DoubleAttr_ObjVal);
		// writing results of different variables
		vector<double> x; // Vector for storing decision variable output 
		x.push_back(0); // Initialize the decision Variable vector
		objOpt = 0;
		//Power Generation
		int arrayInd = 1;
		x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
		Pg = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator Power iterate
		objOpt += (genSolver.getQuadCost())*(Pg)*(Pg)+(genSolver.getLinCost())*(Pg)+(genSolver.getNLCost());
		++arrayInd;
		// Internal node voltage phase angle variables
		x.push_back((decvar[arrayInd]).get(GRB_DoubleAttr_X));
		Thetag = ((decvar[arrayInd]).get(GRB_DoubleAttr_X)); // get the Generator voltage angle iterate	
		connNodegPtr->powerangleMessage( Pg, v, Thetag ); // passes to node object the corresponding iterates of power, angle and v
		delete modelGenQP; 
	}
} // function gpowerangleMessage ends

double Generator::genPower() //const // function genPower begins
{
	return Pg; // returns the Pg iterate
} // function genPower ends

double Generator::objectiveGen() // function objectiveGen begins
{
	return genSolver.getObj(); //returns the evaluated objective for CVXGEN
} // function objectiveGen ends

double Generator::objectiveGenGUROBI() // function objectiveGen begins
{
	return objOpt; //returns the evaluated objective for GUROBI 
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
double Generator::getPMax(){return genSolver.getPMax();}
double Generator::getPMin(){return genSolver.getPMin();}
double Generator::getQuadCost(){return genSolver.getQuadCost();}
double Generator::getLinCost(){return genSolver.getLinCost();}
double Generator::getNLCost(){return genSolver.getNLCost();}
