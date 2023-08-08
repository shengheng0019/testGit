#pragma once

#include <iostream>
#include <scip/scipdefplugins.h>
#include <vector>
#include <string>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"


namespace FPPARA
{
	const double alpha = 1.0;
	const double alpha_quot = 0.9;
	const int maxIter1 = 10000;
	const int maxit1wi = 70;
	const double epInt = 1E-5;
	const double alphadist = 0.005;
	const int trsrnd = 4;
	const int minChangeNum = 100;
	const double trsmin = 0.001;
}


class FP
{
public:

	std::vector<int> intVars;		//index of all integer variables in "std::vector<SCIP_VAR*> vars;"
	std::vector<bool> isBin;			//to indicate whether a variable in intVars is binary

	std::vector<double> relaxedIntVars;		//fractional value
	std::vector<double> roundedIntVars;		//integer value
	std::vector<double> prevRoundedIntVars;		//integer value
	std::vector<double> bestpoint;			//best integer value

	
	int RUNFP();
	void SetIntegerValue(std::unordered_map<std::string, double> patialSol);
	std::map<std::string, double> OutSol();		//output the solution map

    double GetTempObjVal();					//get the original objective of the model
	FP(std::string fileName);
	~FP();

private:
	void setLpObj();						//update the objective of the LP model
	void Round();
	void Solve(int& found, int stage);		
	double getObjVal();
	void Restart();
	

	SCIP* scip;
	double consTerm = 0;						//the constant term in the model
	int changed = 0;							//the number of rounded variables whose value are changed compared to the last iteration
	int stage = 0;								//value: 1/2/3
	int restarts = 0;							//the number of times the Restart() is called

	void Flip();							//flit the value of the variables randomly
	bool isCoincident();					//if the relaxed binary variables have the same values as the rounded ones
	bool isCoincidentS2();
	void RecoverIntVars();					//recover all the integer variables from the lp relaxation to get a MIP model
	void SetNewObjS2();
	void IntroduceAux(int intVarIdx);					//dynamically introduce auxilary variables for the model


	bool CheckFeas();						//check if the solution is feasible

	std::vector<double> intVarsLB;				//to record the lb of all variables
	std::vector<double> intVarsUB;				//to record the ub of all variables
	std::vector<SCIP_VAR*> vars;			//all original variables
	std::vector<SCIP_CONS*> conss;			//all original constraints
	std::vector<SCIP_VAR*> auxVars;			//auxilary variables
	std::vector<SCIP_CONS*> auxLBConss;		//auxilary constraints
	std::vector<SCIP_CONS*> auxUBConss;		//auxilary constraints

	std::vector<double> oriObj;				//original model objective coefficients
	double offSet = 0;

	//parameters
	double alpha;							
	int maxIter1;
	double epInt;
	int maxit1wi;
	double alpha_quot;
	int trsrnd;
	int minChangeNum;
	double trsmin;
	

};


