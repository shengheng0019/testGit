#pragma once
#include "Constraints.h"
#include"Vars.h"
#include"Solution.h"
#include<iostream>
#include<string>
#include "scip/scipdefplugins.h"
#include "scip/type_var.h"
#include "scip/type_scip.h"
#include "scip/scip.h"
#include "Eigen/Dense"
#include "Objective.h"
#include "ColumnPool.h"
#include "scip/struct_scip.h"
#include "heuristicTools.h"
#include "spdlog/spdlog.h"

//using namespace Eigen;

//void BuildSCIPModel(SCIP* scipModel, SCIP_VAR** scipVars, SCIP_CONS** conss, Constraints* constraints, int numConss, Objective obj);

void GetDualSolution(SCIP* scipModel, SCIP_CONS** scipConss, int scipConsNum, Solution& sol); 

Eigen::VectorXd Sol2Vec(Solution& sol);
Eigen::VectorXd Vector2Vec(const vector<double>& sol);

class Lagrangean;

class SubProb
{
public:
	Objective oriObj;                   //@wqy:using by diving
	Constraints* conss;
	int numConss;						//num of conss  //@wqy:using by diving
	unordered_map<int, Vars> varDic;	//key: varid, value: Vars
	unordered_map<int, double> tempObjVal;	//to store the obj of subProb, key:varInd, value:objcoef
	int varSize;						 //@wqy:using by diving
	SubProb(Constraints* constraints, int numConss, const Objective& obj, Constraints* couplingCons, int couplingConsNum, Vars* vars, int varNum);
	SubProb(vector<Constraints>& constraints, const Objective& obj, vector<Constraints>& couplingCons, const vector<Vars>& vars);
	SubProb(vector<Constraints>& constraints, const Objective& obj, vector<Constraints>& couplingCons, const vector<Vars>& vars, bool flag);//@wqy diving using this constructor
	SubProb() = default;
	~SubProb();
	SubProb& operator = (const SubProb& sp);
	void BuildSCIPModel();				//build the scip subProblem
	void BuildSCIPModelNew(FIX Fixed_Info, int subprobNo);                         //@wqy:using by diving
	void BuildSCIPModelNew(FIX fixedInfo, int subprobNo, const std::vector<double>& duals_part);//@wqy:using by diving
	void UpdateTempObj(Solution dualSol);	//
	void UpdateTempObj(double* dualSol, int solSize);
	void UpdateTempObj(const Eigen::VectorXd& dualSol);
	void UpdateTempObjNew(const vector<double>& dualSol);                           //@wqy:using by diving
	void updateSCIPObj(vector<double>objCeff);                                      //@wqy:using by diving
	bool GetSolution(); 
	bool GetSolutionNew(Solution& columnGen);                                      //@wqy:using by diving
	Solution tempSol;					//to store the current sol of the subProb

	ColumnPool colPool;
	vector<int> subOrder2Id;

	friend class Lagrangean;
private:
	SCIP* scipModel;                                                              //@wqy:using by diving
	SCIP_VAR** scipVars;                                                          //@wqy:using by diving
	SCIP_CONS** scipConss;                                                        //@wqy:using by diving
	Eigen::MatrixXd A;				//coupling constraints Matrix                 //@wqy:using by diving
	Eigen::VectorXd c;				//original objective coeficient               //@wqy:using by diving
	void updateSCIPObj();			//update the obj coef of scip model without reconstructing the scip model

	//@wqy add0621, for a clearer order in subproblem variable
	std::vector<SCIP_VAR*> scipVarsVec_;  //vector of scipvars for the troublesome order, instead of "scipVars"     //@wqy:using by diving
	std::vector<SCIP_CONS*>scipConssVec_; //vector of scipcons for the troublesome order, instead of "scipConss"	//@wqy:using by diving
	std::vector<Vars> varsVec_;           //vector of vars, instead of "varDic".									//@wqy:using by diving
	std::vector<Constraints> consVec_;    //vector of constraints, instead of "conss".								//@wqy:using by diving
};

