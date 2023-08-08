#pragma once
#include "scip/scipdefplugins.h"
#include "scip/scip.h"
#include "scip/struct_scip.h"
#include"Constraints.h"
#include"Objective.h"
#include"ColumnPool.h"
#include"heuristicTools.h"
#include"Eigen/Dense"
#include "spdlog/spdlog.h"



struct RMPSOL
{
	//here maybe a tuple
	vector<int> subprob_no;
	vector<double> varValue;
	vector<std::string> varName;
};

class MasterProb
{
public:
	MasterProb() = default;
	~MasterProb() = default;
	MasterProb& operator += (const MasterProb& mp);
	void BuildSCIPModel(Objective obj, vector<Constraints> couplingCons, int couplingConsNum, int subprobNum, FIX fixedInfo);
	void BuildSCIPModelBinary(Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX> Fixed_Info);
	int SolveDual(vector<double>& duals);
	int SolveDualCplex(vector<double>& duals);
	int SolveObj(double& obj);
	bool AddColumn(const Solution& sol, int subprob_no);
	void Solve(RMPSOL& sol, double& objval, bool root_flag);
	void free();

	vector<ColumnPool> columnpool;
	int col_num;
	int sp_num;
	int CheckSet1(int varGlobalIndex) const; //@wqy0626 no use
	int CheckSet1(string varGlobalName) const; //@wqy0626 in use
	int CheckSet1(std::unordered_map<int, double> varValueFixed, std::unordered_map<int, bool> varInt) const; //@wqy0626 no use
	

private:
	SCIP* RMPscipModel;
	SCIP_VAR** RMPscipVars;
	int RMPscipVars_num;
	SCIP_CONS** RMPscipConss;
	int RMPscipConss_num;
	Eigen::MatrixXd A;				
	Eigen::MatrixXd ConsCoefMatrix; 
};



