#pragma once
#include "ColumnPool.h"
#include "MasterProb.h"
#include "SubProbManagement.h"
#include <map>
#include "GeneralBlock.h"




class Node
{
public:
	Node() = default;
	~Node() = default;
	Node(Objective obj_master, vector<Constraints> couplingCons, int couplingCons_num, vector<Constraints> constraints,
		int constraints_num, vector<ColumnPool> columnpool, int subprob_num, std::unordered_map<int, bool>varIntTemp, int nodeNo)
		:obj_(obj_master), couplingCons_(couplingCons), couplingConsNum_(couplingCons_num), constraints_(constraints),
		constraintsNum_(constraints_num), columnpool_(columnpool), subprobNum_(subprob_num), varInt_(varIntTemp), index_(nodeNo) {
		lpRelaxObj_ = -1;
		fixedSpNo_ = -1;
		masterProblem_.columnpool = columnpool;
		masterProblem_.sp_num = subprob_num;
		for (int i = 0; i < subprob_num; i++) {
			fixedInfo_.spIsFixed.push_back(false);
		}
	};
	//Node(const Node& node); //@wqy no used
	//Node& operator = (const Node& node); //@wqy no used

	//inline function
	RMPSOL LpRelaxSol() const{ return lpRelaxSol_; }
	FIX FixedInfo()const { return fixedInfo_; }
	int SubprobNum()const { return subprobNum_; }
	MasterProb MasterProblem()const { return masterProblem_; }
	std::unordered_map<int, bool>VarInt()const { return varInt_; }
	int CouplingConsNum()const { return couplingConsNum_; }
	SubProb SubProblem(int spNo)const { return subProblem_[spNo]; }
	std::vector<Column> GetMasterProblemColumns(int spNo)const { return masterProblem_.columnpool[spNo].Columns(); }

	//get some value function
	std::vector<int> GetMasterColNum()const;

	//set some value function
	void SetFixedInfo(const FIX& fixInfo) { fixedInfo_ = fixInfo; }
	void SetMasterColumnool(int spNo, const ColumnPool& columnPool) { 
		masterProblem_.columnpool[spNo] = columnPool; 
	}; //@wqy0620 the assignment of columns need to check, need deep copy
	void SetLpRelaxObj(double lpRelaxObj) { lpRelaxObj_ = lpRelaxObj; }
	void SetLpRelaxSol(RMPSOL lpRelaxSol) { lpRelaxSol_ = lpRelaxSol; }

	//functional
	bool CheckFixFeasible(int varGlobalIndex)const; //check wheather feasible after fixing
	bool CheckFixFeasible(string varGlobalName) const; //this is in use @wqy0625  
	bool CheckFixFeasible(std::unordered_map<int, double> varValueFixed) const; //@wqy0625 in use, only fix integer var and check
 	void InitSubProblem(std::map<int, General_Block> block);
	void BuildSubProblemModel(int subprob_no);
	void BuildSubProblemModelNew(int subprob_no, const std::vector<double>& duals_part);
	void UpdateSubProbObjCoef(int subprob_no, const std::vector<double>& duals_part);
	bool AddMasterProblemColumn(const Solution& sol, int subprob_no);
	void BuildMasterProblemSCIPModel();

private:
	int index_;
	RMPSOL lpRelaxSol_;
	double lpRelaxObj_;
	int fixedSpNo_;
	Column fixedColumn_;
	FIX fixedInfo_;
	int subprobNum_; //the number of subproblem
	MasterProb masterProblem_;
	vector<SubProb> subProblem_;
	Objective obj_;
	vector<Constraints> couplingCons_;
	int couplingConsNum_;
	vector<Constraints> constraints_;
	int constraintsNum_;
	vector<ColumnPool> columnpool_;
	std::unordered_map<int, bool>varInt_; //record variable is integer or not, int->var.id; bool->yes/no

};



