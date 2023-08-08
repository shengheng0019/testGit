#include "Node.h"



//Node& Node::operator = (const Node& node)//0530@wqy deep copy
//{
//	if (this == &node)
//		return *this;
//
//	this->columnpool = new ColumnPool[node.subprob_num];
//	for (int i = 0; i < node.subprob_num; i++)
//	{
//		this->columnpool[i].var_num = node.columnpool[i].var_num;
//		this->columnpool[i].vars_id = node.columnpool[i].vars_id;
//		this->columnpool[i].columns = node.columnpool[i].columns;
//	}
//	this->constraints = new Constraints[node.constraints_num];
//	for (int i = 0; i < node.constraints_num; i++)
//	{
//		this->constraints[i] = node.constraints[i];
//	}
//	this->constraints_num = node.constraints_num;
//	this->couplingCons = new Constraints[node.couplingCons_num];
//	for (int i = 0; i < node.couplingCons_num; i++)
//	{
//		this->couplingCons[i] = node.couplingCons[i];
//	}
//	this->couplingCons_num = node.couplingCons_num;
//	this->fixed_column = node.fixed_column;
//	this->Fixed_Info = node.Fixed_Info;
//	this->fixed_sp_no = node.fixed_sp_no;
//	this->LPrelax_obj = node.LPrelax_obj;
//	this->LPrelax_sol = node.LPrelax_sol;
//	this->MasterProblem += node.MasterProblem;
//	this->index = node.index;
//	this->obj = node.obj;
//	this->SubProblem = new SubProb[node.subprob_num];
//	for (int i = 0; i < node.subprob_num; i++)
//	{
//		this->SubProblem[i] = node.SubProblem[i];
//	}
//	this->subprob_num = node.subprob_num;
//	this->varInt = node.varInt;
//
//	return *this;
//
//}
//
//Node::Node(const Node& node)
//{
//	//columnpool = node.columnpool;
//	this->columnpool = new ColumnPool[node.subprob_num];
//	for (int i = 0; i < node.subprob_num; i++)
//	{
//		this->columnpool[i].var_num = node.columnpool[i].var_num;
//		this->columnpool[i].vars_id = node.columnpool[i].vars_id;
//		this->columnpool[i].columns = node.columnpool[i].columns;
//	}
//
//	//constraints = node.constraints;
//	this->constraints = new Constraints[node.constraints_num];
//	for (int i = 0; i < node.constraints_num; i++)
//	{
//		this->constraints[i] = node.constraints[i];
//	}
//	this->constraints_num = node.constraints_num;
//	//couplingCons = node.couplingCons;
//	this->couplingCons = new Constraints[node.couplingCons_num];
//	for (int i = 0; i < node.couplingCons_num; i++)
//	{
//		this->couplingCons[i] = node.couplingCons[i];
//	}
//	this->couplingCons_num = node.couplingCons_num;
//	this->fixed_column = node.fixed_column;
//	this->Fixed_Info = node.Fixed_Info;
//	this->fixed_sp_no = node.fixed_sp_no;
//	this->LPrelax_obj = node.LPrelax_obj;
//	this->LPrelax_sol = node.LPrelax_sol;
//	this->MasterProblem += node.MasterProblem;//���
//	//this->MasterProblem.sp_num = node.MasterProblem.sp_num;
//	//this->MasterProblem.col_num = node.MasterProblem.col_num;
//	//this->MasterProblem.columnpool = new ColumnPool[node.subprob_num];
//	//for (int i = 0; i < node.subprob_num; i++) 
//	//{
//	//	this->MasterProblem.columnpool[i].var_num = node.MasterProblem.columnpool[i].var_num;
//	//	this->MasterProblem.columnpool[i].vars_id = node.MasterProblem.columnpool[i].vars_id;
//	//	this->MasterProblem.columnpool[i].columns = node.MasterProblem.columnpool[i].columns;
//	//}
//	this->index = node.index;
//	this->obj = node.obj;
//	//SubProblem = node.SubProblem;
//	this->SubProblem = new SubProb[node.subprob_num];
//	for (int i = 0; i < node.subprob_num; i++)
//	{
//		this->SubProblem[i] = node.SubProblem[i];
//	}
//	this->subprob_num = node.subprob_num;
//}

std::vector<int> Node::GetMasterColNum()const
{
	vector<int>retNum;
	for (int i = 0; i < subprobNum_; i++){
		retNum.push_back(masterProblem_.columnpool[i].Columns().size());
	}
	return retNum;
}

bool Node::CheckFixFeasible(int varGlobalIndex) const
{
	return masterProblem_.CheckSet1(varGlobalIndex);
}

bool Node::CheckFixFeasible(string varGlobalName) const
{
	return masterProblem_.CheckSet1(varGlobalName);
}

bool Node::CheckFixFeasible(std::unordered_map<int, double> varValueFixed) const
{
	//return masterProblem_.CheckSet1(varValueFixed, varInt_);
	return 0;
}

void Node::BuildMasterProblemSCIPModel()
{
	masterProblem_.BuildSCIPModel(obj_, couplingCons_, couplingConsNum_, subprobNum_, fixedInfo_);
}

void Node::InitSubProblem(std::map<int, General_Block> block)
{
	for (int i = 0; i < subprobNum_; i++){
		subProblem_.push_back(SubProb(block[i].bCons, block[i].bobj, block[subprobNum_].bCons, block[i].bVars, true));
	}
}

void  Node::BuildSubProblemModel(int subprob_no)
{
	subProblem_[subprob_no].BuildSCIPModelNew(fixedInfo_, subprob_no);
}

void  Node::UpdateSubProbObjCoef(int subprob_no, const std::vector<double>& duals_part)
{
	subProblem_[subprob_no].UpdateTempObjNew(duals_part);
}

void Node::BuildSubProblemModelNew(int subprob_no, const std::vector<double>& duals_part)
{
	subProblem_[subprob_no].BuildSCIPModelNew(fixedInfo_, subprob_no, duals_part);
}

bool Node::AddMasterProblemColumn(const Solution& sol, int subprob_no)
{
	return masterProblem_.AddColumn(sol, subprob_no);
}