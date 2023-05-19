#include "Node.h"


Node::Node()
{
	
}

Node::~Node()
{
}

Node::Node(ColumnPool* columnpool)
{
	this->columnpool = columnpool;
}

//释放内存
void Node::free() 
{
	delete[] couplingCons;
	delete[] constraints;
	delete[] columnpool;
	delete[] SubProblem;

	MasterProblem.free();
}

Node::Node(Objective obj, Constraints* couplingCons, int couplingCons_num,
	Constraints* constraints, int constraints_num, ColumnPool* columnpool, int subprob_num)
{
	this->obj = obj;
	this->couplingCons = couplingCons;
	this->couplingCons_num = couplingCons_num;
	this->constraints = constraints;
	this->constraints_num = constraints_num;
	this->columnpool = columnpool;
	this->subprob_num = subprob_num;

	this->index = 0;


	this->MasterProblem.columnpool = columnpool;//初始时，直接赋值，之后需要判断才能进主问题列池
	this->MasterProblem.sp_num = subprob_num;
	//this->MasterProblem.col_num = col_num;
}


Node& Node::operator = (const Node& node)
{
	if (this == &node)
		return *this;
	this->columnpool = node.columnpool;
	this->constraints = node.constraints;
	this->constraints_num = node.constraints_num;
	this->couplingCons = node.couplingCons;
	this->couplingCons_num = node.couplingCons_num;
	this->fixed_column = node.fixed_column;
	this->Fixed_Info = node.Fixed_Info;
	this->fixed_sp_no = node.fixed_sp_no;
	this->LPrelax_obj = node.LPrelax_obj;
	this->LPrelax_sol = node.LPrelax_sol;
	this->MasterProblem = node.MasterProblem;
	this->index = node.index;
	this->obj = node.obj;
	this->SubProblem = node.SubProblem;
	this->subprob_num = node.subprob_num;


	return *this;

}

//深拷贝构造函数
Node::Node(const Node& node)
{
	//columnpool = node.columnpool;
	this->columnpool = new ColumnPool[node.subprob_num];
	for (int i = 0; i < node.subprob_num; i++)
	{
		this->columnpool[i].var_num = node.columnpool[i].var_num;
		this->columnpool[i].vars_id = node.columnpool[i].vars_id;
		this->columnpool[i].columns = node.columnpool[i].columns;
	}
	
	//constraints = node.constraints;
	this->constraints = new Constraints[node.constraints_num];
	for (int i = 0; i < node.constraints_num; i++)
	{
		this->constraints[i] = node.constraints[i];
	}
	this->constraints_num = node.constraints_num;
	//couplingCons = node.couplingCons;
	this->couplingCons = new Constraints[node.couplingCons_num];
	for (int i = 0; i < node.couplingCons_num; i++)
	{
		this->couplingCons[i] = node.couplingCons[i];
	}
	this->couplingCons_num = node.couplingCons_num;
	this->fixed_column = node.fixed_column;
	this->Fixed_Info = node.Fixed_Info;
	this->fixed_sp_no = node.fixed_sp_no;
	this->LPrelax_obj = node.LPrelax_obj;
	this->LPrelax_sol = node.LPrelax_sol;
	this->MasterProblem += node.MasterProblem;//深拷贝
	//this->MasterProblem.sp_num = node.MasterProblem.sp_num;
	//this->MasterProblem.col_num = node.MasterProblem.col_num;
	//this->MasterProblem.columnpool = new ColumnPool[node.subprob_num];
	//for (int i = 0; i < node.subprob_num; i++) 
	//{
	//	this->MasterProblem.columnpool[i].var_num = node.MasterProblem.columnpool[i].var_num;
	//	this->MasterProblem.columnpool[i].vars_id = node.MasterProblem.columnpool[i].vars_id;
	//	this->MasterProblem.columnpool[i].columns = node.MasterProblem.columnpool[i].columns;
	//}
	this->index = node.index;
	this->obj = node.obj;
	//SubProblem = node.SubProblem;
	this->SubProblem = new SubProb[node.subprob_num];
	for (int i = 0; i < node.subprob_num; i++)
	{
		this->SubProblem[i] = node.SubProblem[i];
	}
	this->subprob_num = node.subprob_num;
}