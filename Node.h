#pragma once
#include "ColumnPool.h"
#include "MasterProb.h"
#include "SubProbManagement.h"



class Node
{
public:

	int index; //当前节点的序号

	RMPSOL LPrelax_sol;//lp松弛解
	double LPrelax_obj;//lp松弛目标值
	int fixed_sp_no;//固定的列属于哪个子问题
	Column fixed_column;//本次被固定的列
	vector<FIX> Fixed_Info;//固定的信息（禁忌表）

	int subprob_num;//子问题数量

	MasterProb MasterProblem; //主问题
	SubProb* SubProblem; //子问题

	//定义一些约束类，目的是可以直接根据node构建模型，递归调用时直接修改node的约束类就可以
	Objective obj;
	Constraints* couplingCons;
	int couplingCons_num;
	Constraints* constraints;
	int constraints_num;
	ColumnPool* columnpool;


	Node();
	Node(ColumnPool* columnpool);
	Node(Objective obj, Constraints* couplingCons, int couplingCons_num,
		Constraints* constraints, int constraints_num, ColumnPool* columnpool, int subprob_num);
	Node(const Node& node);
	void free();
	~Node();
	Node& operator = (const Node& node);

private:


};


