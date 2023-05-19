#pragma once
#include <scip/scipdefplugins.h>
#include"Constraints.h"
#include"Objective.h"
#include"ColumnPool.h"
#include<Eigen/dense>

struct FIX
{
	int sp_no;
	vector<double>col_val;
	vector<int>col_var_no;//列中每一个值对应的变量编号
	int col_no;//在这个子问题列中的序号。（能用的原因是上下两代的node.MasterProb列池没变，序号可以接着用）
};

struct RMPSOL
{
	vector<int> subprob_no;//当前解属于哪一个子问题
	vector<double> varValue;//解值
};

class MasterProb
{
public:
	MasterProb();
	~MasterProb();
	MasterProb& operator += (const MasterProb& mp);//深拷贝
	void BuildSCIPModel(Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX> Fixed_Info);
	int Solve_Dual(Solution& duals);
	void Solve_Obj(double& obj);
	void AddColumn(Solution sol, int subprob_no);
	void Solve(RMPSOL& sol, double& objval);
	void Test_feasble2(bool& feasible, Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX>& Fixed_Info);
	void free();
;	//这里的列池与节点类里的列池不同，这里的是加入到主问题中的列，节点类中的列池经过筛选加入此，拉格朗日生成的列放入节点类的列池
	ColumnPool* columnpool;
	int col_num;
	int sp_num;

private:
	SCIP* RMPscipModel;
	SCIP_VAR** RMPscipVars;
	int RMPscipVars_num;
	SCIP_CONS** RMPscipConss;
	int RMPscipConss_num;
	Eigen::MatrixXd A;				//coupling constraints Matrix
	Eigen::MatrixXd ConsCoefMatrix; //约束矩阵
};


