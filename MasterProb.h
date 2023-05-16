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
	vector<int>col_var_no;//����ÿһ��ֵ��Ӧ�ı������
	int col_no;//��������������е���š������õ�ԭ��������������node.MasterProb�г�û�䣬��ſ��Խ����ã�
};

struct RMPSOL
{
	vector<int> subprob_no;//��ǰ��������һ��������
	vector<double> varValue;//��ֵ
};

class MasterProb
{
public:
	MasterProb();
	~MasterProb();
	MasterProb& operator += (const MasterProb& mp);//���
	void BuildSCIPModel(Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX> Fixed_Info);
	int Solve_Dual(Solution& duals);
	void Solve_Obj(double& obj);
	void AddColumn(Solution sol, int subprob_no);
	void Solve(RMPSOL& sol, double& objval);
	void Test_feasble2(bool& feasible, Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX>& Fixed_Info);
	void free();
;	//������г���ڵ�������гز�ͬ��������Ǽ��뵽�������е��У��ڵ����е��гؾ���ɸѡ����ˣ������������ɵ��з���ڵ�����г�
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
	Eigen::MatrixXd ConsCoefMatrix; //Լ������
};


