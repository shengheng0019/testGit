#pragma once
#include "Constraints.h"
#include"Vars.h"
#include"Solution.h"
#include<iostream>
#include<string>
#include <scip/scipdefplugins.h>
#include<Eigen/dense>
#include"Objective.h"
#include"ColumnPool.h"


using namespace Eigen;

//void BuildSCIPModel(SCIP* scipModel, SCIP_VAR** scipVars, SCIP_CONS** conss, Constraints* constraints, int numConss, Objective obj);

void GetDualSolution(SCIP* scipModel, SCIP_CONS** scipConss, int scipConsNum, Solution& sol);

VectorXd Sol2Vec(Solution& sol);

class Lagrangean;

class SubProb
{
public:
	Objective oriObj;
	Constraints* conss;
	int numConss;						//约束数量
	unordered_map<int, Vars> varDic;	//子问题包含的全部变量
	unordered_map<int, double> tempObjVal;	//当前目标函数系数，由对偶值计算出来,size=varSize,含0元素
	int varSize;						//子问题变量数
	SubProb(Constraints* constraints, int numConss, const Objective& obj, Constraints* couplingCons, int couplingConsNum, Vars* vars, int varNum);
	SubProb(vector<Constraints>& constraints, const Objective& obj, vector<Constraints>& couplingCons, const vector<Vars>& vars);
	//SubProb(Constraints* constraints, int numConss, Objective obj, Constraints* couplingCons, int couplingConsNum);
	SubProb();
	~SubProb();
	SubProb& operator = (const SubProb& sp);//深拷贝
	void BuildSCIPModel();				//构建子问题
	void UpdateTempObj(Solution dualSol);	//根据dual/lagrangian multiplyer更新tempObjVal
	void UpdateTempObj(double* dualSol, int solSize);//同上
	void UpdateTempObj(const VectorXd& dualSol);	//同上

	bool GetSolution();
	Solution tempSol;					//子问题当前解

	ColumnPool colPool;     //子问题列池
	vector<int> subOrder2Id;

	friend class Lagrangean;
private:
	SCIP* scipModel;
	SCIP_VAR** scipVars;
	SCIP_CONS** scipConss;
	Eigen::MatrixXd A;				//coupling constraints Matrix
	Eigen::VectorXd c;				//original objective coeficient
	void updateSCIPObj();			//update the obj coef of scip model without reconstructing the scip model

};

