#pragma once
#include"SubProbManagement.h"
//#include"GeneralBlock.h"
//#include"FeasibilityJump.h"
#include"ColumnPool.h"
#include<map>
#include<vector>
#include"primalBlock.hpp"

using namespace std;
using namespace Eigen;

void InitDualValwithLp(SCIP* oriModel, double* dualVals, int dualSize, double& objVal, string* consNames);	//求原始模型

//class Blocks
//{
//public:
//	Constraints* conss;
//	int numConss;
//
//	Vars* vars;
//	int numVar;
//
//	Objective obj;
//
//	Blocks() {};
//};


class Lagrangean
{
public:
	SubProb* subProbs;
	int subProbNum;
	Solution masterSol;
	Solution optMasterSol;
	ColumnPool* columnpool;

	VectorXd multiplyer;

	//初始化子问题，初始化乘子
	//Lagrangean(Blocks* blocks_, int numBlock_, double PI = 2.0, int MAXNOIMPNUM = 30, double S = 0.2, int KMAX = 10, double INITIALT = 0.1);
	Lagrangean(map<int, General_Block>& blocks_, SCIP* scip_, vector<SCIP_CONS*> scipcons, vector<SCIP_VAR*> scipvars, double PI = 2.0, int MAXNOIMPNUM = 30, double S = 0.2, int KMAX = 10, double INITIALT = 0.1);
	~Lagrangean();

	void RUN(int maxIter);
	bool InitCG(ColumnPool*& columnpool);

private:
	int iterNum;
	double T;
	double t;
	double stepSize;
	double UB = INFINITY;
	double LB = -INFINITY;

	//Blocks* block;
	MatrixXd CouplingA;
	VectorXd couplingRhs;
	PROP* couplingProp;

	SCIP* scip;
	vector<SCIP_CONS*> scipcons;
	vector<SCIP_VAR*> scipvars;

	//Parameters
	const double PI;
	const int MAXNOIMPNUM;
	const double S;
	const int KMAX;
	double INITIALT;
	vector<int> order2Id;

	bool solveSubProbs(double& val);
	void InitMultiplyer(string* consNames);
	double SearchBestT();
	//bool getFeawithFJ(double& fjVal);
	VectorXd Sol2Vec_new(Solution masterSol);

};