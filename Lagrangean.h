#pragma once
#include"SubProbManagement.h"
#include"GeneralBlock.h"
#include"ColumnPool.h"
#include "feasibilityjump.hh"
#include<map>
#include<vector>
#include<fstream>
#include "SolJson.hpp"
#include<chrono>
#include <filesystem>
#include <fstream>

using namespace std;
using namespace Eigen;

void InitDualValwithLp(SCIP* oriModel, double* dualVals, int dualSize, double& objVal, string* consNames);	//求原始模型

double getUB(string fileName);

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

	VectorXd multiplyer;
	double UB = INFINITY;
	double LB = -INFINITY;

	//初始化子问题，初始化乘子
	//Lagrangean(Blocks* blocks_, int numBlock_, double PI = 2.0, int MAXNOIMPNUM = 30, double S = 0.2, int KMAX = 10, double INITIALT = 0.1);
	Lagrangean(map<int, General_Block> & blocks_, SCIP* scip_, vector<SCIP_CONS*> scipcons, vector<SCIP_VAR*> scipvars, unordered_map<string, int>& Var_Ni, double UB_, std::string fileName_, std::chrono::system_clock::time_point start_global_, double PI = 0.001, int MAXNOIMPNUM = 15, double S = 0.2, int KMAX = 10, double INITIALT = 0.1);
	~Lagrangean();

	void RUN(int maxIter);
	bool InitCG();

private:
	int iterNum;
	double T;
	double t;
	double stepSize;
	

	//Blocks* block;
	MatrixXd CouplingA;
	VectorXd couplingRhs;
	PROP* couplingProp;

	SCIP* scip;
	vector<SCIP_CONS*> scipcons;
	vector<SCIP_VAR*> scipvars;
	unordered_map<string, int> Var_Ni_;
	string fileName;
	std::chrono::system_clock::time_point start_global;

	//Parameters
	const double PI;
	const int MAXNOIMPNUM;
	const double S;
	const int KMAX;
	double INITIALT;

	vector<int> order2Id;

	bool solveSubProbs(double& val);
	void InitMultiplyer(const vector<string>& consNames);
	double SearchBestT();
	bool getFeawithFJ(double& fjVal, VectorXd& vec);
	VectorXd Sol2Vec_new(Solution masterSol);
	
};