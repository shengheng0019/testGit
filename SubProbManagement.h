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
	int numConss;						//Լ������
	unordered_map<int, Vars> varDic;	//�����������ȫ������
	unordered_map<int, double> tempObjVal;	//��ǰĿ�꺯��ϵ�����ɶ�żֵ�������,size=varSize,��0Ԫ��
	int varSize;						//�����������
	SubProb(Constraints* constraints, int numConss, const Objective& obj, Constraints* couplingCons, int couplingConsNum, Vars* vars, int varNum);
	SubProb(vector<Constraints>& constraints, const Objective& obj, vector<Constraints>& couplingCons, const vector<Vars>& vars);
	//SubProb(Constraints* constraints, int numConss, Objective obj, Constraints* couplingCons, int couplingConsNum);
	SubProb();
	~SubProb();
	SubProb& operator = (const SubProb& sp);//���
	void BuildSCIPModel();				//����������
	void UpdateTempObj(Solution dualSol);	//����dual/lagrangian multiplyer����tempObjVal
	void UpdateTempObj(double* dualSol, int solSize);//ͬ��
	void UpdateTempObj(const VectorXd& dualSol);	//ͬ��

	bool GetSolution();
	Solution tempSol;					//�����⵱ǰ��

	ColumnPool colPool;     //�������г�
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

