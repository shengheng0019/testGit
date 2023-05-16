#include"SubProbManagement.h"
#include<string>
#include<unordered_map>

using namespace Eigen;

SubProb& SubProb::operator = (const SubProb& sp)//深拷贝
{
	if (this == &sp)
		return *this;
	this->A = sp.A;
	this->c = sp.c;
	this->colPool = sp.colPool;
	this->conss = new Constraints[sp.numConss];
	for (int i = 0; i < sp.numConss; i++)
	{
		this->conss[i] = sp.conss[i];
	}
	this->numConss = sp.numConss;
	this->oriObj = sp.oriObj;
	this->tempObjVal = sp.tempObjVal;
	this->tempSol = sp.tempSol;
	this->varDic = sp.varDic;
	this->varSize = sp.varSize;
	this->scipModel = new SCIP();
	this->scipModel = sp.scipModel;
	this->scipVars = new SCIP_VAR * [sp.varSize];
	for (int i = 0; i < sp.varSize; i++)
	{
		this->scipVars[i] = sp.scipVars[i];
	}
	this->subOrder2Id = sp.subOrder2Id;
}

bool SubProb::GetSolution()
{
	
	SCIPsetIntParam(scipModel, "display/verblevel", 0);
	if (!SCIPsolve(scipModel) == SCIP_RETCODE::SCIP_OKAY)
	{
		return false;
	}
	SCIPwriteOrigProblem(scipModel, "D://sub.lp", NULL, false);
	//write the result to the tempSol
	int tempVarCounter = 0;

	//wqy改0406，scipVars的顺序和varValue的顺序不一定相同，顺序易变，用名称来定
	//for (unordered_map<int, double>::iterator iter = tempSol.varValue.begin(); iter != tempSol.varValue.end(); iter++)
	//{
	//	iter->second = SCIPgetSolVal(scipModel, SCIPgetBestSol(scipModel), scipVars[tempVarCounter]);
	//	//cout <<"iter->first:" << iter->first << endl;
	//	//cout <<"scipVars[tempVarCounter]->name"<< scipVars[tempVarCounter]->name << endl;
	//	tempVarCounter++;
	//}
	//for (int i = 0; i < tempSol.varValue.size(); i++)
	//{
	//	string name = scipVars[i]->name;
	//	name.erase(0, 1);//把x删掉，留下数字id
	//	tempSol.varValue[stoi(name)] = SCIPgetSolVal(scipModel, SCIPgetBestSol(scipModel), scipVars[i]);
	//	if (tempSol.varValue.find(stoi(name)) == tempSol.varValue.end())
	//	{
	//		cout << "error!!!" << endl;
	//	}
	//}

	for (auto iter = varDic.begin(); iter != varDic.end(); iter++)
	{
		for (int i = 0; i < varDic.size(); i++)
		{
			string scipVarName = scipVars[i]->name;
			if (scipVarName == iter->second.Name)
			{
				tempSol.varValue[iter->first] = SCIPgetSolVal(scipModel, SCIPgetBestSol(scipModel), scipVars[i]);
				break;
			}
			if (i == varDic.size() - 1)
				cout << "Not Found!!" << endl;
		}
	}


	//tempSol.objValue = SCIPgetLPObjval(scipModel);
	tempSol.objValue = SCIPgetSolOrigObj(scipModel, SCIPgetBestSol(scipModel));
	tempSol.solSize = tempSol.varValue.size();
	//释放
	for (int i = 0; i < varDic.size(); i++)
	{
		//cout << scipVars[i]->name << endl;
		SCIPreleaseVar(scipModel, &scipVars[i]);
	}
	for (size_t i = 0; i < numConss; i++)
	{
		SCIPreleaseCons(scipModel, &scipConss[i]);
	}
	
	SCIPfree(&scipModel);

	return true;
}

void GetDualSolution(SCIP* scipModel, SCIP_CONS** scipConss, int scipConsNum, Solution& sol)
{
	SCIPsolve(scipModel);
	sol.solSize = scipConsNum;
	for (int i = 0; i < scipConsNum; i++)
	{
		sol.varValue[i] = SCIPgetDualsolLinear(scipModel, scipConss[i]);
	}
}

VectorXd Sol2Vec(Solution& sol)
{
	VectorXd vec{ sol.varValue.size() };
	int i = 0;
	for (unordered_map<int, double>::iterator iter = sol.varValue.begin(); iter != sol.varValue.end(); iter++)
	{
		vec[i++] = iter->second;
	}
	return vec;
}

SubProb::SubProb()
{
	
}
SubProb::~SubProb()
{
	/*if (conss != NULL)
		delete[] conss;
	for (int i = 0; i < varSize; i++)
	{
		SCIPreleaseVar(scipModel, &scipVars[i]);
		delete[]scipVars[i];
	}
	delete[] scipVars;

	for (int i = 0; i < numConss; i++)
	{
		SCIPreleaseCons(scipModel, &scipConss[i]);
		delete[]scipConss[i];
	}
	delete[] scipConss;

	SCIPfree(&scipModel);
	delete scipModel;*/
}

//SubProb::~SubProb()
//{
//
//}

//注意：需要实现约束、目标、变量类型的拷贝构造函数，尤其注意unordered_map类型
SubProb::SubProb(Constraints* constraints, int numConss, const Objective& obj, Constraints* couplingCons, int couplingConsNum, Vars* vars, int varNum)
	: conss{ nullptr }, numConss{ numConss }, oriObj{ obj }
{
	//copy the vars
	for (int i = 0; i < varNum; i++)
	{
		varDic.insert(pair<int, Vars>(vars[i].Id, vars[i]));
		tempObjVal.insert(pair<int, double>(vars[i].Id, 0.0));
		tempSol.varValue.insert(pair<int, double>(vars[i].Id, 0.0));
	}
	for (unordered_map<int, Vars>::iterator iter = varDic.begin(); iter != varDic.end(); iter++)
	{
		subOrder2Id.push_back(iter->first);
	}
	varSize = varDic.size();
	//copy the constraints
	conss = reinterpret_cast<Constraints*>(operator new(numConss * sizeof(Constraints)));
	for (int i = 0; i < numConss; i++)
	{
		new(&conss[i])Constraints{ constraints[i] };
	}
	//copy the objective
	for (unordered_map<int, double>::iterator iter = oriObj.coef.begin(); iter != oriObj.coef.end(); iter++)
	{
		tempObjVal[iter->first] = iter->second;
	}

	//c and A remain unchanged after this
	c = Eigen::VectorXd::Zero(varSize);
	A = Eigen::MatrixXd::Zero(couplingConsNum, varSize);
	int scipVarsIdx = 0;
	for (unordered_map<int, Vars>::iterator iter = varDic.begin(); iter != varDic.end(); iter++)
	{
		c(scipVarsIdx) = tempObjVal[iter->first];
		for (int i = 0; i < couplingConsNum; i++)
		{
			unordered_map<Vars, double>::iterator subIter = couplingCons[i].exprDic.find(iter->second);
			if (subIter == couplingCons[i].exprDic.end())
				continue;
			A(scipVarsIdx, i) = subIter->second;
		}
		scipVarsIdx++;
	}
	//SCIPsetMessagehdlrQuiet(scipModel, NULL);

}

SubProb::SubProb(vector<Constraints>& constraints, const Objective& obj, vector<Constraints>& couplingCons, const vector<Vars>& vars)
	: conss{ nullptr }, numConss{ (int)constraints.size() }, oriObj{ obj }
{
	//copy the vars
	for (int i = 0; i < vars.size(); i++)
	{
		varDic.insert(pair<int, Vars>(vars[i].Id, vars[i]));
		tempObjVal.insert(pair<int, double>(vars[i].Id, 0.0));
		tempSol.varValue.insert(pair<int, double>(vars[i].Id, 0.0));
		//subOrder2Id.push_back(vars[i].Id);
	}
	for (unordered_map<int, Vars>::iterator iter = varDic.begin(); iter != varDic.end(); iter++)
	{
		subOrder2Id.push_back(iter->first);
	}
	varSize = varDic.size();
	//copy the constraints
	conss = new Constraints[constraints.size()]{};
	for (int i = 0; i < constraints.size(); i++)
	{
		conss[i] = Constraints{ constraints[i] };
	}

	//copy the objective
	for (unordered_map<int, double>::iterator iter = oriObj.coef.begin(); iter != oriObj.coef.end(); iter++)
	{
		tempObjVal[iter->first] = iter->second;
	}

	//c and A remain unchanged after this
	c = Eigen::VectorXd::Zero(varSize);
	A = Eigen::MatrixXd::Zero(couplingCons.size(), varSize);

	for (int ii = 0; ii < subOrder2Id.size(); ii++)
	{
		int varId = subOrder2Id[ii];
		c(ii) = tempObjVal[varId];
		for (int i = 0; i < couplingCons.size(); i++)
		{
			unordered_map<Vars, double>::iterator subIter = couplingCons[i].exprDic.find(varDic[varId]);
			if (subIter == couplingCons[i].exprDic.end())
				continue;
			A(i, ii) = subIter->second;
		}
	}
}

void SubProb::BuildSCIPModel()
{
	SCIPcreate(&scipModel);
	SCIPincludeDefaultPlugins(scipModel);

	SCIPcreateProbBasic(scipModel, "subProb");
	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
	if (oriObj.direction == Direction::max)
	{
		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
	}
	SCIPsetObjsense(scipModel, objSense);

	unordered_map<int, SCIP_VAR*> scipVarMap;
	scipVars = new SCIP_VAR * [varSize];
	int scipVarsIdx = 0;
	for (int ii = 0; ii < subOrder2Id.size(); ii++)
	{
		int varId = subOrder2Id[ii];
		SCIP_VAR* tempVar;
		VarType varType = varDic[varId].Type;
		SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
		if (varType == VarType::Bool)
		{
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
		}
		else if (varType == VarType::Int)
		{
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_INTEGER;
		}
		else if (varType == VarType::Num)
		{
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
		}
		SCIPcreateVarBasic(scipModel, &tempVar, varDic[varId].Name.c_str(), varDic[varId].Lb, varDic[varId].Ub, tempObjVal[varId], scipVarType);
		SCIPaddVar(scipModel, tempVar);

		scipVars[ii] = tempVar;
		scipVarMap.insert(make_pair(varId, tempVar));
	}

	scipConss = new SCIP_CONS * [numConss];
	for (int i = 0; i < numConss; i++)
	{
		Constraints consTemp = conss[i];
		unordered_map<Vars, double>::iterator iter = consTemp.exprDic.begin();
		int tempVarsSize = consTemp.exprDic.size();
		SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];
		int idx = 0;
		while (iter != consTemp.exprDic.end())
		{
			tempVars[idx] = scipVarMap[iter->first.Id];
			tempConsCoefs[idx++] = iter->second;
			iter++;
		}

		double tempLb = consTemp.rhs, tempUb = SCIPinfinity(scipModel);
		if (consTemp.Prop == PROP::eq)
		{
			tempLb = consTemp.rhs;
			tempUb = consTemp.rhs;
		}
		else if (consTemp.Prop == PROP::leq)
		{
			tempUb = consTemp.rhs;
			tempLb = -SCIPinfinity(scipModel);
		}
		string consName = "cons" + std::to_string(i);
		SCIPcreateConsBasicLinear(scipModel, &scipConss[i], consName.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
		SCIPaddCons(scipModel, scipConss[i]);
		SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scipModel), 1);
		//delete[] tempVars;
		//delete[] tempConsCoefs;
		//SCIPwriteOrigProblem(scipModel, "D://sub1.lp", NULL, false);
	}

}

void SubProb::UpdateTempObj(Solution dualSol)
{
	SCIPwriteOrigProblem(scipModel, "D://sub2.lp", NULL, false);
	//to calculate: c - mu * A
	Eigen::VectorXd mu = Sol2Vec(dualSol);
	Eigen::VectorXd cBar = c - A.transpose() * mu;

	//cout
	/*cout << "c:" << endl;
	cout << c;
	cout << endl;
	cout << "A:" << endl;
	cout << A;
	cout << endl;
	cout << "mu: " << endl;
	cout << mu;
	cout << endl;
	cout << "子问题目标函数系数：" << endl;
	cout << cBar;*/

	//update tempObjVal with cBar
	int i = 0;
	for (unordered_map<int, double>::iterator iter = tempObjVal.begin(); iter != tempObjVal.end(); iter++)
	{
		iter->second = cBar(i++);
	}

	updateSCIPObj();
	SCIPwriteOrigProblem(scipModel, "D://sub3.lp", NULL, false);
}

void SubProb::UpdateTempObj(double* dualSol, int solSize)
{
	//to calculate: c - mu * A
	Eigen::VectorXd mu = Eigen::Map<Eigen::VectorXd>(dualSol, solSize);
	Eigen::VectorXd cBar = c - A.transpose() * mu;
	//update tempObjVal with cBar
	int i = 0;
	for (unordered_map<int, double>::iterator iter = tempObjVal.begin(); iter != tempObjVal.end(); iter++)
	{
		iter->second = cBar(i++);
	}
	updateSCIPObj();
}

void SubProb::UpdateTempObj(const VectorXd& dualSol)
{
	//to calculate: c - mu * A
	Eigen::VectorXd cBar = c - A.transpose() * dualSol;
	//update tempObjVal with cBar
	//int i = 0;
	//for (unordered_map<int, double>::iterator iter = tempObjVal.begin(); iter != tempObjVal.end(); iter++)
	//{
	//	iter->second = cBar(i++);
	//}
	for (int ii = 0; ii < subOrder2Id.size(); ii++)
	{
		tempObjVal[subOrder2Id[ii]] = cBar(ii);
	}
	//updateSCIPObj();
}


void SubProb::updateSCIPObj()
{
	int varCounter = 0;
	for (unordered_map<int, double>::iterator iter = tempObjVal.begin(); iter != tempObjVal.end(); iter++)
	{
		SCIPchgVarObj(scipModel, scipVars[varCounter++], iter->second);
	}
}

