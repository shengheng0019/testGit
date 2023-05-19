#include "Lagrangean.h"
using namespace std;
using namespace Eigen;

ColumnPool columnPool;


void InitDualValwithLp(SCIP* scip, double* dualVals, int dualSize, double& objVal, string* consNames)
{
	//int nrows = SCIPgetNConss(scip);
	int ncols = SCIPgetNVars(scip);
	SCIP_CONS** conss = SCIPgetConss(scip);
	SCIP_VAR** vars = SCIPgetVars(scip);

	unsigned mark;
	for (int i = 0; i < ncols; ++i)
	{
		SCIPchgVarType(scip, vars[i], SCIP_VARTYPE_CONTINUOUS, &mark);
	}
	SCIPwriteOrigProblem(scip, "subLP.lp", NULL, false);
	SCIPsolve(scip);
	/*for (int i = 0; i < dualSize; ++i)
	{
		for (int j = 0; j < SCIPgetNConss(scip); j++)
		{
			if (SCIPconsGetName(conss[j]) == consNames[i])
			{
				SCIPgetDualSolVal(scip, conss[i], &dualVals[i], &mark);
				break;
			}
		}
	}*/
	objVal = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
	delete[] conss;
	delete[] vars;
}
Lagrangean::~Lagrangean()
{
	free(subProbs);
}

//Lagrangean::Lagrangean(Blocks* blocks_, int numBlock_, double PI_, int MAXNOIMPNUM_, double S_, int KMAX_, double INITIALT_)
//	: multiplyer{ blocks_[numBlock_ - 1].numConss }, PI{ PI_ }, MAXNOIMPNUM{ MAXNOIMPNUM_ }, S{ S_ }, KMAX{ KMAX_ }, INITIALT{ INITIALT_ }, couplingRhs{ blocks_[numBlock_ - 1].numConss }
//{
//	subProbNum = numBlock_ - 1;
//	//Init Subproblems
//	//subProbs = new SubProb[subProbNum];
//	subProbs = reinterpret_cast<SubProb*>(operator new(subProbNum * sizeof(SubProb)));
//	for (int i = 0; i < subProbNum; i++)
//	{
//		new(&subProbs[i]) SubProb{ blocks_[i].conss, blocks_[i].numConss, blocks_[i].obj, blocks_[subProbNum].conss, blocks_[subProbNum].numConss, blocks_[i].vars, blocks_[i].numVar };
//		for (unordered_map<int, double>::iterator iter = subProbs[i].tempSol.varValue.begin(); iter != subProbs[i].tempSol.varValue.end(); iter++)
//		{
//			masterSol.varValue.insert(pair<int, double>(iter->first, 0.0));
//		}
//	}	
//
//	//rhs of the coupling constraints
//	couplingProp = new PROP[multiplyer.size()];
//	for (int i = 0; i < blocks_[subProbNum].numConss; i++)
//	{
//		couplingRhs[i] = blocks_[subProbNum].conss[i].rhs;
//		couplingProp[i] = blocks_[subProbNum].conss[i].Prop;
//	}
//	InitMultiplyer();
//}

Lagrangean::Lagrangean(map<int, General_Block>& blocks_, SCIP* scip_, vector<SCIP_CONS*> scipcons, vector<SCIP_VAR*> scipvars, double PI_, int MAXNOIMPNUM_, double S_, int KMAX_, double INITIALT_)
	: multiplyer{ blocks_[blocks_.size() - 1].bCons.size() }, scip{ scip_ }, PI{ PI_ }, MAXNOIMPNUM{ MAXNOIMPNUM_ }, S{ S_ }, KMAX{ KMAX_ }, INITIALT{ INITIALT_ },
	couplingRhs{ blocks_[blocks_.size() - 1].bCons.size() }, scipcons{ scipcons }, scipvars{ scipvars }
{
	subProbNum = blocks_.size() - 1;
	//Init Subproblems
	//subProbs = new SubProb[subProbNum];
	subProbs = reinterpret_cast<SubProb*>(operator new(subProbNum * sizeof(SubProb)));
	for (int i = 0; i < subProbNum; i++)
	{
		new(&subProbs[i]) SubProb{ blocks_[i].bCons, blocks_[i].bobj, blocks_[subProbNum].bCons, blocks_[i].bVars };
		//delete &subProbs[i];
		for (unordered_map<int, double>::iterator iter = subProbs[i].tempSol.varValue.begin(); iter != subProbs[i].tempSol.varValue.end(); iter++)
		{
			masterSol.varValue.insert(pair<int, double>(iter->first, 0.0));
		}
	}
	for (int i = 0; i < blocks_[subProbNum].bVars.size(); i++)
	{
		masterSol.varValue.insert(pair<int, double>(blocks_[subProbNum].bVars[i].Id, 0.0));
	}

	for (auto iter = masterSol.varValue.begin(); iter != masterSol.varValue.end(); iter++)
	{
		order2Id.push_back(iter->first);
	}

	//rhs of the coupling constraints
	couplingProp = new PROP[multiplyer.size()];
	for (int i = 0; i < blocks_[subProbNum].bCons.size(); i++)
	{
		couplingRhs[i] = blocks_[subProbNum].bCons[i].rhs;
		couplingProp[i] = blocks_[subProbNum].bCons[i].Prop;
	}

	//init A
	CouplingA = MatrixXd::Zero(multiplyer.size(), masterSol.varValue.size());	//resize(multiplyer.size(), masterSol.varValue.size());
	int AColCount = 0;
	for (auto iter = masterSol.varValue.begin(); iter != masterSol.varValue.end(); iter++)
	{
		for (int i = 0; i < subProbNum; i++)
		{
			bool isFound = false;
			int subProbACount = 0;
			for (auto iter_ = subProbs[i].varDic.begin(); iter_ != subProbs[i].varDic.end(); iter_++)
			{
				if (iter_->first == iter->first)
				{
					CouplingA.col(AColCount++) = subProbs[i].A.col(subProbACount);
					isFound = true;
				}
				subProbACount++;
			}
			if (isFound)
				break;
		}
	}

	for (int ii = 0; ii < order2Id.size(); ii++)
	{
		for (int i = 0; i < subProbNum; i++)
		{
			bool isFound = false;
			for (int j = 0; j < subProbs[i].subOrder2Id.size(); j++)
			{
				if (subProbs[i].subOrder2Id[j] == order2Id[ii])
				{
					CouplingA.col(ii) = subProbs[i].A.col(j);
					isFound = true;
				}
			}
			if (isFound)
				break;
			if (i == subProbNum - 1)
				cout << "Not Found!" << endl;
		}
	}

	/*cout << "couplingA: " << endl;
	for (int i = 0; i < CouplingA.rows(); i++)
	{
		for (int j = 0; j < CouplingA.cols(); j++)
			cout << CouplingA(i, j) << " ";
		cout << endl;
	}*/

	//init multiplyer with lp relaxation
	string* consNames = new string[multiplyer.size()];
	InitMultiplyer(consNames);
	delete[] consNames;
}

VectorXd Lagrangean::Sol2Vec_new(Solution sol)
{
	VectorXd vec{ sol.varValue.size() };
	for (int i = 0; i < order2Id.size(); i++)
		vec[i] = sol.varValue[order2Id[i]];
	return vec;

}

void Lagrangean::RUN(int maxIter)
{
	/*for (int i = 0; i < subProbNum; i++)
	{
		subProbs[i].BuildSCIPModel();
	}*/

	double pi = PI;
	int noImpCounter = 0;
	int maxNoImpNum = MAXNOIMPNUM;
	int iterCounter = 0;
	VectorXd subgradient;
	VectorXd prevMul;
	double prevVal = -INFINITY;
	while (pi > 0.00000 && UB - LB >= 1 && iterCounter < maxIter)
	{
		//double lagVal = SearchBestT();
		double lagVal = 0;
		t = 1.0;
		T = 1.0;

		solveSubProbs(lagVal);
		lagVal += multiplyer.dot(couplingRhs);
		while (lagVal < -1E10)
		{
			multiplyer = multiplyer / 2;
			cout << endl;
			solveSubProbs(lagVal);
			lagVal += multiplyer.dot(couplingRhs);
		}

		prevVal = lagVal;
		cout << "iter: " << iterCounter << " Lagrangean dual objective value: " << lagVal << "  best value: " << LB;

		if (lagVal > LB)
		{
			LB = lagVal;
			noImpCounter = 0;
		}
		else
		{
			noImpCounter++;
			if (noImpCounter == maxNoImpNum)
			{
				pi *= 0.5;
				noImpCounter = 0;
			}
		}

		/*double feasibleFJVal = 0;
		if (iterCounter > 1000)
		{
			if (getFeawithFJ(feasibleFJVal, multiplyer))
			{
				if (feasibleFJVal < UB)
				{
					UB = feasibleFJVal;
				}
			}
		}*/

		subgradient = couplingRhs - CouplingA * Sol2Vec_new(masterSol);
		cout << " penalty term: " << subgradient.dot(multiplyer) << endl;

		//update step size
		stepSize = pi * (UB - LB) / subgradient.dot(subgradient);
		//if (stepSize < 0.000001)		//because subgradient.norm()==0.0
		//	break;
		prevMul = multiplyer;
		for (int i = 0; i < multiplyer.size(); i++)
		{
			if (couplingProp[i] == PROP::geq)
				multiplyer[i] = max(subgradient[i] * stepSize + multiplyer[i] * T, 0.0);
			else if (couplingProp[i] == PROP::leq)
				multiplyer[i] = min(subgradient[i] * stepSize + multiplyer[i] * T, 0.0);
			else
				multiplyer[i] = subgradient[i] * stepSize + multiplyer[i] * T;
		}


		iterCounter++;
	}
	//wqy加，把每一个子问题的列池汇总到lag.columnpool，方便在外面访问
	columnpool = new ColumnPool[subProbNum];
	for (size_t i = 0; i < subProbNum; i++)
	{
		columnpool[i] = subProbs[i].colPool;
	}
}


bool Lagrangean::InitCG(ColumnPool*& columnpool)
{
	columnpool = new ColumnPool[subProbNum];
	for (int i = 0; i < subProbNum; i++)
	{
		subProbs[i].BuildSCIPModel();
		subProbs[i].UpdateTempObj(multiplyer);
		if (!subProbs[i].GetSolution())
		{
			return false;
		}
		for (unordered_map<int, double>::iterator iter = subProbs[i].tempSol.varValue.begin(); iter != subProbs[i].tempSol.varValue.end(); iter++)
		{
			masterSol.varValue[iter->first] = iter->second;
		}
		columnpool[i].AddCol(subProbs[i].tempSol);
		for (auto iter = subProbs[i].varDic.begin(); iter != subProbs[i].varDic.end(); iter++)
		{
			columnpool[i].vars_id.push_back(iter->first);
		}
		columnpool[i].var_num = columnpool[i].vars_id.size();

	}
	return true;
}

bool Lagrangean::solveSubProbs(double& val)
{
	//solve subproblems in order
	val = 0;
	for (int i = 0; i < subProbNum; i++)
	{
		subProbs[i].UpdateTempObj(multiplyer * t);
		subProbs[i].BuildSCIPModel();
		if (!subProbs[i].GetSolution())
		{
			return false;
		}
		for (unordered_map<int, double>::iterator iter = subProbs[i].tempSol.varValue.begin(); iter != subProbs[i].tempSol.varValue.end(); iter++)
		{
			masterSol.varValue[iter->first] = iter->second;
		}
		subProbs[i].colPool.AddCol(subProbs[i].tempSol);
		val += subProbs[i].tempSol.objValue;
		if (subProbs[i].colPool.vars_id.size() == 0)
		{
			for (auto iter = subProbs[i].varDic.begin(); iter != subProbs[i].varDic.end(); iter++)
			{
				subProbs[i].colPool.vars_id.push_back(iter->first);
			}
			subProbs[i].colPool.var_num = subProbs[i].colPool.vars_id.size();
		}
	}
	//cout << "solution: " << masterSol << endl;
	return true;
}

double Lagrangean::SearchBestT()
{
	double tMinus = -1, tPlus = INFINITY;
	double zMinus = -INFINITY, zPlus = INFINITY, s = S;
	int kMax = KMAX;
	double lagDualObj, z = -INFINITY;
	t = 0.0;
	for (int k = 0; k < kMax; k++)
	{
		if (!solveSubProbs(lagDualObj))
		{
			cout << "no sol found!" << endl;
		}
		cout << "value: " << lagDualObj << endl;
		if (lagDualObj > z)
		{
			z = lagDualObj;
			T = t;
			optMasterSol = masterSol;
			//get subgradient corresponding to t
			double subgradient;
			for (int i = 0; i < multiplyer.size(); i++)
			{
				cout << multiplyer[i] << " ";
			}
			cout << endl;
			VectorXd tempVec = couplingRhs - CouplingA * Sol2Vec(masterSol);
			for (int i = 0; i < tempVec.size(); i++)
			{
				cout << tempVec[i] << " ";
			}
			cout << endl;
			subgradient = multiplyer.dot(couplingRhs - CouplingA * Sol2Vec(masterSol));
			if (fabs(subgradient) < 0.001)
			{
				break;
			}
			else if (subgradient > 0)
			{
				tMinus = t;
				zMinus = z;
				if (tPlus == INFINITY)	//if tPlus undefined
					t += s;
				else
				{
					//try to improve T
					t = (zMinus * tMinus + zPlus * tPlus) / (zMinus + zPlus);
					if (!solveSubProbs(lagDualObj))
					{
						cout << "no sol found!" << endl;
					}
					if (lagDualObj > z)
					{
						T = t;
						z = lagDualObj;
						optMasterSol = masterSol;
					}
					break;
				}


			}
			else
			{
				tPlus = t;
				zPlus = z;
				if (tMinus == -1)
					t -= s;
				else
				{
					t = (zMinus * tMinus + zPlus * tPlus) / (zMinus + zPlus);
					if (!solveSubProbs(lagDualObj))
					{
						cout << "no sol found!" << endl;
					}
					if (lagDualObj > z)
					{
						T = t;
						z = lagDualObj;
						optMasterSol = masterSol;
					}
					break;
				}
			}

		}
		else
		{
			t -= (s / 2);
			if (!solveSubProbs(lagDualObj))
			{
				cout << "no sol found!" << endl;
			}
			if (lagDualObj > z)
			{
				T = t;
				z = lagDualObj;
				optMasterSol = masterSol;
			}
			break;
		}
	}

	return z;

}

void Lagrangean::InitMultiplyer(string* consNames)
{

	double* tempMultiplyer = new double[multiplyer.size()];
	double lpRlxVal = 0;
	InitDualValwithLp(scip, tempMultiplyer, multiplyer.size(), lpRlxVal, consNames);
	multiplyer = Eigen::Map<Eigen::VectorXd>(tempMultiplyer, multiplyer.size());
	//初始化UB
	UB = lpRlxVal;
	//UB = 7;
	multiplyer = VectorXd::Zero(multiplyer.size());
}

//bool Lagrangean::getFeawithFJ(double& fjVal)
//{
//	Problem prob = FJ_for_LD(scip, scipvars, scipcons, optMasterSol);
//	if (prob.flag)
//	{
//		fjVal = prob.incumbentObjective;
//		return true;
//	}
//	else
//	{
//		return false;
//	}
//}

