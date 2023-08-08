#include "Lagrangean.h"
#include "LocalBranch.hpp"
using namespace std;
using namespace Eigen;

//ColumnPool columnPool;

double getUB(string fileName)
{
	SCIP* scip_ = NULL;
	SCIPcreate(&scip_);
	SCIPincludeDefaultPlugins(scip_);
	// ����һЩ������
	SCIPreadProb(scip_, fileName.c_str(), NULL);
	//SCIPsetIntParam(scip_, "presolving/maxrounds", 0);
	SCIPsetIntParam(scip_, "limits/solutions", 1);
	if (!SCIPsolve(scip_) == SCIP_RETCODE::SCIP_OKAY)
	{
		return false;
	}
	double UB_ = SCIPgetSolOrigObj(scip_, SCIPgetBestSol(scip_));	//init UB
	std::cout << "UB: " << UB_ << std::endl;
	return UB_;
}


void InitDualValwithLp(SCIP* scip_, double* dualVals, int dualSize, double& objVal, vector<string> consNames)
{
	SCIP* scip;
	SCIPcreate(&scip);
	//SCIPcreateProbBasic(scip, "subProb");
	//SCIPcopyOrigProb(scip_, scip, NULL, NULL, "copyProb");
	const char* suf = "";
	SCIPcopyOrig(scip_, scip, NULL, NULL, suf, true, false, false, NULL);

	//int nrows = SCIPgetNConss(scip);
	SCIP_VAR** vars_ori = SCIPgetVars(scip);
	SCIP_CONS** cons_ori = SCIPgetConss(scip);
	int nvars = SCIPgetNVars(scip);
	int nconss = SCIPgetNConss(scip);
	vector<SCIP_VAR*> Name_Var;
	vector<SCIP_CONS*> consVec;
	for (int i = 0; i < nvars; ++i)
		Name_Var.push_back(vars_ori[i]);
	for (int i = 0; i < nconss; i++)
		consVec.push_back(cons_ori[i]);
	for (int i = 0; i < nvars; i++)
	{
		auto var = Name_Var[i];
		SCIP_VARTYPE const varType_ = SCIPvarGetType(var);
		if (varType_ == SCIP_VARTYPE_CONTINUOUS)
			continue;

		SCIP_Bool infeasible;
		SCIPchgVarType(scip, var, SCIP_VARTYPE_CONTINUOUS, &infeasible);
	}
	SCIPwriteOrigProblem(scip, "LpRelaxation.lp", NULL, false);
	SCIPsolve(scip);
	objVal = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
	std::cout << "lp obj: " << objVal << endl;
	for (int i = 0; i < consNames.size(); i++)
	{
		for (int j = 0; j < nconss; j++)
		{
			if (consNames[i] == SCIPconsGetName(consVec[j]))
			{
				SCIP_Bool infeasible;
				SCIPgetDualSolVal(scip, consVec[j], &dualVals[i], &infeasible);
			}
		}
	}

	std::cout << "get lp realxation and dual value" << endl;
}

Lagrangean::~Lagrangean()
{
	std::free( subProbs);
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

Lagrangean::Lagrangean(map<int, General_Block>& blocks_,SCIP* scip_, vector<SCIP_CONS*> scipcons, vector<SCIP_VAR*> scipvars, unordered_map<string, int>& Var_Ni, double UB_, std::string fileName_, std::chrono::system_clock::time_point start_global_, double PI_, int MAXNOIMPNUM_, double S_, int KMAX_, double INITIALT_)
	: multiplyer{ blocks_[blocks_.size() - 1].bCons.size() }, scip{ scip_ }, UB{ UB_ }, fileName{ fileName_ }, start_global{ start_global_ }, PI{ PI_ }, MAXNOIMPNUM{ MAXNOIMPNUM_ }, S{ S_ }, KMAX{ KMAX_ }, INITIALT{ INITIALT_ },
	couplingRhs{ blocks_[blocks_.size() - 1].bCons.size() }, scipcons{ scipcons }, scipvars{ scipvars }, Var_Ni_{Var_Ni}
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


	vector<string> consNames;
	//rhs of the coupling constraints
	couplingProp = new PROP[multiplyer.size()];
	for (int i = 0; i < blocks_[subProbNum].bCons.size(); i++)
	{
		couplingRhs[i] = blocks_[subProbNum].bCons[i].rhs;
		couplingProp[i] = blocks_[subProbNum].bCons[i].Prop;
		consNames.push_back(blocks_[subProbNum].bCons[i].Name);
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
				std::cout << "Not Found!" << endl;
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
	InitMultiplyer(consNames);
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

	vector<bool> isRise;

	/// reverse Var_Ni_ ,copyright sy
	std::unordered_map<int, std::string> Var_In;
	for (const auto& pair : Var_Ni_)
	{
		Var_In[pair.second] = pair.first;
	}
    std::vector<writeJson::Soljson> soljsons;
    nlohmann::json j_array = nlohmann::json::array();  // ����һ���յ�Json����
    std::ofstream o((fileName + ".json").c_str(), std::ofstream::trunc);
	while (pi > 0.00003 && UB - LB >= 1 && iterCounter < maxIter)
	{
		//double lagVal = SearchBestT();
		/*if(iterCounter % 10 == 0)
			double scale = SearchBestT();*/
		double lagVal = 0;
		t = 1.0;
		T = 1.0;

		solveSubProbs(lagVal);
		lagVal += multiplyer.dot(couplingRhs);

		int numHalf;
		for (numHalf = 0; numHalf < 5 && lagVal < -1E10; numHalf++)
		{
			multiplyer = multiplyer / 2;
			solveSubProbs(lagVal);
			lagVal += multiplyer.dot(couplingRhs);
		}
		if (numHalf == 5)	//unbounded problems that cannot be solved currently
			return;

		//to judge the whole trend
		if (prevVal < lagVal)
		{
			isRise.push_back(true);
			if (isRise.size() > 5)
			{
				isRise.erase(isRise.begin());
				int riseNum = 0;
				for (int i = 0; i < isRise.size(); i++)
					if (isRise[i])
						riseNum++;
				if (riseNum >= 4 && pi < PI * 0.5)
					pi *= 1.1;
			}
		}
		else
		{
			isRise.push_back(false);
			if (isRise.size() > 5)
			{
				isRise.erase(isRise.begin());
				int noRiseNum = 0;
				for (int i = 0; i < isRise.size(); i++)
					if (!isRise[i])
						noRiseNum++;
				if (noRiseNum == 3)
					pi *= 0.9;
			}
		}


		prevVal = lagVal;
		std::cout << "iter: " << iterCounter << " Lagrangean dual objective value: " << lagVal << "  best value: " << LB;
		
		writeJson::Soljson tmpSol;
		auto end_global_ = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_global_ - start_global).count();
		//tmpSol.source = "L " + std::to_string(iterCounter) + " " + std::to_string(duration);
		tmpSol.obj = lagVal;
		for (int i = 0; i < masterSol.varValue.size(); ++i)
		{
			std::string tmpname = Var_In[i];
			tmpSol.solution.emplace(tmpname, masterSol.varValue[i]);
		}

//		nlohmann::json j;
//		tmpSol.to_json(j);
//		o << j << std::endl;
        soljsons.push_back(tmpSol);

		/// write lagrangean solution to json file
		/// copyright :  sy

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
				if(pi > 0.000001)
					pi *= 0.5;
				noImpCounter = 0;
			}
		}

		//update UB
		double feasibleFJVal = 0;
		if (iterCounter == 300)
		{
			if (getFeawithFJ(feasibleFJVal, multiplyer))
			{
				if (feasibleFJVal < UB)
				{
					UB = feasibleFJVal;
				}
			}
		}

		subgradient = couplingRhs - CouplingA * Sol2Vec_new(masterSol);

		//check if the solution is feasible for coupling constraints
		//bool isFeasible = true;
		int infeasibleNum = 0;
		int infeasibleEq = 0;
		for (int i = 0; i < subgradient.size(); i++)
		{
			if (couplingProp[i] == PROP::geq && subgradient[i] > 0)
			{
				infeasibleNum++;
			}
			else if (couplingProp[i] == PROP::leq && subgradient[i] < 0)
			{
				infeasibleNum++;
			}
			else if(couplingProp[i] == PROP::eq && std::fabs(subgradient[i]) > 0.00001 )
			{
				infeasibleEq++;
			}
		}
		
		std::cout << "infeasible <=/>= num = " << infeasibleNum << " infeasible == num = " << infeasibleEq << " " ;

		if (infeasibleNum == 0 && infeasibleEq == 0)
		{
			writeJson::Soljson tmpSol;
			auto end_global = std::chrono::high_resolution_clock::now();
			//auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_global - start_global).count();
			//tmpSol.source = "F" + std::to_string(iterCounter) + " " + std::to_string(duration);
			tmpSol.obj = lagVal;
			for (int i = 0; i < masterSol.varValue.size(); ++i)
			{
				std::string tmpname = Var_In[i];
				tmpSol.solution.emplace(tmpname, masterSol.varValue[i]);
			}

//			nlohmann::json j;
//			tmpSol.to_json(j);
//			o << j << std::endl;
            soljsons.push_back(tmpSol);
            for (auto& soljson : soljsons) {
                nlohmann::json j;
                soljson.to_json(j);  // ��Soljson����ת��Ϊjson
                j_array.push_back(j);  // ��json������ӵ�json������
            }

            o << j_array << std::endl;  // ������json����д���ļ�
			std::cout << "feasible solution found " << std::endl;
			UB = lagVal;
			return;
		}

		double penaltyTerm = subgradient.dot(multiplyer);
		std::cout << " penalty term: " << penaltyTerm << endl;
		if (fabs(penaltyTerm) < 0.1)			//value 1 is an arbitary value to judge if the penalty term is small enough
		{
            for (auto& soljson : soljsons) {
                nlohmann::json j;
                soljson.to_json(j);  // ��Soljson����ת��Ϊjson
                j_array.push_back(j);  // ��json������ӵ�json������
            }
            o << j_array << std::endl;  // ������json����д���ļ�
			std::cout << "terminate because it is close to optimal lagrangean dual" << endl;
			return;
		}

		//update step size
		stepSize = pi * (UB - LB) / subgradient.dot(subgradient);
		//if (stepSize < 0.000001)		//because subgradient.norm()==0.0

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


    ///
    ///  repair a solution generated by lagrange , by sy
    ///
//    std::cout<<"repair a solution generated by lagrange "<<std::endl;
//    std::unordered_map<std::string, double> lagrange_sol;
//    for (int i = 0; i < masterSol.varValue.size(); ++i)
//    {
//        std::string tmpname = Var_In[i];
//        lagrange_sol.emplace(tmpname, masterSol.varValue[i]);
//    }
//    auto [repair_obj,repair_sol] = LocalBranching::RepairwithModel(scip,10,lagrange_sol);
//    if (!repair_sol.empty())
//    {
//        writeJson::Soljson tmpSol;
//        tmpSol.source = "R ";
//        tmpSol.obj = repair_obj;
//        std::map<std::string,double> tmpmap_sol(repair_sol.begin(),repair_sol.end());
//        tmpSol.solution = tmpmap_sol;
//        soljsons.emplace_back(tmpSol);
//    }
    for (auto& soljson : soljsons) {
        nlohmann::json j;
        soljson.to_json(j);  // ��Soljson����ת��Ϊjson
        j_array.push_back(j);  // ��json������ӵ�json������
    }
    o << j_array << std::endl;  // ������json����д���ļ�
}

bool Lagrangean::InitCG()
{
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
		//columnPool.AddCol(subProbs[i].tempSol);
	}
	return true;
}

bool Lagrangean::solveSubProbs(double & val)
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
	}
	//cout << "solution: " << masterSol << endl;
	return true;
}


double Lagrangean::SearchBestT()
{
	double tMinus = -1, tPlus = INFINITY;
	double zMinus=-INFINITY, zPlus=INFINITY, s = S;
	int kMax = KMAX;
	double lagDualObj, z = -INFINITY;
	double vStar;
	t = 0;
	T = 0;
	if (!solveSubProbs(vStar))
	{
		std::cout << "no sol found!" << endl;
	}
	for(int k = 0; k < kMax; k++)
	{
		t += s;
		if (!solveSubProbs(lagDualObj))
		{
			std::cout << "no sol found!" << endl;
		}
		std::cout << "value: " << lagDualObj << endl;
		if (lagDualObj > vStar)
		{
			vStar = lagDualObj;
			T = t;
			optMasterSol = masterSol;
			//get subgradient corresponding to t
			double subgradient;
			VectorXd tempVec = couplingRhs - CouplingA * Sol2Vec(masterSol);
			subgradient = multiplyer.dot(couplingRhs - CouplingA * Sol2Vec(masterSol));
			if (subgradient < 0)
			{
				t -= (s / 2);
				if (!solveSubProbs(lagDualObj))
				{
					std::cout << "no sol found!" << endl;
				}
				if (lagDualObj > vStar)
					T = t;
				return lagDualObj;
			}
			else
			{
				;
			}

		}
		else
		{
			t -= s;
			s /= 2;
		}
	}
	if (T < 0.0000001)
		T = s;
	return vStar;

}

void Lagrangean::InitMultiplyer(const vector<string>& consNames)
{

	double* tempMultiplyer = new double[multiplyer.size()];
	double lpRlxVal;
	InitDualValwithLp(scip, tempMultiplyer, multiplyer.size(), lpRlxVal, consNames);
	multiplyer = Eigen::Map<Eigen::VectorXd>(tempMultiplyer, multiplyer.size());
	//��ʼ��UB
	//UB = 760;
	//multiplyer = VectorXd::Zero(multiplyer.size());
}

bool Lagrangean::getFeawithFJ(double& fjVal, VectorXd& vec)
{
	//double* mul = new double[scipcons.size()]{1};
	//for (int i = 0; i < scipcons.size(); i++)
	//{
	//	mul[i] = 1;
	//}
	//for (int i = 0; i < vec.size(); i++)
	//	mul[i] = 1;
	//double* Sol = new double[masterSol.varValue.size()];
	//int solCounter = 0;
	///*for (unordered_map<int, double>::iterator iter = masterSol.varValue.begin(); iter != masterSol.varValue.end(); iter++)
	//{
	//	Sol[solCounter++] = iter->second;
	//}*/
	////����˳���ձ���ID���㿪ʼ
	//for (int i = 0; i < masterSol.varValue.size(); i++)
	//{
	//	Sol[i] = masterSol.varValue[i];
	//}
	//MySolution mySol = UseFeasibilityJump(scip, scipvars, scipcons, mul, Sol,Var_Ni_);
	///*cout << "prob.violatedConstraints.size(): " << prob.violatedConstraints.size() << endl;*/
	//if (mySol.bestObjective != std::numeric_limits<double>::infinity())
	//{
	//	fjVal = mySol.bestObjective;
	//	std::cout << "FJ SOL: " << fjVal << endl;
	//	for (int i = 0; i < mySol.bestIncumbentAssignment.size(); i++)
	//	{
	//		std::cout << mySol.bestIncumbentAssignment[i] << " ";
	//	}
	//	std::cout << endl;
	//	return true;
	//}
	//else
	//{
	//	return false;
	//}
	return true;
}

