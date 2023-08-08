#include"SubProbManagement.h"
#include<string>
#include<unordered_map>

//using namespace Eigen;

//TODO
SubProb& SubProb::operator = (const SubProb& sp)
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
	return *this;
}

bool SubProb::GetSolution()
{

	if (!SCIPsolve(scipModel) == SCIP_RETCODE::SCIP_OKAY)
	{
		return false;
	}
	SCIPwriteOrigProblem(scipModel, "sub.lp", NULL, false);

	for (auto iter = varDic.begin(); iter != varDic.end(); iter++)
	{
		for (int i = 0; i < varDic.size(); i++)
		{
			string scipVarName = SCIPvarGetName(scipVars[i]);
			if (scipVarName == iter->second.Name)
			{
				tempSol.varValue[iter->first] = SCIPgetSolVal(scipModel, SCIPgetBestSol(scipModel), scipVars[i]);
				break;
			}
			if(i == varDic.size()-1)
				cout << "Not Found!!" << endl;
		}
	}
	

	//tempSol.objValue = SCIPgetLPObjval(scipModel);
	tempSol.objValue = SCIPgetSolOrigObj(scipModel, SCIPgetBestSol(scipModel));
	for (int i = 0; i < varSize; i++)
		SCIPreleaseVar(scipModel, &scipVars[i]);
	for (int i = 0; i < numConss; i++)
		SCIPreleaseCons(scipModel, &scipConss[i]);
	SCIPfree(&scipModel);
	return true;
}

bool SubProb::GetSolutionNew(Solution& columnGen)
{
	//SCIPwriteOrigProblem(scipModel, "D://SubproblemBeforeSolve.lp", nullptr, 0);
	SCIPsetIntParam(scipModel, "display/verblevel", 0);
	SCIPsolve(scipModel);
	if (SCIPgetStatus(scipModel) != SCIP_STATUS_OPTIMAL) {
		cout << "SCIPStatus: "<< SCIPgetStatus(scipModel) << endl;
		spdlog::error("There is some wrong when solving the subproblem!");
		return false;
	}
	
	columnGen.objValue = SCIPgetSolOrigObj(scipModel, SCIPgetBestSol(scipModel));
	columnGen.solSize = varSize;

	//@wqy vector remeber to clear().;
	//@wqy0625 the order of varsVec_ and scipVarVec_ are the same, it may not need to compare the name.
	for (int i = 0; i < varSize; i++) {
		bool flag(false);
		string scipVarName = SCIPvarGetName(scipVarsVec_[i]);
		for (int j = 0; j < varSize; j++) {
			if (scipVarName == varsVec_[j].Name || scipVarName == "t_" + varsVec_[j].Name) {
				columnGen.varValue.insert(make_pair(varsVec_[j].Id, SCIPgetSolVal(scipModel, SCIPgetBestSol(scipModel), scipVarsVec_[i])));
				flag = true;
				break;
			}
		}
		if (!flag) {
			spdlog::error("There is some wrong in converting the subproblem solution into column!");
		}
	}

	//free
	for (int i = 0; i < varSize; i++){
		SCIPreleaseVar(scipModel, &scipVarsVec_[i]);
	}
	for (int i = 0; i < numConss; i++){
		SCIPreleaseCons(scipModel, &scipConssVec_[i]);
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

Eigen::VectorXd Sol2Vec(Solution& sol)
{
	Eigen::VectorXd vec{ sol.varValue.size() };
	int i = 0;
	for (unordered_map<int, double>::iterator iter = sol.varValue.begin(); iter != sol.varValue.end(); iter++)
	{
		vec[i++] = iter->second;
	}
	
	return vec;
}

SubProb::~SubProb()
{
    //@wqy delete0701 the constructor is different
/*	if(conss != NULL)
		delete[] conss;
	for (int i = 0; i < varSize; i++)
	{
		SCIPreleaseVar(scipModel, &scipVars[i]);
		delete []scipVars[i];
	}
	delete[] scipVars;

	for (int i = 0; i < numConss; i++)
	{
		SCIPreleaseCons(scipModel, &scipConss[i]);
		delete []scipConss[i];
	}
	delete[] scipConss;
	SCIPfree(&scipModel);
	delete scipModel;*/
}

Eigen::VectorXd Vector2Vec(const vector<double>& sol)
{
	Eigen::VectorXd vec{ sol.size() };
	for (int i = 0; i< sol.size(); i++ )
	{
		vec[i] = sol[i];
	}
	return vec;
}


SubProb::SubProb(Constraints* constraints, int numConss, const Objective& obj, Constraints* couplingCons, int couplingConsNum, Vars* vars, int varNum)
	: conss{nullptr}, numConss{numConss}, oriObj{obj}
{
	/*for (int i = 0; i < numConss; i++)
	{
		Constraints consTemp = conss[i];
		unordered_map<Vars, double>::iterator iter = consTemp.exprDic.begin();
		while (iter != consTemp.exprDic.end())
		{
			varDic.insert(pair<int, Vars>(iter->first.Id, iter->first));
			tempObjVal.insert(pair<int, double>(iter->first.Id, 0.0));
			iter++;
		}
	}*/
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
	for (unordered_map<int, double>::iterator iter= oriObj.coef.begin(); iter != oriObj.coef.end(); iter++)
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
	/*for (unordered_map<int, Vars>::iterator iter = varDic.begin(); iter != varDic.end(); iter++)
	{
		c(scipVarsIdx) = tempObjVal[iter->first];
		for (int i = 0; i < couplingCons.size(); i++)
		{
			unordered_map<Vars, double>::iterator subIter = couplingCons[i].exprDic.find(iter->second);
			if (subIter == couplingCons[i].exprDic.end())
				continue;
			A(i, scipVarsIdx) = subIter->second;
		}
		scipVarsIdx++;
	}*/

	/*cout << "A: " << endl;
	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < A.cols(); j++)
			cout << A(i, j) << " ";
		cout << endl;
	}*/

}


//@wqy add new constractor for subproblem, the paramter of this function add a bool in the end
//only init the variables which diving using in this constructor
SubProb::SubProb(vector<Constraints>& constraints, const Objective& obj, vector<Constraints>& couplingCons, const vector<Vars>& vars, bool flag)
	: conss{ nullptr }, numConss{ (int)constraints.size() }, oriObj{ obj },varSize{(int)vars.size()}, consVec_{ constraints }
{
	//@wqy here only init varsVec_,consVec_, c and A.
	//init varsVec_
	//sort the varsVec_ by the vars.id
	std::vector<Vars> VarsVecTemp;
	for (int i = 0; i < vars.size(); i++) {
		VarsVecTemp.push_back(vars[i]);
	}	
	std::vector<bool>selectFlag(vars.size(), false);
	for (int j = 0; j < VarsVecTemp.size(); j++) {
		int record = -1;
		Vars varMinTemp; //init large id
		bool assignFlag(true);
		for (int i = 0; i < VarsVecTemp.size(); i++) {
			if (selectFlag[i]) { continue; }
			if (assignFlag) { 
				varMinTemp = VarsVecTemp[i];
				record = i;
				assignFlag = false;
			}
			if (VarsVecTemp[i].Id < varMinTemp.Id) {
				varMinTemp = VarsVecTemp[i];
				record = i;
			}
		}
		selectFlag[record] = true;
		varsVec_.push_back(varMinTemp);
	}
	//c and A remain unchanged after this
	c = Eigen::VectorXd::Zero(varSize);   // according to the order of vars id
	A = Eigen::MatrixXd::Zero(couplingCons.size(), varSize);  //according to the order of vars id and the order of the constraints

	for (int i = 0; i < varSize; i++) {
		if (obj.coef.find(varsVec_[i].Id) == obj.coef.end()) { spdlog::error("Obj coef not found!"); }
		c(i) = obj.coef.find(varsVec_[i].Id)->second;

		for (int j = 0; j < couplingCons.size(); j++)
		{
			auto subIter = couplingCons[j].exprDic.find(varsVec_[i]);
			if (subIter == couplingCons[j].exprDic.end()) { continue; }
			A(j, i) = subIter->second;
		}
	}
}

void SubProb::BuildSCIPModel()
{
	SCIPcreate(&scipModel);
	SCIPincludeDefaultPlugins(scipModel);

	SCIPcreateProbBasic(scipModel, "subProb");
	//SCIPsetIntParam(scipModel, "presolving/maxrounds", 0);
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
		SCIPcreateConsBasicLinear(scipModel, &scipConss[i], consTemp.Name.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
		SCIPaddCons(scipModel, scipConss[i]);

		SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scipModel), 1);
		/*delete[] tempVars;
		delete[] tempConsCoefs;*/
	}

}

//fix integer variable value based on fixedinformation
void SubProb::BuildSCIPModelNew(FIX fixedInfo, int subprobNo)
{
	SCIPcreate(&scipModel);
	SCIPincludeDefaultPlugins(scipModel);
	SCIPcreateProbBasic(scipModel, "subProb");

	//set min or max
	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
	if (oriObj.direction == Direction::max) {
		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
	}
	SCIPsetObjsense(scipModel, objSense);

	//set variables
	unordered_map<int, SCIP_VAR*> scipVarMap;  //using in set constraints
	for (int i = 0; i < varSize; i++) {
		SCIP_VAR* tempVar;   //need to be added to the scip model
		VarType varType = varsVec_[i].Type;
		SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
		if (varType == VarType::Bool) {
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
		}
		else if (varType == VarType::Int) {
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_INTEGER;
		}
		else if (varType == VarType::Num) {
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
		}

		//@wqy subproblem set continous for test
		//scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;

		//@wqy fix the integer variable of the column
		bool fixedFlag = false;
		int fixedValue = 0;
		if (fixedInfo.spIsFixed[subprobNo]) {
			auto fixCol(fixedInfo.fixedColumnValueList[subprobNo]);
			fixedFlag = true;
			for (const auto& fixColVar : fixCol) {
				if (fixColVar.first == varsVec_[i].Id) {
					fixedValue = fixColVar.second;
					break;
				}
			}
		}

		//if the variable is integer and fixed, change its lower and upper bound to the same value
		//@wqy 0621 obj coef set to 1.0 when building the subproblem model, and it will be set in later function.
		if (fixedFlag && scipVarType != SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS) {
			SCIPcreateVarBasic(scipModel, &tempVar, varsVec_[i].Name.c_str(), fixedValue, fixedValue, 1.0, scipVarType);
		}else {
			SCIPcreateVarBasic(scipModel, &tempVar, varsVec_[i].Name.c_str(), varsVec_[i].Lb, varsVec_[i].Ub, 1.0, scipVarType);
		}
		SCIPaddVar(scipModel, tempVar);

		scipVarMap.insert(make_pair(varsVec_[i].Id, tempVar));
		scipVarsVec_.push_back(tempVar); //@wqy 0621 add

		//delete pointer
		//delete[] tempVar;

	}

	//set constraints
	SCIP_CONS* tempCons; //need to be added to the scip model
	for (int i = 0; i < numConss; i++){
		Constraints consTemp = consVec_[i];
		unordered_map<Vars, double>::iterator iter = consTemp.exprDic.begin();
		int tempVarsSize = consTemp.exprDic.size();
		SCIP_VAR** tempVars = new SCIP_VAR* [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];
		int idx = 0;
		while (iter != consTemp.exprDic.end()) {
			tempVars[idx] = scipVarMap[iter->first.Id];
			tempConsCoefs[idx++] = iter->second;
			iter++;
		}

		double tempLb = consTemp.rhs;
		double tempUb = SCIPinfinity(scipModel);
		if (consTemp.Prop == PROP::eq) {
			tempLb = consTemp.rhs;
			tempUb = consTemp.rhs;
		}
		if (consTemp.Prop == PROP::leq) {
			tempUb = consTemp.rhs;
			tempLb = -SCIPinfinity(scipModel);
		}
		string consName = "cons" + std::to_string(i);
		SCIPcreateConsBasicLinear(scipModel, &tempCons, consName.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
		SCIPaddCons(scipModel, tempCons);
		SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scipModel), 1);
		scipConssVec_.push_back(tempCons);

		//delete pointer (@wqy it maybe wrong)
		delete[] tempVars;
		//delete[] tempCons;
		delete[] tempConsCoefs;
	}

	//SCIPwriteOrigProblem(scipModel, "D://sub1.lp", nullptr, 0);
}

void SubProb::BuildSCIPModelNew(FIX fixedInfo, int subprobNo, const std::vector<double>& duals_part)
{
	SCIPcreate(&scipModel);
	SCIPincludeDefaultPlugins(scipModel);
	SCIPcreateProbBasic(scipModel, "subProb");

	//these two vector will be "push_back"coss later ,and the "BuildSCIPModelNew" function is used more than once.They need to be clear.
	//the other two vector (varsVec_ and consVec_) are initialled in the constructor.
	scipVarsVec_.clear();
	scipConssVec_.clear();

	//calculate the objective coef --c - mu * A
	Eigen::VectorXd mu = Vector2Vec(duals_part);
	Eigen::VectorXd cBar = c - A.transpose() * mu;

	//update tempObjVal with cBar
	vector<double> objCeff;
	for (int i = 0; i < varSize; i++) {
		objCeff.push_back(cBar(i));
        //std::cout<<objCeff[i]<<'\t';
	}
    //std::cout<<std::endl;

	//set min or max
	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
	if (oriObj.direction == Direction::max) {
		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
	}
	SCIPsetObjsense(scipModel, objSense);

	//set variables
	unordered_map<int, SCIP_VAR*> scipVarMap;  //using in set constraints
	for (int i = 0; i < varSize; i++) {
		SCIP_VAR* tempVar;   //need to be added to the scip model
		VarType varType = varsVec_[i].Type;
		SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
		if (varType == VarType::Bool) {
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
		}
		else if (varType == VarType::Int) {
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_INTEGER;
		}
		else if (varType == VarType::Num) {
			scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
		}

		//@wqy subproblem set continous for test
		//scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;

		//@wqy fix the integer variable of the column
		bool fixedFlag = false;
		int fixedValue = 0;
		if (fixedInfo.spIsFixed[subprobNo]) {
			//find the fixed column
			std::vector<std::pair<int, double>> fixCol;
			for (int i = 0; i < fixedInfo.spNoList.size(); i++) {
				if (fixedInfo.spNoList[i] == subprobNo) {
					(fixCol = fixedInfo.fixedColumnValueList[i]);
				}
			}
			fixedFlag = true;
			for (const auto& fixColVar : fixCol) {
				if (fixColVar.first == varsVec_[i].Id) {
					fixedValue = fixColVar.second;
					break;
				}
			}
		}

		//if the variable is integer and fixed, change its lower and upper bound to the same value
	
		if (fixedFlag && scipVarType != SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS) {
			SCIPcreateVarBasic(scipModel, &tempVar, varsVec_[i].Name.c_str(), fixedValue, fixedValue, objCeff[i], scipVarType);
		}
		else {
			SCIPcreateVarBasic(scipModel, &tempVar, varsVec_[i].Name.c_str(), varsVec_[i].Lb, varsVec_[i].Ub, objCeff[i], scipVarType);
		}
		SCIPaddVar(scipModel, tempVar);

		scipVarMap.insert(make_pair(varsVec_[i].Id, tempVar));
		scipVarsVec_.push_back(tempVar); //@wqy 0621 add

		//delete pointer
		//delete[] tempVar;
	}

	//set constraints
	SCIP_CONS* tempCons; //need to be added to the scip model
	for (int i = 0; i < numConss; i++) {
		Constraints consTemp = consVec_[i];
		unordered_map<Vars, double>::iterator iter = consTemp.exprDic.begin();
		int tempVarsSize = consTemp.exprDic.size();
		SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];
		int idx = 0;
		while (iter != consTemp.exprDic.end()) {
			tempVars[idx] = scipVarMap[iter->first.Id];
			tempConsCoefs[idx] = iter->second;
			idx++;
			iter++;
		}

		double tempLb = consTemp.rhs;
		double tempUb = SCIPinfinity(scipModel);
		if (consTemp.Prop == PROP::eq) {
			tempLb = consTemp.rhs;
			tempUb = consTemp.rhs;
		}
		if (consTemp.Prop == PROP::leq) {
			tempUb = consTemp.rhs;
			tempLb = -SCIPinfinity(scipModel);
		}
		string consName = "cons" + std::to_string(i);
		SCIPcreateConsBasicLinear(scipModel, &tempCons, consName.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
		SCIPaddCons(scipModel, tempCons);
		SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scipModel), 1);
		scipConssVec_.push_back(tempCons);

		//delete pointer (@wqy it maybe wrong)
		delete[] tempVars;
		//delete[] tempCons;
		delete[] tempConsCoefs;
	}

	//SCIPwriteOrigProblem(scipModel, "D://sub1.lp", nullptr, 0);
}
void SubProb::UpdateTempObj(Solution dualSol)
{
	//to calculate: c - mu * A
	Eigen::VectorXd mu = Sol2Vec(dualSol);
	Eigen::VectorXd cBar = c - A.transpose() * mu;
	//update tempObjVal with cBar
	int i = 0;
	for (unordered_map<int, double>::iterator iter = tempObjVal.begin(); iter != tempObjVal.end(); iter++)
	{
		iter->second = cBar(i++);
	}
	updateSCIPObj();
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

void SubProb::UpdateTempObj(const Eigen::VectorXd& dualSol)
{
	//to calculate: c - mu * A
	Eigen::VectorXd cBar = c - A.transpose() * dualSol;
	//update tempObjVal with cBar
	//int i = 0;
	/*for (unordered_map<int, double>::iterator iter = tempObjVal.begin(); iter != tempObjVal.end(); iter++)
	{
		iter->second = cBar(i);
		i++;
	}*/
	for (int ii = 0; ii < subOrder2Id.size(); ii++)
	{
		tempObjVal[subOrder2Id[ii]] = cBar(ii);
		//SCIPchgVarObj(scipModel, scipVars[ii], cBar(ii));
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

void SubProb::updateSCIPObj(vector<double>objCeff)
{
	for (int i = 0; i < varSize; i++) {
		SCIPchgVarObj(scipModel, scipVarsVec_[i], objCeff[i]);
	}
}

void SubProb::UpdateTempObjNew(const vector<double>& dualSol)
{
	//SCIPwriteOrigProblem(scipModel, "D://sub2.lp", NULL, false);

	//to calculate: c - mu * A
	Eigen::VectorXd mu = Vector2Vec(dualSol);
	Eigen::VectorXd cBar = c - A.transpose() * mu;

	//update tempObjVal with cBar
	vector<double> objCeff;
	for (int i = 0; i < varSize; i++) {
		objCeff.push_back(cBar(i));
	}

	updateSCIPObj(objCeff);

	//SCIPwriteOrigProblem(scipModel, "D://sub3.lp", NULL, false);
}

