#include "utility.h"
#include"DW.h"



//Comparing with 2, the function using scip get many feasible solution. It may not get enough feasible solutions.
int GenerateInitColumns3(std::vector<ColumnPool>& columnpoolInput, std::map<int, General_Block> block, int genColNum, SCIP* scip, double& initScipFeasObj)
{
	int spNum = block.size() - 1;


	//0605 Here a feasible solution of original problem is obtained and added to the columnpool.
	//SCIPsetRealParam(scip, "limits/time", 10);
    SCIPsetIntParam(scip, "display/verblevel", 0);
	SCIPsetIntParam(scip, "limits/solutions", 1);
	SCIPsolve(scip);
	initScipFeasObj = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
	cout << "The feasible solution obj of primal problem is: " << initScipFeasObj << endl;

	vector<pair<string, double>>colNameString;
	SCIP_VAR** varsOri = SCIPgetVars(scip);
	int nVars = SCIPgetNVars(scip);

	for (int varNo = 0; varNo < nVars; varNo++){
		double val = SCIPgetSolVal(scip, SCIPgetBestSol(scip), varsOri[varNo]);
		string name = SCIPvarGetName(varsOri[varNo]);
		colNameString.push_back(make_pair(name, val));
	}

	//maybe it can be deleted later, because it just a check
	if (!DW::Check_All_Cons(block, colNameString)) {
		cout << "The solution found for primal problem is infeasible!" << endl;
	}

	for (int subproblem_no = 0; subproblem_no < spNum; subproblem_no++){
		//assign local variable
		int varsSize = block[subproblem_no].bVars.size();
		int consSize = block[subproblem_no].bCons.size();
		ColumnPool pool4NowSub;
		SCIP* scipModel;
		
		//build the model
		SCIPcreate(&scipModel);
		SCIPincludeDefaultPlugins(scipModel);
		SCIPcreateProbBasic(scipModel, "InitsubProb");
		//SCIPsetIntParam(scipModel, "limits/solutions", 10);
        SCIPsetIntParam(scipModel, "display/verblevel", 0);
		SCIPsetRealParam(scipModel, "limits/time", 30);         // set 30s for scip solve feasible for subproblem
		SCIPsetIntParam(scipModel, "presolving/maxrounds", 0);
		unordered_map<int, SCIP_VAR*> scipVarMap;
		unordered_map<int, Vars> varDic;
		vector<SCIP_VAR*> scipVars;
		SCIP_CONS** scipConss;

		//set dirction
		SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
		if (block[subproblem_no].bobj.direction == Direction::max){
			objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
		}
		SCIPsetObjsense(scipModel, objSense);

		//add variable
		for (int i = 0; i < varsSize; i++){
			Vars varNow = block[subproblem_no].bVars[i];
			pool4NowSub.vars_id.push_back(varNow.Id);
			SCIP_VAR* tempVar;
			VarType varType = varNow.Type;
			SCIP_VARTYPE scipVarType;
			if (varType == VarType::Bool){
				scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
			}else if (varType == VarType::Int){
				scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_INTEGER;
			}else if (varType == VarType::Num){
				scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
			}

			string varName = varNow.Name;
			double lbDeal = varNow.Lb;
			double ubDeal = varNow.Ub;
			if (varNow.Lb < -10000000){
				if (varNow.Ub < -1000){
					lbDeal = varNow.Ub - 1000;
				}else{
					lbDeal = -1000;
				}
			}
			if (varNow.Ub > 1000000){
				if (varNow.Lb > 1000){
					ubDeal = varNow.Lb + 1000;
				}else{
					ubDeal = 1000;
				}
			}
			SCIPcreateVarBasic(scipModel, &tempVar, varName.c_str(), lbDeal, ubDeal, block[subproblem_no].bobj.coef[varNow.Id], scipVarType);
			SCIPaddVar(scipModel, tempVar);
			scipVarMap.insert(make_pair(varNow.Id, tempVar));
			varDic.insert(pair<int, Vars>(varNow.Id, varNow));
			scipVars.push_back(tempVar);
		}

		//add constraints
		scipConss = new SCIP_CONS* [consSize];
		for (int i = 0; i < consSize; i++){
			Constraints consTemp = block[subproblem_no].bCons[i];
			unordered_map<Vars, double>::iterator iter = consTemp.exprDic.begin();
			int tempVarsSize = consTemp.exprDic.size();
			SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
			double* tempConsCoefs = new double[tempVarsSize];
			int idx = 0;
			while (iter != consTemp.exprDic.end()){
				tempVars[idx] = scipVarMap[iter->first.Id];
				tempConsCoefs[idx++] = iter->second;
				iter++;
			}

			double tempLb = consTemp.rhs; 
			double tempUb = SCIPinfinity(scipModel);
			if (consTemp.Prop == PROP::eq){
				tempLb = consTemp.rhs;
				tempUb = consTemp.rhs;
			}else if (consTemp.Prop == PROP::leq){
				tempUb = consTemp.rhs;
				tempLb = -SCIPinfinity(scipModel);
			}
			string consName = "cons" + std::to_string(i);
			SCIPcreateConsBasicLinear(scipModel, &scipConss[i], consTemp.Name.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
			SCIPaddCons(scipModel, scipConss[i]);

			delete[] tempVars;
			delete[] tempConsCoefs;
		}

		//SCIPwriteOrigProblem(scipModel, "C://gen_col.lp", nullptr, 0);

		SCIPsolve(scipModel);

		//get more than one feasible solutions
		int nSols = SCIPgetNSols(scipModel);
		SCIP_SOL** sols = SCIPgetSols(scipModel);

		for (int solNo = 0; solNo < genColNum && solNo < nSols; solNo++){
			Column columnTemp;
			int nvars = SCIPgetNVars(scipModel);
			for (int j = 0; j < nvars; j++){
				double val = SCIPgetSolVal(scipModel, sols[solNo], scipVars[j]);
				string name = SCIPvarGetName(scipVars[j]);
				val = chgInt(val); //0523@wqy modify
				bool isFlag(false);
				for (const auto& varVar : varDic){
					if (name == varVar.second.Name || name == "t_" + varVar.second.Name) {
						columnTemp.varValue.insert(pair<int, double>(varVar.second.Id, val));
						isFlag = 1;
						break;
					}
				}
				if (!isFlag) { 
					cout << "Error1! varValue of column is missing!" << endl;  
					return 0;
				}
			}
			columnTemp.solSize = nvars;
			pool4NowSub.PushBackColumn(columnTemp);
		}
		pool4NowSub.var_num = varsSize;

		//release
		for (int i = 0; i < varDic.size(); i++){
			SCIPreleaseVar(scipModel, &scipVars[i]);
		}
		for (int i = 0; i < consSize; i++){
			SCIPreleaseCons(scipModel, &scipConss[i]);
		}
		SCIPfree(&scipModel);

		//@wqy0605 Split the solution and add it to the columnpool of each subproblem
		Column columnCorrect;
		int nvars = SCIPgetNVars(scip);
		int countFlag = 0;
		for (int j = 0; j < nvars; j++) {
			double val = SCIPgetSolVal(scip, SCIPgetBestSol(scip), varsOri[j]);
			string name = SCIPvarGetName(varsOri[j]);
			val = chgInt(val); //0523@wqy modify			
			for (const auto& varVar : varDic) {
				if (name == varVar.second.Name || name == "t_" + varVar.second.Name) {
					columnCorrect.varValue.insert(pair<int, double>(varVar.second.Id, val));
					countFlag++;
					break;
				}
			}
		}
		if (countFlag != varsSize) {
			cout << "Error2! varValue of column is missing!" << endl;
			return 0;
		}
		columnCorrect.solSize = varsSize;
		pool4NowSub.PushBackColumn(columnCorrect);
		columnpoolInput.push_back(pool4NowSub);
	}
	
 	return 1;
}