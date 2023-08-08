#include"MasterProb.h"
#include"heuristicTools.h"
#include "ilcplex/ilocplex.h"  @wqy coplex1
//#include<cstdlib>
#include <iostream>

//using namespace std;



void MasterProb::free()
{
	/*
		for (int i = 0; i < RMPscipVars_num; i++)
		{
			if (SCIPvarGetName(RMPscipVars[i]) != NULL  && RMPscipVars[i] != NULL)
			{
				SCIPreleaseVar(RMPscipModel, &RMPscipVars[i]);
			}
		}
		RMPscipVars_num = 0;
		for (int i = 0; i < RMPscipConss_num; i++)
		{
			if (RMPscipConss[i] != NULL &&  SCIPvarGetName(RMPscipConss[i]) != NULL)
			{
				SCIPreleaseCons(RMPscipModel, &RMPscipConss[i]);
			}
		}

		SCIPfree(&RMPscipModel);*/
}

//deep copy
MasterProb& MasterProb::operator += (const MasterProb& mp)
{
	if (this == &mp)
		return *this;

	this->A = mp.A;

	this->columnpool = mp.columnpool;

	this->col_num = mp.col_num;
	this->ConsCoefMatrix = mp.ConsCoefMatrix;
	this->RMPscipConss = new SCIP_CONS * [mp.RMPscipConss_num];
	for (int i = 0; i < mp.RMPscipConss_num; i++)
	{
		this->RMPscipConss[i] = mp.RMPscipConss[i];
	}
	this->RMPscipConss_num = mp.RMPscipConss_num;
	this->RMPscipModel = new SCIP();
	this->RMPscipModel = mp.RMPscipModel;
	this->RMPscipVars = new SCIP_VAR * [mp.RMPscipVars_num];
	for (int i = 0; i < mp.RMPscipVars_num; i++)
	{
		this->RMPscipVars[i] = mp.RMPscipVars[i];
	}
	this->RMPscipVars_num = mp.RMPscipVars_num;
	this->sp_num = mp.sp_num;


	return *this;

}


//void MasterProb::BuildSCIPModelBinary(Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX> Fixed_Info)
//{
//	SCIPcreate(&RMPscipModel);
//	SCIPincludeDefaultPlugins(RMPscipModel);
//	SCIPsetIntParam(RMPscipModel, "presolving/maxrounds", 0);
//
//	SCIPcreateProbBasic(RMPscipModel, "R_MasterProb");
//
//	//set max or min
//	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
//	if (obj.direction == Direction::max)
//	{
//		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
//	}
//	SCIPsetObjsense(RMPscipModel, objSense);
//
//	int masterColNum = 0;
//	for (int i = 0; i < subprob_num; i++)
//	{
//		masterColNum += columnpool[i].columns.size();
//	}
//
//	//add variable
//	unordered_map<int, SCIP_VAR*> RMPscipVarMap;
//	RMPscipVars_num = masterColNum;
//	RMPscipVars = new SCIP_VAR * [RMPscipVars_num];
//	int scipVarsIdx = 0;
//	int col_num_all = 0;//the number of all columns
//	bool cout_flag = false;
//
//	ConsCoefMatrix.resize(0, 0);
//	//cout << "111ConsCoefMatrix" << endl;
//	//cout << ConsCoefMatrix;
//	//cout << endl;
//
//	for (int i = 0; i < subprob_num; i++)
//	{
//		//assign column matrix
//		Eigen::MatrixXd colmatrix;
//		colmatrix = Eigen::MatrixXd::Zero(columnpool[i].var_num, columnpool[i].columns.size());
//		int column_idx = 0;
//		//sort vars_id
//		sort(columnpool[i].vars_id.begin(), columnpool[i].vars_id.end());
//
//		for (auto iter = columnpool[i].columns.begin(); iter != columnpool[i].columns.end(); ++iter)
//		{
//			for (int j = 0; j < columnpool[i].var_num; j++)
//			{
//				colmatrix(j, column_idx) = iter->varValue.find(columnpool[i].vars_id[j])->second;
//			}
//			column_idx++;
//		}
//
//		//Take the submatrix of the constraint matrix A of the original problem corresponding to this subproblem
//		Eigen::MatrixXd A_part;
//		A_part = Eigen::MatrixXd::Zero(couplingConsNum, columnpool[i].var_num);
//		for (int t = 0; t < couplingConsNum; t++)
//		{
//			for (int j = 0; j < columnpool[i].var_num; j++)
//			{
//				if (couplingCons[t].idDic.find(columnpool[i].vars_id[j]) == couplingCons[t].idDic.end())
//				{
//					A_part(t, j) = 0;
//				}
//				else
//				{
//					A_part(t, j) = couplingCons[t].idDic.find(columnpool[i].vars_id[j])->second;
//				}
//			}
//		}
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "A_part" << endl;
//			cout << A_part;
//			cout << endl;
//			cout << "colmatrix" << endl;
//			cout << colmatrix;
//			cout << endl;
//		}
//
//
//		//matrix multiple, calculate the constraint coefficient of the master problem
//		Eigen::MatrixXd RMPcoef_matrix = A_part * colmatrix;
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "RMPcoef_matrix" << endl;
//			cout << RMPcoef_matrix;
//			cout << endl;
//		}
//
//
//		Eigen::MatrixXd Coef_temp(RMPcoef_matrix.rows(), ConsCoefMatrix.cols() + RMPcoef_matrix.cols());
//
//		Coef_temp << ConsCoefMatrix, RMPcoef_matrix;
//
//		ConsCoefMatrix = Coef_temp;
//
//		//Take the objective function coefficients of the original problem corresponding to the part of this subproblem
//		Eigen::VectorXd obj_part;
//		obj_part = Eigen::VectorXd::Zero(columnpool[i].var_num);
//		for (int t = 0; t < columnpool[i].var_num; t++)
//		{
//			obj_part(t) = obj.coef[columnpool[i].vars_id[t]];
//		}
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "obj_part" << endl;
//			cout << obj_part;
//			cout << endl;
//
//
//			cout << "obj_part.transpose()" << endl;
//			cout << obj_part.transpose();
//			cout << endl;
//
//
//			cout << "colmatrix" << endl;
//			cout << colmatrix;
//			cout << endl;
//
//			cout << "ConsCoefMatrix:" << endl;
//			cout << ConsCoefMatrix;
//			cout << endl;
//		}
//
//
//		//calculate the objective function coefficients of the master problem
//		Eigen::VectorXd RMPobjcoef_vector = obj_part.transpose() * colmatrix;
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "RMPobjcoef_vector" << endl;
//			cout << RMPobjcoef_vector;
//			cout << endl;
//		}
//
//
//
//		//add SCIP variable
//		for (int j = 0; j < columnpool[i].columns.size(); j++)
//		{
//			SCIP_VAR* tempVar;
//			//set variable type
//			SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY; //@wqy0608 Only here is changed compared to BuildModel function
//			const char* var_name;
//			string varname = "_lambda_" + to_string(i) + "_" + to_string(j);
//			var_name = varname.c_str();
//			SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, 1, RMPobjcoef_vector(j), scipVarType);
//			//SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, SCIPinfinity(RMPscipModel), RMPobjcoef_vector(j), scipVarType); //@wqytest0527
//			//SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, -SCIPinfinity(RMPscipModel), SCIPinfinity(RMPscipModel), RMPobjcoef_vector(j), scipVarType);//@wqy0515test
//			SCIPaddVar(RMPscipModel, tempVar);
//
//
//			RMPscipVars[scipVarsIdx++] = tempVar;
//			RMPscipVarMap.insert(make_pair(col_num_all, tempVar));
//			col_num_all++;
//		}
//
//	}
//
//	//add SCIP constraints
//	RMPscipConss_num = couplingConsNum + subprob_num;
//	RMPscipConss = new SCIP_CONS * [RMPscipConss_num];
//
//	//cout << "ConsCoefMatrix:" << endl;
//	//cout << ConsCoefMatrix;
//	//cout << endl;
//
//	for (int i = 0; i < couplingConsNum; i++)
//	{
//		Constraints consTemp = couplingCons[i];
//
//		int tempVarsSize = col_num_all; //number of master problem variables
//		//SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
//		double* tempConsCoefs = new double[tempVarsSize];
//		int idx = 0;
//		while (idx != col_num_all)
//		{
//			tempConsCoefs[idx] = ConsCoefMatrix(i, idx);
//			idx++;
//		}
//
//		//cout
//		//cout << "constraint " << i << " coef" << endl;
//		//for (int kk = 0; kk < tempVarsSize; kk++)
//		//{
//		//	cout << tempConsCoefs[kk] << " ";
//		//}
//		//cout << endl;
//
//		const char* con_name;
//		string conname = "z_cons" + to_string(i);
//		con_name = conname.c_str();
//		double rhs, lhs;
//		if (couplingCons[i].Prop == 0)
//		{
//			lhs = couplingCons[i].rhs;
//			rhs = SCIPinfinity(RMPscipModel);
//		}
//		else if (couplingCons[i].Prop == 1)
//		{
//			lhs = -SCIPinfinity(RMPscipModel);
//			rhs = couplingCons[i].rhs;
//		}
//		else
//		{
//			lhs = couplingCons[i].rhs;
//			rhs = couplingCons[i].rhs;
//		}
//
//		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], con_name,
//			tempVarsSize, RMPscipVars, tempConsCoefs, lhs, rhs);
//		SCIPaddCons(RMPscipModel, RMPscipConss[i]);
//	}
//	int record_idx = 0;
//	for (int i = couplingConsNum; i < couplingConsNum + subprob_num; i++)
//	{
//		int begin_idx, end_idx;//define start and end index
//
//		int tempVarsSize = columnpool[i - couplingConsNum].columns.size();
//		SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
//		double* tempConsCoefs = new double[tempVarsSize];
//
//		begin_idx = record_idx;
//		end_idx = begin_idx + tempVarsSize;
//		record_idx = end_idx;
//		int idx = 0;
//		for (auto iter = RMPscipVarMap.begin(); iter != RMPscipVarMap.end(); ++iter)
//		{
//			if (iter->first >= begin_idx && iter->first < end_idx)
//			{
//				tempVars[idx] = iter->second;
//				//cout << iter->second->name << endl;;
//				tempConsCoefs[idx] = 1;
//				idx++;
//			}
//		}
//
//		const char* con_name;
//		string conname = "z1_cons" + to_string(i);
//		con_name = conname.c_str();
//		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], con_name,
//			tempVarsSize, tempVars, tempConsCoefs, 1, 1);
//		SCIPaddCons(RMPscipModel, RMPscipConss[i]);
//
//
//	}
//
//
//	SCIPwriteOrigProblem(RMPscipModel, "D://RMP_init.lp", nullptr, 0);
//	SCIPsolve(RMPscipModel);
//
//	if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_OPTIMAL) {
//		cout << "Obj:" << SCIPgetSolOrigObj(RMPscipModel, SCIPgetBestSol(RMPscipModel)) << endl;
//	}
//	else if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_INFEASIBLE) {
//		cout << "Binary RMP is infeasible!" << endl;
//	}
//	else {
//		cout << "Binary RMP is solved abnormally!" << endl;
//	}
//
//
//
//}

//@wqy 0625save
//void MasterProb::BuildSCIPModel(Objective obj, vector<Constraints> couplingCons, int couplingConsNum, int subprobNum, FIX fixedInfo)
//{
//	SCIPcreate(&RMPscipModel);
//	SCIPincludeDefaultPlugins(RMPscipModel);
//	SCIPsetIntParam(RMPscipModel, "presolving/maxrounds", 0);
//	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
//	SCIPcreateProbBasic(RMPscipModel, "R_MasterProb");
//
//	//set max or min
//	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
//	if (obj.direction == Direction::max)
//	{
//		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
//	}
//	SCIPsetObjsense(RMPscipModel, objSense);
//
//	int masterColNum = 0;
//	for (int i = 0; i < subprobNum; i++){
//		masterColNum += columnpool[i].Columns().size();
//	}
//
//	//add variable
//	unordered_map<int, SCIP_VAR*> RMPscipVarMap;
//	RMPscipVars_num = masterColNum;  //the number of variable which equals to the number of column in columnpool
//	RMPscipVars = new SCIP_VAR * [RMPscipVars_num];
//	int scipVarsIdx = 0;
//	int col_num_all = 0;//the number of all columns
//	bool cout_flag = 0;
//
//	ConsCoefMatrix.resize(0, 0);
//
//	for (int i = 0; i < subprobNum; i++){
//		//assign column matrix
//		Eigen::MatrixXd colmatrix;
//		colmatrix = Eigen::MatrixXd::Zero(columnpool[i].var_num, columnpool[i].Columns().size());
//		int column_idx = 0;
//
//		//@wqy Before make the colmatrix, it will sort the value of the column by value's id here, so the value in the column can be unordered.
//		//sort vars_id
//		std::sort(columnpool[i].vars_id.begin(), columnpool[i].vars_id.end());
//
//		for (const auto& colIter : columnpool[i].columns){
//			for (int j = 0; j < columnpool[i].var_num; j++){
//				if (colIter.varValue.find(columnpool[i].vars_id[j]) == colIter.varValue.end()){
//					spdlog::error("There is some wrong in assigning the column matrix!");
//				}
//				colmatrix(j, column_idx) = colIter.varValue.find(columnpool[i].vars_id[j])->second;
//			}
//			column_idx++;
//		}
//
//		//Take the submatrix from the constraint matrix A of the original problem corresponding to this subproblem
//		Eigen::MatrixXd A_part;
//		A_part = Eigen::MatrixXd::Zero(couplingConsNum, columnpool[i].var_num);
//		for (int t = 0; t < couplingConsNum; t++){
//			for (int j = 0; j < columnpool[i].var_num; j++){
//				if (couplingCons[t].idDic.find(columnpool[i].vars_id[j]) == couplingCons[t].idDic.end()){
//					A_part(t, j) = 0;
//				}else{
//					A_part(t, j) = couplingCons[t].idDic.find(columnpool[i].vars_id[j])->second;
//				}
//			}
//		}
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "A_part" << endl;
//			cout << A_part;
//			cout << endl;
//			cout << "colmatrix" << endl;
//			cout << colmatrix;
//			cout << endl;
//		}
//
//		//matrix multiple, calculate the constraint coefficient of the master problem
//		Eigen::MatrixXd RMPcoef_matrix = A_part * colmatrix;
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "RMPcoef_matrix" << endl;
//			cout << RMPcoef_matrix;
//			cout << endl;
//		}
//
//		// @wqy modify for gcc
//		if (i == 0) {
//			ConsCoefMatrix = RMPcoef_matrix;
//		}else {
//			Eigen::MatrixXd Coef_temp(RMPcoef_matrix.rows(), ConsCoefMatrix.cols() + RMPcoef_matrix.cols());
//			Coef_temp << ConsCoefMatrix, RMPcoef_matrix;
//			ConsCoefMatrix = Coef_temp;
//		}
//
//		//Take the objective function coefficients of the original problem corresponding to the part of this subproblem
//		Eigen::VectorXd obj_part;
//		obj_part = Eigen::VectorXd::Zero(columnpool[i].var_num);
//		for (int t = 0; t < columnpool[i].var_num; t++){
//			obj_part(t) = obj.coef[columnpool[i].vars_id[t]];
//		}
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "obj_part" << endl;
//			cout << obj_part;
//			cout << endl;
//
//
//			cout << "obj_part.transpose()" << endl;
//			cout << obj_part.transpose();
//			cout << endl;
//
//
//			cout << "colmatrix" << endl;
//			cout << colmatrix;
//			cout << endl;
//
//			cout << "ConsCoefMatrix:" << endl;
//			cout << ConsCoefMatrix;
//			cout << endl;
//		}
//
//
//		//calculate the objective function coefficients of the master problem
//		Eigen::VectorXd RMPobjcoef_vector = obj_part.transpose() * colmatrix;
//
//		//cout
//		if (cout_flag)
//		{
//			cout << "RMPobjcoef_vector" << endl;
//			cout << RMPobjcoef_vector;
//			cout << endl;
//		}
//
//		//add SCIP variable
//		for (int j = 0; j < columnpool[i].columns.size(); j++){
//			SCIP_VAR* tempVar;
//			//set variable type
//			SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
//			const char* var_name;
//			string varname = "_lambda_" + to_string(i) + "_" + to_string(j);
//			var_name = varname.c_str();
//			//SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, 1, RMPobjcoef_vector(j), scipVarType);
//			SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, SCIPinfinity(RMPscipModel), RMPobjcoef_vector(j), scipVarType); //@wqytest0527 (0,+inf)
//			SCIPaddVar(RMPscipModel, tempVar);
//
//			RMPscipVars[scipVarsIdx++] = tempVar;
//			RMPscipVarMap.insert(make_pair(col_num_all, tempVar));
//			col_num_all++;
//		}
//
//	}
//
//	//add SCIP constraints
//	RMPscipConss_num = couplingConsNum + subprobNum;
//	RMPscipConss = new SCIP_CONS* [RMPscipConss_num];
//
//	for (int i = 0; i < couplingConsNum; i++){
//		Constraints consTemp = couplingCons[i];
//
//		int tempVarsSize = col_num_all; //number of master problem variables
//		double* tempConsCoefs = new double[tempVarsSize];
//		int idx = 0;
//		while (idx != col_num_all){
//			tempConsCoefs[idx] = ConsCoefMatrix(i, idx);
//			idx++;
//		}
//
//		string conName = "z_cons" + to_string(i);
//		double rhs, lhs;
//		if (couplingCons[i].Prop == 0){
//			lhs = couplingCons[i].rhs;
//			rhs = SCIPinfinity(RMPscipModel);
//		}else if (couplingCons[i].Prop == 1){
//			lhs = -SCIPinfinity(RMPscipModel);
//			rhs = couplingCons[i].rhs;
//		}else{
//			lhs = couplingCons[i].rhs;
//			rhs = couplingCons[i].rhs;
//		}
//
//		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], conName.c_str(),
//			tempVarsSize, RMPscipVars, tempConsCoefs, lhs, rhs);
//		//@wqy test 0620
//		if (SCIPaddCons(RMPscipModel, RMPscipConss[i]) != SCIP_OKAY) {
//			spdlog::error("There is some wrong in add constraint to Master SCIP Model!");
//		}
//		
//	}
//	int record_idx = 0;
//	for (int i = couplingConsNum; i < couplingConsNum + subprobNum; i++)
//	{
//		int begin_idx;
//		int end_idx;//define start and end index
//
//		int tempVarsSize = columnpool[i - couplingConsNum].columns.size();
//		SCIP_VAR** tempVars = new SCIP_VAR* [tempVarsSize];
//		double* tempConsCoefs = new double[tempVarsSize];
//
//		begin_idx = record_idx;
//		end_idx = begin_idx + tempVarsSize;
//		record_idx = end_idx;
//
//		for (int idx = begin_idx; idx < end_idx; idx++) {
//			if (RMPscipVarMap.find(idx) == RMPscipVarMap.end()) {
//				spdlog::error("There is some wrong in build constraints of master problem!");
//			}else {
//				tempVars[idx - begin_idx] = RMPscipVarMap.find(idx)->second;
//				tempConsCoefs[idx - begin_idx] = 1;
//			}
//		}
//
//		string conName = "z1_cons" + to_string(i);
//		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], conName.c_str(),
//			tempVarsSize, tempVars, tempConsCoefs, 1, 1);
//		SCIPaddCons(RMPscipModel, RMPscipConss[i]);
//	}
//	SCIPwriteOrigProblem(RMPscipModel, "D://RMP_init.lp", nullptr, 0);
//}

void MasterProb::BuildSCIPModel(Objective obj, vector<Constraints> couplingCons, int couplingConsNum, int subprobNum, FIX fixedInfo)
{
	SCIPcreate(&RMPscipModel);
	SCIPincludeDefaultPlugins(RMPscipModel);
	SCIPsetIntParam(RMPscipModel, "presolving/maxrounds", 0);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPcreateProbBasic(RMPscipModel, "R_MasterProb");

	//set max or min
	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
	if (obj.direction == Direction::max)
	{
		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
	}
	SCIPsetObjsense(RMPscipModel, objSense);

	int masterColNum = 0;
	for (int i = 0; i < subprobNum; i++) {
		masterColNum += columnpool[i].Columns().size();
	}

	//add variable
	unordered_map<int, SCIP_VAR*> RMPscipVarMap;
	RMPscipVars_num = masterColNum;  //the number of variable which equals to the number of column in columnpool
	RMPscipVars = new SCIP_VAR * [RMPscipVars_num];
	int scipVarsIdx = 0;
	int col_num_all = 0;//the number of all columns
	bool cout_flag = 0;

	ConsCoefMatrix.resize(0, 0);

	for (int i = 0; i < subprobNum; i++) {
		//assign column matrix
		Eigen::MatrixXd colmatrix;
		colmatrix = Eigen::MatrixXd::Zero(columnpool[i].var_num, columnpool[i].Columns().size());

		//@wqy Before make the colmatrix, it will sort the value of the column by value's id here, so the value in the column can be unordered.
		//sort vars_id
		std::sort(columnpool[i].vars_id.begin(), columnpool[i].vars_id.end());

		for (int colNo = 0; colNo < columnpool[i].Columns().size();colNo++) {
			for (int j = 0; j < columnpool[i].var_num; j++) {
				/*if (columnpool[i].Columns()[colNo].varValue.find(columnpool[i].vars_id[j]) == columnpool[i].Columns()[colNo].varValue.end()) {
					spdlog::error("There is some wrong in assigning the column matrix!");
				}*/			
				//colmatrix(j, colNo) = columnpool[i].FindColumnValue(colNo, columnpool[i].vars_id[j]); 
				//colmatrix(j, colNo) = columnpool[i].Columns()[colNo].varValue.find(columnpool[i].vars_id[j])->second;
				//colmatrix(j, colNo) = columnpool[i].Columns()[colNo].varValue[columnpool[i].vars_id[j]];
				colmatrix(j, colNo) = columnpool[i].columns_[colNo].varValue[columnpool[i].vars_id[j]];
			}
		}

		//Take the submatrix from the constraint matrix A of the original problem corresponding to this subproblem
		Eigen::MatrixXd A_part;
		A_part = Eigen::MatrixXd::Zero(couplingConsNum, columnpool[i].var_num);
		for (int t = 0; t < couplingConsNum; t++) {
			for (int j = 0; j < columnpool[i].var_num; j++) {
				if (couplingCons[t].idDic.find(columnpool[i].vars_id[j]) == couplingCons[t].idDic.end()) {
					A_part(t, j) = 0;
				}
				else {
					A_part(t, j) = couplingCons[t].idDic.find(columnpool[i].vars_id[j])->second;
				}
			}
		}

		//cout
		if (cout_flag)
		{
			cout << "A_part" << endl;
			cout << A_part;
			cout << endl;
			cout << "colmatrix" << endl;
			cout << colmatrix;
			cout << endl;
		}

		//matrix multiple, calculate the constraint coefficient of the master problem
		Eigen::MatrixXd RMPcoef_matrix = A_part * colmatrix;

		//cout
		if (cout_flag)
		{
			cout << "RMPcoef_matrix" << endl;
			cout << RMPcoef_matrix;
			cout << endl;
		}

		// @wqy modify for gcc
		if (i == 0) {
			ConsCoefMatrix = RMPcoef_matrix;
		}
		else {
			Eigen::MatrixXd Coef_temp(RMPcoef_matrix.rows(), ConsCoefMatrix.cols() + RMPcoef_matrix.cols());
			Coef_temp << ConsCoefMatrix, RMPcoef_matrix;
			ConsCoefMatrix = Coef_temp;
		}

		//Take the objective function coefficients of the original problem corresponding to the part of this subproblem
		Eigen::VectorXd obj_part;
		obj_part = Eigen::VectorXd::Zero(columnpool[i].var_num);
		for (int t = 0; t < columnpool[i].var_num; t++) {
			obj_part(t) = obj.coef[columnpool[i].vars_id[t]];
		}

		//cout
		if (cout_flag)
		{
			cout << "obj_part" << endl;
			cout << obj_part;
			cout << endl;


			cout << "obj_part.transpose()" << endl;
			cout << obj_part.transpose();
			cout << endl;


			cout << "colmatrix" << endl;
			cout << colmatrix;
			cout << endl;

			cout << "ConsCoefMatrix:" << endl;
			cout << ConsCoefMatrix;
			cout << endl;
		}


		//calculate the objective function coefficients of the master problem
		Eigen::VectorXd RMPobjcoef_vector = obj_part.transpose() * colmatrix;

		//cout
		if (cout_flag)
		{
			cout << "RMPobjcoef_vector" << endl;
			cout << RMPobjcoef_vector;
			cout << endl;
		}

		//add SCIP variable
		for (int j = 0; j < columnpool[i].Columns().size(); j++) {
			SCIP_VAR* tempVar;
			//set variable type
			SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
			const char* var_name;
			string varname = "_lambda_" + to_string(i) + "_" + to_string(j);
			var_name = varname.c_str();
			//SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, 1, RMPobjcoef_vector(j), scipVarType);
			SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, SCIPinfinity(RMPscipModel), RMPobjcoef_vector(j), scipVarType); //@wqytest0527 (0,+inf)
			SCIPaddVar(RMPscipModel, tempVar);

			RMPscipVars[scipVarsIdx++] = tempVar;
			RMPscipVarMap.insert(make_pair(col_num_all, tempVar));
			col_num_all++;
		}

	}

	//add SCIP constraints
	RMPscipConss_num = couplingConsNum + subprobNum;
	RMPscipConss = new SCIP_CONS * [RMPscipConss_num];

	for (int i = 0; i < couplingConsNum; i++) {
		Constraints consTemp = couplingCons[i];

		int tempVarsSize = col_num_all; //number of master problem variables
		double* tempConsCoefs = new double[tempVarsSize];
		int idx = 0;
		while (idx != col_num_all) {
			tempConsCoefs[idx] = ConsCoefMatrix(i, idx);
			idx++;
		}

		string conName = "z_cons" + to_string(i);
		double rhs, lhs;
		if (couplingCons[i].Prop == PROP::geq) {
			lhs = couplingCons[i].rhs;
			rhs = SCIPinfinity(RMPscipModel);
		}
		else if (couplingCons[i].Prop == PROP::leq) {
			lhs = -SCIPinfinity(RMPscipModel);
			rhs = couplingCons[i].rhs;
		}
		else {
			lhs = couplingCons[i].rhs;
			rhs = couplingCons[i].rhs;
		}

		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], conName.c_str(),
			tempVarsSize, RMPscipVars, tempConsCoefs, lhs, rhs);
		//@wqy test 0620
		if (SCIPaddCons(RMPscipModel, RMPscipConss[i]) != SCIP_OKAY) {
			spdlog::error("There is some wrong in add constraint to Master SCIP Model!");
		}

	}
	int record_idx = 0;
	for (int i = couplingConsNum; i < couplingConsNum + subprobNum; i++)
	{
		int begin_idx;
		int end_idx;//define start and end index

		int tempVarsSize = columnpool[i - couplingConsNum].Columns().size();
		SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];

		begin_idx = record_idx;
		end_idx = begin_idx + tempVarsSize;
		record_idx = end_idx;

		for (int idx = begin_idx; idx < end_idx; idx++) {
			if (RMPscipVarMap.find(idx) == RMPscipVarMap.end()) {
				spdlog::error("There is some wrong in build constraints of master problem!");
			}
			else {
				tempVars[idx - begin_idx] = RMPscipVarMap.find(idx)->second;
				tempConsCoefs[idx - begin_idx] = 1;
			}
		}

		string conName = "z1_cons" + to_string(i);
		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], conName.c_str(),
			tempVarsSize, tempVars, tempConsCoefs, 1, 1);
		SCIPaddCons(RMPscipModel, RMPscipConss[i]);
	}
	//SCIPwriteOrigProblem(RMPscipModel, "D://RMP_init.lp", nullptr, 0);
}


int MasterProb::SolveDual(vector<double>& duals)
{
	//SCIPwriteOrigProblem(RMPscipModel, "D://testRMP.lp", nullptr, 0);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPsolve(RMPscipModel);

	if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_OPTIMAL){
		for (int i = 0; i < RMPscipConss_num; i++){
			unsigned mark;
			double var;
			SCIPgetDualSolVal(RMPscipModel, RMPscipConss[i], &var, &mark);
			//@wqy 0523 test
			var = chgInt(var);
			duals.push_back(var); //the order of duals is the same with that of coupling constraints
		}
	}else{
		spdlog::error("There is some wrong when solving the dual value of the masterproblem!");
		return 0;
	}
    //@wqy0627 add for avoiding the stopping
    for(int i = 0; i < duals.size(); i++){
        if(duals[i] > 1e15){return -1;}
    }
	return 1;
}


//int MasterProb::SolveDualCplex(vector<double>& duals)
//{
//	SCIPwriteOrigProblem(RMPscipModel, "D://testRMP.lp", nullptr, 0);
//
//	IloEnv env;
//	IloModel model(env);
//	IloCplex cplex(env);
//	IloNumArray dualss;
//	IloNumVarArray vars;
//	IloRangeArray cons(env);
//	IloObjective obj;
//	cplex.importModel(model, "D://testRMP.lp", obj, vars, cons);
//	cplex.extract(model);
//	std::ostream nullStream(nullptr);
//	cplex.setOut(nullStream);
//	cplex.solve();
//	if (cplex.getStatus() != IloAlgorithm::Optimal) {
//		spdlog::error("No optimal solution found.");
//		return 0;
//	}
//	double var = 0;
//	for (int i = 0; i < cons.getSize(); i++){
//		var = cplex.getDual(cons[i]);
//		var = chgInt(var);
//		duals.push_back(var);
//	}
//
//	cplex.end();
//	env.end();
//	return 1;
//}

int MasterProb::SolveDualCplex(vector<double>& duals)
{
	SCIPwriteOrigProblem(RMPscipModel, "testRMP.lp", nullptr, 0);

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(env);
	IloNumArray dualss;
	IloNumVarArray vars;
	IloRangeArray cons(env);
	IloObjective obj;
	cplex.importModel(model, "testRMP.lp", obj, vars, cons);
	cplex.extract(model);
	std::ostream nullStream(nullptr);
	cplex.setOut(nullStream);
	cplex.solve();
	if (cplex.getStatus() != IloAlgorithm::Optimal) {
		spdlog::error("No optimal solution found.");
		return 0;
	}
	double var = 0;
	for (int i = 0; i < cons.getSize(); i++){
		var = cplex.getDual(cons[i]);
		var = chgInt(var);
		duals.push_back(var);
	}

	cplex.end();
	env.end();
	return 1;
}


//@wqy 0627save
//bool MasterProb::AddColumn(const Solution& sol, int subprob_no)
//{
//	Column column;
//	column.solSize = sol.solSize;
//	column.varValue = sol.varValue;
//	column.objValue = sol.objValue;
//	
//	//@wqy 0626 judge whether to repeat
//	for (int i = 0; i < columnpool[subprob_no].Columns().size(); i++) {
//		bool flagDiff(false);
//		for (const auto& var : columnpool[subprob_no].columns_[i].varValue) {
//			if (sol.varValue.find(var.first) == sol.varValue.end()) {
//				spdlog::error("there is some wrong in MasterProb's AddColumn function!");
//			}
//			if (abs(var.second - sol.varValue.find(var.first)->second) > EPSFORINT) {
//				flagDiff = true;
//				break;
//			}
//		}
//		if (!flagDiff) {
//			spdlog::error("Repeat column appeared! The column that is being added to the subproblem {0} is equal to the column {1} in the columnpool!", subprob_no, i);
//			return 0;
//		}
//	}
//	columnpool[subprob_no].AddColumn(column);
//	return 1;
//}

bool MasterProb::AddColumn(const Solution& sol, int subprob_no)
{
	Column column;
	column.solSize = sol.solSize;
	column.varValue = sol.varValue;
	column.objValue = sol.objValue;
	if (columnpool[subprob_no].AddCol(column)) {
		columnpool[subprob_no].AddColumn(column);
		return 1;
	}else {
		spdlog::error("Repeat column appeared! The column that is being added to the subproblem {0} is repeat!", subprob_no);
		return 0;
	}

}

void MasterProb::Solve(RMPSOL& sol, double& objval, bool root_flag)
{
	//SCIPwriteOrigProblem(RMPscipModel, "D://test11.lp", "cip", FALSE);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);  //scip quiet
	SCIPsolve(RMPscipModel);
	int sp_no = 0;
	int add_num = 0;
	for (int i = 0; i < RMPscipVars_num; i++){
		if (i >= columnpool[sp_no].Columns().size() + add_num) {
			add_num += columnpool[sp_no].Columns().size();
			sp_no++;
		}
		sol.subprob_no.push_back(sp_no);
		sol.varValue.push_back(chgInt(SCIPgetSolVal(RMPscipModel, SCIPgetBestSol(RMPscipModel), RMPscipVars[i])));
		sol.varName.push_back(SCIPvarGetName(RMPscipVars[i]));
	}
	objval = SCIPgetSolOrigObj(RMPscipModel, SCIPgetBestSol(RMPscipModel));
}

int MasterProb::SolveObj(double& obj)
{
	//SCIPwriteOrigProblem(RMPscipModel, "D://lastRMP.lp", nullptr, 0);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPsolve(RMPscipModel);
	obj = SCIPgetSolOrigObj(RMPscipModel, SCIPgetBestSol(RMPscipModel));

	//SCIP* scipModelTemp = nullptr;
	//SCIPcreate(&scipModelTemp);
	//SCIPincludeDefaultPlugins(scipModelTemp);
	//SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	//SCIPreadProb(scipModelTemp, "D://lastRMP.lp", NULL);
	//SCIPsolve(scipModelTemp);
	//cout << "recordTest -> SCIPStatus: " << SCIPgetStatus(scipModelTemp) << endl;


	if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_OPTIMAL){
		
		return 1;
	}else{
		cout << "SolveObj -> SCIPStatus: " << SCIPgetStatus(RMPscipModel) << endl;
	}
	return 0;

}


int MasterProb::CheckSet1(int varGlobalIndex) const
{
	//SCIPwriteOrigProblem(RMPscipModel, "D://beforeCheck.lp", nullptr, 0);
	SCIPfreeTransform(RMPscipModel); 
	SCIPchgVarLb(RMPscipModel, RMPscipVars[varGlobalIndex], 1.0);
	SCIPchgVarUb(RMPscipModel, RMPscipVars[varGlobalIndex], 1.0);
	//SCIPwriteOrigProblem(RMPscipModel, "D://afterCheck.lp.lp", nullptr, 0);
	SCIPsolve(RMPscipModel);

	bool feasible(true);
	if (SCIPgetStatus(RMPscipModel) != SCIP_STATUS_OPTIMAL){
		feasible = false;
	}
	cout <<"CheckSet1 -> SCIPStatus"<< SCIPgetStatus(RMPscipModel) << endl;
	SCIPfreeTransform(RMPscipModel);
	SCIPchgVarLb(RMPscipModel, RMPscipVars[varGlobalIndex], 0);
	SCIPchgVarUb(RMPscipModel, RMPscipVars[varGlobalIndex], SCIPinfinity(RMPscipModel));
	return feasible ? 1 : 0;
}

int MasterProb::CheckSet1(string varGlobalName) const
{
	bool flag(false);
	for (int i = 0; i < RMPscipVars_num; i++) {
		if (SCIPvarGetName(RMPscipVars[i]) != varGlobalName && SCIPvarGetName(RMPscipVars[i]) != "t_" + varGlobalName) {
			continue;
		}
		flag = true;
		//SCIPwriteOrigProblem(RMPscipModel, "D://beforeCheck.lp", nullptr, 0);
		SCIPfreeTransform(RMPscipModel);
		SCIPchgVarLb(RMPscipModel, RMPscipVars[i], 1.0);
		SCIPchgVarUb(RMPscipModel, RMPscipVars[i], 1.0);
		//SCIPwriteOrigProblem(RMPscipModel, "D://afterCheck.lp.lp", nullptr, 0);
		SCIPsolve(RMPscipModel);

		bool feasible(true);
		if (SCIPgetStatus(RMPscipModel) != SCIP_STATUS_OPTIMAL) {
			feasible = false;
		}
		//cout << "CheckSet1 -> SCIPStatus" << SCIPgetStatus(RMPscipModel) << endl;
		SCIPfreeTransform(RMPscipModel);
		SCIPchgVarLb(RMPscipModel, RMPscipVars[i], 0);
		SCIPchgVarUb(RMPscipModel, RMPscipVars[i], SCIPinfinity(RMPscipModel));
		return feasible ? 1 : 0;
	}
	
	if (!flag) {
		spdlog::error("there is some wrong in CheckSet1 function!");
	}
	return 0;
}
