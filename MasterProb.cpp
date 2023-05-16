#include"MasterProb.h"

//using namespace std;

MasterProb::MasterProb()
{

}

void MasterProb::free()
{

	for (size_t i = 0; i < RMPscipVars_num; i++)
	{
		if (RMPscipVars[i]->name != NULL && RMPscipModel->mem != NULL && RMPscipVars[i] != NULL)
		{
			SCIPreleaseVar(RMPscipModel, &RMPscipVars[i]);
		}	
	}
	RMPscipVars_num = 0;
	for (size_t i = 0; i < RMPscipConss_num; i++)
	{
		if (RMPscipConss[i] != NULL && RMPscipModel->mem != NULL && RMPscipConss[i]->name != NULL)
		{
			SCIPreleaseCons(RMPscipModel, &RMPscipConss[i]);
		}		
	} 

	SCIPfree(&RMPscipModel);
}

//���
MasterProb& MasterProb::operator += (const MasterProb& mp)
{
	if (this == &mp)
		return *this;

	this->A = mp.A;
	this->columnpool = new ColumnPool[mp.sp_num];
	for (int i = 0; i < mp.sp_num; i++)
	{
		this->columnpool[i] = mp.columnpool[i];
	}	
	this->col_num = mp.col_num;
	this->ConsCoefMatrix = mp.ConsCoefMatrix;
	this->RMPscipConss = new SCIP_CONS*[mp.RMPscipConss_num];
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

//MasterProb::MasterProb(Objective obj, Constraints* couplingCons, int couplingConsNum)
//{
//	A = Eigen::MatrixXd::Zero(couplingConsNum, obj.coef.size());
//	for (size_t i = 0; i < couplingConsNum; i++)
//	{
//		for (size_t j = 0; j < obj.coef.size(); j++)
//		{
//			A(i, j) = couplingCons[i].Coef[j];
//		}
//	}
//}

MasterProb::~MasterProb()
{
}


void MasterProb::BuildSCIPModel(Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX> Fixed_Info)
{
	SCIPcreate(&RMPscipModel);
	SCIPincludeDefaultPlugins(RMPscipModel);
	SCIPsetIntParam(RMPscipModel, "presolving/maxrounds", 0);

	SCIPcreateProbBasic(RMPscipModel, "R_MasterProb");

#pragma region COUT�г��е���
	//for (size_t i = 0; i < subprob_num; i++)
	//{
	//	cout << "������ " << i << " :" << endl;
	//	for (auto iter = columnpool[i].columns.begin(); iter != columnpool[i].columns.end(); iter++)
	//	{
	//		for (auto var_iter = iter->varValue.begin(); var_iter != iter->varValue.end(); var_iter++)
	//		{
	//			cout << var_iter->first << " : " << var_iter->second << endl;
	//		}
	//		cout << "%%%%%%%%%" << endl;
	//	}
	//}

#pragma endregion



	//�������or��С������
	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
	if (obj.direction == Direction::max)
	{
		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
	}
	SCIPsetObjsense(RMPscipModel, objSense);

	int col_num = 0;
	for (size_t i = 0; i < subprob_num; i++)
	{
		col_num += columnpool[i].columns.size();
	}

	//��ӱ���
	unordered_map<int, SCIP_VAR*> RMPscipVarMap;
	RMPscipVars_num = col_num;//���������������г����е�������
	RMPscipVars = new SCIP_VAR * [RMPscipVars_num];
	int scipVarsIdx = 0;
	int col_num_all = 0;//��������
	bool cout_flag = 0;

	ConsCoefMatrix.resize(0,0);
	//cout << "111ConsCoefMatrix" << endl;
	//cout << ConsCoefMatrix;
	//cout << endl;

	for (size_t i = 0; i < subprob_num; i++)
	{
		//��ֵ�г��е���ϵ������
		Eigen::MatrixXd colmatrix;
		colmatrix = Eigen::MatrixXd::Zero(columnpool[i].var_num, columnpool[i].columns.size());
		int column_idx = 0;
		//����vars_id
		sort(columnpool[i].vars_id.begin(), columnpool[i].vars_id.end());

		for (auto iter = columnpool[i].columns.begin(); iter != columnpool[i].columns.end(); ++iter)
		{
			for (size_t j = 0; j < columnpool[i].var_num; j++)
			{
				colmatrix(j, column_idx) = iter->varValue.find(columnpool[i].vars_id[j])->second;
			}
			column_idx++;
		}

		//ȡԭ����Լ������A�Ķ�Ӧ����������Ӿ���
		Eigen::MatrixXd A_part;
		A_part = Eigen::MatrixXd::Zero(couplingConsNum, columnpool[i].var_num);
		for (size_t t = 0; t < couplingConsNum; t++)
		{
			for (size_t j = 0; j < columnpool[i].var_num; j++)
			{
				if (couplingCons[t].idDic.find(columnpool[i].vars_id[j]) == couplingCons[t].idDic.end())//wqy��0403���п����Ҳ���
				{
					A_part(t, j) = 0;
				}
				else
				{
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


		//������ˣ�����������Լ��ϵ��
		Eigen::MatrixXd RMPcoef_matrix = A_part * colmatrix;

		//cout
		if (cout_flag)
		{
			cout << "RMPcoef_matrix" << endl;
			cout << RMPcoef_matrix;
			cout << endl;
		}


		Eigen::MatrixXd Coef_temp(RMPcoef_matrix.rows(), ConsCoefMatrix.cols() + RMPcoef_matrix.cols());

		Coef_temp << ConsCoefMatrix, RMPcoef_matrix;

		ConsCoefMatrix = Coef_temp;

		//ȡԭ����Ŀ�꺯��ϵ����Ӧ��������Ĳ���
		Eigen::VectorXd obj_part;
		obj_part = Eigen::VectorXd::Zero(columnpool[i].var_num);
		for (size_t t = 0; t < columnpool[i].var_num; t++)
		{
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


		//����������Ŀ�꺯��ϵ��
		Eigen::VectorXd RMPobjcoef_vector = obj_part.transpose() * colmatrix;

		//cout
		if (cout_flag)
		{
			cout << "RMPobjcoef_vector" << endl;
			cout << RMPobjcoef_vector;
			cout << endl;
		}


		//���SCIP����
		for (size_t j = 0; j < columnpool[i].columns.size(); j++)
		{
			SCIP_VAR* tempVar;
			//�����������
			SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
			const char* var_name;
			string varname = "_lambda_" + to_string(i) + "_" + to_string(j);
			var_name = varname.c_str();
			SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, 1, RMPobjcoef_vector(j), scipVarType);
			//SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, -SCIPinfinity(RMPscipModel), SCIPinfinity(RMPscipModel), RMPobjcoef_vector(j), scipVarType);//wqy0515test
			SCIPaddVar(RMPscipModel, tempVar);


			RMPscipVars[scipVarsIdx++] = tempVar;
			RMPscipVarMap.insert(make_pair(col_num_all, tempVar));
			col_num_all++;
		}

	}

	//���SCIPԼ��

	RMPscipConss_num = couplingConsNum + subprob_num;
	RMPscipConss = new SCIP_CONS * [RMPscipConss_num];

	//cout << "��ConsCoefMatrix:" << endl;
	//cout << ConsCoefMatrix;
	//cout << endl;

	for (int i = 0; i < couplingConsNum; i++)
	{
		Constraints consTemp = couplingCons[i];

		int tempVarsSize = col_num_all; //�������������
		//SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];
		int idx = 0;
		while (idx != col_num_all)
		{
			tempConsCoefs[idx] = ConsCoefMatrix(i, idx);
			idx++;
		}

		//cout
		//cout << "��" << i << "��Լ��ϵ����" << endl;
		//for (size_t kk = 0; kk < tempVarsSize; kk++)
		//{
		//	cout << tempConsCoefs[kk] << " ";
		//}
		//cout << endl;

		const char* con_name;
		string conname = "z_cons" + to_string(i);
		con_name = conname.c_str();
		double rhs, lhs;
		if (couplingCons[i].Prop == 0)
		{
			lhs = couplingCons[i].rhs;
			rhs = SCIPinfinity(RMPscipModel);
		}
		else if (couplingCons[i].Prop == 1)
		{
			lhs = -SCIPinfinity(RMPscipModel);
			rhs = couplingCons[i].rhs;
		}
		else
		{
			lhs = couplingCons[i].rhs;
			rhs = couplingCons[i].rhs;
		}

		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], con_name,
			tempVarsSize, RMPscipVars, tempConsCoefs, lhs, rhs);
		SCIPaddCons(RMPscipModel, RMPscipConss[i]);
	}
	int record_idx = 0;
	for (int i = couplingConsNum; i < couplingConsNum + subprob_num; i++)
	{
		int col_no;
		int begin_idx, end_idx;//������ֹ���

		int tempVarsSize = columnpool[i - couplingConsNum].columns.size();
		SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];

		begin_idx = record_idx;
		end_idx = begin_idx + tempVarsSize;
		record_idx = end_idx;
		int idx = 0;
		for (auto iter = RMPscipVarMap.begin(); iter != RMPscipVarMap.end(); ++iter)
		{
			if (iter->first >= begin_idx && iter->first < end_idx)
			{
				tempVars[idx] = iter->second;
				//cout << iter->second->name << endl;;
				tempConsCoefs[idx] = 1;
				idx++;
			}
		}

		const char* con_name;
		string conname = "1_cons" + to_string(i);
		con_name = conname.c_str();
		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], con_name,
			tempVarsSize, tempVars, tempConsCoefs, 1, 1);
		SCIPaddCons(RMPscipModel, RMPscipConss[i]);


	}


	SCIPwriteOrigProblem(RMPscipModel, "D://RMP_init.lp", "cip", FALSE);
	//SCIPsolve(scipModel);


}

int MasterProb::Solve_Dual(Solution& duals)
{

	SCIPwriteOrigProblem(RMPscipModel, "D://testRMP.lp", "cip", FALSE);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPsolve(RMPscipModel);

	if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_OPTIMAL)
	{
		duals.solSize = RMPscipConss_num;

		//��ż�����
		for (int i = 0; i < RMPscipConss_num; i++)
		{
			unsigned mark;
			double var;
			SCIPgetDualSolVal(RMPscipModel, RMPscipConss[i], &var, &mark);
			//cout << RMPscipConss[i]->name << ":" << var << endl;

			duals.varValue.insert(pair<int, double>(i, var)); //�������int��Լ����˳�򣬺��治����keyȥȡvalue,ֻ�õ�˳��

		}
	}
	else if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_INFEASIBLE)
	{
		printf("�������޿��н⣡��\n");
		return 0;
	}
	else
	{
		cout << to_string(SCIPgetStatus(RMPscipModel)) << endl;
		printf("����������쳣����\n");
		return 0;
	}
	

	//cout << "���Ž⣺" << endl;
    //cout << SCIPgetSolOrigObj(RMPscipModel, SCIPgetBestSol(RMPscipModel)) << endl;

	//ԭ�����
	//for (int i = 0; i < RMPscipVars_num; i++)
	//{
	//	duals.push_back(SCIPgetSolVal(RMPscipModel, SCIPgetBestSol(RMPscipModel),RMPscipVars[i]));
	//}

	//Solve��Ϳ��԰���Դ�ͷŵ��ˣ�Ȼ�������½�������ģ��
	//--�ͷ���Դ
	//for (size_t i = 0; i < RMPscipVars_num; i++)
	//{
	//	SCIPreleaseVar(RMPscipModel, &RMPscipVars[i]);
	//}
	//for (size_t i = 0; i < RMPscipConss_num; i++)
	//{
	//	SCIPreleaseCons(RMPscipModel, &RMPscipConss[i]);
	//}
	//SCIPfree(&RMPscipModel);
	return 1;


}


void MasterProb::AddColumn(Solution sol, int subprob_no)
{
	Column column;
	column.solSize = sol.solSize;
	column.varValue = sol.varValue;
	column.objValue = sol.objValue;
	columnpool[subprob_no].AddCol(column);

}

void MasterProb::Solve(RMPSOL& sol, double& objval)
{
	SCIPwriteOrigProblem(RMPscipModel, "D://test11.lp", "cip", FALSE);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPsolve(RMPscipModel);
	int sp_no = 0;
	int add_num = 0;
	for (int i = 0; i < RMPscipVars_num; i++)
	{
		sol.varValue.push_back(SCIPgetSolVal(RMPscipModel, SCIPgetBestSol(RMPscipModel), RMPscipVars[i]));
	}
	for (size_t i = 0; i < RMPscipVars_num; i++)
	{
		if (i >= columnpool[sp_no].columns.size() + add_num)
		{
			add_num += columnpool[sp_no].columns.size();
			sp_no++;
		}
		sol.subprob_no.push_back(sp_no);
	}
	objval = SCIPgetSolOrigObj(RMPscipModel, SCIPgetBestSol(RMPscipModel));
}

void MasterProb::Solve_Obj(double& obj)
{
	SCIPwriteOrigProblem(RMPscipModel, "D://lastRMP.lp", "cip", FALSE);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPsolve(RMPscipModel);
	obj = SCIPgetSolOrigObj(RMPscipModel, SCIPgetBestSol(RMPscipModel));
	//cout
	//for (size_t i = 0; i < RMPscipVars_num; i++)
	//{
	//	cout << RMPscipVars[i]->name<<":"<< SCIPgetSolVal(RMPscipModel, SCIPgetBestSol(RMPscipModel), RMPscipVars[i]) << endl;
	//}
}

//��BuildRMP��ֻ࣬�Ƕ�һ������Ƿ����
void MasterProb::Test_feasble2(bool& feasible, Objective obj, Constraints* couplingCons, int couplingConsNum, int subprob_num, vector<FIX>& Fixed_Info)
{
	SCIPcreate(&RMPscipModel);
	SCIPincludeDefaultPlugins(RMPscipModel);
	SCIPsetIntParam(RMPscipModel, "presolving/maxrounds", 0);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	


	SCIPcreateProbBasic(RMPscipModel, "R_MasterProb");

#pragma region COUT�г��е���
	//for (size_t i = 0; i < subprob_num; i++)
	//{
	//	cout << "������ " << i << " :" << endl;
	//	for (auto iter = columnpool[i].columns.begin(); iter != columnpool[i].columns.end(); iter++)
	//	{
	//		for (auto var_iter = iter->varValue.begin(); var_iter != iter->varValue.end(); var_iter++)
	//		{
	//			cout << var_iter->first << " : " << var_iter->second << endl;
	//		}
	//		cout << "%%%%%%%%%" << endl;
	//	}
	//}

#pragma endregion


	//�������or��С������
	SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
	if (obj.direction == Direction::max)
	{
		objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
	}
	SCIPsetObjsense(RMPscipModel, objSense);

	int col_num = 0;
	for (size_t i = 0; i < subprob_num; i++)
	{
		col_num += columnpool[i].columns.size();
	}

	//��ӱ���
	unordered_map<int, SCIP_VAR*> RMPscipVarMap;
	RMPscipVars_num = col_num;//���������������г����е�������
	RMPscipVars = new SCIP_VAR * [RMPscipVars_num];
	int scipVarsIdx = 0;
	int col_num_all = 0;//��������
	bool cout_flag = 0;

	ConsCoefMatrix.resize(0,0);//���°����������ֵΪ�գ���Ȼһֱ�ۼ�

	for (size_t i = 0; i < subprob_num; i++)
	{
		//��ֵ�г��е���ϵ������
		Eigen::MatrixXd colmatrix;
		colmatrix = Eigen::MatrixXd::Zero(columnpool[i].var_num, columnpool[i].columns.size());
		int column_idx = 0;
		//����vars_id
		sort(columnpool[i].vars_id.begin(), columnpool[i].vars_id.end());

		for (auto iter = columnpool[i].columns.begin(); iter != columnpool[i].columns.end(); ++iter)
		{
			for (size_t j = 0; j < columnpool[i].var_num; j++)
			{
				colmatrix(j, column_idx) = iter->varValue.find(columnpool[i].vars_id[j])->second;
			}
			column_idx++;
		}

		//ȡԭ����Լ������A�Ķ�Ӧ����������Ӿ���
		Eigen::MatrixXd A_part;
		A_part = Eigen::MatrixXd::Zero(couplingConsNum, columnpool[i].var_num);
		for (size_t t = 0; t < couplingConsNum; t++)
		{
			for (size_t j = 0; j < columnpool[i].var_num; j++)
			{
				if (couplingCons[t].idDic.find(columnpool[i].vars_id[j]) == couplingCons[t].idDic.end())//wqy��0403���п����Ҳ���
				{
					A_part(t, j) = 0;
				}
				else
				{
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


		//������ˣ�����������Լ��ϵ��
		Eigen::MatrixXd RMPcoef_matrix = A_part * colmatrix;

		//cout
		if (cout_flag)
		{
			cout << "RMPcoef_matrix" << endl;
			cout << RMPcoef_matrix;
			cout << endl;
		}


		Eigen::MatrixXd Coef_temp(RMPcoef_matrix.rows(), ConsCoefMatrix.cols() + RMPcoef_matrix.cols());

		Coef_temp << ConsCoefMatrix, RMPcoef_matrix;

		ConsCoefMatrix = Coef_temp;

		//ȡԭ����Ŀ�꺯��ϵ����Ӧ��������Ĳ���
		Eigen::VectorXd obj_part;
		obj_part = Eigen::VectorXd::Zero(columnpool[i].var_num);
		for (size_t t = 0; t < columnpool[i].var_num; t++)
		{
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
		}


		//����������Ŀ�꺯��ϵ��
		Eigen::VectorXd RMPobjcoef_vector = obj_part.transpose() * colmatrix;

		//cout
		if (cout_flag)
		{
			cout << "RMPobjcoef_vector" << endl;
			cout << RMPobjcoef_vector;
			cout << endl;
		}



		//���SCIP����
		for (size_t j = 0; j < columnpool[i].columns.size(); j++)
		{
			SCIP_VAR* tempVar;
			//�����������
			SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
			const char* var_name;
			string varname = "_lambda_" + to_string(i) + "_" + to_string(j);
			var_name = varname.c_str();
			SCIPcreateVarBasic(RMPscipModel, &tempVar, var_name, 0, 1, RMPobjcoef_vector(j), scipVarType);
			SCIPaddVar(RMPscipModel, tempVar);


			RMPscipVars[scipVarsIdx++] = tempVar;
			RMPscipVarMap.insert(make_pair(col_num_all, tempVar));
			col_num_all++;
		}

	}

	//���SCIPԼ��

	RMPscipConss_num = couplingConsNum + subprob_num;
	RMPscipConss = new SCIP_CONS * [RMPscipConss_num];

	for (int i = 0; i < couplingConsNum; i++)
	{
		Constraints consTemp = couplingCons[i];

		int tempVarsSize = col_num_all; //�������������
		//SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];
		int idx = 0;
		while (idx != col_num_all)
		{
			tempConsCoefs[idx] = ConsCoefMatrix(i, idx);
			idx++;
		}

		//cout
		//cout << "��" << i << "��Լ��ϵ����" << endl;
		//for (size_t kk = 0; kk < tempVarsSize; kk++)
		//{
		//	cout << tempConsCoefs[kk] << " ";
		//}
		//cout << endl;

		const char* con_name;
		string conname = "z_cons" + to_string(i);
		con_name = conname.c_str();
		double rhs, lhs;
		if (couplingCons[i].Prop == 0)
		{
			lhs = couplingCons[i].rhs;
			rhs = SCIPinfinity(RMPscipModel);
		}
		else if (couplingCons[i].Prop == 1)
		{
			lhs = -SCIPinfinity(RMPscipModel);
			rhs = couplingCons[i].rhs;
		}
		else
		{
			lhs = couplingCons[i].rhs;
			rhs = couplingCons[i].rhs;
		}

		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], con_name,
			tempVarsSize, RMPscipVars, tempConsCoefs, lhs, rhs);
		SCIPaddCons(RMPscipModel, RMPscipConss[i]);
	}
	int record_idx = 0;
	for (int i = couplingConsNum; i < couplingConsNum + subprob_num; i++)
	{
		int col_no;
		int begin_idx, end_idx;//������ֹ���

		int tempVarsSize = columnpool[i - couplingConsNum].columns.size();
		SCIP_VAR** tempVars = new SCIP_VAR * [tempVarsSize];
		double* tempConsCoefs = new double[tempVarsSize];

		begin_idx = record_idx;
		end_idx = begin_idx + tempVarsSize;
		record_idx = end_idx;
		int idx = 0;
		for (auto iter = RMPscipVarMap.begin(); iter != RMPscipVarMap.end(); ++iter)
		{
			if (iter->first >= begin_idx && iter->first < end_idx)
			{
				tempVars[idx] = iter->second;
				//cout << iter->second->name << endl;
				tempConsCoefs[idx] = 1;
				idx++;
			}
		}

		const char* con_name;
		string conname = "1_cons" + to_string(i);
		con_name = conname.c_str();
		SCIPcreateConsBasicLinear(RMPscipModel, &RMPscipConss[i], con_name,
			tempVarsSize, tempVars, tempConsCoefs, 1, 1);
		SCIPaddCons(RMPscipModel, RMPscipConss[i]);


	}


	
	SCIPwriteOrigProblem(RMPscipModel, "D://fixed.lp", "cip", FALSE);
	SCIPsetIntParam(RMPscipModel, "display/verblevel", 0);
	SCIPsolve(RMPscipModel);

	if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_OPTIMAL)
	{
		//printf("The model has a feasible solution!\n");
		feasible = true;
	}
	else if (SCIPgetStatus(RMPscipModel) == SCIP_STATUS_INFEASIBLE)
	{
		//printf("The model is infeasible.\n");
		feasible = false;
	}
	else
	{
		cout << to_string(SCIPgetStatus(RMPscipModel)) << endl;
		feasible = false;
	}

}