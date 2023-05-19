//#include "general.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include "linear.hpp"
#include"DW.h"
#include"Constraints.h"
#include"Objective.h"
#include"Lagrangean.h"
#include "Dection.hpp"
//#include "ThreadPool.h"




using namespace std;

static const int THREAD_NUMS_MAX = 7;
static const int SUB_THRESHOLD = 1;

//ThreadPool g_Pool(THREAD_NUMS_MAX);
void Input_InitColumns(ColumnPool*&columnpool_Input, std::map<int, General_Block> block);
void Generate_InitColumns(ColumnPool*& columnpool_Input, std::map<int, General_Block> block,int col_num);
void Generate_InitColumns_2(ColumnPool*& columnpool_Input, std::map<int, General_Block> block, int col_num,SCIP*scip);


int main()
{
	clock_t time_1 = clock();
	//=================超图划分开始=========================//
	string Folder = "E:\\MipSolver\\collectionunzip";
    //vector<pair<string, string>> NameList;
	auto Linear = FindLinear("D:\\result2.csv");

	ofstream outfile("resul3.csv", std::ios::app);
	//遍历NameList

	string tempname = "D:\\a1c1s1.mps";
	//tempname = "D:\\test50bin100con.lp";
	//tempname = "D:\\test50bin100con3.lp";
	//tempname = "D:\\test50bin100con3_conti.lp";
	//tempname = "D:\\blend2.mps";
	//tempname = "G:\\collectionunzip\\f2gap40400.mps";//加重复列
	//tempname = "G:\\collectionunzip\\f2gap401600.mps";//算半天不知道算啥呢
	//tempname = "G:\\collectionunzip\\f2gap201600.mps"; //同上
	//tempname = "D:\\primalmodel_conti.lp";
	//tempname = "D:\\mipdata\\bab1.mps";//报错 sp.scipVars[i]为空
	//tempname = "D:\\mipdata\\berlin.mps";//RMP无可行解，求最优太慢
	//tempname = "D:\\mipdata\\neos-1112787.mps";
	//tempname = "D:\\mipdata\\beavma.mps";
	tempname = "D:\\blend2_10000.lp";


	
	auto start = std::chrono::high_resolution_clock::now();
	int block_num = 6;
	auto [block, scipmodel, scipvars, scipcons, _ori_vars_x, _ori_vars_z] = Dection(tempname, block_num);
	//SCIPfree(&scipmodel);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	//=================超图划分结束=========================//
	int mode = 2;//1为lag，2为自生成 ------------------------------------------ mode -----------
	clock_t time_2 = clock();
	for (size_t i = 0; i <= block_num; i++)
	{
		cout<<block[i].bVars.size()<<","<<block[i].bCons.size()<<endl;
	}
	Objective master_obj;
	for (size_t i = 0; i < block_num; i++)
	{
		for (auto v : block[i].bobj.coef)
		{
			master_obj.coef[v.first] = v.second;
		}
	}
	cout << "master obj size: " << master_obj.coef.size() << endl;
	//=================拉格朗日算法开始=========================//

	
	//=================拉格朗日算法结束=========================//
    clock_t time_3 = clock();
	//cout
	int cout_flagg = 0;
	if (cout_flagg)
	{
		for (auto iter = block.begin(); iter != block.end(); iter++)
		{
			cout << "BLOCK " << iter->first << endl;
			cout << endl;
			cout << "Vars.size(): " << iter->second.bVars.size() << endl;
			cout << "Vars.Id" << endl;

			for (auto Var_iter = iter->second.bVars.begin(); Var_iter != iter->second.bVars.end(); Var_iter++)
			{
				cout << Var_iter->Name << "  ";
			}

			cout << endl;
			cout << "Cons.size(): " << iter->second.bCons.size() << endl;
			cout << "Cons.Id" << endl;

			for (auto Con_iter = iter->second.bCons.begin(); Con_iter != iter->second.bCons.end(); Con_iter++)
			{
				cout << Con_iter->Name << "  ";
			}
			cout << endl;

		}
	}
	
	//=================Diving算法开始=========================//
	//自定义初始列，测试用
	ColumnPool* columnpool_Input = new ColumnPool[block.size() - 1];
	//Input_InitColumns(columnpool_Input, block);
	if (mode == 2)
	{
		int col_numm = 4;
		Generate_InitColumns_2(columnpool_Input, block, col_numm, scipmodel);
		DW::Diving(columnpool_Input,block,master_obj);
	}

	if (mode == 1)
	{
		Lagrangean lag{ block, scipmodel, scipcons, scipvars };
		lag.RUN(1000);
		DW::Diving(lag.columnpool, block, master_obj);
	}
	
	//=================Diving算法结束=========================//
	clock_t time_4 = clock();

	//时间统计
	cout << "超图划分用时：" << endl;
	cout << (time_2 - time_1)/CLOCKS_PER_SEC << "秒" << endl;
	cout << "拉格朗日用时：" << endl;
	cout << (time_3 - time_2) / CLOCKS_PER_SEC << "秒" << endl;
	cout << "Diving用时：" << endl;
	cout << (time_4 - time_3) / CLOCKS_PER_SEC << "秒" << endl;
}

void Input_InitColumns(ColumnPool*& columnpool_Input, std::map<int, General_Block> block)
{
	columnpool_Input[0].var_num = 25;
	vector<int>id0;
	for (size_t i = 0; i < block[0].bVars.size(); i++)
	{
		id0.push_back(block[0].bVars[i].Id);
	}
	columnpool_Input[0].vars_id = id0;
#pragma region column1
	Column column1;
	column1.solSize = 25;
	column1.varValue.insert(pair<int, double>(0, 0));
	column1.varValue.insert(pair<int, double>(3, 2));
	column1.varValue.insert(pair<int, double>(4, 10));
	column1.varValue.insert(pair<int, double>(5, 1));
	column1.varValue.insert(pair<int, double>(7, 0));
	column1.varValue.insert(pair<int, double>(8, 10));
	column1.varValue.insert(pair<int, double>(10, 0));
	column1.varValue.insert(pair<int, double>(12, 0));
	column1.varValue.insert(pair<int, double>(13, 2));
	column1.varValue.insert(pair<int, double>(14, 0));
	column1.varValue.insert(pair<int, double>(15, 0));
	column1.varValue.insert(pair<int, double>(16, 10));
	column1.varValue.insert(pair<int, double>(17, 2));
	column1.varValue.insert(pair<int, double>(18, 0));
	column1.varValue.insert(pair<int, double>(19, 0));
	column1.varValue.insert(pair<int, double>(21, 10));
	column1.varValue.insert(pair<int, double>(22, 0));
	column1.varValue.insert(pair<int, double>(25, 0));
	column1.varValue.insert(pair<int, double>(26, 0));
	column1.varValue.insert(pair<int, double>(31, 0));
	column1.varValue.insert(pair<int, double>(37, 0));
	column1.varValue.insert(pair<int, double>(40, 10));
	column1.varValue.insert(pair<int, double>(42, 1));
	column1.varValue.insert(pair<int, double>(44, 0));
	column1.varValue.insert(pair<int, double>(49, 1));
	columnpool_Input[0].columns.insert(column1);
#pragma endregion

#pragma region column2
	Column column2;
	column2.solSize = 25;
	column2.varValue.insert(pair<int, double>(0, 0));
	column2.varValue.insert(pair<int, double>(3, 1));
	column2.varValue.insert(pair<int, double>(4, 10));
	column2.varValue.insert(pair<int, double>(5, 0));
	column2.varValue.insert(pair<int, double>(7, 0));
	column2.varValue.insert(pair<int, double>(8, 10));
	column2.varValue.insert(pair<int, double>(10, 0));
	column2.varValue.insert(pair<int, double>(12, 0));
	column2.varValue.insert(pair<int, double>(13, 2));
	column2.varValue.insert(pair<int, double>(14, 0));
	column2.varValue.insert(pair<int, double>(15, 0));
	column2.varValue.insert(pair<int, double>(16, 10));
	column2.varValue.insert(pair<int, double>(17, 2));
	column2.varValue.insert(pair<int, double>(18, 0));
	column2.varValue.insert(pair<int, double>(19, 0));
	column2.varValue.insert(pair<int, double>(21, 0));
	column2.varValue.insert(pair<int, double>(22, 1));
	column2.varValue.insert(pair<int, double>(25, 0));
	column2.varValue.insert(pair<int, double>(26, 0));
	column2.varValue.insert(pair<int, double>(31, 0));
	column2.varValue.insert(pair<int, double>(37, 0));
	column2.varValue.insert(pair<int, double>(40, 10));
	column2.varValue.insert(pair<int, double>(42, 1));
	column2.varValue.insert(pair<int, double>(44, 0));
	column2.varValue.insert(pair<int, double>(49, 1));
	columnpool_Input[0].columns.insert(column2);
#pragma endregion

#pragma region column3
	Column column3;
	column3.solSize = 25;
	column3.varValue.insert(pair<int, double>(0, 0));
	column3.varValue.insert(pair<int, double>(3, 0));
	column3.varValue.insert(pair<int, double>(4, 10));
	column3.varValue.insert(pair<int, double>(5, 10));
	column3.varValue.insert(pair<int, double>(7, 0));
	column3.varValue.insert(pair<int, double>(8, 2));
	column3.varValue.insert(pair<int, double>(10, 0));
	column3.varValue.insert(pair<int, double>(12, 0));
	column3.varValue.insert(pair<int, double>(13, 10));
	column3.varValue.insert(pair<int, double>(14, 0));
	column3.varValue.insert(pair<int, double>(15, 0));
	column3.varValue.insert(pair<int, double>(16, 0));
	column3.varValue.insert(pair<int, double>(17, 2));
	column3.varValue.insert(pair<int, double>(18, 0));
	column3.varValue.insert(pair<int, double>(19, 0));
	column3.varValue.insert(pair<int, double>(21, 10));
	column3.varValue.insert(pair<int, double>(22, 0));
	column3.varValue.insert(pair<int, double>(25, 0));
	column3.varValue.insert(pair<int, double>(26, 0));
	column3.varValue.insert(pair<int, double>(31, 0));
	column3.varValue.insert(pair<int, double>(37, 0));
	column3.varValue.insert(pair<int, double>(40, 1));
	column3.varValue.insert(pair<int, double>(42, 1));
	column3.varValue.insert(pair<int, double>(44, 0));
	column3.varValue.insert(pair<int, double>(49, 1));
	columnpool_Input[0].columns.insert(column3);
#pragma endregion
	
#pragma region column4
	Column column4;
	column4.solSize = 25;
	column4.varValue.insert(pair<int, double>(0, 0));
	column4.varValue.insert(pair<int, double>(3, 1));
	column4.varValue.insert(pair<int, double>(4, 10));
	column4.varValue.insert(pair<int, double>(5, 0));
	column4.varValue.insert(pair<int, double>(7, 0));
	column4.varValue.insert(pair<int, double>(8, 10));
	column4.varValue.insert(pair<int, double>(10, 0));
	column4.varValue.insert(pair<int, double>(12, 0));
	column4.varValue.insert(pair<int, double>(13, 2));
	column4.varValue.insert(pair<int, double>(14, 0));
	column4.varValue.insert(pair<int, double>(15, 0));
	column4.varValue.insert(pair<int, double>(16, 10));
	column4.varValue.insert(pair<int, double>(17, 2));
	column4.varValue.insert(pair<int, double>(18, 0));
	column4.varValue.insert(pair<int, double>(19, 0));
	column4.varValue.insert(pair<int, double>(21, 10));
	column4.varValue.insert(pair<int, double>(22, 1));
	column4.varValue.insert(pair<int, double>(25, 0));
	column4.varValue.insert(pair<int, double>(26, 0));
	column4.varValue.insert(pair<int, double>(31, 0));
	column4.varValue.insert(pair<int, double>(37, 0));
	column4.varValue.insert(pair<int, double>(40, 1));
	column4.varValue.insert(pair<int, double>(42, 1));
	column4.varValue.insert(pair<int, double>(44, 0));
	column4.varValue.insert(pair<int, double>(49, 1));
	columnpool_Input[0].columns.insert(column4);
#pragma endregion

#pragma region column5
	Column column5;
	column5.solSize = 25;
	column5.varValue.insert(pair<int, double>(0, 0));
	column5.varValue.insert(pair<int, double>(3, 1));
	column5.varValue.insert(pair<int, double>(4, 1));
	column5.varValue.insert(pair<int, double>(5, 0));
	column5.varValue.insert(pair<int, double>(7, 0));
	column5.varValue.insert(pair<int, double>(8, 2));
	column5.varValue.insert(pair<int, double>(10, 0));
	column5.varValue.insert(pair<int, double>(12, 0));
	column5.varValue.insert(pair<int, double>(13, 1));
	column5.varValue.insert(pair<int, double>(14, 1));
	column5.varValue.insert(pair<int, double>(15, 0));
	column5.varValue.insert(pair<int, double>(16, 0));
	column5.varValue.insert(pair<int, double>(17, 2));
	column5.varValue.insert(pair<int, double>(18, 0));
	column5.varValue.insert(pair<int, double>(19, 0));
	column5.varValue.insert(pair<int, double>(21, 0));
	column5.varValue.insert(pair<int, double>(22, 1));
	column5.varValue.insert(pair<int, double>(25, 0));
	column5.varValue.insert(pair<int, double>(26, 1));
	column5.varValue.insert(pair<int, double>(31, 0));
	column5.varValue.insert(pair<int, double>(37, 0));
	column5.varValue.insert(pair<int, double>(40, 1));
	column5.varValue.insert(pair<int, double>(42, 1));
	column5.varValue.insert(pair<int, double>(44, 0));
	column5.varValue.insert(pair<int, double>(49, 1));
	columnpool_Input[0].columns.insert(column5);
#pragma endregion

	columnpool_Input[1].var_num = 25;
	vector<int>id1;
	for (size_t i = 0; i < block[1].bVars.size(); i++)
	{
		id1.push_back(block[1].bVars[i].Id);
	}
	columnpool_Input[1].vars_id = id1;
#pragma region column6
	Column column6;
	column6.solSize = 25;
	column6.varValue.insert(pair<int, double>(1, 0));
	column6.varValue.insert(pair<int, double>(2, 0));
	column6.varValue.insert(pair<int, double>(6, 3));
	column6.varValue.insert(pair<int, double>(9, 0));
	column6.varValue.insert(pair<int, double>(11, 0));
	column6.varValue.insert(pair<int, double>(20, 0));
	column6.varValue.insert(pair<int, double>(23, 0));
	column6.varValue.insert(pair<int, double>(24, 1));
	column6.varValue.insert(pair<int, double>(27, 0));
	column6.varValue.insert(pair<int, double>(28, 0));
	column6.varValue.insert(pair<int, double>(29, 0));
	column6.varValue.insert(pair<int, double>(30, 0));
	column6.varValue.insert(pair<int, double>(32, 1));
	column6.varValue.insert(pair<int, double>(33, 0));
	column6.varValue.insert(pair<int, double>(34, 0));
	column6.varValue.insert(pair<int, double>(35, 1));
	column6.varValue.insert(pair<int, double>(36, 0));
	column6.varValue.insert(pair<int, double>(38, 0));
	column6.varValue.insert(pair<int, double>(39, 0));
	column6.varValue.insert(pair<int, double>(41, 1));
	column6.varValue.insert(pair<int, double>(43, 0));
	column6.varValue.insert(pair<int, double>(45, 3));
	column6.varValue.insert(pair<int, double>(46, 0));
	column6.varValue.insert(pair<int, double>(47, 0));
	column6.varValue.insert(pair<int, double>(48, 0));
	columnpool_Input[1].columns.insert(column6);
#pragma endregion

#pragma region column7
	Column column7;
	column7.solSize = 25;
	column7.varValue.insert(pair<int, double>(1, 0));
	column7.varValue.insert(pair<int, double>(2, 0));
	column7.varValue.insert(pair<int, double>(6, 3));
	column7.varValue.insert(pair<int, double>(9, 0));
	column7.varValue.insert(pair<int, double>(11, 0));
	column7.varValue.insert(pair<int, double>(20, 0));
	column7.varValue.insert(pair<int, double>(23, 0));
	column7.varValue.insert(pair<int, double>(24, 1));
	column7.varValue.insert(pair<int, double>(27, 10));
	column7.varValue.insert(pair<int, double>(28, 0));
	column7.varValue.insert(pair<int, double>(29, 0));
	column7.varValue.insert(pair<int, double>(30, 0));
	column7.varValue.insert(pair<int, double>(32, 1));
	column7.varValue.insert(pair<int, double>(33, 0));
	column7.varValue.insert(pair<int, double>(34, 10));
	column7.varValue.insert(pair<int, double>(35, 1));
	column7.varValue.insert(pair<int, double>(36, 10));
	column7.varValue.insert(pair<int, double>(38, 10));
	column7.varValue.insert(pair<int, double>(39, 0));
	column7.varValue.insert(pair<int, double>(41, 1));
	column7.varValue.insert(pair<int, double>(43, 0));
	column7.varValue.insert(pair<int, double>(45, 3));
	column7.varValue.insert(pair<int, double>(46, 0));
	column7.varValue.insert(pair<int, double>(47, 10));
	column7.varValue.insert(pair<int, double>(48, 10));
	columnpool_Input[1].columns.insert(column7);
#pragma endregion

#pragma region column8
	Column column8;
	column8.solSize = 25;
	column8.varValue.insert(pair<int, double>(1, 0));
	column8.varValue.insert(pair<int, double>(2, 10));
	column8.varValue.insert(pair<int, double>(6, 3));
	column8.varValue.insert(pair<int, double>(9, 0));
	column8.varValue.insert(pair<int, double>(11, 0));
	column8.varValue.insert(pair<int, double>(20, 0));
	column8.varValue.insert(pair<int, double>(23, 0));
	column8.varValue.insert(pair<int, double>(24, 1));
	column8.varValue.insert(pair<int, double>(27, 0));
	column8.varValue.insert(pair<int, double>(28, 10));
	column8.varValue.insert(pair<int, double>(29, 0));
	column8.varValue.insert(pair<int, double>(30, 0));
	column8.varValue.insert(pair<int, double>(32, 0));
	column8.varValue.insert(pair<int, double>(33, 0));
	column8.varValue.insert(pair<int, double>(34, 0));
	column8.varValue.insert(pair<int, double>(35, 0));
	column8.varValue.insert(pair<int, double>(36, 0));
	column8.varValue.insert(pair<int, double>(38, 0));
	column8.varValue.insert(pair<int, double>(39, 0));
	column8.varValue.insert(pair<int, double>(41, 1));
	column8.varValue.insert(pair<int, double>(43, 0));
	column8.varValue.insert(pair<int, double>(45, 1));
	column8.varValue.insert(pair<int, double>(46, 0));
	column8.varValue.insert(pair<int, double>(47, 0));
	column8.varValue.insert(pair<int, double>(48, 10));
	columnpool_Input[1].columns.insert(column8);
#pragma endregion

#pragma region column9
	Column column9;
	column9.solSize = 25;
	column9.varValue.insert(pair<int, double>(1, 0));
	column9.varValue.insert(pair<int, double>(2, 0));
	column9.varValue.insert(pair<int, double>(6, 3));
	column9.varValue.insert(pair<int, double>(9, 0));
	column9.varValue.insert(pair<int, double>(11, 0));
	column9.varValue.insert(pair<int, double>(20, 0));
	column9.varValue.insert(pair<int, double>(23, 0));
	column9.varValue.insert(pair<int, double>(24, 0));
	column9.varValue.insert(pair<int, double>(27, 10));
	column9.varValue.insert(pair<int, double>(28, 10));
	column9.varValue.insert(pair<int, double>(29, 0));
	column9.varValue.insert(pair<int, double>(30, 2));
	column9.varValue.insert(pair<int, double>(32, 1));
	column9.varValue.insert(pair<int, double>(33, 0));
	column9.varValue.insert(pair<int, double>(34, 10));
	column9.varValue.insert(pair<int, double>(35, 0));
	column9.varValue.insert(pair<int, double>(36, 0));
	column9.varValue.insert(pair<int, double>(38, 0));
	column9.varValue.insert(pair<int, double>(39, 0));
	column9.varValue.insert(pair<int, double>(41, 1));
	column9.varValue.insert(pair<int, double>(43, 0));
	column9.varValue.insert(pair<int, double>(45, 1));
	column9.varValue.insert(pair<int, double>(46, 0));
	column9.varValue.insert(pair<int, double>(47, 0));
	column9.varValue.insert(pair<int, double>(48, 10));
	columnpool_Input[1].columns.insert(column9);
#pragma endregion

#pragma region column10
	Column column10;
	column10.solSize = 25;
	column10.varValue.insert(pair<int, double>(1, 0));
	column10.varValue.insert(pair<int, double>(2, 0));
	column10.varValue.insert(pair<int, double>(6, 3));
	column10.varValue.insert(pair<int, double>(9, 10));
	column10.varValue.insert(pair<int, double>(11, 0));
	column10.varValue.insert(pair<int, double>(20, 0));
	column10.varValue.insert(pair<int, double>(23, 0));
	column10.varValue.insert(pair<int, double>(24, 0));
	column10.varValue.insert(pair<int, double>(27, 0));
	column10.varValue.insert(pair<int, double>(28, 10));
	column10.varValue.insert(pair<int, double>(29, 0));
	column10.varValue.insert(pair<int, double>(30, 10));
	column10.varValue.insert(pair<int, double>(32, 1));
	column10.varValue.insert(pair<int, double>(33, 0));
	column10.varValue.insert(pair<int, double>(34, 0));
	column10.varValue.insert(pair<int, double>(35, 0));
	column10.varValue.insert(pair<int, double>(36, 0));
	column10.varValue.insert(pair<int, double>(38, 0));
	column10.varValue.insert(pair<int, double>(39, 0));
	column10.varValue.insert(pair<int, double>(41, 1));
	column10.varValue.insert(pair<int, double>(43, 0));
	column10.varValue.insert(pair<int, double>(45, 1));
	column10.varValue.insert(pair<int, double>(46, 0));
	column10.varValue.insert(pair<int, double>(47, 10));
	column10.varValue.insert(pair<int, double>(48, 10));
	columnpool_Input[1].columns.insert(column10);
#pragma endregion


}

void Generate_InitColumns(ColumnPool*& columnpool_Input, std::map<int, General_Block> block, int col_num)
{
	SCIP* scipModel;
	SCIP_VAR** scipVars;
	SCIP_CONS** scipConss;

	int sp_num = block.size() - 1;
	for (int subproblem_no = 0; subproblem_no < sp_num; subproblem_no++)
	{
		//创建模型
		SCIPcreate(&scipModel);
		SCIPincludeDefaultPlugins(scipModel);
		SCIPcreateProbBasic(scipModel, "InitsubProb");
		SCIPsetIntParam(scipModel, "limits/solutions", 10);
		SCIPsetIntParam(scipModel, "misc/multistart/maxstarts", 10);
		SCIPsetIntParam(scipModel, "presolving/maxrounds", 0);
		unordered_map<int, SCIP_VAR*> scipVarMap;
		unordered_map<int, Vars> varDic;
		SCIP_VAR** scipVars;
		scipVars = new SCIP_VAR * [block[subproblem_no].bVars.size()];
		
		//设置目标
		SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
		if (block[subproblem_no].bobj.direction == Direction::max)
		{
			objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
		}
		SCIPsetObjsense(scipModel, objSense);
		//添加变量
		for (int i = 0; i < block[subproblem_no].bVars.size(); i++)
		{
			Vars var_now = block[subproblem_no].bVars[i];
			columnpool_Input[subproblem_no].vars_id.push_back(var_now.Id);
			SCIP_VAR* tempVar;
			VarType varType = var_now.Type;
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
			string VarName = var_now.Name;
			SCIPcreateVarBasic(scipModel, &tempVar, VarName.c_str(), var_now.Lb, var_now.Ub, block[subproblem_no].bobj.coef[var_now.Id], scipVarType);//wqy可能有误1，目标函数系数
			SCIPaddVar(scipModel, tempVar);
			scipVarMap.insert(make_pair(var_now.Id, tempVar));
			varDic.insert(pair<int, Vars>(var_now.Id, var_now));
			scipVars[i] = tempVar;
		}
		//添加约束
		scipConss = new SCIP_CONS * [block[subproblem_no].bCons.size()];
		for (int i = 0; i < block[subproblem_no].bCons.size(); i++)
		{
			Constraints consTemp = block[subproblem_no].bCons[i];
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
			SCIPcreateConsBasicLinear(scipModel, &scipConss[i], consTemp.Name.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
			SCIPaddCons(scipModel, scipConss[i]);

			delete[] tempVars;
			delete[] tempConsCoefs;	
		}

		SCIPwriteOrigProblem(scipModel, "D://gen_col.lp", NULL, false);
		SCIPsolve(scipModel);

		//赋值列池
		columnpool_Input[subproblem_no].var_num = block[subproblem_no].bVars.size();
		

		//获得多个可行解
		int nsols = SCIPgetNSols(scipModel);
		SCIP_SOL** sols = SCIPgetSols(scipModel);
		for (int i = 0; i < nsols && i < col_num; ++i)
		{
			Column column_temp;
			

			SCIP_SOL* sol = sols[i];
			double obj = SCIPgetSolOrigObj(scipModel, sol);
			//printf("Solution %d: objective = %f\n", i + 1, obj);
			//SCIP_VAR** vars = SCIPgetVars(scipModel);
			int nvars = SCIPgetNVars(scipModel);
			for (int j = 0; j < nvars; ++j)
			{
				SCIP_Real val = SCIPgetSolVal(scipModel, sol, scipVars[j]);
				//printf("  %s = %f\n", SCIPvarGetName(scipVars[j]), val);
				
				for (auto iter = varDic.begin(); iter != varDic.end(); iter++)
				{
					string name = scipVars[j]->name;
					if (name == iter->second.Name)
					{
						column_temp.varValue.insert(pair<int, double>(iter->second.Id, val)); 
						break;
					}
					if (++iter == varDic.end())
					{
						cout << "Error! Not Found!!" << endl;
					}
					iter--;
					
				}	
			}
			column_temp.solSize = nvars;

			columnpool_Input[subproblem_no].columns.insert(column_temp);
		}

		//SCIPfreeSols(scipModel, &sols, nsols);
		//SCIPfree(scipModel);
	}
}

void Generate_InitColumns_2(ColumnPool*& columnpool_Input, std::map<int, General_Block> block, int col_num,  SCIP* scip)

{
	SCIP* scipModel;
	SCIP_VAR** scipVars;
	SCIP_CONS** scipConss;

	int sp_num = block.size() - 1;
	for (int subproblem_no = 0; subproblem_no < sp_num; subproblem_no++)
	{
		for (int column_num = 0; column_num < col_num; column_num++)
		{
			//创建模型
			SCIPcreate(&scipModel);
			SCIPincludeDefaultPlugins(scipModel);
			SCIPcreateProbBasic(scipModel, "InitsubProb");
			//SCIPsetIntParam(scipModel, "misc/allowcontinue", 1);
			SCIPsetIntParam(scipModel, "limits/solutions", 100);
			SCIPsetBoolParam(scipModel, "constraints/countsols/collect", TRUE);
			//SCIPsetIntParam(scipModel, "misc/multistart/maxstarts", 10);
			SCIPsetIntParam(scipModel, "presolving/maxrounds", 0);
			unordered_map<int, SCIP_VAR*> scipVarMap;
			unordered_map<int, Vars> varDic;
			SCIP_VAR** scipVars;
			scipVars = new SCIP_VAR * [block[subproblem_no].bVars.size()];

			//设置目标
			SCIP_OBJSENSE objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MINIMIZE;
			if (block[subproblem_no].bobj.direction == Direction::max)
			{
				objSense = SCIP_OBJSENSE::SCIP_OBJSENSE_MAXIMIZE;
			}
			SCIPsetObjsense(scipModel, objSense);
			//添加变量
			for (int i = 0; i < block[subproblem_no].bVars.size(); i++)
			{
				Vars var_now = block[subproblem_no].bVars[i];
				if (column_num == 0)
				{
					columnpool_Input[subproblem_no].vars_id.push_back(var_now.Id);
				}
				SCIP_VAR* tempVar;
				VarType varType = var_now.Type;
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
				string VarName = var_now.Name;
				int lb_deal = var_now.Lb, ub_deal = var_now.Ub;
				if (var_now.Lb <  -10000000)
				{
					if (var_now.Ub < -1000)
					{
						lb_deal = var_now.Ub - 1000;
					}
					else
					{
						lb_deal = -1000;
					}

				}
				if (var_now.Ub > 1000000)
				{
					if (var_now.Lb > 1000)
					{
						ub_deal = var_now.Lb + 1000;
					}
					else
					{
						ub_deal = 1000;
					}
					
				}
				//SCIPcreateVarBasic(scipModel, &tempVar, VarName.c_str(), var_now.Lb, var_now.Ub, block[subproblem_no].bobj.coef[var_now.Id], scipVarType);//wqy可能有误1，目标函数系数
				SCIPcreateVarBasic(scipModel, &tempVar, VarName.c_str(), lb_deal, ub_deal, rand() % 21 - 10, scipVarType);
				SCIPaddVar(scipModel, tempVar);
				scipVarMap.insert(make_pair(var_now.Id, tempVar));
				varDic.insert(pair<int, Vars>(var_now.Id, var_now));
				scipVars[i] = tempVar;
			}
			//添加约束
			scipConss = new SCIP_CONS * [block[subproblem_no].bCons.size()];
			for (int i = 0; i < block[subproblem_no].bCons.size(); i++)
			{
				Constraints consTemp = block[subproblem_no].bCons[i];
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
				SCIPcreateConsBasicLinear(scipModel, &scipConss[i], consTemp.Name.c_str(), tempVarsSize, tempVars, tempConsCoefs, tempLb, tempUb);
				SCIPaddCons(scipModel, scipConss[i]);

				delete[] tempVars;
				delete[] tempConsCoefs;
			}

			SCIPwriteOrigProblem(scipModel, "D://gen_col.lp", NULL, false);

			SCIPsolve(scipModel);
			Column column_temp;

			//double obj = SCIPgetSolOrigObj(scipModel, sol);
			//printf("Solution %d: objective = %f\n", i + 1, obj);
			//SCIP_VAR** vars = SCIPgetVars(scipModel);
			int nvars = SCIPgetNVars(scipModel);
			for (int j = 0; j < nvars; ++j)
			{
				double val = SCIPgetSolVal(scipModel, SCIPgetBestSol(scipModel), scipVars[j]);
				//printf("  %s = %f\n", SCIPvarGetName(scipVars[j]), val);

				for (auto iter = varDic.begin(); iter != varDic.end(); iter++)
				{
					string name = scipVars[j]->name;
					if (name == iter->second.Name)
					{
						column_temp.varValue.insert(pair<int, double>(iter->second.Id, val));
						break;
					}
					if (++iter == varDic.end())
					{
						cout << "Error! Not Found!!" << endl;
					}
					iter--;

				}
			}
			column_temp.solSize = nvars;
			cout << columnpool_Input[subproblem_no].columns.size() << endl;
			columnpool_Input[subproblem_no].columns.insert(column_temp);
			cout << columnpool_Input[subproblem_no].columns.size() << endl;

			//赋值列池
			columnpool_Input[subproblem_no].var_num = block[subproblem_no].bVars.size();

			for (int i = 0; i < varDic.size(); i++)
			{
				SCIPreleaseVar(scipModel, &scipVars[i]);
			}
			for (size_t i = 0; i < block[subproblem_no].bCons.size(); i++)
			{
				SCIPreleaseCons(scipModel, &scipConss[i]);
			}
			SCIPfree(&scipModel);

			//新增0516 加入正常列测试 begin
			int test_flag = 1;
			if (test_flag)
			{
				if (column_num == 0)//只加一遍
				{
					Column column_correct;
					SCIPsolve(scip);
					int ncols = SCIPgetNVars(scip);
					SCIP_VAR** vars_ori = SCIPgetVars(scip);


					for (auto iter = varDic.begin(); iter != varDic.end(); iter++)
					{
						//cout << "t_ + iter->second.Name:" << "t_" + iter->second.Name << endl;
						for (int j = 0; j < ncols; j++)
						{
							string name = vars_ori[j]->name;
							//cout << "name:" << name << endl;

							if (name != "t_" + iter->second.Name)
							{
								if (j == ncols - 1)
								{
									cout << "Error2! Not Found!!" << endl;
								}
								continue;
							}
							double val = SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars_ori[j]);
							column_correct.varValue.insert(pair<int, double>(iter->second.Id, val));

							break;
						}


					}

					column_correct.solSize = ncols;
					cout << columnpool_Input[subproblem_no].columns.size() << endl;
					columnpool_Input[subproblem_no].columns.insert(column_correct);
					cout << columnpool_Input[subproblem_no].columns.size() << endl;

					//double obj = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
					////cout
					//for (size_t i = 0; i < ncols; i++)
					//{
					//	cout << vars_ori[i]->name << ":" << SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars_ori[i]) << endl;
					//}
				}
			}
			
			//新增0516 加入正常列测试 end
		}
		
		
    }

	



}