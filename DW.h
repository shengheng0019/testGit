#pragma once
#include"Node.h"
#include<map>
#include"primalBlock.hpp"
#include "ThreadPool.h"
#include <windows.h>
#include <iostream>
#include <psapi.h>
#pragma comment(lib,"psapi.lib")
#define EPS 1e-4

//using namespace std;
class DW
{
public:
	DW();
	~DW();

	static void LoadModel(Node& node, std::map<int, General_Block> block);
	static void UpdateModel(Node& node, std::map<int, General_Block>block);
	static void Diving(ColumnPool* columnpool, std::map<int,General_Block> block, Objective master_obj);
	static void Get_InitColumns(ColumnPool*& columnpool);
	static int ColumnGeneration(Node& node);
	static int Fixed_Column(Node& node,int no);
	static bool Test_feasible(Node node);
	static bool Check_All_Cons(std::map<int, General_Block>block, vector<double>col_val);
	static bool Check_Part_Cons(std::map<int, General_Block>block, unordered_map<int, double>var_value, int block_no);
	static int Get_Tabu_list(Node node, vector<vector<unordered_map<int, double>>>& tabu_list);

private:


};


