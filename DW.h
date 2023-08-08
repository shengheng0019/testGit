#ifndef _DW_
#define _DW_

#include"Node.h"
#include "SolJson.hpp"
//#include "ThreadPool.h"
#include <iostream>
#include "fstream"
#include "chrono"
#include<stack>
//#include "scip/scipdefplugins.h"
/*#include <psapi.h>
#pragma comment(lib,"psapi.lib")*/
#define EPS 1e-4

//using namespace std;
class DW
{
public:
	DW() = default;
	~DW() = default;

	static void LoadModel(Node& node, std::map<int, General_Block> block);
	static void UpdateModel(Node& node);
	static void Diving_new(const vector<ColumnPool>& columnpool, std::map<int, General_Block> block, Objective master_obj, SCIP* scipmodel,std::map<std::string, writeJson::Soljson>& all_solsmap, std::string part_name, std::chrono::system_clock::time_point start);
	static int ColumnGeneration(Node& node, bool root_flag);
	static int ColumnGenerationNew(Node& node, bool root_flag);
	static bool Check_All_Cons(std::map<int, General_Block>block, vector<double>col_val);
	static bool Check_All_Cons(std::map<int, General_Block>block, vector<pair<string, double>>colNameString);
	static bool Check_Part_Cons(std::map<int, General_Block>block, unordered_map<int, double>var_value, int block_no);

private:


};

//function declaration
int FindFixedInformation(const Node& node, std::string columnName, FIX& fixedInfoCurrent);
void SortNodeColumn(const Node& node, std::vector<std::tuple<string, int, double>>& columnNameIndexTransolution);
/*void UpdateNode(const Node &savedNode, Node &node, const std::vector<FIX> &fixedInfos);*/
void UpdateNodeNew(Node& node, const FIX& fixedInfos);
void UpdateSavedNodeColumnpol(Node& savedNode, const Node& node);
void CalculateFixedInfo(std::stack<FIX>& fixedInformationList, const Node& node, int maxDiscrepancy);
int FindFeasSol(const Node& node, SCIP* scipmodel, std::map<int, General_Block> block, double& obj); //find feasible solution based on all of the integer variable fixed


#endif // !_DW_




