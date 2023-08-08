#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
//double eps = 1e-8;

pair< std::unordered_map<Vars, double>, double> MIR(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
{
	//double cont_coefficient = 0;
	//for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	//{
	//	if (static_cast<int>(iter->first.Type) == 1)
	//	{
	//		cont_coefficient = iter->second;//获取连续变量的系数
	//		break;
	//	}
	//}
	double eps = 1e-8;
	unordered_map<Vars, double> exprDic_add = exprDic;//不改变量，只改系数作为新约束
	double f0 = rhs - floor(rhs);
	vector<double> fi;
	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		if (static_cast<int>(iter->first.Type) != 1)
		{
			fi.push_back(iter->second - floor(iter->second));	//fi向量记录整数变量的分数差
		}
	}
	
	int fi_count = 0;

	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		if (static_cast<int>(iter->first.Type) == 1)
		{
			if (iter->second < 0)
			{
				iter->second = (1 / (1 - f0)) * iter->second;//获取连续变量的系数
			}
			else
			{
				iter->second = 0;
			}
			
		}
		else
		{
			iter->second = (floor(iter->second) + (((fi[fi_count] - f0) >= 0) ? (fi[fi_count] - f0) : 0) / (1 - f0));//(floor(iter->second) + (((fi[fi_count] - f0) >= 0) ? (fi[fi_count] - f0) : 0) / (1 - f0))

			fi_count += 1;
		}
	}
	//rhs = floor(rhs);
	double rhs_temp = floor(rhs);


	for (auto it = exprDic_add.begin(); it != exprDic_add.end(); )
	{
		if (it->second == 0)
		{ 
			exprDic_add.erase(it++);
		}
		else 
		{ 
			++it;
		}
	}

	return make_pair(exprDic_add, rhs_temp);
}



