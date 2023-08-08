#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

pair< std::unordered_map<Vars, double>, double> GMI(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
{
	unordered_map<Vars, double> exprDic_add = exprDic;//不改变量，只改系数作为新约束
	double f0 = rhs - floor(rhs);
	vector<double> fi;
	//vector<double> gi;
	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		if (static_cast<int>(iter->first.Type) != 1)
		{
			fi.push_back(iter->second - floor(iter->second));	//fi向量记录整数变量的分数差
		}
		//if (static_cast<int>(iter->first.Type) == 1)
		//{
		//	gi.push_back(iter->second - floor(iter->second));	//gi向量记录连续变量的分数差
		//}
	}

	int fi_count = 0;
	//int gi_count = 0;

	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		if (static_cast<int>(iter->first.Type) != 1)
		{
			if (fi[fi_count] <= f0)
			{
				iter->second = fi[fi_count] / f0;
			}
			else
			{
				iter->second = (1 - fi[fi_count]) / (1 - f0);
			}

			fi_count += 1;
			
			

		}
		else
		{
			if (iter->second >= 0)
			{
				iter->second = iter->second / f0;//获取连续变量的系数
			}
			else
			{
				iter->second = -iter->second / (1 - f0);
			}



			/*if (gi[gi_count] >= 0)
			{
				iter->second = gi[gi_count] / f0;
			}
			else
			{
				iter->second = -gi[gi_count] / (1 - f0);
			}
			gi_count += 1;*/
		}
	}
	double rhs_temp = 1;

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
	
	//在generate_cut返回后，改成了大于等于1
	return make_pair(exprDic_add, rhs_temp);
}

