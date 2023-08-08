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
	//		cont_coefficient = iter->second;//��ȡ����������ϵ��
	//		break;
	//	}
	//}
	double eps = 1e-8;
	unordered_map<Vars, double> exprDic_add = exprDic;//���ı�����ֻ��ϵ����Ϊ��Լ��
	double f0 = rhs - floor(rhs);
	vector<double> fi;
	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		if (static_cast<int>(iter->first.Type) != 1)
		{
			fi.push_back(iter->second - floor(iter->second));	//fi������¼���������ķ�����
		}
	}
	
	int fi_count = 0;

	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		if (static_cast<int>(iter->first.Type) == 1)
		{
			if (iter->second < 0)
			{
				iter->second = (1 / (1 - f0)) * iter->second;//��ȡ����������ϵ��
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



