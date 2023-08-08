#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Lifting.hpp"
#include "MIR.hpp"
#include "GMI.hpp"

using namespace std;

double eps = 1e-8;

//int DynamicProgramming(vector<int> weight, vector<int> value, int bagWeight)
//{
//	// 初始化
//	vector<int> dp(bagWeight + 1, 0);
//	for (int i = 0; i < weight.size(); i++) 
//	{ // 遍历物品
//		for (int j = bagWeight; j >= weight[i]; j--) 
//		{ // 遍历背包容量
//			dp[j] = max(dp[j], dp[j - weight[i]] + value[i]);
//		}
//	}
//	//cout << dp[bagWeight] << endl;
//
//	//这地方需要释放vector数组?
//	vector<int>().swap(dp);
//	return dp[bagWeight];
//
//}
//
//vector<int> find_cover(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
//{
//	vector<int> minimal_cover;//返回变量Id
//	unordered_map<int, double> cover_set;	//<exprDic.int(Id), idDic.double>
//	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
//	{
//		cover_set.emplace(iter->first.Id, iter->second);
//	}
//	vector<pair<int, double>> order_vector;
//	for (auto it : cover_set) 
//	{
//		order_vector.push_back(it);
//	}
//	sort(order_vector.begin(), order_vector.end(), [](pair<int, double>& a, pair<int, double>& b) {return a.second > b.second; });
//	//此时order_vector为排序后的集合
//	double efficient_sum = 0;
//	int cover_num = 0;	//覆盖集中包含变量的个数
//	int over_num = 0;	//记录系数大于右端项的个数 正常来说=0
//	int start_number = 8;//从第几个元素开始取最小覆盖集
//	for (int i = start_number; i < order_vector.size(); i++)
//	{
//		if (order_vector[i].second > rhs)//系数大于右端项，直接略过并生成当前变量=0的约束？？？？？？？？？？？？？？
//		{
//			//此处需要添加生成的当前变量=0的约束！！！！！！！！！！！！！！1
//			over_num += 1;
//			continue;
//		}
//		else
//		{
//			cover_num += 1;
//			efficient_sum += order_vector[i].second;
//			if (efficient_sum > rhs)
//			{
//				for (int j = 0; j < cover_num; j++)
//				{
//					minimal_cover.push_back(order_vector[i - cover_num + j + 1].first);
//				}
//
//				break;
//			}
//		}
//
//		
//	}
//	
//	return minimal_cover;
//}
//
//
//pair< std::unordered_map<Vars, double>, double> Lifting(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
//{
//	//若01变量系数<0，则将其取正，并在右端项+同样的数
//	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
//	{
//		if (iter->second < 0)
//		{
//			iter->second = -iter->second;
//			rhs += iter->second;
//		}
//	}
//
//	vector<int> minimal_cover;
//	minimal_cover = find_cover(exprDic, idDic, rhs);	//cover内包含变量Id的vector
//	//如果minimla cover之后不止一个，那么将会用二维vector进行存储，下边算法中应该变成两层for循环
//
//	vector<int> out_cover;								//cover内不包含变量Id的vector
//	vector<int> alpha;
//	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
//	{
//		if (!(std::find(minimal_cover.begin(), minimal_cover.end(), iter->first.Id) != minimal_cover.end()))
//		{
//			out_cover.push_back(iter->first.Id);
//		}
//	}
//
//
//
//	for (int t = 0; t < out_cover.size(); t++)//若find_cover函数中的over_num>0则这地方需要修改
//	{
//		
//
//		vector<int> weight;
//		vector<int> value;
//
//		Vars a;
//		double b = 0;
//		
//		for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
//		{
//			if (iter->first.Id == out_cover[t])
//			{
//				a = iter->first;
//				break;
//			}
//		}
//		b = exprDic[a];
//
//		for (int k = 0; k < alpha.size(); k++)
//		{
//			for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
//			{
//				if (iter->first.Id == out_cover[k])
//				{
//					weight.push_back(iter->second);
//					break;
//				}
//			}
//			value.push_back(alpha[k]);
//		}
//
//		for (int i = t; i < t + minimal_cover.size(); i++)
//		{
//			for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
//			{
//				if (iter->first.Id == minimal_cover[i - t])
//				{
//					weight.push_back(iter->second);
//					break;
//				}
//			}
//			value.push_back(1);
//		}
//
//		double bagWeight = rhs - b;//是否存在类型转换的问题？
//		int z = DynamicProgramming(weight, value, bagWeight);
//		alpha.push_back(minimal_cover.size() - 1 - z);
//
//		//delete weight;
//
//	}
//
//	cout << "alpha vector:";
//	for (int t = 0; t < out_cover.size(); t++)
//	{
//		cout << alpha[t] << "  ";
//	}
//	cout << endl;
//
//	unordered_map<Vars, double> exprDic_add = exprDic;//不改变量，只改系数作为新约束
//	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
//	{
//
//		for (int i = 0; i < minimal_cover.size(); i++)
//		{
//			if (iter->first.Id == minimal_cover[i])
//			{
//				iter->second = 1;
//				break;
//			}
//		}
//		for (int i = 0; i < out_cover.size(); i++)
//		{
//			if (iter->first.Id == out_cover[i])
//			{
//				iter->second = alpha[i];
//				break;
//			}
//		}
//
//	}
//
//
//	return make_pair(exprDic_add, minimal_cover.size());
//}


inline std::map<int, General_Block> GenerateCut(std::map<int, General_Block>& GeneralBlocks)
{
	int CoverCut_num = 0;
	int MIRCut_num = 0;
	int GMICut_num= 0;
	/*std::map<int, General_Block> cut_block;
	cut_block = GeneralBlocks;*/
	cout << "GenerateCut--Blocks1:" << endl;
	for (auto& block : GeneralBlocks)
	{
		int bCons_origin_size = block.second.bCons.size();
		cout << "index:" << block.first << endl;
		cout << "block.cons.exprDic:" << endl;
		cout << "block.cons.bCons_origin_size:" << bCons_origin_size << endl;
		for (int index = 0; index < bCons_origin_size; index++)
		{
			int int_num = 0;
			int cont_num = 0;
			int binary_num = 0;

			int int_num_positive = 0;
			int cont_num_positive = 0;
			int binary_num_positive = 0;

			int cont_negative = 0;
			bool rhs_float = 0;
			for (auto iter = block.second.bCons[index].exprDic.begin(); iter != block.second.bCons[index].exprDic.end(); ++iter)
			{
				/*cout << iter->first.Id << " ";
				cout << iter->first.Name << " ";
				cout << static_cast<int>(iter->first.Type) << " ";*/
				if (static_cast<int>(iter->first.Type) == 0)
				{
					int_num += 1;
					if (iter->first.Lb > -eps)//iter->first.Lb >= 0
					{
						int_num_positive += 1;
					}

				}
				if (static_cast<int>(iter->first.Type) == 1)
				{
					cont_num += 1;
					if (iter->first.Lb > -eps)//iter->first.Lb >= 0
					{
						cont_num_positive += 1;
					}
				}
				if (static_cast<int>(iter->first.Type) == 2)
				{
					binary_num += 1;
					if (iter->first.Lb > -eps)//iter->first.Lb >= 0
					{
						binary_num_positive += 1;
					}
				}

			}

			if (block.second.bCons[index].rhs - floor(block.second.bCons[index].rhs) != 0)
			{
				rhs_float = 1;
			}
			

			if (binary_num == block.second.bCons[index].exprDic.size() && (block.second.bCons[index].Prop) == PROP::leq && binary_num >= 8)//这地方正常是等于！！！！！！！！！！！！
			{
				pair <std::unordered_map<Vars, double>, double> exprDic_CoverCut;
				exprDic_CoverCut = Lifting(block.second.bCons[index].exprDic, block.second.bCons[index].idDic, block.second.bCons[index].rhs);
				Constraints bCons_add;
				bCons_add = block.second.bCons[index];
				bCons_add.exprDic = exprDic_CoverCut.first;
				bCons_add.rhs = exprDic_CoverCut.second;
				bCons_add.Prop = PROP::leq;
				block.second.bCons.push_back(bCons_add);
				//block.second.bCons[index].exprDic.insert(exprDic_CoverCut);
				CoverCut_num += 1;
			}

			if ((int_num + binary_num) >= 1 && block.second.bCons[index].Prop == PROP::leq && rhs_float && int_num == int_num_positive &&
				cont_num == cont_num_positive && binary_num == binary_num_positive && (int_num + binary_num) >= 2)
			{
				pair <std::unordered_map<Vars, double>, double> exprDic_MIRCut;
				exprDic_MIRCut = MIR(block.second.bCons[index].exprDic, block.second.bCons[index].idDic, block.second.bCons[index].rhs);
				Constraints bCons_add;
				bCons_add = block.second.bCons[index];
				bCons_add.exprDic = exprDic_MIRCut.first;
				bCons_add.rhs = exprDic_MIRCut.second;
				bCons_add.Prop = PROP::leq;
				block.second.bCons.push_back(bCons_add);
				MIRCut_num += 1;

			}
			if ((int_num + binary_num) >= 1 && block.second.bCons[index].Prop == PROP::eq && rhs_float && int_num == int_num_positive &&
				cont_num == cont_num_positive && binary_num == binary_num_positive)
			{
				pair <std::unordered_map<Vars, double>, double> exprDic_GMICut;
				exprDic_GMICut = GMI(block.second.bCons[index].exprDic, block.second.bCons[index].idDic, block.second.bCons[index].rhs);
				Constraints bCons_add;
				bCons_add = block.second.bCons[index];
				bCons_add.exprDic = exprDic_GMICut.first;
				bCons_add.rhs = exprDic_GMICut.second;
				bCons_add.Prop = PROP::geq;
				block.second.bCons.push_back(bCons_add);
				GMICut_num += 1;
			}


			//cout << endl;


		}
		//cout << endl;
	}

	cout << "产生lifting cut的个数为：" << CoverCut_num << endl;
	cout << "产生MIR cut的个数为：" << MIRCut_num << endl;
	cout << "产生GMI cut的个数为：" << GMICut_num << endl;


	return GeneralBlocks;

}


