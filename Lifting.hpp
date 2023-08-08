#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

int DynamicProgramming(vector<double> weight, vector<double> value, double bagWeight)
{
	cout << "bagWeight=" << bagWeight << endl;
	////初始化
	//vector<int> dp(bagWeight + 1, 0);
	//for (int i = 0; i < weight.size(); i++) 
	//{ // 遍历物品
	//	for (int j = bagWeight; j >= weight[i]; j--) 
	//	{ // 遍历背包容量
	//		dp[j] = max(dp[j], dp[j - weight[i]] + value[i]);
	//	}
	//}
	
	int pack_size = weight.size();

	// 二维数组
	vector<vector<double>> dp(weight.size(), vector<double>(bagWeight + 1, 0));

	// 初始化
	for (int j = weight[0]; j <= bagWeight; j++) {
		dp[0][j] = value[0];
	}

	// weight数组的大小 就是物品个数
	for (int i = 1; i < weight.size(); i++) { // 遍历物品
		for (int j = 0; j <= bagWeight; j++) { // 遍历背包容量
			if (j < weight[i]) dp[i][j] = dp[i - 1][j];
			else dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - weight[i]] + value[i]);

		}
	}

	cout << dp[weight.size() - 1][bagWeight] << endl;

	//这地方需要释放vector数组?
	//vector<double>().swap(dp);
	return dp[weight.size() - 1][bagWeight];

}

vector<int> find_cover(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
{
	vector<int> minimal_cover;//返回变量Id
	unordered_map<int, double> cover_set;	//<exprDic.int(Id), idDic.double>
	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	{
		cover_set.emplace(iter->first.Id, iter->second);
	}
	vector<pair<int, double>> order_vector;
	for (auto it : cover_set) 
	{
		order_vector.push_back(it);
	}
	sort(order_vector.begin(), order_vector.end(), [](pair<int, double>& a, pair<int, double>& b) {return a.second > b.second; });
	//此时order_vector为排序后的集合
	double efficient_sum = 0;
	int cover_num = 0;	//覆盖集中包含变量的个数
	int over_num = 0;	//记录系数大于右端项的个数 正常来说=0
	int start_number = 0;//从第几个元素开始取最小覆盖集
	for (int i = start_number; i < order_vector.size(); i++)
	{
		if (order_vector[i].second > rhs)//系数大于右端项，直接略过并生成当前变量=0的约束？？？？？？？？？？？？？？
		{
			//此处需要添加生成的当前变量=0的约束！！！！！！！！！！！！！！1
			over_num += 1;
			continue;
		}
		else
		{
			cover_num += 1;
			efficient_sum += order_vector[i].second;
			if (efficient_sum > rhs)
			{
				for (int j = 0; j < cover_num; j++)
				{
					minimal_cover.push_back(order_vector[i - cover_num + j + 1].first);
				}

				break;
			}
		}

		
	}
	
	return minimal_cover;
}


pair< std::unordered_map<Vars, double>, double> Lifting(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
{
	//若01变量系数<0，则将其取正，并在右端项+同样的数
	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	{
		if (iter->second < 0)
		{
			iter->second = -iter->second;
			rhs += iter->second;
		}
	}

	vector<int> minimal_cover;
	minimal_cover = find_cover(exprDic, idDic, rhs);	//cover内包含变量Id的vector
	//如果minimla cover之后不止一个，那么将会用二维vector进行存储，下边算法中应该变成两层for循环


	//20230711注释
	//vector<int> out_cover;								//cover内不包含变量Id的vector
	//vector<int> alpha;
	//for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	//{
	//	if (!(std::find(minimal_cover.begin(), minimal_cover.end(), iter->first.Id) != minimal_cover.end()))
	//	{
	//		out_cover.push_back(iter->first.Id);
	//	}
	//}


	//for (int t = 0; t < out_cover.size(); t++)//若find_cover函数中的over_num>0则这地方需要修改
	//{
	//	

	//	vector<double> weight;
	//	vector<double> value;

	//	Vars a;
	//	double b = 0;
	//	
	//	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	//	{
	//		if (iter->first.Id == out_cover[t])
	//		{
	//			a = iter->first;
	//			break;
	//		}
	//	}
	//	b = exprDic[a];

	//	for (int k = 0; k < alpha.size(); k++)
	//	{
	//		for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	//		{
	//			if (iter->first.Id == out_cover[k])
	//			{
	//				weight.push_back(iter->second);
	//				break;
	//			}
	//		}
	//		value.push_back(alpha[k]);
	//	}

	//	for (int i = t; i < t + minimal_cover.size(); i++)
	//	{
	//		for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	//		{
	//			if (iter->first.Id == minimal_cover[i - t])
	//			{
	//				weight.push_back(iter->second);
	//				break;
	//			}
	//		}
	//		value.push_back(1);
	//	}

	//	double bagWeight = rhs - b;//是否存在类型转换的问题？
	//	int z = DynamicProgramming(weight, value, bagWeight);
	//	alpha.push_back(minimal_cover.size() - 1 - z);

	//	//delete weight;

	//}

	//cout << "alpha vector:";
	//for (int t = 0; t < out_cover.size(); t++)
	//{
	//	cout << alpha[t] << "  ";
	//}
	//cout << endl;




	unordered_map<Vars, double> exprDic_add = exprDic;//不改变量，只改系数作为新约束
	for (auto iter = exprDic_add.begin(); iter != exprDic_add.end(); ++iter)
	{
		bool cover_flag = false;
		for (int i = 0; i < minimal_cover.size(); i++)
		{
			if (iter->first.Id == minimal_cover[i])
			{
				iter->second = 1;
				cover_flag = true;
				break;
			}
		}

		if (!cover_flag)//当前变量非cover内元素
		{
			iter->second = 0;
		}

		//20230711注释
		/*for (int i = 0; i < out_cover.size(); i++)
		{
			if (iter->first.Id == out_cover[i])
			{
				iter->second = alpha[i];
				break;
			}
		}*/

	}

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


	return make_pair(exprDic_add, minimal_cover.size() - 1);
}



