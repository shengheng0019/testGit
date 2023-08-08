#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

int DynamicProgramming(vector<double> weight, vector<double> value, double bagWeight)
{
	cout << "bagWeight=" << bagWeight << endl;
	////��ʼ��
	//vector<int> dp(bagWeight + 1, 0);
	//for (int i = 0; i < weight.size(); i++) 
	//{ // ������Ʒ
	//	for (int j = bagWeight; j >= weight[i]; j--) 
	//	{ // ������������
	//		dp[j] = max(dp[j], dp[j - weight[i]] + value[i]);
	//	}
	//}
	
	int pack_size = weight.size();

	// ��ά����
	vector<vector<double>> dp(weight.size(), vector<double>(bagWeight + 1, 0));

	// ��ʼ��
	for (int j = weight[0]; j <= bagWeight; j++) {
		dp[0][j] = value[0];
	}

	// weight����Ĵ�С ������Ʒ����
	for (int i = 1; i < weight.size(); i++) { // ������Ʒ
		for (int j = 0; j <= bagWeight; j++) { // ������������
			if (j < weight[i]) dp[i][j] = dp[i - 1][j];
			else dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - weight[i]] + value[i]);

		}
	}

	cout << dp[weight.size() - 1][bagWeight] << endl;

	//��ط���Ҫ�ͷ�vector����?
	//vector<double>().swap(dp);
	return dp[weight.size() - 1][bagWeight];

}

vector<int> find_cover(std::unordered_map<Vars, double> exprDic, unordered_map<int, double> idDic, double rhs)
{
	vector<int> minimal_cover;//���ر���Id
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
	//��ʱorder_vectorΪ�����ļ���
	double efficient_sum = 0;
	int cover_num = 0;	//���Ǽ��а��������ĸ���
	int over_num = 0;	//��¼ϵ�������Ҷ���ĸ��� ������˵=0
	int start_number = 0;//�ӵڼ���Ԫ�ؿ�ʼȡ��С���Ǽ�
	for (int i = start_number; i < order_vector.size(); i++)
	{
		if (order_vector[i].second > rhs)//ϵ�������Ҷ��ֱ���Թ������ɵ�ǰ����=0��Լ������������������������������
		{
			//�˴���Ҫ������ɵĵ�ǰ����=0��Լ������������������������������1
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
	//��01����ϵ��<0������ȡ���������Ҷ���+ͬ������
	for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	{
		if (iter->second < 0)
		{
			iter->second = -iter->second;
			rhs += iter->second;
		}
	}

	vector<int> minimal_cover;
	minimal_cover = find_cover(exprDic, idDic, rhs);	//cover�ڰ�������Id��vector
	//���minimla cover֮��ֹһ������ô�����ö�άvector���д洢���±��㷨��Ӧ�ñ������forѭ��


	//20230711ע��
	//vector<int> out_cover;								//cover�ڲ���������Id��vector
	//vector<int> alpha;
	//for (auto iter = exprDic.begin(); iter != exprDic.end(); ++iter)
	//{
	//	if (!(std::find(minimal_cover.begin(), minimal_cover.end(), iter->first.Id) != minimal_cover.end()))
	//	{
	//		out_cover.push_back(iter->first.Id);
	//	}
	//}


	//for (int t = 0; t < out_cover.size(); t++)//��find_cover�����е�over_num>0����ط���Ҫ�޸�
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

	//	double bagWeight = rhs - b;//�Ƿ��������ת�������⣿
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




	unordered_map<Vars, double> exprDic_add = exprDic;//���ı�����ֻ��ϵ����Ϊ��Լ��
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

		if (!cover_flag)//��ǰ������cover��Ԫ��
		{
			iter->second = 0;
		}

		//20230711ע��
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



