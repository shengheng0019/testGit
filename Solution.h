#pragma once
#include<vector>
#include<unordered_map>
#include<iostream>

using namespace std;

//ostream& operator << (ostream& out, Solution& sol);

class Solution
{
public:
	unordered_map<int, double> varValue;
	int solSize;
	double objValue;
	//Solution(int solSize);
	//~Solution();
	Solution(const Solution& sol);
	Solution operator=(const Solution& sol);
	Solution();
	friend ostream& operator << (ostream& out, Solution& sol);
};

//ostream& operator << (ostream& out, Solution& sol)
//{
//	for (auto iter = sol.varValue.begin(); iter != sol.varValue.begin(); iter++)
//	{
//		out << iter->second;
//	}
//	return out;
//}

