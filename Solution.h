#pragma once
#include<vector>
#include<unordered_map>

class Solution
{
public:
	std::unordered_map<int, double> varValue;
	int solSize;
	double objValue;
	Solution(const Solution& sol);
	Solution operator=(const Solution& sol);
	Solution();
};