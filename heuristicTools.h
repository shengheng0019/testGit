#pragma once
#include <cmath>
#include<iostream>
#include<vector>
#include"ColumnPool.h"



constexpr double EPSFORINT = 1e-6;

double chgInt(double inaccuracy);

struct FIX //one FIX corresponding to one node
{
	std::vector<std::vector<std::pair<int, double>>>fixedColumnValueList; //list of fixed columns
	std::vector<int>spNoList; //list of fix columns' spno, correspending to "fixedColumnValueList"
	std::vector<bool>spIsFixed;   //size equal to the number of subproblem, 1 represent this subproblem is fixed, otherwise 0.
	std::vector<ColumnPool> recordColumnpool; //record all the columnpool at that node, and the subproblem corresponding to the newest fixed column only leave one column
};