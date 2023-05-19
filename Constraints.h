#pragma once
#include"Vars.h"
#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>

using namespace std;

enum PROP
{
	geq = 0,
	leq = 1,
	eq = 2
};

template<>
struct hash<Vars>
{
	size_t operator()(const Vars& vars) const noexcept
	{
		return hash<int>()(vars.Id);
	}
};



class Constraints
{
public:
	//vector<double> Coef;
	//double Coefficient;
	//Vars* ContainedVars;
	unordered_map<Vars, double> exprDic;
	unordered_map<int, double> idDic;//int->var.id
	//int Id;
	//double Lb;
	string Name;
	PROP Prop;
	//double Ub;
	double rhs;

	Constraints& operator = (const Constraints& cons);
};
