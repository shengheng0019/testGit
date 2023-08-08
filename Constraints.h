#pragma once
#include"Vars.h"
#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>

using namespace std;

enum class PROP
{
	geq=0,
	leq=1,
	eq=2
};

template<>
struct std::hash<Vars>
{
	size_t operator()(const Vars& vars) const noexcept;
};


class Constraints
{
public:
	vector<double> Coef;
	double Coefficient;
	//Vars* ContainedVars;
	unordered_map<Vars, double> exprDic;
	unordered_map<int, double> idDic;
	int Id;
	//double Lb;
	string Name;
	PROP Prop;
	//double Ub;
	double rhs;	//right hand side

	Constraints(const Constraints& cons_);
	Constraints& operator = (const Constraints& cons);
	Constraints();
	//~Constraints();
};



