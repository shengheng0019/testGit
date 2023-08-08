#pragma once
#include"Vars.h"
#include<iostream>
#include<string>
#include<unordered_map>
#include<vector>

enum class Direction
{
	min=0,
	max=1
};

class Objective
{
public: 
	std::unordered_map<int, double> coef;	//key: vars.Id, value: coef of each var
	Vars* containedVars;
	Direction direction;
	int id;
	double lb;
	double ub;

	Objective(const Objective& obj_);
	Objective();
	//~Objective();
};

