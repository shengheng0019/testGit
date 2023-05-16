#pragma once
#include<iostream>
#include<string>

using namespace std;

enum class VarType
{
	Int = 0,
	Num = 1,
	Bool = 2
};

class Vars
{
public:
	int Id;
	double Lb;
	string Name;
	VarType Type;
	double Ub;
	Vars& operator = (const Vars& var);

private:
};

bool operator==(const Vars& var1, const Vars& var2);
//{
//	if (var1.Id == var2.Id)
//		return true;
//	return false;
//}
