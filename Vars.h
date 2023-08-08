#pragma once
#include<iostream>
#include<string>



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
    std::string Name;
    VarType Type;
    double Ub;
    Vars(const Vars& var_);
    Vars();
    Vars operator=(const Vars& var1);

private:
};

bool operator==(const Vars& var1, const Vars& var2);
