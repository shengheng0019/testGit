#pragma once
#include<vector>
#include"Vars.h"
#include"Constraints.h"
#include"Objective.h"

class General_Block
{
public:
    std::vector<Constraints> bCons;
    std::vector<Vars> bVars;
    Objective bobj;
    General_Block();
};

