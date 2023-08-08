#pragma once
#include<iostream>
#include<queue>
#include<unordered_set>
#include <map>
#include<vector>
#include"Constraints.h"
#include"Vars.h"
#include "GeneralBlock.h"
#include <scip/scipdefplugins.h>
#include<chrono>

#define INF 99999999


void SetVarsBounds(std::unordered_map<int, Vars>& vars, const std::vector<Constraints>& conss);

std::pair<std::unordered_map<int, Vars>, std::vector<Constraints>> getAllVarsandCons(std::map<int, General_Block>& blocks);
