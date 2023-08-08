#pragma once
#include"ColumnPool.h"
#include<Map>
#include"GeneralBlock.h"
#include "scip/scipdefplugins.h"


int GenerateInitColumns3(std::vector<ColumnPool>& columnpool_Input, std::map<int, General_Block> block, int col_num, SCIP* scip, double& initScipFeasObj);