#pragma once
#include <iostream>
#include <scip/scipdefplugins.h>
#include <vector>
#include <cmath>
#include "ADMMFunc.h"
#include "primalBlock.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"


void UseADMM(SCIP* scipmodel, std::vector<SCIP_VAR*> scipvars, std::vector<DUALBlock> dualblock, std::string tempname, std::shared_ptr<spdlog::logger> spd,std::chrono::system_clock::time_point start, std::unordered_map<std::string, int> In_Var_Ni, std::map<std::string,writeJson::Soljson>& all_sol) {

		int num_vars = SCIPgetNVars(scipmodel);
		
		SCIP_VAR** varss = SCIPgetVars(scipmodel);
		
		std::map<std::string, double> init_value;
		for (int i = 0; i < num_vars; i++)
		{
			std::string name = SCIPvarGetName(varss[i]);
			init_value.insert(std::pair<std::string, double>(name, 0));
		}
		SCIPfreeTransform(scipmodel);
		double feasible_obj = 0;
		double offset = SCIPgetOrigObjoffset(scipmodel);

		NewIterUpdate(dualblock, init_value, feasible_obj, scipvars, offset, tempname, spd,start, In_Var_Ni, all_sol);
}
