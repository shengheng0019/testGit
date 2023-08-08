#pragma once
#include "ADMMFunc.h"
#include "DUALBlock.h"
#include "Partition.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"



void UseADMM(SCIP* scipmodel, std::vector<SCIP_VAR*> scipvars, std::vector<DUALBlock> dualblock, std::string tempname, std::shared_ptr<spdlog::logger> spd,std::chrono::system_clock::time_point start, std::unordered_map<std::string, int> In_Var_Ni, std::map<std::string,writeJson::Soljson>& all_sol);
