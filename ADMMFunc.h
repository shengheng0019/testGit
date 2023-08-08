#pragma once
#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <scip/scipdefplugins.h>
#include "DUALBlock.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "SolJson.hpp"



std::unordered_map<std::string, double> get_obj_value_vector(std::vector<SCIP_VAR*> ori_vars);

std::vector<std::unordered_map<std::string, double>> get_con_coef_vector(SCIP* ori_scip, std::vector<SCIP_CONS*> ori_cons);

std::vector<std::vector<double>> get_left_right_side_vector(SCIP* scip, std::vector<SCIP_CONS*> cons);

void TryFeasXdot(int block_num, std::vector<std::map<std::string, double>>& iter_z_vars, std::vector<std::map<std::string, double>>& final_x_vars,
    std::vector<DUALBlock> dual_blocks, std::vector<std::map<std::string, double>> multiplier, double rho,
    std::map<std::string, double> ori_vars, std::vector<double>& iter_obj, int& all_feas_flag, int& count_p);

void TryFeasFP(std::unordered_map<std::string, double> Sol, std::string tempname, int& feas_flag);

void NewIterUpdate(std::vector<DUALBlock> dual_blocks, std::map<std::string, double> init_value, double feasible_obj, std::vector<SCIP_VAR*> scipvars, double offset, std::string tempname, std::shared_ptr<spdlog::logger> spd,std::chrono::system_clock::time_point start, std::unordered_map<std::string, int> In_Var_Ni, std::map<std::string,writeJson::Soljson>& all_sol);

