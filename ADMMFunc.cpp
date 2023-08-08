#pragma once

#include <fstream>
#include "ADMMFunc.h"
#include "DUALBlock.h"
#include "FP.h"
#include "feasibilityjump.hh"
#include "RINS.hpp"

#define SET_NAME_x(name) x##name
#define SET_NAME_z(name) z##name
#define SET_NAME_ori_z(name) ori_z##name
#define SET_NAME_cons(name) cons##name

std::map<std::string, int> name_to_index;

std::vector<double> sum1_record;
std::vector<double> sum2_record;
std::vector<double> rho_record;

extern std::string Folder;

//Obtain the coefficient of the objective function corresponding to the variable
std::unordered_map<std::string, double> get_obj_value_vector(std::vector<SCIP_VAR*> ori_vars)
{
    std::unordered_map<std::string, double> obj_values;
    for (int i = 0; i < ori_vars.size(); i++)
    {
        obj_values.insert(std::pair<std::string, double>(SCIPvarGetName(ori_vars[i]), SCIPvarGetObj(ori_vars[i])));
    }
    return obj_values;
}

//Obtain the coefficients in the constraint matrix corresponding to the variables
std::vector<std::unordered_map<std::string, double>> get_con_coef_vector(SCIP* ori_scip, std::vector<SCIP_CONS*> ori_cons)
{
    std::vector<std::unordered_map<std::string, double>> cons_coefs;
    for (int j = 0; j < ori_cons.size(); j++)
    {
        std::unordered_map<std::string, double> con_coefs;
        double* coef = SCIPgetValsLinear(ori_scip, ori_cons[j]);
        int nnonz = SCIPgetNVarsLinear(ori_scip, ori_cons[j]);
        SCIP_VAR** vars = SCIPgetVarsLinear(ori_scip, ori_cons[j]);

        for (int i = 0; i < nnonz; i++)
        {
            con_coefs.insert(std::pair<std::string, double>(SCIPvarGetName(vars[i]), coef[i]));
        }
        cons_coefs.push_back(con_coefs);
    }
    return cons_coefs;
}

//Gets the left - and right-end entries corresponding to the constraint
std::vector<std::vector<double>> get_left_right_side_vector(SCIP* ori_scip, std::vector<SCIP_CONS*> ori_cons)
{
    std::vector<std::vector<double>> left_right_side;

    unsigned int success = 1;

    for (int j = 0; j < ori_cons.size(); j++)
    {
        std::vector<double> left_right_side_inner;
        double left_side = SCIPconsGetLhs(ori_scip, ori_cons[j], &success);
        double right_side = SCIPconsGetRhs(ori_scip, ori_cons[j], &success);
        left_right_side_inner.push_back(left_side);
        left_right_side_inner.push_back(right_side);
        left_right_side.push_back(left_right_side_inner);
    }
    return left_right_side;
}

//A subproblem model for splitting dually coupled variables by block is established
void BuildProb_Xdot(std::vector<SCIP_VAR*> x_vars, std::vector<SCIP_VAR*> z_vars, std::vector<SCIP_CONS*> block_cons, std::map<std::string,double>  multiplier,
                    double rho, std::map<std::string, double> ori_vars, std::map<std::string,double>& iter_z_vars, std::map<std::string, double>& final_x_vars, double& obj)
{
    std::unordered_map<std::string, double> coefx = get_obj_value_vector(x_vars);
    std::unordered_map<std::string, double> coefz = get_obj_value_vector(z_vars);

    SCIP* _scip1;
    SCIPcreate(&_scip1);

    // include default plugins
    SCIPincludeDefaultPlugins(_scip1);

    // creates empty problem and initializes all solving data structures
    SCIPcreateProbBasic(_scip1, "test");

    //sets objective sense of problem
    SCIPsetObjsense(_scip1, SCIP_OBJSENSE_MINIMIZE);

    // changes the value of an existing SCIP_Real parameter
    /*SCIPsetRealParam(_scip1, "limits/gap", 0.00001);*/
    SCIPsetIntParam(_scip1, "limits/solutions", 1);
    SCIPsetIntParam(_scip1, "presolving/maxrounds", 0);
    SCIPsetMessagehdlrQuiet(_scip1, 1);

    std::vector<SCIP_VAR*> vars_x;
    for (int i = 0; i < x_vars.size(); i++)
    {
        SCIP_VARTYPE type = SCIPvarGetType(x_vars[i]);
        double lb = SCIPvarGetLbGlobal(x_vars[i]);
        double ub = SCIPvarGetUbGlobal(x_vars[i]);
        int temp_i = SCIPvarGetIndex(x_vars[i]);
        SCIP_VAR* SET_NAME_x(temp_i) = nullptr;
        std::string name = (SCIPvarGetName(x_vars[i]));
        double coef_x_i = coefx[name];
        SCIPcreateVarBasic(_scip1, &SET_NAME_x(temp_i), name.c_str(), lb, ub, coef_x_i, type);
        /*SCIPcreateVarBasic(_scip1, &SET_NAME_x(temp_i), name.c_str(), lb, ub, coef_x_i, SCIP_VARTYPE_CONTINUOUS);*/
        SCIPaddVar(_scip1, SET_NAME_x(temp_i));
        vars_x.push_back(SET_NAME_x(temp_i));
    }
    std::vector<SCIP_VAR*> vars_z;
    for (int i = 0; i < z_vars.size(); i++)
    {
        SCIP_VARTYPE type = SCIPvarGetType(z_vars[i]);
        double lb = SCIPvarGetLbGlobal(z_vars[i]);
        double ub = SCIPvarGetUbGlobal(z_vars[i]);
        int temp_i = SCIPvarGetIndex(z_vars[i]);
        SCIP_VAR* SET_NAME_z(temp_i) = nullptr;
        std::string name = (SCIPvarGetName(z_vars[i]));
        double coef_z_i = coefz[SCIPvarGetName(z_vars[i])];
        coef_z_i += multiplier[name] + rho / 2 - rho * ori_vars[name];
        //Check whether the coefficient is negative and the upper bound of the variable is infinite
        if (coef_z_i < 0 && ub == SCIPinfinity(_scip1))
        {
            coef_z_i = 0;
        }
        SCIPcreateVarBasic(_scip1, &SET_NAME_z(temp_i), name.c_str(), lb, ub, coef_z_i, type);
        /*SCIPcreateVarBasic(_scip1, &SET_NAME_z(temp_i), name.c_str(), lb, ub, coef_z_i, SCIP_VARTYPE_CONTINUOUS);*/
        SCIPaddVar(_scip1, SET_NAME_z(temp_i));
        vars_z.push_back(SET_NAME_z(temp_i));
    }
    std::vector<std::unordered_map<std::string, double>> ori_cons;
    ori_cons = get_con_coef_vector(_scip1, block_cons);
    std::vector<std::vector<double>> ori_left_right_side;
    ori_left_right_side = get_left_right_side_vector(_scip1, block_cons);
    std::vector<SCIP_CONS*> cons;
    for (int i = 0; i < ori_cons.size(); i++)
    {
        SCIP_CONS* SET_NAME_cons(i) = nullptr;
        std::string name = "1cons";
        std::string name2 = std::to_string(i);
        name = name + name2;
        SCIPcreateConsBasicLinear(_scip1, &SET_NAME_cons(i), name.c_str(), 0, NULL, NULL, ori_left_right_side[i][0], ori_left_right_side[i][1]);
        for (int j = 0; j < vars_x.size(); j++)
        {
            std::string name_real = SCIPvarGetName(vars_x[j]);
            SCIPaddCoefLinear(_scip1, SET_NAME_cons(i), vars_x[j], ori_cons[i][name_real]);
        }
        for (int j = 0; j < vars_z.size(); j++)
        {
            std::string name_real = SCIPvarGetName(vars_z[j]);
            SCIPaddCoefLinear(_scip1, SET_NAME_cons(i), vars_z[j], ori_cons[i][name_real]);
        }
        SCIPaddCons(_scip1, SET_NAME_cons(i));
        cons.push_back(SET_NAME_cons(i));
    }

    SCIPsolve(_scip1);

    for (int j = 0; j < vars_z.size(); j++)
    {
        std::string name_real = SCIPvarGetName(vars_z[j]);
        iter_z_vars[name_real] = SCIPgetSolVal(_scip1, SCIPgetBestSol(_scip1), vars_z[j]);
    }
    for (int j = 0; j < vars_x.size(); j++)
    {
        std::string name_real = SCIPvarGetName(vars_x[j]);
        final_x_vars[name_real] = SCIPgetSolVal(_scip1, SCIPgetBestSol(_scip1), vars_x[j]);
    }

    SCIP_SOL* sol;
    sol = SCIPgetBestSol(_scip1);
    SCIP_Real objval = SCIPgetSolOrigObj(_scip1, sol);
    obj = objval;

    for (int i = 0; i < vars_x.size(); i++)
    {
        SCIPreleaseVar(_scip1, &vars_x[i]);
    }
    for (int i = 0; i < vars_z.size(); i++)
    {
        SCIPreleaseVar(_scip1, &vars_z[i]);
    }
    for (int i = 0; i < cons.size(); i++)
    {
        SCIPreleaseCons(_scip1, &cons[i]);
    }

    SCIPfree(&_scip1);
}

//Try to fix the dual coupling variable value in each block to solve the feasible solution of each block
void TryFeas_Xdot(std::vector<SCIP_VAR*> x_vars, std::vector<SCIP_VAR*> z_vars, std::vector<SCIP_CONS*> block_cons, std::map<std::string, double>  multiplier,
                  double rho, std::map<std::string, double> ori_vars, std::map<std::string, double>& iter_z_vars, std::map<std::string, double>& final_x_vars, double& obj, int& flag)
{
    std::unordered_map<std::string, double> coefx = get_obj_value_vector(x_vars);
    std::unordered_map<std::string, double> coefz = get_obj_value_vector(z_vars);

    SCIP* _scip1;
    SCIPcreate(&_scip1);

    // include default plugins
    SCIPincludeDefaultPlugins(_scip1);

    // creates empty problem and initializes all solving data structures
    SCIPcreateProbBasic(_scip1, "test");

    //sets objective sense of problem
    SCIPsetObjsense(_scip1, SCIP_OBJSENSE_MINIMIZE);

    // changes the value of an existing SCIP_Real parameter
    /*SCIPsetRealParam(_scip1, "limits/gap", 0.00001);*/
    SCIPsetRealParam(_scip1, "limits/time", 200);
    SCIPsetIntParam(_scip1, "presolving/maxrounds", 0);
    SCIPsetMessagehdlrQuiet(_scip1, 1);

    std::vector<SCIP_VAR*> vars_x;
    for (int i = 0; i < x_vars.size(); i++)
    {
        SCIP_VARTYPE type = SCIPvarGetType(x_vars[i]);
        double lb = SCIPvarGetLbGlobal(x_vars[i]);
        double ub = SCIPvarGetUbGlobal(x_vars[i]);
        int temp_i = SCIPvarGetIndex(x_vars[i]);
        SCIP_VAR* SET_NAME_x(temp_i) = nullptr;
        std::string name = (SCIPvarGetName(x_vars[i]));
        double coef_x_i = coefx[name];
        SCIPcreateVarBasic(_scip1, &SET_NAME_x(temp_i), name.c_str(), lb, ub, coef_x_i, type);
        SCIPaddVar(_scip1, SET_NAME_x(temp_i));
        vars_x.push_back(SET_NAME_x(temp_i));
    }
    std::vector<SCIP_VAR*> vars_ori;
    for (int i = 0; i < z_vars.size(); i++)
    {
        SCIP_VARTYPE type = SCIPvarGetType(z_vars[i]);
        double lb = SCIPvarGetLbGlobal(z_vars[i]);
        double ub = SCIPvarGetUbGlobal(z_vars[i]);
        int temp_i = SCIPvarGetIndex(z_vars[i]);
        SCIP_VAR* SET_NAME_ori_z(temp_i) = nullptr;
        std::string name = (SCIPvarGetName(z_vars[i]));
        std::string name_real = SCIPvarGetName(z_vars[i]);
        double ori_z_vars = ori_vars[name_real];
        if (type != SCIP_VARTYPE_CONTINUOUS)
        {
            ori_z_vars = std::round(ori_z_vars);
        }
        double coef_z_i = coefz[name_real];
        SCIPcreateVarBasic(_scip1, &SET_NAME_ori_z(temp_i), name.c_str(), ori_z_vars, ori_z_vars, coef_z_i, SCIP_VARTYPE_CONTINUOUS);
        SCIPaddVar(_scip1, SET_NAME_ori_z(temp_i));
        vars_ori.push_back(SET_NAME_ori_z(temp_i));
    }
    std::vector<std::unordered_map<std::string, double>> ori_cons;
    ori_cons = get_con_coef_vector(_scip1, block_cons);
    std::vector<std::vector<double>> ori_left_right_side;
    ori_left_right_side = get_left_right_side_vector(_scip1, block_cons);
    std::vector<SCIP_CONS*> cons;
    for (int i = 0; i < ori_cons.size(); i++)
    {
        SCIP_CONS* SET_NAME_cons(i) = nullptr;
        std::string name = "1cons";
        std::string name2 = std::to_string(i);
        name = name + name2;
        SCIPcreateConsBasicLinear(_scip1, &SET_NAME_cons(i), name.c_str(), 0, NULL, NULL, ori_left_right_side[i][0], ori_left_right_side[i][1]);
        for (int j = 0; j < vars_x.size(); j++)
        {
            std::string name_real = SCIPvarGetName(vars_x[j]);
            SCIPaddCoefLinear(_scip1, SET_NAME_cons(i), vars_x[j], ori_cons[i][name_real]);
        }
        for (int j = 0; j < vars_ori.size(); j++)
        {
            std::string name_real = SCIPvarGetName(vars_ori[j]);
            SCIPaddCoefLinear(_scip1, SET_NAME_cons(i), vars_ori[j], ori_cons[i][name_real]);
        }
        SCIPaddCons(_scip1, SET_NAME_cons(i));
        cons.push_back(SET_NAME_cons(i));
    }

    SCIPsolve(_scip1);
    if (SCIPgetStatus(_scip1) == SCIP_STATUS_OPTIMAL)
        flag = 1;

    for (int j = 0; j < vars_ori.size(); j++)
    {
        std::string name_real = SCIPvarGetName(vars_ori[j]);
        iter_z_vars[name_real] = SCIPgetSolVal(_scip1, SCIPgetBestSol(_scip1), vars_ori[j]);
    }
    for (int j = 0; j < vars_x.size(); j++)
    {
        std::string name_real = SCIPvarGetName(vars_x[j]);
        final_x_vars[name_real] = SCIPgetSolVal(_scip1, SCIPgetBestSol(_scip1), vars_x[j]);
    }

    SCIP_SOL* sol;
    sol = SCIPgetBestSol(_scip1);
    SCIP_Real objval = SCIPgetSolOrigObj(_scip1, sol);
    obj = objval;

    for (int i = 0; i < vars_x.size(); i++)
    {
        SCIPreleaseVar(_scip1, &vars_x[i]);
    }
    for (int i = 0; i < vars_ori.size(); i++)
    {
        SCIPreleaseVar(_scip1, &vars_ori[i]);
    }
    for (int i = 0; i < cons.size(); i++)
    {
        SCIPreleaseCons(_scip1, &cons[i]);
    }

    SCIPfree(&_scip1);
}

//Update the Lagrange multiplier
void UpdateMultiplier(std::vector<std::map<std::string,double>>&  multiplier, double rho, std::vector<SCIP_VAR*> z_vars, std::vector<std::map<std::string,double>> iter_z_vars, std::map<std::string, double> ori_vars)
{
    for (int i = 0; i < iter_z_vars.size(); i++)
    {
        for (int j = 0; j < z_vars.size(); j++)
        {
            std::string name_real = SCIPvarGetName(z_vars[j]);
            multiplier[i][name_real] = multiplier[i][name_real] - rho * (ori_vars[name_real] - iter_z_vars[i][name_real]);
        }
    }
}

//The solution generated by the last iteration is brought into the original problem to obtain the target value
double ValidCorrect(std::vector<DUALBlock> dual_blocks, std::map<std::string, double> ori_vars, std::vector<std::map<std::string, double>> final_x_vars,std::vector<std::map<std::string, double>> iter_z_vars, double offset, std::string tempname)
{
    double final_obj = 0;
    std::vector<double> final_z_obj;
    std::unordered_map<std::string, double> Sol;
    for (int k = 0; k < dual_blocks.size(); k++)
    {
        std::vector<SCIP_VAR*> x_vars = dual_blocks[k].Independent_vars;
        std::vector<SCIP_VAR*> z_vars = dual_blocks[k].Coupling_vars;
        std::vector<SCIP_CONS*> block_cons = dual_blocks[k].Cons;
        std::unordered_map<std::string, double> coefx = get_obj_value_vector(x_vars);
        std::unordered_map<std::string, double> coefz = get_obj_value_vector(z_vars);

        for (int i = 0; i < x_vars.size(); i++)
        {
            std::string name = (SCIPvarGetName(x_vars[i]));
            Sol.insert(std::pair<std::string, double>(name, final_x_vars[k][name]));
            double coef_x_i = coefx[name];
            final_obj += final_x_vars[k][name] * coef_x_i;
        }
        double temp_z_obj_record = 0;
        for (int i = 0; i < z_vars.size(); i++)
        {
            std::string name = (SCIPvarGetName(z_vars[i]));
            if (SCIPvarGetType(z_vars[i]) == SCIP_VARTYPE_BINARY)
            {
                ori_vars[name] = std::round(ori_vars[name]);
                if (ori_vars[name] > 1)
                    ori_vars[name] = 1;
                if (ori_vars[name] < 0)
                    ori_vars[name] = 0;
            }
            else if (SCIPvarGetType(z_vars[i]) == SCIP_VARTYPE_INTEGER)
            {
                ori_vars[name] = std::round(ori_vars[name]);
            }
            Sol.insert(std::pair<std::string, double>(name, ori_vars[name]));
            double coef_z_i = coefz[name];
            temp_z_obj_record += ori_vars[name] * coef_z_i;
        }
        final_z_obj.push_back(temp_z_obj_record);
    }
    FP fpsolver(tempname);
    fpsolver.SetIntegerValue(Sol);
    int feas_flag = fpsolver.RUNFP();
    if (feas_flag == 1)
    {
        std::cout << "ADMM and FP Find a Feasible Solution!!!" << std::endl;
        return 0;
    }

    final_obj += final_z_obj[0];
    std::cout << "Final Obj is :" << final_obj << std::endl;
    return final_obj;
}

//Determine whether the algorithm meets the convergence condition
bool Convergence(double& rho, std::vector<std::map<std::string, double>> iter_z_vars, std::map<std::string, double> last_ori_vars, std::map<std::string, double> ori_vars, std::vector<SCIP_VAR*> z_vars, int count, double feasible_obj,
                 std::vector<double> block_obj, std::map<std::string, double> init_value)
{
    int flag1 = 0;
    int flag2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum_for_step = 0;
    for (int j = 0; j < z_vars.size(); j++)
    {
        double temp_value1 = 0;
        for (int i = 0; i < iter_z_vars.size(); i++)
        {
            temp_value1 += pow(ori_vars[SCIPvarGetName(z_vars[j])] - iter_z_vars[i][SCIPvarGetName(z_vars[j])], 2);
        }
        sum1 += temp_value1;
    }
    sum_for_step = sum1;
    sum1 = sqrt(sum1);
    if (sum1 <= 0.05)
    {
        flag1 = 1;
    }
    for (int j = 0; j < z_vars.size(); j++)
    {
        double temp_value1 = abs(ori_vars[SCIPvarGetName(z_vars[j])] - last_ori_vars[SCIPvarGetName(z_vars[j])]);
        sum2 += temp_value1;
    }
    sum2 = sum2 * rho;
    if (sum2 <= 0.005)
    {
        flag2 = 1;
    }
    double miu = 20;
    if (sum1 > miu * sum2)
    {
        rho = rho * 2;
    }
    else if (sum2 > miu * sum1)
    {
        rho = rho / 1.5;
    }
    /*double lamda = 0.1;
    double iter_obj_sum = CountLRValue(block_obj, z_vars, init_value);
    rho = lamda * (feasible_obj - iter_obj_sum) / sum_for_step;*/

    sum1_record.push_back(sum1);
    sum2_record.push_back(sum2);
    rho_record.push_back(rho);
    if (count == 200)
        int tong = 0;
    return flag1 && flag2 || count > 10;
}

void TryFeasXdot(int block_num, std::vector<std::map<std::string, double>>& iter_z_vars, std::vector<std::map<std::string, double>>& final_x_vars,
                 std::vector<DUALBlock> dual_blocks, std::vector<std::map<std::string, double>> multiplier, double rho,
                 std::map<std::string, double> ori_vars, std::vector<double>& iter_obj, int& all_feas_flag, int& count_p)
{
    int i = block_num;
    std::map<std::string, double> iter_z_vars_temp = iter_z_vars[i];
    std::map<std::string, double> final_x_vars_temp = final_x_vars[i];
    std::vector<SCIP_VAR*> x_vars = dual_blocks[i].Independent_vars;
    std::vector<SCIP_VAR*> z_vars = dual_blocks[i].Coupling_vars;
    std::vector<SCIP_CONS*> block_cons = dual_blocks[i].Cons;

    double obj = 0;
    int flag = 0;
    TryFeas_Xdot(x_vars, z_vars, block_cons, multiplier[i], rho, ori_vars, iter_z_vars[i], final_x_vars[i], obj, flag);
    if (flag == 1)
    {
        iter_obj[i] = obj;
    }
    else
    {
        all_feas_flag = 0;
        iter_z_vars[i] = iter_z_vars_temp;
        final_x_vars[i] = final_x_vars_temp;
    }
    count_p++;
}

void TryFeasFP(std::unordered_map<std::string, double> Sol, std::string tempname, int& feas_flag, std::chrono::system_clock::time_point start, std::map<std::string,writeJson::Soljson>& all_sol)
{
    FP fpsolver(tempname);
    fpsolver.SetIntegerValue(Sol);
    feas_flag = fpsolver.RUNFP();
    std::map<std::string, double> fpSol = fpsolver.OutSol();

    // save FP result
    std::string part_name = tempname.substr(Folder.size() +1,tempname.size()-Folder.size()-5);
    std::ofstream o((part_name + "ADMMFP.json").c_str(), std::ofstream::trunc);
    writeJson::Soljson tmpSol;
    tmpSol.source = "FP";
    tmpSol.obj = fpsolver.GetTempObjVal();
    tmpSol.duration = double (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count()) / 1000;
    nlohmann::json j;
    tmpSol.to_json(j);
    o << j << std::endl;

    std::unordered_map<std::string, double> fpSol2(fpSol.begin(), fpSol.end());
    double tmp_obj = fpsolver.GetTempObjVal();
    auto Rins_sol = RINS::RINSrun(tempname, fpSol2, tmp_obj, 10, 1000, 0.001);
    if (!Rins_sol.sol_.empty())
    {
        double rins_sol_obj = Rins_sol.obj_;
        std::cout << "rins sol obj is:" << rins_sol_obj << std::endl;

        // save RINS result
        std::ofstream o_r((part_name + "ADMMRINS.json").c_str(), std::ofstream::trunc);
        writeJson::Soljson tmpSol2;
        tmpSol2.source = "RINS";
        tmpSol2.obj = rins_sol_obj;
        tmpSol2.duration = double (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count()) / 1000;
        nlohmann::json j2;
        tmpSol2.to_json(j2);
        o_r << j2 << std::endl;

        all_sol[part_name + "_ADMM_FP"] = tmpSol;
        all_sol[part_name + "_ADMM_RINS"] = tmpSol2;
    }

}


//ADMM algorithm
void NewIterUpdate(std::vector<DUALBlock> dual_blocks, std::map<std::string, double> init_value, double feasible_obj, std::vector<SCIP_VAR*> scipvars, double offset, std::string tempname, std::shared_ptr<spdlog::logger> spd,std::chrono::system_clock::time_point start, std::unordered_map<std::string, int> In_Var_Ni, std::map<std::string,writeJson::Soljson>& all_sol)
{
    bool usePresolved = 0;
    size_t maxTotalSolutions = 1;

    //Initializes the Lagrange multiplier and step size
    std::vector<std::map<std::string, double>> multiplier(dual_blocks.size());
    for (int j = 0; j < multiplier.size(); j++)
    {
        std::map<std::string, double> temp_multiplier;
        for (int i = 0; i < dual_blocks[0].Coupling_vars.size(); i++)
        {
            std::string name_real = SCIPvarGetName(dual_blocks[0].Coupling_vars[i]);
            temp_multiplier.insert(std::pair<std::string, double>(name_real, 0.0));
        }
        multiplier[j] = temp_multiplier;
    }
    double rho = 50;
    //Record the solution for each block coupled variable
    std::vector<std::map<std::string, double>> iter_z_vars(dual_blocks.size());
    for (int i = 0; i < iter_z_vars.size(); i++)
    {
        std::map<std::string, double> temp_z;
        for (int j = 0; j < dual_blocks[0].Coupling_vars.size(); j++)
        {
            std::string name_real = SCIPvarGetName(dual_blocks[0].Coupling_vars[j]);
            temp_z.insert(std::pair<std::string, double>(name_real, 0.0));
        }
        iter_z_vars[i] = temp_z;
    }
    //Record the solution for each block corner block variable (last iteration verification solution used)
    std::vector<std::map<std::string, double>> final_x_vars(dual_blocks.size());
    for (int i = 0; i < final_x_vars.size(); i++)
    {
        std::map<std::string, double> temp_final_x_vars;
        for (int j = 0; j < dual_blocks[i].Independent_vars.size(); j++)
        {
            std::string name_real = SCIPvarGetName(dual_blocks[i].Independent_vars[j]);
            temp_final_x_vars.insert(std::pair<std::string, double>(name_real, 0.0));
        }
        final_x_vars[i] = temp_final_x_vars;
    }
    //Record the solution for the original coupling variable
    std::vector<std::map<std::string, double>> every_ori_vars;
    std::vector<std::vector<std::map<std::string, double>>> every_iter_multiplier;
    every_iter_multiplier.push_back(multiplier);
    std::vector<std::vector<double>> every_iter_obj;
    std::vector<double> every_rho;
    every_rho.push_back(rho);
    std::vector<double> iter_obj(dual_blocks.size());
    std::map<std::string, double> ori_vars;
    for (int i = 0; i < dual_blocks[0].Coupling_vars.size(); i++)
    {
        std::string name_real = SCIPvarGetName(dual_blocks[0].Coupling_vars[i]);
        ori_vars.insert(std::pair<std::string, double>(name_real, init_value[name_real]));
    }
    every_ori_vars.push_back(ori_vars);

    int count = 0;
    //Collect the target value of the feasible solution
    std::vector<double> feas_objs;
    //Collect the number of constraint violations
    std::vector<int> vio_cons_num_now;
    //Collect historical iteration information
    std::vector<std::map<std::string, double>> ori_vars_iter;
    std::vector<std::vector<std::map<std::string, double>>>  final_x_vars_iter;
    //Judge the generation that converges best in history
    double best_converg_iter = 1e10;
    int best_converg_iter_num = 0;
    while (1)
    {
        std::map<std::string, double> last_ori_vars = ori_vars;
        if (count != 0)
        {
            std::unordered_map<std::string, double> Sol;
            int all_feas_flag = 1;
            int count_p = 0;

            for (int i = 0; i < dual_blocks.size(); i++)
            {
                std::map<std::string, double> iter_z_vars_temp = iter_z_vars[i];
                std::map<std::string, double> final_x_vars_temp = final_x_vars[i];
                std::vector<SCIP_VAR*> x_vars = dual_blocks[i].Independent_vars;
                std::vector<SCIP_VAR*> z_vars = dual_blocks[i].Coupling_vars;
                std::vector<SCIP_CONS*> block_cons = dual_blocks[i].Cons;
                for (auto x : x_vars)
                {
                    std::string name = SCIPvarGetName(x);
                    Sol.insert(std::pair<std::string, double>(name, final_x_vars_temp[name]));
                }
            }
            for (auto z : dual_blocks[0].Coupling_vars)
            {
                std::string name = SCIPvarGetName(z);

                if (SCIPvarGetType(z) == SCIP_VARTYPE_BINARY)
                {
                    ori_vars[name] = std::round(ori_vars[name]);
                    if (ori_vars[name] > 1)
                        ori_vars[name] = 1;
                    if (ori_vars[name] < 0)
                        ori_vars[name] = 0;
                }
                else if (SCIPvarGetType(z) == SCIP_VARTYPE_INTEGER)
                {
                    ori_vars[name] = std::round(ori_vars[name]);
                }
                Sol.insert(std::pair<std::string, double>(name, ori_vars[name]));
            }
            std::unordered_map<int, double> in_values;
            for (int zz = 0; zz < scipvars.size(); zz++)
            {
                string name = SCIPvarGetName(scipvars[zz]);
                in_values.insert(std::pair<int, double>(In_Var_Ni[name], Sol[name]));
            }
            const char* filename = tempname.c_str();
            int feas_flag = 0;
            int FJ_flag = 0;

            /*for (int i = 0; i < dual_blocks.size(); i++)
            {
                TryFeasXdot(i, iter_z_vars, final_x_vars, dual_blocks, multiplier, rho, ori_vars, iter_obj, all_feas_flag, count_p);
            }*/
            TryFeasFP(Sol, tempname, feas_flag, start, all_sol);

            FJ_flag = Use_FJ(In_Var_Ni, in_values, filename);

            if (FJ_flag == 1)
            {
                std::cout << "ADMM and FJ find a Feasible Solution!!!" << std::endl;
                return;
            }
            if (feas_flag == 1)
            {
                std::cout << "ADMM and FP find a Feasible Solution!!!" << std::endl;
                return;
            }
            //        if (all_feas_flag == 1)
            //        {
            //            //At present, the fixed value of the coupling variable is feasible in all blocks, and the global feasible solution is obtained
            //            ori_vars_iter.push_back(ori_vars);
            //            final_x_vars_iter.push_back(final_x_vars);
            //            Convergence(rho, iter_z_vars, last_ori_vars, ori_vars, dual_blocks[0].Coupling_vars, count, feasible_obj, iter_obj, init_value);
            //std::cout << "ADMM Find a Feasible Solution!!!" << std::endl;
            //            break;
            //        }
        }
        for (int i = 0; i < dual_blocks.size(); i++)
        {
            std::vector<SCIP_VAR*> x_vars = dual_blocks[i].Independent_vars;
            std::vector<SCIP_VAR*> z_vars = dual_blocks[i].Coupling_vars;
            std::vector<SCIP_CONS*> block_cons = dual_blocks[i].Cons;
            double obj = 0;
            //Check whether the coupled variable has an outlier, and take its lower bound if it does
            for (int zz = 0; zz < dual_blocks[0].Coupling_vars.size(); zz++)
            {
                string name = SCIPvarGetName(dual_blocks[0].Coupling_vars[zz]);
                double temp_lb = SCIPvarGetLbGlobal(dual_blocks[0].Coupling_vars[zz]);
                if (ori_vars[name] > 1e10)
                    ori_vars[name] = temp_lb;
            }
            BuildProb_Xdot(x_vars, z_vars, block_cons, multiplier[i], rho, ori_vars, iter_z_vars[i], final_x_vars[i], obj);
            iter_obj[i] = obj;
        }
        every_iter_obj.push_back(iter_obj);
        {
            //Fix a solution to the ADMM's z variable with any heuristic
            if (count % 150 == 0 && count != 0)
            {

            }
            else
            {
                //Instead of trying to fix z with a heuristic, update z directly with the average of the three blocks
                for (int p = 0; p < dual_blocks[0].Coupling_vars.size(); p++)
                {
                    double temp_value = 0;
                    double temp_value2 = 0;
                    std::string name_real = SCIPvarGetName(dual_blocks[0].Coupling_vars[p]);
                    for (int i = 0; i < dual_blocks.size(); i++)
                    {
                        temp_value += iter_z_vars[i][name_real];
                        temp_value2 += multiplier[i][name_real];
                    }
                    temp_value2 = temp_value2 / (dual_blocks.size() * rho);
                    ori_vars[name_real] = temp_value / dual_blocks.size();
                    if (SCIPvarGetType(dual_blocks[0].Coupling_vars[p]) != SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(dual_blocks[0].Coupling_vars[p]) != SCIP_VARTYPE_INTEGER)
                    {
                        ori_vars[name_real] = std::round(ori_vars[name_real]);
                        if (ori_vars[name_real] > 1)
                            ori_vars[name_real] = 1;
                        if (ori_vars[name_real] < 0)
                            ori_vars[name_real] = 0;
                    }
                }
            }
            ori_vars_iter.push_back(ori_vars);
            final_x_vars_iter.push_back(final_x_vars);
        }
        if (Convergence(rho, iter_z_vars, last_ori_vars, ori_vars, dual_blocks[0].Coupling_vars, count, feasible_obj, iter_obj, init_value))
            break;
        UpdateMultiplier(multiplier, rho, dual_blocks[0].Coupling_vars, iter_z_vars, ori_vars);
        every_iter_multiplier.push_back(multiplier);
        count++;
        every_rho.push_back(rho);
    }

    /*cout << "Convergence changes under the primal criteria:" << endl;
    for (int p = 0; p < sum1_record.size(); p++)
    {
        cout << sum1_record[p] << endl;
    }
    cout << "Convergence changes under dual criteria:" << endl;
    for (int p = 0; p < sum2_record.size(); p++)
    {
        cout << sum2_record[p] << endl;
    }
    cout << "Step values collected during iteration:" << endl;
    for (int p = 0; p < rho_record.size(); p++)
    {
        cout << rho_record[p] << endl;
    }
    cout << "Feasible solution target values collected during iteration:" << endl;
    for (int p = 0; p < feas_objs.size(); p++)
    {
        cout << feas_objs[p] << endl;
    }
    cout << "The number of constraint violations collected during iteration:" << endl;
    for (int p = 0; p < vio_cons_num_now.size(); p++)
    {
        cout << vio_cons_num_now[p] << endl;
    }*/
    /*for (int i = 0; i < sum1_record.size(); i++)
    {
        if (sum1_record[i] < best_converg_iter)
        {
            best_converg_iter = sum1_record[i];
            best_converg_iter_num = i;
        }
    }
    ValidCorrect(dual_blocks, ori_vars_iter[best_converg_iter_num], final_x_vars_iter[best_converg_iter_num], iter_z_vars, offset, tempname);*/
    std::cout << "ADMM algorithm ends!" << std::endl;
}



