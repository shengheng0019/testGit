#include "Dection.hpp"
#include"Lagrangean.h"
#include"SetVarsBounds.h"
#include "FP.h"
#include "ADMM.h"
#include "DW.h"
#include "utility.h"
#include "spdlog/spdlog.h"
#include <spdlog/sinks/basic_file_sink.h>
#include "RINS.hpp"

extern std::map<std::string,writeJson::Soljson> all_solsmap;;

void LagrangeanFunc(std::map<int, General_Block> block_, SCIP* scipmodel, vector<SCIP_VAR*> scipvars, vector<SCIP_CONS*> scipcons, std::unordered_map<std::string, int> Var_Ni,std::string tempname,std::string part_name ,std::chrono::system_clock::time_point start, std::shared_ptr<spdlog::logger> spd)
{
    //        double UB = getUB(tempname);
    double UB = 1024;
    //auto [vars, conss] = getAllVarsandCons(block_);
    //SetVarsBounds(vars, conss);
    ////update blocks
    //for (int blkIdx = 0; blkIdx < block_.size(); blkIdx++)
    //{
    //    for (int varIdx = 0; varIdx < block_[blkIdx].bVars.size(); varIdx++)
    //    {
    //        block_[blkIdx].bVars[varIdx].Lb = vars[block_[blkIdx].bVars[varIdx].Id].Lb;
    //        block_[blkIdx].bVars[varIdx].Ub = vars[block_[blkIdx].bVars[varIdx].Id].Ub;
    //    }
    //}
    Lagrangean lag{ block_, scipmodel, scipcons, scipvars, Var_Ni, UB, part_name, start };
    lag.RUN(1);
    std::unordered_map<string, double> lagSol;
    for (auto iter = Var_Ni.begin(); iter != Var_Ni.end(); iter++)
    {
        lagSol.insert(std::pair<std::string, double>(iter->first, lag.masterSol.varValue[iter->second]));
    }

    Use_FJ(Var_Ni, lag.masterSol.varValue, tempname.c_str());

    FP fpSolver(tempname);
    fpSolver.SetIntegerValue(lagSol);
    fpSolver.RUNFP();
    std::map<std::string, double> fpSol = fpSolver.OutSol();
    std::ofstream o((part_name + "FP.json").c_str(), std::ofstream::trunc);

    nlohmann::json j_array = nlohmann::json::array();
    writeJson::Soljson tmpSol;
    tmpSol.source = "FP";
    tmpSol.obj = fpSolver.GetTempObjVal();
    tmpSol.solution = fpSol;
    tmpSol.duration = double (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count()) / 1000;
    all_solsmap["Lag_1"] = tmpSol;
    nlohmann::json j;
    tmpSol.to_json(j);

    //j_array.emplace_back(j);


    std::unordered_map<std::string,double> fpSol2(fpSol.begin(),fpSol.end());
    double tmp_obj = tmpSol.obj;
    auto Rins_sol = RINS::RINSrun(tempname,fpSol2,tmp_obj,10,1000,0.001);
    std::cout<<"RINS -> start  the obj before rins is "<<tmpSol.obj<<std::endl;
    if (Rins_sol.sol_.empty())
    {
        std::cout<<"RINS -> failed"<<std::endl;
    }
    else
    {
        std::ofstream o_rins((part_name + "FPRINS.json").c_str(), std::ofstream::trunc);
        std::cout<<"RINS -> success  the obj after rins is "<<Rins_sol.obj_<<std::endl;
        writeJson::Soljson tmpSol2;
        tmpSol2.source = "RINS";
        tmpSol2.obj = Rins_sol.obj_;
        std::map<std::string,double> tmpSol2Sol(Rins_sol.sol_.begin(),Rins_sol.sol_.end());
        tmpSol2.solution = tmpSol2Sol;
        tmpSol2.duration = double (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count()) / 1000;
        all_solsmap["Lag_2"] = tmpSol2;
        nlohmann::json j2;
        tmpSol2.to_json(j2);
        o_rins << j2 << std::endl;
        //j_array.emplace_back(j2);
    }
    o << j << std::endl;
}

//void DWFunc(std::map<int, General_Block>block_, SCIP*scipmodel, Objective master_obj)
//{
//    std::cout << "DWFunc1" << std::endl;
//    vector<ColumnPool> columnpoolInput;
//    int genColNum = 20;
//    double initScipFeasObj = 0;
//
//    if (GenerateInitColumns3(columnpoolInput, block_, genColNum, scipmodel, initScipFeasObj))
//    {
//        clock_t time_9 = clock();
//        DW::Diving_new(columnpoolInput, block_, master_obj, scipmodel, all_solsmap);
//        clock_t time_10 = clock();
//        cout << "SCIP init feasible obj: " << initScipFeasObj << endl;
//        cout << "Diving using: " << endl;
//        cout << (time_10 - time_9) / CLOCKS_PER_SEC << "s" << endl;
//    }
//    else
//    {
//        cout << "Diving do not start" << endl;
//    }
//
//}

void DWFunc(std::map<int, General_Block>block_, SCIP*scipmodel, Objective master_obj,std::map<std::string,writeJson::Soljson>& all_solsmap, std::string part_name, std::chrono::system_clock::time_point start)
{
    std::cout << "DWFunc1" << std::endl;
    vector<ColumnPool> columnpoolInput;
    int genColNum = 20;
    double initScipFeasObj = 0;

    if (GenerateInitColumns3(columnpoolInput, block_, genColNum, scipmodel, initScipFeasObj))
    {
        clock_t time_9 = clock();
        DW::Diving_new(columnpoolInput, block_, master_obj, scipmodel, all_solsmap, part_name, start);
        clock_t time_10 = clock();
        cout << "SCIP init feasible obj: " << initScipFeasObj << endl;
        cout << "Diving using: " << endl;
        cout << (time_10 - time_9) / CLOCKS_PER_SEC << "s" << endl;
    }
    else
    {
        cout << "Diving do not start" << endl;
    }

}

void ADMMFunc(SCIP* scipmodel, std::vector<SCIP_VAR*> scipvars, std::vector<DUALBlock> dualblock, std::string tempname, std::shared_ptr<spdlog::logger> spd,std::chrono::system_clock::time_point start, std::unordered_map<std::string, int> In_Var_Ni, std::map<std::string,writeJson::Soljson>& all_sol)
{
    std::cout << "ADMMFunc" << std::endl;
	std::cout << "ADMM start!!!" << std::endl;
	UseADMM(scipmodel, scipvars, dualblock, tempname, spd,start, In_Var_Ni, all_sol);
	std::cout << "ADMM End!!!" << std::endl;
}