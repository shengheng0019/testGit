#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <map>
#include <set>
#include <iostream>
#include <set>
#include <vector>
#include <cstdlib>
#include <queue>
#include <unordered_map>
#include <random>
#include <unordered_set>




namespace RINS
{

    inline std::unordered_map<std::string,SCIP_VARTYPE> getallVarType(SCIP* scip)
    {
        std::unordered_map<std::string,SCIP_VARTYPE> Rvar_type;
        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        int nvars = SCIPgetNVars(scip);
        std::vector<SCIP_VAR*> fj_vars(vars_ori, vars_ori + nvars);
        for (int i = 0; i < nvars; ++i)
        {
            auto var = fj_vars[i];
            SCIP_VARTYPE const varType_ = SCIPvarGetType(var);
            std::string name = SCIPvarGetName(var);
            Rvar_type[name] = varType_;
        }
        return Rvar_type;
    }

    inline std::unordered_set<std::string> getIntVars(SCIP* scip)
    {
        std::unordered_set<std::string> int_vars;
        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        int nvars = SCIPgetNVars(scip);
        std::vector<SCIP_VAR*> fj_vars(vars_ori, vars_ori + nvars);
        for (int i = 0; i < nvars; ++i)
        {
            auto var = fj_vars[i];
            SCIP_VARTYPE const varType_ = SCIPvarGetType(var);
            std::string name = SCIPvarGetName(var);
            if (varType_ != SCIP_VARTYPE_CONTINUOUS)
                int_vars.insert(name);
        }
        return int_vars;
    }


    inline void relaxModel(SCIP* scip)
    {
//        SCIPfreeTransform(scip);
        std::cout<<"relax stage : "<<SCIPgetStage(scip)<<std::endl;
        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        int nvars = SCIPgetNVars(scip);
        std::vector<SCIP_VAR*> fj_vars(vars_ori, vars_ori + nvars);
        for (int i = 0; i < nvars; ++i)
        {
            auto var = fj_vars[i];
            SCIP_VARTYPE const varType_ = SCIPvarGetType(var);
            if (varType_ == SCIP_VARTYPE_CONTINUOUS)
                continue;

            SCIP_Bool infeasible;
            SCIPchgVarType(scip, var, SCIP_VARTYPE_CONTINUOUS, &infeasible);
        }

    }

    inline void convertModel(SCIP* scip,std::unordered_map<std::string,SCIP_VARTYPE> var_type)
    {
        std::cout<<"convert stage : "<<SCIPgetStage(scip)<<std::endl;
        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        int nvars = SCIPgetNVars(scip);
        std::vector<SCIP_VAR*> fj_vars(vars_ori, vars_ori + nvars);

        for (int i = 0; i < nvars; ++i)
        {
            auto var = fj_vars[i];
            auto type = var_type[SCIPvarGetName(var)];
            if (type == SCIP_VARTYPE_CONTINUOUS)
            {
                continue;
            }
            else
            {
                SCIP_Bool infeasible;
                SCIPchgVarType(scip, var, type, &infeasible);
            }
        }
    }

    inline std::unordered_map<std::string,int> getfixedVars(std::unordered_map<std::string,double> incumbent_sol,std::unordered_map<std::string,double> relax_sol,const std::unordered_set<std::string>& int_vars)
    {
        std::cout<<"getfixedVars stage : "<<incumbent_sol.size()<<"---"<<relax_sol.size()<<"---"<<int_vars.size()<<std::endl;
        std::unordered_map<std::string,int> fixed_vars;
        for (auto const& var : int_vars)
        {
            if (std::abs(incumbent_sol[var] - relax_sol[var]) < 1e-5)
            {
                fixed_vars[var] = std::round(incumbent_sol[var]);
            }
        }
        return fixed_vars;
    }

    //Fix the variables that have the same values in the incumbent and in the current continuous relaxation
    inline void fixedModel(SCIP* scip,const std::unordered_map<std::string,int>& fixed_vars)
    {
//        std::cout<<"fixed begin : "<< SCIPgetStage(scip)<<std::endl;
        for (auto const& var : fixed_vars)
        {
            SCIP_VAR* var_ = SCIPfindVar(scip, var.first.c_str());
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIPfixVar(scip, var_,var.second,&infeasible,&fixed);
        }
    }

    //Set an objective cutoff based on the objective value of the current incumbent
    //Solve a sub-MIP on the remaining variables
    inline void setObjCutoff(SCIP* scip,double obj)
    {
        SCIP_Real cutoff = obj - 1;
        std::cout<<"setObjCutoff value  : "<<cutoff<<std::endl;
        SCIPsetObjlimit(scip, cutoff);
    }

    inline std::unordered_map<std::string,double> solverelaxmodel(SCIP* scip,std::unordered_map<std::string,SCIP_VARTYPE> var_type)
    {
        std::cout<<"solverelaxmodel stage : "<<SCIPgetStage(scip)<<std::endl;
        SCIPsolve(scip);
        SCIP_SOL* sol = SCIPgetBestSol(scip);
        std::unordered_map<std::string,double> sol_;

        for (auto tmpvar : var_type)
        {
            auto var = SCIPfindVar(scip,tmpvar.first.c_str());
            std::string name = SCIPvarGetName(var);
            sol_[name] = SCIPgetSolVal(scip, sol, var);
        }
        SCIPfreeTransform(scip);
        return sol_;
    }

    class RINSResult
    {
    public:
        double obj_;
        std::unordered_map<std::string,double> sol_;
    };

    inline RINSResult RINSrun(std::string filename,std::unordered_map<std::string,double>& incumbent_sol,double& incumbent_obj,double meanTime,
                              long int findnodes,const double limit_gap)
    {
        double begin_obj = incumbent_obj;
//        SCIP* scip = NULL;
//        SCIPcreate(&scip);
//        // 使用 SCIPcopy() 来复制原始问题
//        SCIP_Bool valid;
//        SCIPcopyOrig(_scip, scip, NULL, NULL, "_copy", TRUE, FALSE, TRUE, &valid);
        SCIP* scip = NULL;
        SCIPcreate(&scip);
        SCIPincludeDefaultPlugins(scip);
        //将string类型的filename转换成const char*类型
        // 读入模型的mps文件
        SCIPreadProb(scip, filename.c_str(), NULL);
        SCIPsetMessagehdlrQuiet(scip, true);
//        SCIPsetIntParam(scip, "presolving/maxrounds", 0);
        int nrows = SCIPgetNConss(scip);
        int ncols = SCIPgetNVars(scip);
        std::cout<<"RINS -> load model success "<<"the num of cons is "<<nrows<<"  the num of vars is "<<ncols<<std::endl;

        //getall vars type
        auto vars_type = getallVarType(scip);
        //get all int vars
        auto int_vars = getIntVars(scip);
        std::cout<<"RINS -> get all vars type success,the num of var is : "<<vars_type.size()<<"  the num of int var is : "<<int_vars.size()<<std::endl;
        int flag = 0;
        while (flag<100)
        {
            std::cout<<"==============iter : "<<flag<<"=================="<<std::endl;
            flag++;
            //Relax the model
            relaxModel(scip);
            auto relax_sol = solverelaxmodel(scip,vars_type);
            std::cout<<"RINS -> relax model success "<<std::endl;
            //convert the model
            convertModel(scip,vars_type);
            //Fix the variables that have the same values in the incumbent and in the current continuous relaxation
            auto fixed_vars = getfixedVars(incumbent_sol,relax_sol,int_vars);
            if (fixed_vars.empty() || fixed_vars.size() == int_vars.size())
            {
                std::cout<<"RINS -> fixed model failed "<<std::endl;
                break;
            }
            fixedModel(scip,fixed_vars);
            std::cout<<"RINS -> fixed model success  "<<"total fixed num is : "<<fixed_vars.size()<<std::endl;
            //Set an objective cutoff based on the objective value of the current incumbent
            setObjCutoff(scip,incumbent_obj);

            //Solve a sub-MIP on the remaining variables
            //set the parameters
            //SCIPsetIntParam(scip, "limits/solutions", MIN_VALID_NUM);
//            SCIPsetIntParam(scip, "limits/nodes", findnodes);
            SCIPsetRealParam(scip, "limits/time", meanTime);
            //set gaplimit
//            SCIPsetRealParam(scip, "limits/gap", limit_gap);
            SCIPsolve(scip);
            //get status
            SCIP_STATUS status = SCIPgetStatus(scip);
            int nSols = SCIPgetNSols(scip);
            if (nSols > 0)
            {
                SCIP_SOL* bestSol = SCIPgetBestSol(scip);
                auto tmp_cumobj = SCIPgetSolOrigObj(scip, bestSol);
                if (tmp_cumobj < incumbent_obj)
                {
                    incumbent_obj = tmp_cumobj;
                    std::cout<<"find new solution value is : "<<incumbent_obj<<std::endl;
                    std::unordered_map<std::string,double> sol_;
                    for (auto const& var : vars_type)
                    {
//                    std::cout<<"mip sol : "<<var.first<<std::endl;
                        SCIP_VAR* var_ = SCIPfindVar(scip, var.first.c_str());
                        sol_[var.first] = SCIPgetSolVal(scip, bestSol, var_);
                    }
                    incumbent_sol = sol_;
                }
                else
                {
                    std::cout<<"RINS -> find new solution but not better "<<std::endl;
                    break;
                }
            }
            else if (status == SCIP_STATUS_INFEASIBLE)
            {
                std::cout<<"RINS -> infeasible "<<std::endl;
                break;
            }
            else
            {
                findnodes *= 2;
                meanTime *= 2;
            }
            SCIPfreeTransform(scip);
        }
        RINSResult result;
        if (incumbent_obj < begin_obj)
        {
            result.obj_ = incumbent_obj;
            result.sol_ = incumbent_sol;
        }
        std::cout<<"================RINS -> end ================="<<std::endl;
        return result;
    }

}