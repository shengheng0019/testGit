#pragma once
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

namespace ReadFileLocalBranch
{
	SCIP* readmodel(std::string filename)
	{
        SCIP* scip = NULL;
        SCIPcreate(&scip);
        SCIPincludeDefaultPlugins(scip);
        // ����һЩ������
        SCIPsetIntParam(scip, "presolving/maxrounds", 0);
        SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", 0);
        SCIPsetBoolParam(scip, "constraints/linear/upgrade/indicator", false);
        //��string���͵�filenameת����const char*����
        auto filename1 = filename.c_str();
        // ����ģ�͵�mps�ļ�
        SCIPreadProb(scip, filename1, NULL);
        int nrows = SCIPgetNConss(scip);
        int ncols = SCIPgetNVars(scip);
        std::cout<<"load model success !"<<"the num of cons is "<<nrows<<"  the num of vars is "<<ncols<<std::endl;
        //SCIPsolve(scip);
        return scip;
	}

    double random()
    {
        // ���������
        double r = (double)rand() / RAND_MAX;
        return r;
    }

    int randomInt(int a,int b)
	{
        std::mt19937 gen(1024); // ʹ��÷ɭ��ת�㷨��ʼ��һ�������������
        std::uniform_int_distribution<> dis(a, b); //[a,b]

        // ���ɲ���������
        return dis(gen);
	}

    int randomDouble(int a, int b)
    {
        std::mt19937 gen(1024); // ʹ��÷ɭ��ת�㷨��ʼ��һ�������������
        std::uniform_real_distribution<> dis(a, b);
        // ���ɲ���������
        return dis(gen);
    }

}

namespace LocalBranching
{
	std::unordered_map<std::string, double> getSolution(SCIP* scip)
	{
        std::unordered_map<std::string, double> tmpsolution;
        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        int ncols = SCIPgetNVars(scip);
        for (int i = 0; i < ncols; i++)
		{
			SCIP_VAR* var = vars_ori[i];
			std::string varname = SCIPvarGetName(var);
            auto var_type = SCIPvarGetType(var);
			if (var_type == SCIP_VARTYPE_BINARY)
			{
                auto r = ReadFileLocalBranch::random();
                if (r > 0.5)
                {
	                tmpsolution[varname] = 1;
				}
				else
				{
					tmpsolution[varname] = 0;
                }
			}
			else if(var_type == SCIP_VARTYPE_INTEGER)
			{
				auto lb = SCIPvarGetLbLocal(var);
                auto ub = SCIPvarGetUbLocal(var);
				if (ub - lb > 1000)
				{
                    tmpsolution[varname] = lb;
				}
				else
				{
                    auto value = ReadFileLocalBranch::randomInt(lb, ub);
                    tmpsolution[varname] = value;
				}
			}
			else
			{
				auto lb = SCIPvarGetLbLocal(var);
				auto ub = SCIPvarGetUbLocal(var);
				if (ub - lb > 1000)
				{
                    tmpsolution[varname] = lb;
				}
				else
				{
                    auto value = ReadFileLocalBranch::randomDouble(lb, ub);
                    tmpsolution[varname] = value;
				}
			}
		}
        
        return tmpsolution;
	}

    double computeObj(SCIP* scip,const std::unordered_map<std::string,double>& sol)
	{
        std::cout<<"compute obj begin"<<std::endl;
        double obj = 0.0;
        int nvars = SCIPgetNVars(scip);
        std::cout<<"the num of vars is "<<nvars<<std::endl;
        SCIP_VAR** vars = SCIPgetVars(scip);

        for (int i = 0; i < nvars; i++)
        {
            std::string var_name = SCIPvarGetName(vars[i]);
            auto it = sol.find(var_name.substr(2));

            if (it != sol.end())
            {
                double var_value = it->second;
                double obj_coef = SCIPvarGetObj(vars[i]);
//                std::cout<<"  var_value : "<<var_value<<"  obj_coef : "<<obj_coef<<std::endl;
                obj += obj_coef * var_value;
            }
        }
        std::cout<<"compute obj end : "<<obj<<std::endl;
        return obj;
	}

    
    std::tuple<std::unordered_map<std::string, double>, std::unordered_map<std::string, double>, std::unordered_map<std::string, double>> LocalBranchingFun(SCIP* _scip,const std::unordered_map<std::string, double>& tmpsolution,double neighbor)
	{
        int _ncols = SCIPgetNVars(_scip);
        int _nrows = SCIPgetNConss(_scip);
        std::cout << "the num of vars of model is " << _ncols << "  the num of cons of model is " << _nrows << std::endl;

        std::unordered_set<std::string> binary_fixed;
        std::unordered_set<std::string> binary_free;
        std::unordered_set<std::string> redundant_var;
        std::unordered_map<std::string, SCIP_VAR*> varmap;
        SCIP* scip = NULL;
        SCIPcreate(&scip);
        // ʹ�� SCIPcopy() ������ԭʼ����
        SCIP_Bool valid;
        SCIPcopyOrig(_scip, scip, NULL, NULL, "_copy", TRUE, FALSE, TRUE, &valid);
        // ����������ʱ��Ϊ3��
        SCIP_RETCODE retcode = SCIPsetRealParam(scip, "limits/time", 3.0);
        if (retcode != SCIP_OKAY)
        {
            SCIPerrorMessage("Error setting parameter limits/time\n");
            SCIPprintError(retcode);
        }
        SCIPsetMessagehdlrQuiet(scip, true);
        int ncols = SCIPgetNVars(scip);
        int nrows = SCIPgetNConss(scip);
        std::cout<<"the num of vars of copymodel is "<<ncols<<"  the num of cons of copymodel is "<<nrows<<std::endl;

        // ����һ���µĽ�
        SCIP_SOL* sol;
        SCIPcreateSol(scip, &sol, NULL);

        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        for (int i = 0; i < ncols; ++i)
        {
			SCIP_VAR* var = vars_ori[i];
			std::string varname = SCIPvarGetName(var);
            //std::cout << varname << std::endl;
			varmap[varname] = var;

            //change obj
            SCIP_Real newcoef = 0;
            SCIPchgVarObj(scip,var,newcoef);
        }
        // Ϊÿ����������ֵ
        for (const auto& item : tmpsolution)
        {
            SCIPsetSolVal(scip, sol, varmap.at(item.first), item.second);
            auto var_type = SCIPvarGetType(varmap.at(item.first));
            if (var_type == SCIP_VARTYPE_BINARY)
            {
	            if (item.second == 1)
	            {
		            binary_fixed.insert(item.first);
	            }
	            else
	            {
		            binary_free.insert(item.first);
	            }
            }
        }
        
        // ���ÿ��Լ��

        std::unordered_map<std::string, double> violate_constraints;
        SCIP_CONS** conss = SCIPgetConss(scip);
        for (int i = 0; i < nrows; ++i)
        {
            SCIP_CONS* consdata = conss[i];
            SCIP_RESULT result;
        	SCIPcheckCons(scip,consdata,sol,TRUE,TRUE,TRUE,&result);
			if (result == SCIP_INFEASIBLE)
			{
                double tmpvalue = 0;
                double delta = 0;
                // ���Լ����ϵ������Ӧ�ı���
                int nnonz = SCIPgetNVarsLinear(scip, consdata);
                SCIP_VAR** vars = SCIPgetVarsLinear(scip, consdata);
                SCIP_Real* coefs = SCIPgetValsLinear(scip, consdata);
                // ����һ��Լ���еķ���Ԫ��
                for (int j = 0; j < nnonz; ++j)
                {
                    auto temp_name = SCIPvarGetName(vars[j]);
                    //std::cout<<"coef is "<<coefs[j]<<"  x is  "<< tmpsolution.at(temp_name) <<std::endl;
                    tmpvalue += coefs[j] * tmpsolution.at(temp_name);
                }
                unsigned success;
                auto Lb = SCIPconsGetLhs(scip, consdata, &success);
                auto Ub = SCIPconsGetRhs(scip, consdata, &success);
                std::string consname = SCIPconsGetName(consdata);
                if (Lb==Ub)
                {
                    delta = Ub - tmpvalue;
                    //std::cout << Lb << "  " << tmpvalue << "  " << delta << std::endl;
                }
                else if (Lb > -SCIPinfinity(scip))
                {
	                delta = Lb - tmpvalue;
                    //std::cout << Lb << "  " << tmpvalue << "  " <<delta<< std::endl;
                }
				else
				{
	                delta = Ub - tmpvalue;
                    //std::cout << Ub << "  " << tmpvalue << "  " << delta << std::endl;
				}

                //add variable
                SCIP_VAR* y;
                std::string yname = "y_" + consname;
                SCIPcreateVarBasic(scip, &y,yname.c_str() , 0, 1, 1.0, SCIP_VARTYPE_BINARY);
                SCIPaddVar(scip,y);
                SCIPaddCoefLinear(scip, consdata, y, delta);
                std::cout<<"violate constraints name is  "<< consname<< "delta is " << delta << std::endl;
                violate_constraints[consname] = delta;
                varmap[yname] = y;
                binary_fixed.insert(yname);
                redundant_var.insert(yname);
			}
        }
        std::cout<<"the number of infeasible constraints is  "<<redundant_var.size()<<std::endl;
        if (redundant_var.size() == 0)
        {
            std::unordered_map<std::string, double> tmpmask = { {"mask",1.0} };
            return { tmpmask,tmpsolution,violate_constraints };
        }

        SCIP_CONS* lb_cut;
        SCIPcreateConsBasicLinear(scip,&lb_cut, "lb_cut", 0, NULL, NULL, neighbor-binary_fixed.size()+1, SCIPinfinity(scip));
        //for (auto fixed : binary_free)
        //{
        //	SCIPaddCoefLinear(scip, lb_cut, varmap.at(fixed), 1);
        //}
        for (auto fixed : binary_fixed)
        {
			SCIPaddCoefLinear(scip, lb_cut, varmap.at(fixed), -1);
		}	        
        SCIPaddCons(scip, lb_cut);
        SCIPwriteOrigProblem(scip, "test.lp", NULL, NULL);
        SCIPsolve(scip);
        std::cout << "first solve status : " << SCIPgetStatus(scip) << " first solve obj : " << SCIPgetPrimalbound(scip) << std::endl;
        SCIP_SOL* endsol = SCIPgetBestSol(scip);
        std::unordered_map<std::string, double> solution;
        if (endsol != NULL) // ����Ƿ���ڽ�
        {
            // ��ȡģ�͵����б���
            SCIP_VAR** vars = SCIPgetVars(scip);
            int nvars = SCIPgetNVars(scip);



            // ����ÿ����������ȡ���ڽ��е�ֵ
            for (int i = 0; i < nvars; ++i)
            {
                std::string var_name = SCIPvarGetName(vars[i]);
                var_name = var_name.substr(2);
                if (redundant_var.find(var_name) == redundant_var.end())
                {
                    //std::cout<< var_name << std::endl;
                    double var_value = SCIPgetSolVal(scip, endsol, vars[i]);
                    solution[var_name] = var_value;
                }
            }
        }
        else
        {
            std::cout << "No feasible solution found." << std::endl;
        }
        std::cout << "first solve size " << solution.size() << std::endl;

        SCIPchgLhsLinear(scip,lb_cut, -SCIPinfinity(scip));
        SCIPchgRhsLinear(scip, lb_cut, neighbor - binary_fixed.size());
        SCIPwriteOrigProblem(scip, "test2.lp", NULL, NULL);
		SCIPsolve(scip);
        std::cout<<"second solve status : "<<SCIPgetStatus(scip)<<" second solve obj : "<<SCIPgetPrimalbound(scip)<<std::endl;
        SCIP_SOL* endsol_copy = SCIPgetBestSol(scip);
        std::unordered_map<std::string, double> solution_copy;
        if (endsol_copy != NULL) // ����Ƿ���ڽ�
        {
            // ��ȡģ�͵����б���
            SCIP_VAR** vars = SCIPgetVars(scip);
            int nvars = SCIPgetNVars(scip);
            // ����ÿ����������ȡ���ڽ��е�ֵ
            for (int i = 0; i < nvars; ++i)
            {
                std::string var_name = SCIPvarGetName(vars[i]);
                var_name = var_name.substr(2);
                if (redundant_var.find(var_name) == redundant_var.end())
                {
                    double var_value = SCIPgetSolVal(scip, endsol_copy, vars[i]);
                    solution_copy[var_name] = var_value;
                }
            }
        }
        else
        {
            std::cout << "No feasible solution found." << std::endl;
        }
        std::cout << "second solve size " << solution_copy.size() << std::endl;
        SCIPfree(&scip);
        return { solution,solution_copy,violate_constraints };
	}


    std::tuple<std::unordered_map<std::string, double>, std::unordered_map<std::string, double>> RepairFun(SCIP* _scip, const std::unordered_map<std::string, double>& tmpsolution,double time)
    {
        int _ncols = SCIPgetNVars(_scip);
        int _nrows = SCIPgetNConss(_scip);
        std::cout << "the num of vars of model is " << _ncols << "  the num of cons of model is " << _nrows << std::endl;


        std::unordered_set<std::string> redundant_var;
        std::unordered_map<std::string, SCIP_VAR*> varmap;
        SCIP* scip = NULL;
        SCIPcreate(&scip);
        // ʹ�� SCIPcopy() ������ԭʼ����
        SCIP_Bool valid;
        SCIPcopyOrig(_scip, scip, NULL, NULL, "_copy", TRUE, FALSE, TRUE, &valid);
        // ����������ʱ��Ϊ3��
        SCIP_RETCODE retcode = SCIPsetRealParam(scip, "limits/time", time);
        if (retcode != SCIP_OKAY)
        {
            SCIPerrorMessage("Error setting parameter limits/time\n");
            SCIPprintError(retcode);
        }
        SCIPsetMessagehdlrQuiet(scip, true);
        int ncols = SCIPgetNVars(scip);
        int nrows = SCIPgetNConss(scip);
        std::cout << "the num of vars of copymodel is " << ncols << "  the num of cons of copymodel is " << nrows << std::endl;

        // ����һ���µĽ�
        SCIP_SOL* sol;
        SCIPcreateSol(scip, &sol, NULL);

        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        for (int i = 0; i < ncols; ++i)
        {
            SCIP_VAR* var = vars_ori[i];
            std::string varname = SCIPvarGetName(var);
            //std::cout << varname << std::endl;
            varmap[varname] = var;

            //change obj
            SCIP_Real newcoef = 0;
            SCIPchgVarObj(scip, var, newcoef);
        }
        // Ϊÿ����������ֵ
        for (const auto& item : tmpsolution)
        {
            SCIPsetSolVal(scip, sol, varmap.at(item.first), item.second);
        }

        // ���ÿ��Լ��
        std::unordered_map<std::string, double> violate_constraints;
        SCIP_CONS** conss = SCIPgetConss(scip);
        for (int i = 0; i < nrows; ++i)
        {
            SCIP_CONS* consdata = conss[i];
            SCIP_RESULT result;
            SCIPcheckCons(scip, consdata, sol, TRUE, TRUE, TRUE, &result);
            if (result == SCIP_INFEASIBLE)
            {
                double tmpvalue = 0;
                double delta = 0;
                // ���Լ����ϵ������Ӧ�ı���
                int nnonz = SCIPgetNVarsLinear(scip, consdata);
                SCIP_VAR** vars = SCIPgetVarsLinear(scip, consdata);
                SCIP_Real* coefs = SCIPgetValsLinear(scip, consdata);
                // ����һ��Լ���еķ���Ԫ��
                for (int j = 0; j < nnonz; ++j)
                {
                    auto temp_name = SCIPvarGetName(vars[j]);
                    //std::cout<<"coef is "<<coefs[j]<<"  x is  "<< tmpsolution.at(temp_name) <<std::endl;
                    tmpvalue += coefs[j] * tmpsolution.at(temp_name);
                }
                unsigned success;
                auto Lb = SCIPconsGetLhs(scip, consdata, &success);
                auto Ub = SCIPconsGetRhs(scip, consdata, &success);
                std::string consname = SCIPconsGetName(consdata);
                if (Lb == Ub)
                {
                    delta = Ub - tmpvalue;
                    //std::cout << Lb << "  " << tmpvalue << "  " << delta << std::endl;
                }
                else if (Lb > -SCIPinfinity(scip))
                {
                    delta = Lb - tmpvalue;
                    //std::cout << Lb << "  " << tmpvalue << "  " <<delta<< std::endl;
                }
                else
                {
                    delta = Ub - tmpvalue;
                    //std::cout << Ub << "  " << tmpvalue << "  " << delta << std::endl;
                }

                //add variable
                SCIP_VAR* y;
                std::string yname = "y_" + consname;
                SCIPcreateVarBasic(scip, &y, yname.c_str(), 0, 1, 1.0, SCIP_VARTYPE_BINARY);
                SCIPaddVar(scip, y);
                SCIPaddCoefLinear(scip, consdata, y, delta);
                //std::cout << "violate constraints name is  " << consname << "delta is " << delta << std::endl;
                violate_constraints[consname] = delta;
                varmap[yname] = y;
                redundant_var.insert(yname);
            }
        }
        std::cout << "the number of infeasible constraints is  " << redundant_var.size() << std::endl;
        if (violate_constraints.size() == 0)
        {
            return { tmpsolution ,violate_constraints };
        }

        SCIPwriteOrigProblem(scip, "test.lp", NULL, NULL);
        SCIPsolve(scip);
        std::cout << "first solve status : " << SCIPgetStatus(scip) << " first solve obj : " << SCIPgetPrimalbound(scip) << std::endl;
        SCIP_SOL* endsol = SCIPgetBestSol(scip);
        std::unordered_map<std::string, double> solution;
        if (endsol != NULL) // ����Ƿ���ڽ�
        {
            // ��ȡģ�͵����б���
            SCIP_VAR** vars = SCIPgetVars(scip);
            int nvars = SCIPgetNVars(scip);

            // ����ÿ����������ȡ���ڽ��е�ֵ
            for (int i = 0; i < nvars; ++i)
            {
                std::string var_name = SCIPvarGetName(vars[i]);
                var_name = var_name.substr(2);
                if (redundant_var.find(var_name) == redundant_var.end())
                {
                    //std::cout<< var_name << std::endl;
                    double var_value = SCIPgetSolVal(scip, endsol, vars[i]);
                    solution[var_name] = var_value;
                }
            }
        }
        else
        {
            std::cout << "No feasible solution found." << std::endl;
        }
        //std::cout << "first solve size " << solution.size() << std::endl;


        SCIPfree(&scip);
        return { solution,violate_constraints };
    }



    double Lbrun(std::string filename)
	{
        auto scip = ReadFileLocalBranch::readmodel(filename);
        auto first_sol = LocalBranching::getSolution(scip);

        std::queue<std::unordered_map<std::string, double>> sols_q;
        sols_q.push(first_sol);

        int count = 0;
        std::unordered_map<std::string, double> feasible_sol;
        std::vector<std::unordered_map<std::string,double>> violate_cons_list;
        while (!sols_q.empty())
        {
            std::cout<<"the num of sols is "<<sols_q.size()<<std::endl;
            auto sol = sols_q.front();
            sols_q.pop();
            auto [lb_solution1,lb_solution2,violate_cons] = LocalBranching::LocalBranchingFun(scip, sol, 20);
            if (lb_solution1.size() == 1)
            {
                feasible_sol = lb_solution2;
                break;
            }
            if (lb_solution1.size() != 0)
            {
                sols_q.push(lb_solution1);
            }
            if (lb_solution2.size() != 0)
            {
                sols_q.push(lb_solution2);
            }
            count++;
            violate_cons_list.emplace_back(violate_cons);

            if (violate_cons_list.size() > 1)
            {
                bool maps_are_equal = std::equal(violate_cons_list.rbegin()->begin(), violate_cons_list.rbegin()->end(),
                    (violate_cons_list.rbegin() + 1)->begin());
                if (maps_are_equal)
                {
                    break;
                }
            }

            if (count > 1000)
            {
                break;
            }
        }
        std::cout << "local branching is end " << std::endl;
        if (feasible_sol.size() > 0)
        {
            std::cout << "find a feasible solution " << std::endl;
        }
        else
        {
            std::cout << "can not  find a feasible solution " << std::endl;
        }

        auto obj = LocalBranching::computeObj(scip, feasible_sol);
        std::cout << "feasible obj is :  " << obj <<"  epoch is "<< count<< std::endl;
        return obj;
	}


    double Repairrun(std::string filename,double time)
	{
        auto scip = ReadFileLocalBranch::readmodel(filename);
        auto first_sol = LocalBranching::getSolution(scip);

        std::queue<std::unordered_map<std::string, double>> sols_q;
        sols_q.push(first_sol);

        int count = 0;
        std::unordered_map<std::string, double> feasible_sol;
        std::vector<std::unordered_map<std::string, double>> violate_cons_list;
        while (!sols_q.empty())
        {
            //std::cout << "the num of sols is " << sols_q.size() << std::endl;
            auto sol = sols_q.front();
            sols_q.pop();
            auto [repari_sol, violate_cons] = RepairFun(scip,sol,time);
            //std::cout << "repair sol size is " << repari_sol.size() << std::endl;
            sols_q.push(repari_sol);
            count++;
            violate_cons_list.emplace_back(violate_cons);
            if (violate_cons.size() == 0)
            {
				feasible_sol = repari_sol;
                std::cout << "end reason is : find a feasible solution" << std::endl;
				break;
            }
            if (repari_sol.size() == 0)
            {
                std::cout << "end reason is : solve time is too little , cannot find a feasible solution" << std::endl;
                break;
            }
            if (violate_cons_list.size() > 1)
            {
                bool maps_are_equal = std::equal(violate_cons_list.rbegin()->begin(), violate_cons_list.rbegin()->end(),
                    (violate_cons_list.rbegin() + 1)->begin());
                if (maps_are_equal)
                {
                    std::cout << "end reason is : Endless loop" << std::endl;
                    break;
                }
            }

            if (count > 1000)
            {
                std::cout << "end reason is : The maximum number of iterations is reached" << std::endl;
                break;
            }
            
        }
        std::cout << "local branching is end " << std::endl;
        if (feasible_sol.size() > 0)
        {
            std::cout << "find a feasible solution " << std::endl;
        }
        else
        {
            std::cout << "can not  find a feasible solution " << std::endl;
        }

        auto obj = LocalBranching::computeObj(scip, feasible_sol);
        std::cout << "feasible obj is :  " << obj << "  epoch is " << count << std::endl;
        return obj;
	}

    std::pair<double,std::unordered_map<std::string, double>> RepairwithModel(SCIP* scip,double time,const std::unordered_map<std::string, double>& first_sol)
    {
        std::cout<<"lagrangean sol size is "<<first_sol.size()<<std::endl;
        std::queue<std::unordered_map<std::string, double>> sols_q;
        sols_q.push(first_sol);

        int count = 0;
        std::unordered_map<std::string, double> feasible_sol;
        std::vector<std::unordered_map<std::string, double>> violate_cons_list;
        while (!sols_q.empty())
        {
            //std::cout << "the num of sols is " << sols_q.size() << std::endl;
            auto sol = sols_q.front();
            sols_q.pop();
            auto [repari_sol, violate_cons] = RepairFun(scip,sol,time);
            //std::cout << "repair sol size is " << repari_sol.size() << std::endl;
            sols_q.push(repari_sol);
            count++;
            violate_cons_list.emplace_back(violate_cons);
            if (violate_cons.size() == 0)
            {
                feasible_sol = repari_sol;
                std::cout << "end reason is : find a feasible solution" << std::endl;
                break;
            }
            if (repari_sol.size() == 0)
            {
                std::cout << "end reason is : solve time is too little , cannot find a feasible solution" << std::endl;
                break;
            }
            if (violate_cons_list.size() > 1)
            {
                bool maps_are_equal = std::equal(violate_cons_list.rbegin()->begin(), violate_cons_list.rbegin()->end(),
                                                 (violate_cons_list.rbegin() + 1)->begin());
                if (maps_are_equal)
                {
                    std::cout << "end reason is : Endless loop" << std::endl;
                    break;
                }
            }

            if (count > 1000)
            {
                std::cout << "end reason is : The maximum number of iterations is reached" << std::endl;
                break;
            }

        }
        std::cout << "local branching is end " << std::endl;
        if (feasible_sol.size() > 0)
        {
            std::cout << "find a feasible solution " << std::endl;
        }
        else
        {
            std::cout << "can not  find a feasible solution " << std::endl;
        }

        auto obj = LocalBranching::computeObj(scip, feasible_sol);
        std::cout << "feasible obj is :  " << obj << "  epoch is " << count << std::endl;

        return {obj,feasible_sol};
    }


}
