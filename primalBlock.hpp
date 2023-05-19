#pragma once
#include "Partition.hpp"
#include <map>
#include <set>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "HYPEm.hpp"
#include "Vars.h"
#include "Constraints.h"
#include "Objective.h"

class ADMMreturn
{
public:
    Scip* _Scip;
    vector<SCIP_VAR*> _ori_all_vars;
    vector<SCIP_CONS*> _ori_all_cons;
    vector<SCIP_VAR*> _ori_vars_x;
    vector<SCIP_VAR*> _ori_vars_z;
};

class General_Block
{
public:
    std::vector<Constraints> bCons;
    std::vector<Vars> bVars;
    Objective bobj;
};

inline bool Set_find_in(std::set<int64_t> temp_1,std::set<int64_t> temp_2)
{
    std::set<int64_t> diff;
    set_difference(temp_1.begin(), temp_1.end(), temp_2.begin(), temp_2.end(), inserter(diff, diff.begin()));
    if (diff.empty())
    {
        return true;
    }
    else
    {
	    return false;
    }
}

//temp_1 - temp_2
inline std::set<int64_t> Set_diff(std::set<int64_t> temp_1, std::set<int64_t> temp_2)
{
    std::set<int64_t> diff;
    set_difference(temp_1.begin(), temp_1.end(), temp_2.begin(), temp_2.end(), inserter(diff, diff.begin()));
    return diff;
}


inline std::pair<std::map<int, block>,int> PrimalBlock(std::vector<part::Partition> parts, std::unordered_map<int, std::vector<uint64_t>> matrix, std::unordered_set<int64_t>& couping1)
{
    int num_of_block = 1;
    std::map<int, block> Blocks;
    std::set<int64_t> last_block_cons(couping1.begin(), couping1.end());
    for (int i = 0; i < parts.size(); ++i)
    {
        std::set<int64_t> nodelist(parts[i].getNodes().begin(), parts[i].getNodes().end());
        //cout << "第" << temp_block.getId() << "块:有" << nodelist.size() << "个变量" << endl;
        block temp_b;
        temp_b.Vars = nodelist;
        Blocks.emplace(parts[i].getId(), temp_b);

        std::set<int64_t> getEdges(parts[i].getEdges().begin(), parts[i].getEdges().end());
        auto need_allot = Set_diff(getEdges, last_block_cons);
        if (need_allot.size()!=0)
        {
            num_of_block++;
        }
        Blocks[i].Cons = need_allot;
    }

    block temp_b;
    temp_b.Cons = last_block_cons;
    Blocks.emplace(parts.size(), temp_b);
    //处理最后一个block的节点信息
    for (auto it = Blocks[parts.size()].Cons.begin(); it != Blocks[parts.size()].Cons.end(); ++it)
    {
        auto nodelist = matrix[*it];
        for (auto value : nodelist)
        {
            Blocks[parts.size()].Vars.emplace(value);
        }
    }
    return { Blocks ,num_of_block};
}


inline std::tuple<vector<SCIP_VAR*>,vector<SCIP_VAR*>> DualBlock(std::vector<part::Partition> parts, std::unordered_map<int, std::vector<uint64_t>> matrix_dual, std::unordered_set<int64_t>& couping2, std::unordered_map<int, SCIP_VAR*> Var_map)
{
    int num_of_block = 1;
    std::map<int, block> Blocks;
    std::set<int64_t> last_block_vars(couping2.begin(),couping2.end());
    vector<SCIP_VAR*> ori_vars_x;
    vector<SCIP_VAR*> ori_vars_z;
    for (auto v : Var_map)
    {
	    if (last_block_vars.find(v.first)!= last_block_vars.end())
	    {
            ori_vars_z.emplace_back(v.second);
	    }
	    else
	    {
            ori_vars_x.emplace_back(v.second);
	    }
    }



    //for (int i = 0; i < parts.size(); ++i)
    //{
    //    std::set<int64_t> nodelist(parts[i].getNodes().begin(), parts[i].getNodes().end());
    //    //cout << "第" << temp_block.getId() << "块:有" << nodelist.size() << "个变量" << endl;
    //    block temp_b;
    //    temp_b.Cons = nodelist;
    //    Blocks.emplace(parts[i].getId(), temp_b);
    //    std::set<int64_t> getEdges(parts[i].getEdges().begin(), parts[i].getEdges().end());
    //    auto need_allot = Set_diff(getEdges, last_block_vars);
    //    if (need_allot.size() != 0)
    //    {
    //        num_of_block++;
    //    }
    //    Blocks[i].Vars = need_allot;

    //}
    //block temp_b;
    //temp_b.Vars = last_block_vars;
    //Blocks.emplace(parts.size(), temp_b);

    ////处理最后一个block的节点信息
    //for (auto it = Blocks[parts.size()].Vars.begin(); it != Blocks[parts.size()].Vars.end(); ++it)
    //{
    //    auto nodelist = matrix_dual[*it];
    //    for (auto value : nodelist)
    //    {
    //        Blocks[parts.size()].Cons.emplace(value);
    //    }
    //}


    return { ori_vars_x ,ori_vars_z };
}

inline std::map<int, block> BendersBlock(std::vector<part::Partition> parts, std::vector<std::vector<int64_t>> matrix_continuous, std::vector<std::vector<int64_t>> matrix_int,std::unordered_set<int64_t>& couping2,
								  std::map<int64_t, int64_t> map_continuous_index, std::map<int64_t, int64_t> map_int_index)
{
    std::map<int, block> Blocks;
    std::set<int64_t> last_block_vars(couping2.begin(), couping2.end());

    for (int i = 0; i < parts.size(); ++i)
    {
        std::set<int64_t> nodelist(parts[i].getNodes().begin(), parts[i].getNodes().end());
        //cout << "第" << temp_block.getId() << "块:有" << nodelist.size() << "个变量" << endl;
        block temp_b;
        temp_b.Cons = nodelist;
        Blocks.emplace(parts[i].getId(), temp_b);
        std::set<int64_t> getEdges(parts[i].getEdges().begin(), parts[i].getEdges().end());
        auto need_allot = Set_diff(getEdges, last_block_vars);
        std::set<int64_t> true_vars;
        for (auto v : need_allot)
        {
            true_vars.emplace(map_continuous_index[v]);
        }
        Blocks[i].Vars = true_vars;
        //for (auto edge : need_allot)
        //{
        //    std::set<int64_t> temp_cons(matrix_dual[edge].begin(), matrix_dual[edge].end());
        //    if (Set_find_in(temp_cons, nodelist))
        //    {
        //        Blocks[i].Vars.emplace(edge);
        //    }
        //    else
        //    {
        //        last_block_vars.emplace(edge);
        //    }
        //    have_loading.emplace(edge);
        //}
    }

    std::set<int64_t> true_coupling;
    std::set<int64_t> temp_cons;
    for (auto v :last_block_vars )
    {
	    true_coupling.emplace(map_continuous_index[v]);
        auto nodelist = matrix_continuous[v];
	    for (auto value : nodelist)
	    {
		    temp_cons.emplace(value);
	    }
    }
    for (int i = 0; i < matrix_int.size(); ++i)
    {
        true_coupling.emplace(map_int_index[i]);
        auto nodelist = matrix_int[i];
        for (auto value : nodelist)
        {
	        temp_cons.emplace(value);
        }
    }
    
    block temp_b;
    temp_b.Vars = true_coupling;
    temp_b.Cons = temp_cons;
    Blocks.emplace(parts.size(), temp_b);

    return Blocks;
}



inline std::map<int,General_Block> TransformerBlock(SCIP* _scip,std::map<int,block> origin_block, std::unordered_map<int, SCIP_VAR*> Var_map,std::unordered_map<int, SCIP_CONS*> Cons_map,int num_of_part, std::unordered_map<std::string, int> Var_Ni)
{

	std::map<int, General_Block> Blocks;
    std::unordered_map<int, string> Var_mapping_all;
    General_Block temp_block_withoutcons;
    int part_index = 0;
	for (auto it = origin_block.begin(); it != origin_block.end(); ++it)
	{
		if (it->second.Cons.size() == 0)
		{
            auto vars = it->second.Vars;
            for (auto var : vars)
            {
                Vars temp_var;
                auto varType = SCIPvarGetType(Var_map[var]);
                if (varType == SCIP_VARTYPE_BINARY)
                {
                    temp_var.Type = VarType::Bool;
                }
                else if (varType == SCIP_VARTYPE_INTEGER)
                {
                    temp_var.Type = VarType::Int;
                }
                else if (varType == SCIP_VARTYPE_CONTINUOUS)
                {
                    temp_var.Type = VarType::Num;
                }
                temp_var.Name = SCIPvarGetName(Var_map[var]);
                //temp_var.Id = SCIPvarGetIndex(Var_map[var]);
                temp_var.Id = Var_Ni[temp_var.Name];
                temp_var.Lb = SCIPvarGetLbOriginal(Var_map[var]);
                temp_var.Ub = SCIPvarGetUbOriginal(Var_map[var]);
                temp_block_withoutcons.bVars.emplace_back(temp_var);
                temp_block_withoutcons.bobj.coef[var] = SCIPvarGetObj(Var_map[var]);
            }
			continue;
		}
        std::unordered_map<int,Vars> Var_mapping;
		General_Block temp_block;
		auto vars = it->second.Vars;
		for (auto var : vars)
		{
            Vars temp_var;
            auto varType = SCIPvarGetType(Var_map[var]);
            if (varType == SCIP_VARTYPE_BINARY)
            {
	            temp_var.Type = VarType::Bool;
			}
			else if (varType == SCIP_VARTYPE_INTEGER)
			{
				temp_var.Type = VarType::Int;
			}
			else if (varType == SCIP_VARTYPE_CONTINUOUS)
			{
				temp_var.Type = VarType::Num;
            }
            temp_var.Name = SCIPvarGetName(Var_map[var]);
            //temp_var.Id = SCIPvarGetIndex(Var_map[var]);
            temp_var.Id = Var_Ni[temp_var.Name];
            temp_var.Lb = SCIPvarGetLbOriginal(Var_map[var]);
            temp_var.Ub = SCIPvarGetUbOriginal(Var_map[var]);
            temp_block.bVars.emplace_back(temp_var);
            temp_block.bobj.coef[var] = SCIPvarGetObj(Var_map[var]);
            Var_mapping[temp_var.Id] = temp_var;
		}
        if (it->first < origin_block.size()-1)
        {
            for (auto v : Var_mapping)
            {
                Var_mapping_all[v.first] = v.second.Name;
            }
        }

        unsigned success;
        auto cons= it->second.Cons;
		for (auto c :cons )
		{
            Constraints temp_cons;
            auto Lb = SCIPconsGetLhs(_scip,Cons_map[c],&success);
            auto Ub = SCIPconsGetRhs(_scip,Cons_map[c],&success);
            auto consname = SCIPconsGetName(Cons_map[c]);
            temp_cons.Name = consname;
            //cout << temp_cons.Name << endl;
            if (Lb==Ub)
            {
	            temp_cons.Prop = PROP::eq;
                temp_cons.rhs = Ub;
            }
            else if (Lb > -SCIPinfinity(_scip))
            {
	            temp_cons.Prop = PROP::geq;
				temp_cons.rhs = Lb;
            }
            else
            {
	            temp_cons.Prop = PROP::leq;
                temp_cons.rhs = Ub;
            }
            int nnonz = SCIPgetNVarsLinear(_scip, Cons_map[c]);
            SCIP_VAR** vars = SCIPgetVarsLinear(_scip, Cons_map[c]);
            SCIP_Real* coefs = SCIPgetValsLinear(_scip, Cons_map[c]);
            for (int i = 0; i < nnonz; ++i)
            {
                auto tempname = SCIPvarGetName(vars[i]);
	            temp_cons.exprDic[Var_mapping[Var_Ni[tempname]]] = coefs[i];
                temp_cons.idDic[Var_Ni[tempname]] = coefs[i];
            }
            temp_block.bCons.emplace_back(temp_cons);
		}
        Blocks.emplace(part_index, temp_block);
        part_index++;
	}


    auto endBlock = Blocks[part_index - 1];
	if (temp_block_withoutcons.bVars.size() != 0)
	{
        Blocks.emplace(part_index,endBlock);
        Blocks[part_index-1] = temp_block_withoutcons;
	}


	return Blocks;
}


inline void Check(std::map<int, General_Block> blocks, Objective master_obj)
{
    SCIP* scipModel;
    SCIPcreate(&scipModel);
    SCIPincludeDefaultPlugins(scipModel);
    SCIPcreateProbBasic(scipModel, "subProb");
    SCIPsetObjsense(scipModel, SCIP_OBJSENSE_MINIMIZE);

    int num_of_vars = 0;
    int num_of_cons = 0;
    unordered_map<int, SCIP_VAR*> scipVarMap;
    for (int i = 0; i < blocks.size()-1; ++i)
    {
	    for (auto v :blocks[i].bVars)
	    {
            SCIP_VAR* tempVar;
            VarType varType = v.Type;
            SCIP_VARTYPE scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
            if (varType == VarType::Bool)
            {
                scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_BINARY;
            }
            else if (varType == VarType::Int)
            {
                scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_INTEGER;
            }
            else if (varType == VarType::Num)
            {
                scipVarType = SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS;
            }

            SCIPcreateVarBasic(scipModel, &tempVar, NULL,v.Lb, v.Ub,master_obj.coef[v.Id], scipVarType);
            SCIPaddVar(scipModel, tempVar);
            scipVarMap[v.Id] = tempVar;
	    }
    }
    for (int i = 0; i < blocks.size(); ++i)
    {
	    for (auto c :blocks[i].bCons)
	    {
		    SCIP_CONS* tempCons;
            double tempLb = c.rhs, tempUb = SCIPinfinity(scipModel);
            if (c.Prop == PROP::eq)
            {
                tempLb = c.rhs;
                tempUb = c.rhs;
            }
            else if (c.Prop == PROP::leq)
            {
                tempUb = c.rhs;
                tempLb = -SCIPinfinity(scipModel);
            }

			SCIPcreateConsBasicLinear(scipModel, &tempCons, c.Name.c_str(), 0, NULL, NULL, tempLb, tempUb);
			for (auto v : c.exprDic)
			{
				SCIPaddCoefLinear(scipModel, tempCons, scipVarMap[v.first.Id], v.second);
			}
			SCIPaddCons(scipModel, tempCons);
	    }
    }
    SCIPsolve(scipModel);



}
