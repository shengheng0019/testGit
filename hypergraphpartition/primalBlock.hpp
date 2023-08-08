#pragma once
#include <algorithm>
#include <execution>
#include <iterator>


#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <future>
#include "Vars.h"
#include "Constraints.h"
#include "Objective.h"
#include "DUALblock.h"
#include "GeneralBlock.h"

class block
{
public:
    std::set<uint64_t> Cons;
    std::set<uint64_t> Vars;
    //std::map<int, double> obj;
    block() {}
};



inline bool Set_find_in(std::set<uint64_t> temp_1, std::set<uint64_t> temp_2)
{
    std::set<uint64_t> diff;
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
inline std::set<uint64_t> Set_diff(std::set<uint64_t> temp_1, std::set<uint64_t> temp_2)
{
    std::set<uint64_t> diff;
    set_difference(temp_1.begin(), temp_1.end(), temp_2.begin(), temp_2.end(), inserter(diff, diff.begin()));
    return diff;
}


inline std::map<uint64_t, block> PrimalBlock(std::unordered_map<int, std::set<uint64_t>> parts, std::unordered_map<int, std::vector<uint64_t>> matrix)
{
    std::map<uint64_t, block> Blocks;

    for (int i = 0; i < parts.size(); ++i)
    {
        Blocks[i].Vars = parts[i];
    }
    for (int i = 0; i < matrix.size(); ++i)
    {
	    std::set<uint64_t> nodes(matrix[i].begin(),matrix[i].end());
        bool mark = 0;
        for (int j = 0; j < parts.size(); ++j)
        {
	        if (Set_find_in(nodes,parts[j]))
	        {
		        Blocks[j].Cons.emplace(i);
		        mark = 1;
		        break;
	        }
        }
        if (mark == 0)
        {
            Blocks[parts.size()].Cons.emplace(i);
        }
    }


    //�������һ��block�Ľڵ���Ϣ
    for (auto it = Blocks[parts.size()].Cons.begin(); it != Blocks[parts.size()].Cons.end(); ++it)
    {
	    std::set<uint64_t> nodes(matrix[*it].begin(), matrix[*it].end());
        for (auto value : nodes)
        {
            Blocks[parts.size()].Vars.emplace(value);
        }
    }
    return Blocks;
}




//���̴߳�����Ч�ʲ��� 68ֵ��һ��
inline std::map<uint64_t, block> PrimalBlock2(std::unordered_map<int, std::set<uint64_t>> parts, std::unordered_map<int, std::vector<uint64_t>> matrix)
{
    std::map<uint64_t, block> Blocks;

    for (int i = 0; i < parts.size(); ++i)
    {
        Blocks[i].Vars = parts[i];
    }

    std::for_each(std::execution::par, matrix.begin(), matrix.end(), [&](auto& row) {
        std::set<uint64_t> nodes(row.second.begin(), row.second.end());
        bool mark = 0;
        for (int j = 0; j < parts.size(); ++j)
        {
            if (Set_find_in(nodes, parts[j]))
            {
                Blocks[j].Cons.emplace(row.first);
                mark = 1;
                break;
            }
        }
        if (mark == 0)
        {
            Blocks[parts.size()].Cons.emplace(row.first);
        }
        });

    //�������һ��block�Ľڵ���Ϣ
    for (auto it = Blocks[parts.size()].Cons.begin(); it != Blocks[parts.size()].Cons.end(); ++it)
    {
	    std::set<uint64_t> nodes(matrix[*it].begin(), matrix[*it].end());
        for (auto value : nodes)
        {
            Blocks[parts.size()].Vars.emplace(value);
        }
    }
    return Blocks;
}


//��ӻ����� ��89��ʼ
inline std::pair<std::map<int, block>, int>  PrimalBlock3(std::unordered_map<int, std::set<uint64_t>> parts, std::unordered_map<int, std::vector<uint64_t>> matrix)
{
    std::map<int, block> Blocks;
    std::mutex mtx; // ������

    for (int i = 0; i < parts.size(); ++i)
    {
        Blocks[i].Vars = parts[i];
    }

    std::for_each(std::execution::par, matrix.begin(), matrix.end(), [&](auto& row) {
        std::set<uint64_t> nodes(row.second.begin(), row.second.end());
        bool mark = 0;
        for (int j = 0; j < parts.size(); ++j)
        {
            if (Set_find_in(nodes, parts[j]))
            {
                std::unique_lock<std::mutex> lock(mtx); // ����
                Blocks[j].Cons.emplace(row.first);
                mark = 1;
                lock.unlock(); // ����
                break;
            }
        }
        if (mark == 0)
        {
            std::unique_lock<std::mutex> lock(mtx); // ����
            Blocks[parts.size()].Cons.emplace(row.first);
            lock.unlock(); // ����
        }
        });

    //�������һ��block�Ľڵ���Ϣ
    for (auto it = Blocks[parts.size()].Cons.begin(); it != Blocks[parts.size()].Cons.end(); ++it)
    {
        std::set<uint64_t> nodes(matrix[*it].begin(), matrix[*it].end());
        for (auto value : nodes)
        {
            std::unique_lock<std::mutex> lock(mtx); // ����
            Blocks[parts.size()].Vars.emplace(value);
            lock.unlock(); // ����
        }
    }
    int num_of_block = 0;
	for (auto it = Blocks.begin(); it != Blocks.end(); ++it)
	{
		if (it->second.Cons.size() > 0)
		{
			num_of_block++;
		}
	}

    return {Blocks,num_of_block};
}


inline std::map<int, General_Block> TransformerBlock(SCIP* _scip, std::map<int, block> origin_block, std::unordered_map<int, SCIP_VAR*> Var_map, std::unordered_map<int, SCIP_CONS*> Cons_map, std::unordered_map<std::string, int> Var_Ni)
{

    std::map<int, General_Block> Blocks;
    std::unordered_map<int, std::string> Var_mapping_all;
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
        std::unordered_map<int, Vars> Var_mapping;
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
        if (it->first < origin_block.size() - 1)
        {
            for (auto v : Var_mapping)
            {
                Var_mapping_all[v.first] = v.second.Name;
            }
        }

        unsigned success;
        auto cons = it->second.Cons;
        for (auto c : cons)
        {
            Constraints temp_cons;
            auto Lb = SCIPconsGetLhs(_scip, Cons_map[c], &success);
            auto Ub = SCIPconsGetRhs(_scip, Cons_map[c], &success);
            auto consname = SCIPconsGetName(Cons_map[c]);
            temp_cons.Name = consname;
            //cout << temp_cons.Name << endl;
            if (Lb == Ub)
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
        Blocks.emplace(part_index, endBlock);
        Blocks[part_index - 1] = temp_block_withoutcons;
    }


    return Blocks;
}


inline std::vector<DUALBlock> DualBlock(std::unordered_map<int, std::set<uint64_t>> parts_, std::unordered_map<int, std::vector<uint64_t>> matrix,std::unordered_map<int, std::vector<uint64_t>> dual_matrix, std::unordered_map<int, SCIP_VAR*> Var_map, std::unordered_map<int, SCIP_CONS*> Cons_map)
{
    //parts could be 0 2 3 4 5
    int num_of_block = 0;
    std::unordered_map<int, std::set<uint64_t>> parts;
    for (auto tmppart : parts_)
    {
        parts[num_of_block] = tmppart.second;
        num_of_block++;
    }
    ///
    std::vector<DUALBlock> tempBlocks;
    std::set<uint64_t> coupling_vars;
    std::mutex mtx; // ������

    std::vector<std::set<uint64_t>> Block_all_vars(parts.size());

	for (auto temp_block : parts)
	{
		DUALBlock tmpblock;
		for (auto tmpcon : temp_block.second)
		{
			tmpblock.Cons.emplace_back(Cons_map[tmpcon]);

            auto tmp_con_vars = matrix[tmpcon];
            Block_all_vars[temp_block.first].insert(tmp_con_vars.begin(), tmp_con_vars.end());
		}
        tempBlocks.emplace_back(tmpblock);
	}


    std::for_each(std::execution::par, dual_matrix.begin(), dual_matrix.end(), [&](auto& row) {
        std::set<uint64_t> nodes(row.second.begin(), row.second.end());
        bool mark = 0;
        for (int j = 0; j < parts.size(); ++j)
        {
            if (Set_find_in(nodes, parts[j]))
            {
                std::unique_lock<std::mutex> lock(mtx); // ����
                tempBlocks[j].Independent_vars.emplace_back(Var_map[row.first]);
                mark = 1;
                lock.unlock(); // ����
                break;
            }
        }
        if (mark == 0)
        {
            std::unique_lock<std::mutex> lock(mtx); // ����
            coupling_vars.emplace(row.first);
            lock.unlock(); // ����
        }
        });


    /*for (int i = 0; i < tempBlocks.size(); ++i)
    {
        auto this_block_vars = Block_all_vars[i];
        for (auto tmpvar : this_block_vars)
        {
	        if (coupling_vars.find(tmpvar) != coupling_vars.end())
	        {
	        	tempBlocks[i].Coupling_vars.emplace_back(Var_map[tmpvar]);
	        }
        }
    }*/

	for (int i = 0; i < tempBlocks.size(); i++)
	{
		for (auto tmpcouping : coupling_vars)
		{
			tempBlocks[i].Coupling_vars.emplace_back(Var_map[tmpcouping]);
		}
	}

    int dual_vars = 0;
    int dual_cons = 0;
    for (auto tempdualB : tempBlocks)
    {
        std::cout << "the num of vars in block " << " : " << tempdualB.Independent_vars.size() << "the num of cons in block " << tempdualB.Cons.size() << std::endl;
        dual_vars += tempdualB.Independent_vars.size();
        dual_cons += tempdualB.Cons.size();
    }
    dual_vars += coupling_vars.size();
    std::cout << "the num of vars after dual : " << dual_vars << "-" << Var_map.size() << std::endl;
    std::cout << "the num of cons after dual : " << dual_cons << "-" << Cons_map.size() << std::endl;

    return  tempBlocks;

}