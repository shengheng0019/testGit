#pragma once
#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <map>
#include <set>
#include <iostream>
#include <future>
#include "Minipartplus5.hpp"
#include "primalBlock.hpp"
#include "HYPEm.hpp"
#include <condition_variable>


namespace ParaDetec
{
    inline std::unordered_map<int, std::set<uint64_t>> ParallelDetection(std::unordered_map<int, std::vector<uint64_t>> matrix, std::unordered_map<int, std::vector<uint64_t>> matrix_dual,int k)
    {
        // Start tasks in separate threads.
        std::future<decltype(RunHYPE(matrix, k))> fut1 = std::async(std::launch::async, RunHYPE, matrix, k);
        std::future<decltype(MiniPartStage1(matrix, matrix_dual, k))> fut2 = std::async(std::launch::async, MiniPartStage1, matrix, matrix_dual, k);
        // Get the results.
        auto initial_parts = fut1.get();
        auto [vm, hg] = fut2.get();

        auto parts = MiniPartStage2(vm, hg, k, initial_parts);

        return parts;
    }


    inline std::map<int,block> OnePartition(std::unordered_map<int, std::vector<uint64_t>> matrix, std::unordered_map<int, std::vector<uint64_t>> matrix_dual, int k)
    {
        auto parts = ParaDetec::ParallelDetection(matrix, matrix_dual, k);
        auto [Blocks, num_of_part] = PrimalBlock3(parts, matrix);
        return Blocks;
    }


    class waitingPlan
    {
    public:
        std::map<int, block> waitingBlock;
        double score;
        int num_of_blocks;
        double couping;

        waitingPlan(std::map<int, block> waitingBlock)
        {
            this->waitingBlock = waitingBlock;
            this->score = 0;
            this->num_of_blocks = 0;
            this->couping = 0;
        }

        void computeScore()
        {
            int num_of_blocks = this->waitingBlock.size(); //the end of this is couping
            std::vector<int> consnum_list;
            for (const auto& tmpblock :this->waitingBlock)
            {
                if (tmpblock.second.Cons.size() > 0)
                {
                    consnum_list.emplace_back(tmpblock.second.Cons.size());
                    this->num_of_blocks++;
                }
            }
            double sum = std::accumulate(consnum_list.begin(), consnum_list.end(), 0.0);
            double mean = sum / consnum_list.size();
            for(const auto& tmpint : consnum_list)
            {
                this->score += (tmpint - mean) * (tmpint - mean);
            }
            this->couping = 100.0 * (this->waitingBlock[num_of_blocks-1].Cons.size()) / sum;
        }

        void compute_agg_score(double epsilon)
        {
            int num_of_blocks = this->waitingBlock.size(); //the end of this is couping
            int total_cons = 0;
            for (auto tmp_block : this->waitingBlock)
            {
                if (tmp_block.second.Cons.size() > 0)
                {
                    this->num_of_blocks++;
                    total_cons += tmp_block.second.Cons.size();
                }
            }
            this->couping = 100.0 * (this->waitingBlock[num_of_blocks - 1].Cons.size()) / total_cons;
            this->score = (epsilon * this->couping) / this->num_of_blocks;
        }
    };


    inline std::map<int, block> SelfAdaPartition(std::vector<int> partk_list,std::unordered_map<int, std::vector<uint64_t>> matrix, std::unordered_map<int, std::vector<uint64_t>> matrix_dual)
    {
        std::vector<waitingPlan> waitingPlanList;
        const int mincons_one_part = 50;
        const int minvars_one_part = 50;
        int max_num_part = std::min(matrix.size() / mincons_one_part, matrix_dual.size() / minvars_one_part);
        std::vector<int> partk_list_copy;
        for (const auto& tmpint: partk_list)
        {
            if (tmpint < max_num_part)
            {
                partk_list_copy.emplace_back(tmpint);
            }
        }
        if (partk_list_copy.size() == 0)
        {
            partk_list_copy.emplace_back(2);
        }
        std::cout<<"----partition size is : "<<partk_list_copy.size()<<"--------"<< std::endl;
        // 创建一个future的vector
        std::vector<std::future<std::map<int, block>>> futures;
        for (int k : partk_list_copy)
        {
            futures.push_back(std::async(std::launch::async, &OnePartition, matrix, matrix_dual, k));
        }
        for (auto& future : futures)
        {
            std::map<int, block> result = future.get();
            waitingPlan tmpPlan(result);
            tmpPlan.compute_agg_score(1);
            waitingPlanList.emplace_back(tmpPlan);
        }
        auto minElement = std::min_element(waitingPlanList.begin(), waitingPlanList.end(),
                                           [](const waitingPlan& a, const waitingPlan& b) {
                                               return a.score < b.score;
                                           }
        );

        waitingPlan minScorePlan = *minElement;
        std::cout<< "minScorePlan.coupling is  : " << minScorePlan.couping<<"   minScorePlan.block is  : "<<minScorePlan.num_of_blocks<<"   minScorePlan.score is  "<< minScorePlan.score <<" total plan :  "<< waitingPlanList.size()<< std::endl;
        return minScorePlan.waitingBlock;
    }

}



inline std::tuple<std::map<int, General_Block>, SCIP*, vector<SCIP_VAR*>, vector<SCIP_CONS*>, std::vector<DUALBlock>, std::unordered_map<std::string, int>>  Dection(std::string filename, int k)
{
    vector<int> partk_list = { 2,8,32};
    //filename = "D:\\Codes\\VsProject\\SCIPlinear\\test.lp";
    SCIP* scip = NULL;
    SCIPcreate(&scip);
    SCIPincludeDefaultPlugins(scip);
    // 设置一些求解参数
    SCIPsetIntParam(scip, "presolving/maxrounds", 0);
    //将string类型的filename转换成const char*类型
    auto filename1 = filename.c_str();
    // 读入模型的mps文件
    SCIPreadProb(scip, filename1, NULL);
    std::unordered_map<int, SCIP_VAR*> Var_map;
    std::unordered_map<int, SCIP_CONS*> Cons_map;
    std::unordered_map<std::string, int> Var_Ni;
    //std::unordered_map<std::string, int> Cons_Ni;
    int nrows = SCIPgetNConss(scip);
    int ncols = SCIPgetNVars(scip);
    std::cout << ncols << "," << nrows << std::endl;
    if (ncols > 500000 || nrows > 500000)
    {
        SCIPfree(&scip);
        return std::make_tuple(std::map<int, General_Block>(), scip, vector<SCIP_VAR*>(), vector<SCIP_CONS*>(), std::vector<DUALBlock>(), std::unordered_map<std::string, int>());
    }
    std::unordered_map<int, std::vector<uint64_t>> matrix;
    std::unordered_map<int, std::vector<uint64_t>> matrix_dual;
    SCIP_CONS** conss = SCIPgetConss(scip);
    SCIP_VAR** vars_ori = SCIPgetVars(scip);
    std::unordered_set<std::string> name_of_vars;
    for (int i = 0; i < nrows; ++i)
    {
        // 获得线性约束的数据
        SCIP_CONS* consdata = conss[i];
        //cout << "约束 " << SCIPconsGetName(consdata) << "：" << endl;
        Cons_map[i] = consdata;

        int nnonz = SCIPgetNVarsLinear(scip, consdata);
        SCIP_VAR** vars = SCIPgetVarsLinear(scip, consdata);
        SCIP_Real* coefs = SCIPgetValsLinear(scip, consdata);
        std::vector<uint64_t> temp_v;
        std::set<int> temp_set;

        for (int j = 0; j < nnonz; ++j)
        {
            auto temp_name = SCIPvarGetName(vars[j]);
            if (!name_of_vars.count(temp_name))
            {
                Var_Ni[temp_name] = name_of_vars.size();
                Var_map[name_of_vars.size()] = vars[j];
                name_of_vars.emplace(temp_name);

            }
            //cout << Var_Ni[temp_name] << endl;
            temp_v.emplace_back(Var_Ni[temp_name]);
            matrix_dual[Var_Ni[temp_name]].emplace_back(i);
            //SCIP_VARTYPE_CONTINUOUS
        }
        matrix[i] = temp_v;
    }

    //SCIPwriteOrigProblem(scip,"test.lp",NULL,NULL);
    std::cout << matrix_dual.size() << "," << matrix.size() << std::endl;
    std::cout<<"matrix_varsize: "<< name_of_vars.size()<< std::endl;


    vector<SCIP_CONS*> fj_cons(conss, conss + nrows);
    vector<SCIP_VAR*> fj_vars(vars_ori, vars_ori + ncols);




    //parallel primal and dual
    //std::async
    auto future1 = std::async(std::launch::async, [&]() {
        // 处理 "primal block"
//        auto tmpBlocks = ParaDetec::OnePartition(matrix, matrix_dual, k);
//        return TransformerBlock(scip, tmpBlocks, Var_map, Cons_map, Var_Ni);
        auto Blocks = ParaDetec::SelfAdaPartition(partk_list,matrix,matrix_dual);
        return TransformerBlock(scip, Blocks, Var_map, Cons_map, Var_Ni);
    });

    auto future2 = std::async(std::launch::async, [&]() {

        auto parts_dual = ParaDetec::ParallelDetection(matrix_dual, matrix, k);
        return DualBlock(parts_dual, matrix, matrix_dual, Var_map, Cons_map);
    });


    std::map<int, General_Block> GeneralBlocks;
    std::vector<DUALBlock> Redualblocks;


    GeneralBlocks = future1.get();
    Redualblocks = future2.get();




    int num_col = 0;
    int num_row = 0;
    for (int i = 0; i < GeneralBlocks.size() - 1; ++i)
    {
        num_col += GeneralBlocks[i].bVars.size();
        num_row += GeneralBlocks[i].bCons.size();
    }
    num_row += GeneralBlocks[GeneralBlocks.size() - 1].bCons.size();
    std::cout << "num of vars: " << ncols << "-" << num_col << std::endl;
    std::cout << "num of cons: " << nrows << "-" << num_row << std::endl;
    for (int v = 0; v < GeneralBlocks.size(); ++v)
    {
        std::cout << "the " << v << "block : num of cons :" << GeneralBlocks[v].bCons.size() << " num of vars: " << GeneralBlocks[v].bVars.size() << std::endl;
    }

    return make_tuple(GeneralBlocks, scip, fj_vars, fj_cons, Redualblocks, Var_Ni);
}


