#pragma once

#include "gurobi_c++.h"
#include <map>
#include "primalBlock.hpp"


double random() {
    // 生成随机数
    double r = (double)rand() / RAND_MAX;
    return r;
}


using namespace std;


tuple<map<int, block>, int, map<int, block>,int> Dection(string filename, int k)
{
    vector<vector<int64_t>> matrix;
    
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env,filename);
    //model.getEnv().set(GRB_IntParam_Presolve, 1);
    //model.update();
    //model.presolve();
    std::unordered_map<std::string, int> var_index_map;
    // 遍历所有变量，保存变量名和对应的索引
    for (int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i) {
        GRBVar v = model.getVar(i);
        var_index_map[v.get(GRB_StringAttr_VarName)] = i;
    }
    int64_t min_Cons = model.get(GRB_IntAttr_NumVars);
    int num_cons = model.get(GRB_IntAttr_NumConstrs);

    vector<vector<int64_t>> matrix_dual(model.get(GRB_IntAttr_NumVars));
    //cout << num_cons << "--1111" << endl;
    // 遍历每个约束
    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i) 
    {
        vector<int64_t> temp;
        GRBConstr c = model.getConstr(i);
        //cout << "Constraint " << c.get(GRB_StringAttr_ConstrName) << ":" << endl;
        // 获取约束矩阵的第i行非零元素
        GRBLinExpr expr = model.getRow(c);
        int nzcount = expr.size();
        if (nzcount<min_Cons)
        {
	        min_Cons = nzcount;
        }
        for (int j = 0; j < nzcount; ++j) 
        {
            //cout << expr.getCoeff(j) << " * ";
            //cout << expr.getVar(j).get(GRB_StringAttr_VarName);
            auto var_Id = var_index_map[expr.getVar(j).get(GRB_StringAttr_VarName)];
            temp.emplace_back(var_Id);
            matrix_dual[var_Id].emplace_back(i);
        }
        matrix.emplace_back(temp);
        //cout << endl;
    }
    //model.write("test2.lp");


    vector<part::Partition> parts = RunHYPE(matrix, k);
    map<int, block> Blocks = PrimalBlock(parts,matrix);

    vector<part::Partition> parts2 = RunHYPE(matrix_dual, k);
    map<int, block> Blocks2 = DualBlock(parts2, matrix_dual);

    return make_tuple(Blocks, matrix.size(), Blocks2,matrix_dual.size());
}

