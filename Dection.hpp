#pragma once

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <map>
#include <set>
#include <iostream>
#include <set>

#include "HYPEm.hpp"
#include "primalBlock.hpp"



using namespace std;




std::tuple<std::map<int, General_Block>,SCIP*, vector<SCIP_VAR*>, vector<SCIP_CONS*>, vector<SCIP_VAR*>, vector<SCIP_VAR*>> Dection(string filename,int k)
{
    //filename = "D:\\Codes\\VsProject\\SCIPlinear\\test.lp";
    SCIP* scip = NULL;
    SCIPcreate(&scip);
    SCIPincludeDefaultPlugins(scip);
    // ����һЩ������
    SCIPsetIntParam(scip, "presolving/maxrounds", 0);
    SCIPsetBoolParam(scip, "constraints/linear/upgrade/indicator", false);
    //��string���͵�filenameת����const char*����
    auto filename1 = filename.c_str();
    // ����ģ�͵�mps�ļ�
    SCIPreadProb(scip, filename1, NULL);
    
    SCIPwriteOrigProblem(scip, "D://primalmodel.lp", NULL, NULL);

    ////wqy0517����b
    //SCIPsolve(scip);
    //double obj = SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
    //cout << "obj:" << obj << endl;
    ////wqy0517����e

    std::unordered_map<int, SCIP_VAR*> Var_map;
    std::unordered_map<int, SCIP_CONS*> Cons_map;
    std::unordered_map<std::string, int> Var_Ni;
    //std::unordered_map<std::string, int> Cons_Ni;
    int nrows = SCIPgetNConss(scip);
    int ncols = SCIPgetNVars(scip);
    cout << ncols << "," << nrows << endl;
    if (ncols>10000000 || nrows>10000000)
    {
	    return std::make_tuple(std::map<int, General_Block>(), scip, vector<SCIP_VAR*>(), vector<SCIP_CONS*>(), vector<SCIP_VAR*>(), vector<SCIP_VAR*>());
    }
    std::unordered_map<int, std::vector<uint64_t>> matrix;
    std::unordered_map<int, std::vector<uint64_t>> matrix_dual;
    SCIP_CONS** conss = SCIPgetConss(scip);
    SCIP_VAR** vars_ori = SCIPgetVars(scip);
    std::unordered_set<string> name_of_vars;
    for (int i = 0; i < nrows; ++i)
    {
        // �������Լ��������
        SCIP_CONS* consdata = conss[i];
        //cout << "Լ�� " << conss[i] << "��" << endl;
        Cons_map[i] = consdata;
        // ���Լ����ϵ������Ӧ�ı���
        int nnonz = SCIPgetNVarsLinear(scip, consdata);
        SCIP_VAR** vars = SCIPgetVarsLinear(scip, consdata);
        SCIP_Real* coefs = SCIPgetValsLinear(scip, consdata);
        std::vector<uint64_t> temp_v;
        std::set<int> temp_set;
        // ����һ��Լ���еķ���Ԫ��
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
    cout << matrix_dual.size() << "," << matrix.size() << endl;
    cout << "matrix_varsize: " << name_of_vars.size() << endl;
    

    auto [parts, couping1] = RunHYPE(matrix, k);
    auto [Blocks,num_of_part] = PrimalBlock(parts, matrix, couping1);
    cout << "���ֵĿ���:  " << num_of_part << endl;
    for (auto var_map : Blocks)
    {
        cout << "�� " << var_map.first << "��: Լ������ :" << var_map.second.Cons.size() << " ��������: " << var_map.second.Vars.size() << endl;
    }
    auto GeneralBlocks = TransformerBlock(scip,Blocks,Var_map,Cons_map,num_of_part,Var_Ni);

	vector<SCIP_CONS*> fj_cons(conss, conss + nrows);
	vector<SCIP_VAR*> fj_vars(vars_ori, vars_ori + ncols);
    //����
    int num_col = 0;
    int num_row = 0;
    for (int i = 0; i < GeneralBlocks.size()-1; ++i)
    {
	    num_col += GeneralBlocks[i].bVars.size();
    	num_row += GeneralBlocks[i].bCons.size();
    }
    num_row += GeneralBlocks[GeneralBlocks.size()-1].bCons.size();
    std::cout << "��������: "<<ncols<<"-"<<num_col<<endl;
    std::cout << "Լ������: " << nrows << "-" << num_row << endl;
    for (int v = 0; v < GeneralBlocks.size(); ++v)
    {
        cout << "�� " <<v << "��: Լ������ :" << GeneralBlocks[v].bCons.size() << " ��������: " << GeneralBlocks[v].bVars.size() << endl;
    }


    //auto [parts2, couping2] = RunHYPE(matrix_dual, k);
    //auto [vars_x, vars_z] = DualBlock(parts2,matrix_dual,couping2,Var_map);
    vector<SCIP_VAR*> vars_x;
    vector<SCIP_VAR*> vars_z;


    return make_tuple(GeneralBlocks,scip,fj_vars,fj_cons,vars_x,vars_z);
}
