#include"DW.h"
#include "spdlog/spdlog.h"
#include<stack>

using namespace std;

DW::DW()
{
}

DW::~DW()
{
}

void showMemoryInfo(void)
{
    HANDLE handle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
    cout << "内存使用:" << pmc.WorkingSetSize / 1024 << "KB"
        << "(" << pmc.WorkingSetSize / (1024 * 1024) << "M)" << endl;
}


void DW::UpdateModel(Node& node, std::map<int, General_Block>block)
{
    node.MasterProblem.BuildSCIPModel(node.obj, node.couplingCons, node.couplingCons_num, node.subprob_num, node.Fixed_Info);
    //不需要重新初始化子问题
}

void DW::LoadModel(Node& node, std::map<int, General_Block> block)
{
    node.MasterProblem.BuildSCIPModel(node.obj, node.couplingCons, node.couplingCons_num, node.subprob_num, node.Fixed_Info);
    //初始化子问题
    node.SubProblem = new SubProb[node.subprob_num];
    for (size_t i = 0; i < node.subprob_num; i++)
    {
        node.SubProblem[i] = SubProb::SubProb(block[i].bCons, block[i].bobj, block[node.subprob_num].bCons, block[i].bVars);
    }

}

//列生成,输出当前问题的松弛最优解
int DW::ColumnGeneration(Node& node)
{
    bool breakCG = false;
    int iterrr = 0;
    int repeat_times = 0;
    cout << "======= node " << node.index << " Column generation begin! =======" << endl;
    while (!breakCG)
    {
        //cout << "列生成第" << ++iterrr << "代，开始" << endl;
        iterrr++;
        double objjj = 0;
        node.MasterProblem.Solve_Obj(objjj);
        cout << "第" << iterrr - 1 << "代主问题最优解为： " << objjj << endl;
        Solution duals; //对偶变量

        //求解限制主问题
        if (node.MasterProblem.Solve_Dual(duals) == 0)
        {
            return -1;
        }

        //cout << "对偶变量：" << endl;
        //for (size_t i = 0; i < duals.solSize; i++)
        //{
        //    cout << duals.varValue[i] << endl;//wqy
        //}


        //设置子问题目标函数系数
        for (int subprob_no = 0; subprob_no < node.subprob_num; subprob_no++)
        {
            bool Is_fixed = false;
            for (size_t idx = 0; idx < node.Fixed_Info.size(); idx++)
            {
                //cout << "subprob_no:" << subprob_no << endl;
                //cout << "node.Fixed_Info[idx].sp_no:" << node.Fixed_Info[idx].sp_no << endl;
                if (subprob_no == node.Fixed_Info[idx].sp_no)
                {
                    Is_fixed = true;
                }
            }
            if (Is_fixed)//如果这个子问题的列被固定了，那么将不需要更新和求解这个子问题
            {
                continue;
            }
            //分离出对应子问题的对偶变量
            Solution duals_part;
            int duals_part_num = node.couplingCons_num;//coupling constraints对应
            duals_part.solSize = duals_part_num;
            for (auto iter = duals.varValue.begin(); iter!= duals.varValue.end();iter++)
            {
                duals_part.varValue.insert(pair<int,double>(iter->first,iter->second));
                if (duals_part.varValue.size() == duals_part_num)
                {
                    break;
                }
            }

            node.SubProblem[subprob_no].BuildSCIPModel();
            node.SubProblem[subprob_no].UpdateTempObj(duals_part);
        }


        //求解子问题
        breakCG = true;
        
        for (int subprob_no = 0; subprob_no < node.subprob_num; subprob_no++)
        {
            bool Is_fixed = false;
            for (size_t idx = 0; idx < node.Fixed_Info.size(); idx++)
            {
                //cout << "subprob_no:" << subprob_no << endl;
                //cout << "node.Fixed_Info[idx].sp_no:" << node.Fixed_Info[idx].sp_no << endl;
                if (subprob_no == node.Fixed_Info[idx].sp_no)
                {
                    Is_fixed = true;
                }
            }
            if (Is_fixed)//如果这个子问题的列被固定了，那么将不需要更新和求解这个子问题
            {
                continue;
            }

            //Solution sol;
            //double objval;
            //node.SubProblem[subprob_no].Solve(sol, objval);
            node.SubProblem[subprob_no].GetSolution();

            //cout
            //for (auto iter = node.SubProblem[subprob_no].tempSol.varValue.begin(); iter != node.SubProblem[subprob_no].tempSol.varValue.end(); iter++)
            //{
            //    cout <<"x" << iter->first << ":" << iter->second << endl;
            //}
            //cout << endl;

            //判断reduced_cost
            double reduced_cost = node.SubProblem[subprob_no].tempSol.objValue - duals.varValue[node.couplingCons_num + subprob_no];
            //cout << "reduced_cost:" << reduced_cost << endl;
            //cout << "objValue:" << node.SubProblem[subprob_no].tempSol.objValue << endl;
            //cout << "duals:" << duals.varValue[node.couplingCons_num + subprob_no] << endl;
            if (reduced_cost < -0.0001)
            {
                breakCG = false;//有一个子问题可以加列的话，则不终止
                //cout
                //cout << node.MasterProblem.columnpool[subprob_no].columns.size() << endl;
                int col_num = node.MasterProblem.columnpool[subprob_no].columns.size();
                //加列

                node.MasterProblem.AddColumn(node.SubProblem[subprob_no].tempSol, subprob_no);
                if (node.MasterProblem.columnpool[subprob_no].columns.size() == col_num)
                {
                    spdlog::warn("SubProblem {} Add repetitive column. reduced_cost = {}", subprob_no, reduced_cost);
                    repeat_times++;
                    //cout << "添加重复列！" << endl;
                    if (repeat_times == 50)
                    {
                        spdlog::error("Add repetitive column more than 50 times!");
                        return -2;
                    }
                }

                //cout
                //cout << node.MasterProblem.columnpool[subprob_no].columns.size()<<endl;
                //cout << "子问题" << subprob_no << "加列" << endl;

                ////wqy测试0517 begin
                //cout << "子问题" << subprob_no << "加列：" << endl;

                //for (auto iterr = node.SubProblem[subprob_no].tempSol.varValue.begin(); iterr != node.SubProblem[subprob_no].tempSol.varValue.end(); iterr++)
                //{
                //    cout << iterr->first << " : " << iterr->second << endl;
                //}

                ////wqy测试0517 end

            }
            ////wqy测试0517 begin
            //else
            //{
            //    cout << "子问题" << subprob_no << "未加列：" << endl;
            //}
            ////wqy测试0517 end
        }

        //判断节点中的列池（lagrange加入的列）是否有能加入主问题中的列
        //待加入wqy

        if (!breakCG)
        {
            //根据列池更新（重构）主问题
            node.MasterProblem.BuildSCIPModel(node.obj, node.couplingCons, node.couplingCons_num, node.subprob_num, node.Fixed_Info);
        }
        else
        {
            //输出列生成松弛解
            cout <<"===============Column generation end! Relax solution is：" << objjj<<"================ " << endl;
            
        }

    }

    //求解主问题
    RMPSOL sol;
    double objval;
    node.MasterProblem.Solve(sol, objval);
    node.LPrelax_obj = objval;
    node.LPrelax_sol = sol;


#pragma region COUT列池中的列
    //for (size_t i = 0; i < node.subprob_num; i++)
    //{
    //	cout << "子问题 " << i << " :" << endl;
    //	for (auto iter = node.columnpool[i].columns.begin(); iter != node.columnpool[i].columns.end(); iter++)
    //	{
    //		for (auto var_iter = iter->varValue.begin(); var_iter != iter->varValue.end(); var_iter++)
    //		{
    //			cout << var_iter->first << " : " << var_iter->second << endl;
    //		}
    //		cout << "%%%%%%%%%" << endl;
    //	}
    //}

#pragma endregion

    cout << "MasterProblem.free()前" << endl;
    showMemoryInfo();
    node.MasterProblem.free();
    cout << "MasterProblem.free()后" << endl;
    showMemoryInfo();


}


void DW::Diving(ColumnPool* columnpool, std::map<int, General_Block> block, Objective master_obj)
{
    cout << "===================Diving start!===================" << endl;
#pragma region 检验列是否可行
    for (size_t i = 0; i < block.size() - 1; i++)
    {
        int no = 0;
        for (auto iter = columnpool[i].columns.begin(); iter != columnpool[i].columns.end(); iter++)
        {
            no++;
            if (!Check_Part_Cons(block, iter->varValue, i))
            {
                cout << "Subproblem" << i << "column" << no << endl;
                cout << "Error! Initial columns is infeasible!" << endl;
                break;
            }
            else
            {
                //cout << i << endl;
            }
        }
        cout << "Subproblem " << i << " has " << no << " columns" << endl;
    }
#pragma endregion

    //加入一个列 使得主问题可行
    //for (size_t i = 0; i < block.size() - 1; i++)
    //{
    //    Column column_rhs;
    //    columnpool[i].columns.insert(column_rhs);
    //}


#pragma region 赋值node初始化所需变量
    int CouplingCons_num = block[block.size() - 1].bCons.size();
    Constraints* CouplingCons = new Constraints[CouplingCons_num];
    for (size_t i = 0; i < CouplingCons_num; i++)
    {
        CouplingCons[i] = block[block.size() - 1].bCons[i];
    }
    int subprob_num = block.size() - 1;
    int IndependentCons_num = 0;
    for (size_t i = 0; i < block.size() - 1; i++)
    {
        IndependentCons_num += block[i].bCons.size();
    }
    Constraints* IndependentCons = new Constraints[IndependentCons_num];
    int cons_idx = 0;
    for (size_t i = 0; i < block.size() - 1; i++)
    {
        for (size_t j = 0; j < block[i].bCons.size(); j++)
        {
            IndependentCons[cons_idx] = block[i].bCons[j];
            cons_idx++;
        }
    }
#pragma endregion

    Node node(master_obj, CouplingCons, CouplingCons_num, IndependentCons, IndependentCons_num, columnpool, subprob_num);
    
    //PROCESS_MEMORY_COUNTERS_EX pmc;
    //GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc), sizeof(pmc));
    //std::cout << "当前进程占用的内存为: " << pmc.WorkingSetSize / 1024 << "KB" << std::endl;
    showMemoryInfo();
    bool break_flag = false;//diving循环break标志，1为截止
    bool root_flag = true;//是否为根节点
    std::stack<Node> node_list;
    int maxDiscrepancy = 10;
    int node_no = 0;//节点的序号 
    node.index = node_no;
    node_no++;
    int maxDepth;
    vector<double>final_feasible_obj; //用于存放最后所有的可行解的目标函数值
    vector<vector<double>>final_feasible_solution;//用于存放最后所有的可行解的解值

//=========================Diving循环开始==================================
    while (!break_flag)
    {
#pragma region 加载限制主问题和子问题模型

        if (root_flag)
        {
            //初始化模型
            LoadModel(node, block);
            root_flag = false;
        }
        else
        {
            //根据固定的列修改模型
            UpdateModel(node, block);
        }
#pragma endregion

        

        //列生成
        int CG_Result = ColumnGeneration(node);
        if (CG_Result == -1)
        {
            //cout << "列生成异常结束！限制主问题无可行解" << endl;
            spdlog::error("列生成异常结束！限制主问题无可行解");
            break;
        }
        if (CG_Result == -2)
        {
            //cout << "列生成异常结束！限制主问题无可行解" << endl;
            spdlog::error("列生成中止！原因：加入重复列");
        }
        else
        {
            showMemoryInfo();
            //计算一共有多少变量
            int master_prob_col_num = 0;
            for (size_t i = 0; i < node.subprob_num; i++)
            {
                bool flag_conti = false;
                for (size_t j = 0; j < node.Fixed_Info.size(); j++)
                {
                    if (i == node.Fixed_Info[j].sp_no)
                    {
                        flag_conti = true;
                    }
                }
                if (flag_conti)
                {
                    continue;
                }
                master_prob_col_num += node.MasterProblem.columnpool[i].columns.size();

            }

            bool Is_add_node = false;
            showMemoryInfo();
            cout << "nodelist.size() = " << node_list.size() << endl;
            //for (size_t i = maxDiscrepancy; i > 0; i--)//考虑压进栈的顺序
            for (size_t i = master_prob_col_num; i > 0; i--)
            {

                Node node_temp(node);
#pragma region cout
                //cout << node_temp.columnpool->columns.size() << endl;
               //cout << node.columnpool->columns.size() << endl;
               //cout << node_temp.columnpool[0].columns.size() << endl;
               //cout << node.columnpool[0].columns.size() << endl;
               //cout << node_temp.columnpool[1].columns.size() << endl;
               //cout << node.columnpool[1].columns.size() << endl;
               //

               //cout << "1-node.columnpool[0].var_num = " << node.columnpool[0].var_num << endl;
               //cout << "1-node_temp.columnpool[0].var_num = " << node_temp.columnpool[0].var_num << endl;

               //node_temp.columnpool[0].var_num = 666;
               //cout << "2-node.columnpool[0].var_num = " << node.columnpool[0].var_num << endl;
               //cout << "2-node_temp.columnpool[0].var_num = " << node_temp.columnpool[0].var_num << endl;
#pragma endregion

            //确定固定列, + 修改主问题列池，以及节点列池（可不改）
            //返回值为1代表固定的列在主问题解中取值不为0，才可以固定
                if (Fixed_Column(node_temp, i - 1))
                {
                    if (Test_feasible(node_temp))//可行才加wqy0412，减小内存占用
                    {
                        node_temp.index = node_no;
                        node_no++;
                        Is_add_node = true;
                        //把节点加入节点池
                        node_list.push(node_temp);
                        cout << "Generate node " << node_temp.index << ";  ";
                    }

                }
#pragma region cout
                //cout << node_temp.columnpool->columns.size() << endl;
                //cout << node.columnpool->columns.size() << endl;
                //cout << node_temp.MasterProblem.columnpool[0].columns.size() << endl;
                //cout << node.MasterProblem.columnpool[0].columns.size() << endl;
                //cout << node_temp.MasterProblem.columnpool[1].columns.size() << endl;
                //cout << node.MasterProblem.columnpool[1].columns.size() << endl;
#pragma endregion

            }
            showMemoryInfo();
            cout << "nodelist.size() = " << node_list.size() << endl;
            if (!Is_add_node)
            {
                cout << "此节点未生成新节点！";
            }
            cout << endl;
        }

        

        bool Infea_flag = false;//Diving能否找到可行解，1为找不到
        while (true)
        {
            //生成新的节点
            if (node_list.size() != 0)
            {
                node = node_list.top();
                node_list.pop();
            }


            if (node_list.size() == 0)
            {
                cout << "===================Diving end!===================" << endl;
                break_flag = true;
                if (final_feasible_obj.size() != 0)
                {
                    cout << "Diving所有可行解为：" << endl;
                }               
                for (int k = 0; k < final_feasible_obj.size(); k++)
                {
                    cout << k + 1 << " : " << final_feasible_obj[k] << endl;
                }
                break;
            }
            cout << endl;
            showMemoryInfo();
            cout << "Select node " << node.index<<", fixed "<<node.Fixed_Info.size()<<"columns" << endl;
            if (node.Fixed_Info.size() != node.subprob_num)
            {
                break;
            }
            else
            {
                //cout << "Diving结果为:" << endl;
                //把子问题的列混合
                vector<double>solution_mass;
                vector<int>solution_idx_mass;
                for (size_t i = 0; i < node.subprob_num; i++)
                {
                    for (size_t k = 0; k < node.Fixed_Info[i].col_val.size(); k++)
                    {
                        solution_mass.push_back(node.Fixed_Info[i].col_val[k]);
                        solution_idx_mass.push_back(node.Fixed_Info[i].col_var_no[k]);
                    }
                }
                //按照变量序号排序

                vector<int>Index;
                for (size_t i = 0; i < solution_mass.size(); i++)
                {
                    Index.push_back(i);
                }

                sort(Index.begin(), Index.end(),
                    [&](const int& a, const int& b) {
                    return (solution_idx_mass[a] < solution_idx_mass[b]); });

                vector<double>solution;
                for (size_t i = 0; i < solution_mass.size(); i++)
                {
                    solution.push_back(solution_mass[Index[i]]);
                }

                double obj = 0;
                for (auto iter = master_obj.coef.begin(); iter != master_obj.coef.end(); iter++)
                {
                    obj += iter->second * solution_mass[Index[iter->first]];
                }

                if (Check_All_Cons(block, solution))
                {
                    final_feasible_obj.push_back(obj);
                    final_feasible_solution.push_back(solution);
                    cout << "%%%%%%%%%%% Found a feasible for primal problem! %%%%%%%%% "<<obj<<" %%%%%%" << endl;
                }

            }
            cout << "node.free()前" << endl;
            showMemoryInfo();
            node.free();
            cout << "node.free()后" << endl;
            showMemoryInfo();
        }

        //if (Infea_flag)//找不可行解输出，后接FJ
        //{
        //    if (final_feasible_obj.size() == 0)
        //    {
        //        cout << "找不到可行解！" << endl;
        //    }
        //    
        //    //设置列禁忌表
        //    //Get_Tabu_list(node, tabu_list);
        //    break;
        //}

    }


}

void DW::Get_InitColumns(ColumnPool*& columnpool)
{
    


}


int DW::Fixed_Column(Node& node_temp, int no)
{
    
    FIX fix_info;
    //对解排序
    vector<double> sol_vec;

    for (size_t i = 0; i < node_temp.LPrelax_sol.varValue.size(); i++)
    {
        bool fix_flag = false;
        //固定过的列，在排序中赋值为0，不能继续固定
        for (size_t j = 0; j < node_temp.Fixed_Info.size(); j++)
        {
            if (node_temp.Fixed_Info[j].sp_no == node_temp.LPrelax_sol.subprob_no[i])
            {
                fix_flag = true;
            }
        }

        if (fix_flag)
        {
            sol_vec.push_back(-1);
        }
        else
        {
            sol_vec.push_back(node_temp.LPrelax_sol.varValue[i]);
        }

    }

    vector<int> index(sol_vec.size(), 0);
    for (int i = 0; i != index.size(); i++)
    {
        index[i] = i;
    }
    sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
        return (sol_vec[a] > sol_vec[b]); });


    //选第一个列固定，找到固定的列属于哪个子问题,以及对应列
    int fixed_column_no = index[no];
    if (sol_vec[fixed_column_no] == -1)
    {
        cout << "error!" << endl; 
    }
    else if(sol_vec[fixed_column_no] == 0)
    {
       // return 0;
    }
    int sp_no = 0;
    Column fixed_column;
    for (size_t i = 0; i < node_temp.subprob_num; i++)
    {
        if (int(fixed_column_no - node_temp.MasterProblem.columnpool[i].columns.size()) >= 0)
        {
            fixed_column_no -= node_temp.MasterProblem.columnpool[i].columns.size();
            sp_no++;
        }
        else
        {
            int idexx = 0;
            for (auto iter = node_temp.MasterProblem.columnpool[i].columns.begin(); iter != node_temp.MasterProblem.columnpool[i].columns.end(); ++iter)
            {
                if (idexx == fixed_column_no)
                {
                    fix_info.col_no = idexx;
                    fixed_column.solSize = iter->solSize;
                    fixed_column.varValue = iter->varValue;
                    break;
                }
                idexx++;
            }
            break;
        }
    }


    fix_info.sp_no = sp_no;
    for (auto iter = fixed_column.varValue.begin(); iter != fixed_column.varValue.end(); iter++)
    {
        fix_info.col_val.push_back(iter->second);
        fix_info.col_var_no.push_back(iter->first);
    }
    node_temp.Fixed_Info.push_back(fix_info);

    //按照固定的列，处理列池
    ColumnPool new_columnpool;
    Column fixed_col;
    new_columnpool.vars_id = node_temp.MasterProblem.columnpool[fix_info.sp_no].vars_id;
    new_columnpool.var_num = node_temp.MasterProblem.columnpool[fix_info.sp_no].var_num;
    int index_no = 0;
    for (auto iter = node_temp.MasterProblem.columnpool[fix_info.sp_no].columns.begin(); iter != node_temp.MasterProblem.columnpool[fix_info.sp_no].columns.end(); iter++)
    {
        if (index_no == fix_info.col_no)
        {
            fixed_col.isAdded = iter->isAdded;
            fixed_col.objValue = iter->objValue;
            fixed_col.solSize = iter->solSize;
            fixed_col.varValue = iter->varValue;
            break;
        }
        index_no++;
    }

    new_columnpool.columns.insert(fixed_col);
    node_temp.MasterProblem.columnpool[fix_info.sp_no] = new_columnpool;

    return 1;

}


bool DW::Test_feasible(Node node)
{
    bool feasible = 0;
    node.MasterProblem.Test_feasble2(feasible, node.obj, node.couplingCons, node.couplingCons_num, node.subprob_num, node.Fixed_Info);

    return feasible;
}

bool DW::Check_All_Cons(std::map<int, General_Block>block, vector<double>col_val)
{
    bool Check_flag = true;
    vector<bool>cons_check;
    for (auto iter = block.begin(); iter != block.end(); iter++)
    {
        for (int i = 0; i < iter->second.bCons.size(); i++)
        {
            double lhs = 0;
            for (auto iter_var = iter->second.bCons[i].idDic.begin(); iter_var != iter->second.bCons[i].idDic.end(); iter_var++)
            {
                lhs += iter_var->second * col_val[iter_var->first];
            }
            if (iter->second.bCons[i].Prop == 0)
            {
                if (lhs < iter->second.bCons[i].rhs)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "约束违反！" << endl;
                }
            }
            else if (iter->second.bCons[i].Prop == 1)
            {
                if (lhs > iter->second.bCons[i].rhs)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "约束违反！" << endl;
                }
            }
            else
            {
                if (lhs != iter->second.bCons[i].rhs)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "约束违反！" << endl;
                }
            }
        }
    }
    return Check_flag;
}


bool DW::Check_Part_Cons(std::map<int, General_Block>block, unordered_map<int,double>var_value, int block_no)
{
    bool Check_flag = true;
    vector<bool>cons_check;
    int idx_block_no = 0;
    for (auto iter = block.begin(); iter != block.end(); iter++)
    {
        if (idx_block_no != block_no)
        {
            continue;
        }
        idx_block_no++;

        for (int i = 0; i < iter->second.bCons.size(); i++)
        {
            double lhs = 0;

            for (auto iter_var = iter->second.bCons[i].idDic.begin(); iter_var != iter->second.bCons[i].idDic.end(); iter_var++)
            {
                lhs += iter_var->second * var_value[iter_var->first];

            }
            if (iter->second.bCons[i].Prop == 0)
            {
                if (lhs - iter->second.bCons[i].rhs < -EPS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "约束违反！" << endl;
                    cout << "左端项为：" << lhs << "  右端项为：" << iter->second.bCons[i].rhs << "  符号为：>=" << endl;
                }
            }
            else if (iter->second.bCons[i].Prop == 1)
            {
                if (lhs - iter->second.bCons[i].rhs > EPS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "约束违反！" << endl;
                    cout << "左端项为：" << lhs << "  右端项为：" << iter->second.bCons[i].rhs << "  符号为：<=" << endl;
                }
            }
            else
            {
                if (lhs - iter->second.bCons[i].rhs > EPS || lhs - iter->second.bCons[i].rhs < -EPS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "约束违反！" << endl;
                    cout << "左端项为：" << lhs << "  右端项为：" << iter->second.bCons[i].rhs << "  符号为：=" << endl;
                }
            }
        }
    }
    return Check_flag;
}

int DW::Get_Tabu_list(Node node, vector<vector<unordered_map<int, double>>>& tabu_list)
{

}