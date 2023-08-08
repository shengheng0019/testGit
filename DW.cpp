#include"DW.h"

constexpr double EPSFORREDUCEDCOST = 1e-4;
constexpr double EPSFORCONSTRAINTS = 1e-4;

//using namespace std;


/*void showMemoryInfo(void)
{
    HANDLE handle = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
    cout << "memory using:" << pmc.WorkingSetSize / 1024 << "KB"
        << "(" << pmc.WorkingSetSize / (1024 * 1024) << "M)" << endl;
}*/


void DW::UpdateModel(Node& node)
{
    node.BuildMasterProblemSCIPModel();
    //No need to re-init the subproblem model
}

void DW::LoadModel(Node& node, std::map<int, General_Block> block)
{
    node.BuildMasterProblemSCIPModel();
    //init subproblem
    node.InitSubProblem(block);//not build model here, just initial some value
}

//output: relax solution of current node
int DW::ColumnGeneration(Node& node, bool root_flag)
{
    //some local variable
    bool breakCG = false;
    int iterrr = 0;
    int repeatTimes = 0;

    //cout << "======= Column generation begin! =======" << endl;
    spdlog::info("========== Column generation begin! ==========");
    while (!breakCG){
        iterrr++;
        //@wqy0612test,if column generation over 50 iters, break
        if (iterrr > 50) {
            break;
        }

        double masterObj = 0;
        if (node.MasterProblem().SolveObj(masterObj) == 0){
            return -1;
        }
        //cout << "The optimal obj of gen " << iterrr - 1 << " master problem is�� " << masterObj << endl;
        spdlog::info("The optimal obj of gen {0} master problem is: {1}", iterrr - 1, masterObj);

        vector<double> duals; //dual variable
        //Solve RMP
        node.MasterProblem().SolveDual(duals);

        //set objective function coeff of subproblem
        for (int subprob_no = 0; subprob_no < node.SubprobNum(); subprob_no++){
            //split out the dual variable corresponding to the subproblem
            int duals_part_num = node.CouplingConsNum();//correspond the coupling constraints
            vector<double> duals_part(duals.begin(), duals.begin() + duals_part_num);

            node.BuildSubProblemModel(subprob_no);
            node.UpdateSubProbObjCoef(subprob_no, duals_part);
        }

        //solve the subproblem
        breakCG = true;

        for (int subprob_no = 0; subprob_no < node.SubprobNum(); subprob_no++)
        {
            double subProbObj(-1);
            Solution columnGen;
            node.SubProblem(subprob_no).GetSolutionNew(columnGen);

            //calculate reduced_cost
            double reducedCost = subProbObj - duals[node.CouplingConsNum() + subprob_no];

            if (reducedCost < -EPSFORREDUCEDCOST){
                breakCG = false;//as long as one subproblem can be added column, it is not terminated

                //add column
                if (!node.MasterProblem().AddColumn(columnGen, subprob_no)){
                    spdlog::warn("SubProblem {} Add repetitive column. reduced_cost = {}", subprob_no, reducedCost);
                    repeatTimes++;
                    //cout << "add repetitve columns" << endl;
                    if (repeatTimes == 10){
                        spdlog::error("Add repetitive column more than 10 times!");
                        return -2;
                    }
                }
            }
            ////@wqy test 0517 begin
            //else
            //{
            //    //cout << "subproblem " << subprob_no << "do not add column��" << endl;
            //    //cout << "do not add column reduced_cost:" << reduced_cost << endl;
            //}
            ////@wqy test 0517 end
        }

        //Todo
        //check whether any column generated by lagrange algorithm can be added into the columnpool of the master problem

        if (!breakCG){
            //update master problem based on the new columnpool
            node.BuildMasterProblemSCIPModel();
        }else{
            //output the relax obj
            cout << "===============Column generation end! Relax solution is: " << masterObj << "================ " << endl;
        }

    }

    //solve final master problem
    node.BuildMasterProblemSCIPModel();
    RMPSOL sol;
    double objval;
    node.MasterProblem().Solve(sol, objval, root_flag);
    node.SetLpRelaxObj(objval);
    node.SetLpRelaxSol(sol);

#pragma region COUT the columns in columnpool
    //for (int i = 0; i < node.subprob_num; i++)
    //{
    //	cout << "subproblem " << i << " :" << endl;
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
    return 1;
}

//output: relax solution of current node
//@wqy0622 This function mix the building, updating objcoef and solving of the subproblem.
int DW::ColumnGenerationNew(Node& node, bool root_flag)
{
    //some local variable
    bool breakCG = false;
    int iterrr = 0;
    int repeatTimes = 0;

    //cout << "======= Column generation begin! =======" << endl;
    spdlog::info("========== Column generation begin! ==========");
    while (!breakCG) {
        iterrr++;
        //@wqy0612test,if column generation over 50 iters, break
        if (iterrr > 10) {
            break;
        }

        double masterObj = 0;
        if (node.MasterProblem().SolveObj(masterObj) == 0) {
            return -1;
        }
        //cout << "The optimal obj of gen " << iterrr - 1 << " master problem is " << masterObj << endl;
        spdlog::info("The optimal obj of gen {0} master problem is: {1}", iterrr - 1, masterObj);

        vector<double> duals; //dual variable
        //Solve RMP
        node.MasterProblem().SolveDualCplex(duals);
        //node.MasterProblem().SolveDual(duals);

        //build the subproblem and solve
        breakCG = true;
        for (int subprob_no = 0; subprob_no < node.SubprobNum(); subprob_no++) {
            //build the subproblem
            //split out the dual variable corresponding to the subproblem
            int duals_part_num = node.CouplingConsNum();//correspond the coupling constraints
            vector<double> duals_part(duals.begin(), duals.begin() + duals_part_num);

            node.BuildSubProblemModelNew(subprob_no, duals_part);

            //solve the problem
            Solution columnGen;
            node.SubProblem(subprob_no).GetSolutionNew(columnGen);
            
            //calculate reduced_cost
            double reducedCost = columnGen.objValue - duals[node.CouplingConsNum() + subprob_no];

            if (reducedCost < -EPSFORREDUCEDCOST) {
                breakCG = false;//as long as one subproblem can be added column, it is not terminated
                //add column
                int colNumBefore = node.GetMasterColNum()[subprob_no];
                node.AddMasterProblemColumn(columnGen, subprob_no);
                int colNumAfter = node.GetMasterColNum()[subprob_no];               
                if (colNumAfter == colNumBefore) {
                    spdlog::warn("SubProblem {} Add repetitive column. reduced_cost = {}", subprob_no, reducedCost);
                    repeatTimes++;
                    //cout << "add repetitve columns: " << endl;
                    if (repeatTimes == 10) {
                        spdlog::error("Add repetitive column more than 10 times!");
                        return -2;
                    }
                }
            }
        }

        //Todo
        //check whether any column generated by lagrange algorithm can be added into the columnpool of the master problem

        if (!breakCG) {
            //update master problem based on the new columnpool
            node.BuildMasterProblemSCIPModel();
        }
        else {
            //output the relax obj
            cout << "===============Column generation end! Relax solution is: " << masterObj << "================ " << endl;
        }

    }

    //solve final master problem
    node.BuildMasterProblemSCIPModel();
    RMPSOL sol;
    double objval;
    node.MasterProblem().Solve(sol, objval, root_flag);
    node.SetLpRelaxObj(objval);
    node.SetLpRelaxSol(sol);

#pragma region COUT the columns in columnpool
    //for (int i = 0; i < node.subprob_num; i++)
    //{
    //	cout << "subproblem " << i << " :" << endl;
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
    return 1;
}


void DW::Diving_new(const vector<ColumnPool>& columnpool, std::map<int, General_Block> block, Objective master_obj, SCIP* scipmodel, std::map<std::string,writeJson::Soljson>& all_solsmap, std::string part_name, std::chrono::system_clock::time_point start)
{
    spdlog::info("===================Diving start!===================" );
#pragma region check whether each column in columnpool is feasible
    for (int i = 0; i < block.size() - 1; i++){
        int no = 0;
        for (int j = 0; j < columnpool[i].Columns().size(); j++) {
            no++;
            if (!Check_Part_Cons(block, columnpool[i].Columns()[j].varValue, i)) {
                cout << "Subproblem" << i << "column" << no << endl;
                cout << "Error! Initial columns is infeasible!" << endl;
                //break;
            }else{
                //cout << i << endl;
            }
        }
        spdlog::info("Subproblem {} has {} columns", i, no);
    }
#pragma endregion
    // Init the variable of node
    int nodeNo = 0;//index of node 
    int CouplingCons_num = block[block.size() - 1].bCons.size();
    vector<Constraints> CouplingCons;
    for (int i = 0; i < CouplingCons_num; i++){
        CouplingCons.push_back(block[block.size() - 1].bCons[i]);
    }
    int subprob_num = block.size() - 1;
    int IndependentCons_num = 0;
    for (int i = 0; i < subprob_num; i++){
        IndependentCons_num += block[i].bCons.size();
    }
    vector<Constraints> IndependentCons;
    for (int i = 0; i < block.size() - 1; i++){
        for (const auto& bCon : block[i].bCons){
            IndependentCons.push_back(bCon);
        }
    }
    //preprocess 0613 
    std::unordered_map<int, bool>varIntTemp; //record variable is integer or not, int->var.id; bool->yes/no
    for (int i = 0; i < block.size(); i++){
        for (int j = 0; j < block[i].bVars.size(); j++){
            if (block[i].bVars[j].Type == VarType::Num){
                varIntTemp.insert(make_pair(block[i].bVars[j].Id, false));
            }else{
                varIntTemp.insert(make_pair(block[i].bVars[j].Id, true));
            }
        }
    }

    Node node(master_obj, CouplingCons, CouplingCons_num, IndependentCons, IndependentCons_num, columnpool, subprob_num, varIntTemp, nodeNo);

    //some local variable
    bool divingBreakFlag = false;//flag for judge diving process break
    bool rootFlag = true;//flag for root node
    std::stack<FIX> fixedInformationList;//the set of fixed information, each element in the stack is the fixed info for a node
    int maxDiscrepancy = 20;
    std::vector<double>finalFeasibleObj; //store the obj value of the final feasible solution
    std::vector<std::vector<double>>finalFeasibleSolution;//store the variable value of the final feasible solution
    spdlog::info("2");
    //=========================Diving Algrithm Begin==================================
    while (!divingBreakFlag){

        //load the model of masterproblem and the subproblem 
        if (rootFlag){
            //init model
            LoadModel(node, block);
            rootFlag = false;
        }else{
            //update model corresponding to the columnpool
            UpdateModel(node);
        }

        //column generation
        int CG_Result = ColumnGenerationNew(node, rootFlag);
        if (CG_Result == -1){
            spdlog::error("Column generation exception! There is no feasible solution to the RMP");
        }
        if (CG_Result == -2){
            spdlog::error("Column generation process aborted! Cause: the same column is added");
        }

        CalculateFixedInfo(fixedInformationList, node, maxDiscrepancy);

        //there is no node in nodelist, then break.
        if (fixedInformationList.empty()) { break; }
            
        //if all subproblems is fixed, the node is finish
        while (!fixedInformationList.empty()){
            //get the top of the fixedinformation stack
            FIX fixedInfosTop = fixedInformationList.top();
            fixedInformationList.pop();

            //update(generate) the node based on this node itself and the fixedInformation(change the information of node)
            //0613 @wqy modify: just fix current fixed information instand of all fixed information
            //0620 @wqy: only change the fixInfo and the columnpool of the node in this function
            UpdateNodeNew(node, fixedInfosTop);

            //judge whether all subproblem is fixed in this node
            if (node.FixedInfo().spNoList.size() == node.SubprobNum()) {
                std::vector<double>solution;
                double objFin = 0;
                if (FindFeasSol(node, scipmodel, block, objFin)){
                    finalFeasibleObj.push_back(objFin);
                    spdlog::info("%%%%%%%%%%% Found a feasible for primal problem! %%%%%%%%% {} %%%%%%", objFin);
                    divingBreakFlag = true;//only find one feasible solution!

                    //save DW result
                    std::ofstream o((part_name + "DW.json").c_str(), std::ofstream::trunc);
                    writeJson::Soljson dwSol;
                    dwSol.obj = objFin;
                    dwSol.source = "D";
                    dwSol.duration = double (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count()) / 1000;
                    string solKey = "D1";
                    nlohmann::json j;
                    dwSol.to_json(j);
                    o << j << std::endl;

                    all_solsmap.insert(make_pair(solKey, dwSol));
                    //break;
                }

            }
            else {
                break; //if any subproblem is not fixed, break
            }
        }

    }
    cout << "===================Diving end!===================" << endl;

    //output the result of Diving
    if (finalFeasibleObj.size() != 0){
        cout << "All feasible found by Diving is: " << endl;
    }
    else{
        cout << "Diving algorithm does not find any feasible solution!" << endl;
    }
    for (int k = 0; k < finalFeasibleObj.size(); k++){
        cout << k + 1 << " : " << finalFeasibleObj[k] << endl;
    }

}


bool DW::Check_All_Cons(std::map<int, General_Block>block, std::vector<double>col_val)
{
    bool Check_flag = true;
    std::vector<bool>cons_check;
    for (auto iter = block.begin(); iter != block.end(); iter++)
    {
        for (int i = 0; i < iter->second.bCons.size(); i++)
        {
            double lhs = 0;
            for (auto iter_var = iter->second.bCons[i].idDic.begin(); iter_var != iter->second.bCons[i].idDic.end(); iter_var++)
            {
                lhs += iter_var->second * col_val[iter_var->first];
            }

            //@wqy0523edit using EPS instead of <=
            if (iter->second.bCons[i].Prop == PROP::geq)
            {
                if (lhs - iter->second.bCons[i].rhs < -EPSFORCONSTRAINTS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                }
            }
            else if (iter->second.bCons[i].Prop == PROP::leq)
            {
                if (lhs - iter->second.bCons[i].rhs > EPSFORCONSTRAINTS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                }
            }
            else
            {
                if (abs(lhs - iter->second.bCons[i].rhs) < -EPSFORCONSTRAINTS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                }
            }
        }
    }
    return Check_flag;
}


bool DW::Check_All_Cons(std::map<int, General_Block>block, std::vector<std::pair<string, double>>colNameString)
{
    bool Check_flag = true;
    std::vector<bool>cons_check;
    for (auto iter = block.begin(); iter != block.end(); iter++)
    {
        for (int i = 0; i < iter->second.bCons.size(); i++)
        {
            double lhs = 0;
            for (auto iter_var = iter->second.bCons[i].exprDic.begin(); iter_var != iter->second.bCons[i].exprDic.end(); iter_var++)
            {
                //It is an inefficient way, just for test
                double valNow = 0;
                bool foundFlag = false;
                for (const auto& colPair : colNameString)
                {
                    if (colPair.first == "t_" + iter_var->first.Name) {
                        valNow = colPair.second;
                        foundFlag = true;
                        break;
                    }
                }
                if (!foundFlag) {
                    cout << "Error! Function Check All Cons not found the vars!" << endl;
                }

                lhs += iter_var->second * valNow;
            }

            //@wqy0523edit using EPS instead of <=
            if (iter->second.bCons[i].Prop == PROP::geq)
            {
                if (lhs - iter->second.bCons[i].rhs < -EPSFORCONSTRAINTS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                }
            }
            else if (iter->second.bCons[i].Prop == PROP::leq)
            {
                if (lhs - iter->second.bCons[i].rhs > EPSFORCONSTRAINTS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                }
            }
            else
            {
                if (abs(lhs - iter->second.bCons[i].rhs) < -EPSFORCONSTRAINTS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                }
            }
        }
    }
    return Check_flag;
}

bool DW::Check_Part_Cons(std::map<int, General_Block>block, unordered_map<int, double>var_value, int block_no)
{
    bool Check_flag = true;
    std::vector<bool>cons_check;
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
            if (iter->second.bCons[i].Prop == PROP::geq)
            {
                if (lhs - iter->second.bCons[i].rhs < -EPS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                    cout << "lhs is: " << lhs << "  rhs is: " << iter->second.bCons[i].rhs << "  sign is: >=" << endl;
                }
            }
            else if (iter->second.bCons[i].Prop == PROP::leq)
            {
                if (lhs - iter->second.bCons[i].rhs > EPS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                    cout << "lhs is: " << lhs << "  rhs is: " << iter->second.bCons[i].rhs << "  sign is: <=" << endl;
                }
            }
            else
            {
                if (lhs - iter->second.bCons[i].rhs > EPS || lhs - iter->second.bCons[i].rhs < -EPS)
                {
                    Check_flag = false;
                    cout << iter->second.bCons[i].Name << "constraint violated: " << endl;
                    cout << "lhs is: " << lhs << "  rhs is: " << iter->second.bCons[i].rhs << "  sign is: <=" << endl;
                }
            }
        }
    }
    return Check_flag;
}

//the order of the columns maybe have some problems 0620@wqy
//@wqy save 0625
//int FindFixedInformation(const Node& node, int columnIndex, FIX& fixedInfoCurrent){
//    int varGlobalIndex = columnIndex;
//    fixedInfoCurrent = node.FixedInfo(); //shallow copy
//
//    fixedInfoCurrent.recordColumnpool.clear();
//    for (int k = 0; k < node.SubprobNum(); k++) {
//        fixedInfoCurrent.recordColumnpool.push_back(node.MasterProblem().columnpool[k]);
//    }
//
//    //find the column correspending to the columnIndex
//    int sp_no = 0;
//    for (int i = 0; i < node.SubprobNum(); i++){
//        if (int(columnIndex - node.GetMasterColNum()[i]) >= 0){
//            columnIndex -= node.GetMasterColNum()[i];
//            sp_no++;
//        }else{
//            int idexx = 0;
//            //the problem is here, the order of "unordered_set"columnpool maybe not right 0620 way
//            unordered_map<int, double>testMap;
//            unordered_map<int, double>testMapCopy;
//            for (int i = 0; i < 884; i++) {
//                testMap.insert(std::make_pair(i, i + 2));
//            }
//            testMapCopy = testMap;
//            std::unordered_set<Column> columnsTemp;
//            columnsTemp = node.GetMasterProblemColumns(i);
//            //for (const auto& colIter : node.MasterProblem().columnpool[i].columns){
//            for (const auto& colIter : columnsTemp) {
//                if (idexx == columnIndex){
//                    fixedInfoCurrent.spNoList.push_back(sp_no);
//                    if(fixedInfoCurrent.spIsFixed[sp_no]){
//                        std::cout << "Error! fix column has problem!" << std::endl; //just check, delete later
//                    }
//                    fixedInfoCurrent.spIsFixed[sp_no] = true;
//                    ColumnPool columnpoolTemp;
//                    Column columnTemp;
//                    columnTemp.varValue = colIter.varValue;
//
//                    columnpoolTemp.vars_id = fixedInfoCurrent.recordColumnpool[sp_no].vars_id;
//                    columnpoolTemp.var_num = fixedInfoCurrent.recordColumnpool[sp_no].var_num;
//                    columnpoolTemp.columns.insert(columnTemp);
//                    fixedInfoCurrent.recordColumnpool[sp_no] = columnpoolTemp;  //fixInfo.recordColumnpool have only one column in newest fixed subproblem
//                    //int ittt2(0);
//                    std::vector<std::pair<int, double>> fixColTemp;
//                    for (const auto& varPair : colIter.varValue){
//                        fixColTemp.push_back(varPair);
//                    }
//                    fixedInfoCurrent.fixedColumnValueList.push_back(fixColTemp);
//                    break;
//                }
//                idexx++;
//            }
//            break;
//        }
//    }
//
//    //Check feasibility, 0613 here @wqy edit, only fix integer variable
//    //@wqy need modify 0620, it fixs all of the variables now.
//    return node.CheckFixFeasible(varGlobalIndex) ? 1 : 0;
//}

//@wqy0625 solve the problem of order
int FindFixedInformation(const Node& node, string columnName, FIX& fixedInfoCurrent) {

    if (!node.CheckFixFeasible(columnName)) { return 0; }

    //@wqy 0626 it can use the name of masterproblem variable to find the column, because the order of columnpool is unchanged.
    //find columnValue corresponding to the columnName
    std::unordered_map<int, double> varValueFixed;
    //split the columnName, i = 1 because the first character of columnName is '_'.
    std::string spNoSpilt;
    std::string colNoSpilt;
    int count = 0;
    for (int i = 1; i < columnName.length(); i++) {
        if (columnName[i] == '_') {
            count++;
        }
        else {
            if (count == 1) {
                spNoSpilt.push_back(columnName[i]);
            }
            if (count == 2) {
                colNoSpilt.push_back(columnName[i]);
            }
        }      
    }
    varValueFixed = node.GetMasterProblemColumns(stoi(spNoSpilt))[stoi(colNoSpilt)].varValue;


    //assign the fixedInfoCurrent
    fixedInfoCurrent = node.FixedInfo(); //shallow copy
    fixedInfoCurrent.recordColumnpool.clear();
    for (int k = 0; k < node.SubprobNum(); k++) {
        fixedInfoCurrent.recordColumnpool.push_back(node.MasterProblem().columnpool[k]);
    }

    //set spNoList and spIsFixed
    fixedInfoCurrent.spNoList.push_back(stoi(spNoSpilt));
    fixedInfoCurrent.spIsFixed[stoi(spNoSpilt)] = true;

    //set fixedColumnValueList
    std::vector<std::pair<int, double>> fixColTemp; //the order of varvalue is unimportant
    for (const auto& varPair : varValueFixed) {
        fixColTemp.push_back(varPair);
    }
    fixedInfoCurrent.fixedColumnValueList.push_back(fixColTemp);

    //set recordColumnpool
    //fixInfo.recordColumnpool have only one column in newest fixed subproblem
    ColumnPool columnpoolTemp;
    Column columnTemp;
    columnTemp.varValue = varValueFixed;
    columnpoolTemp.vars_id = fixedInfoCurrent.recordColumnpool[stoi(spNoSpilt)].vars_id;
    columnpoolTemp.var_num = fixedInfoCurrent.recordColumnpool[stoi(spNoSpilt)].var_num;
    columnpoolTemp.AddColumn(columnTemp);
    fixedInfoCurrent.recordColumnpool[stoi(spNoSpilt)] = columnpoolTemp;

    return 1;
}

void SortNodeColumn(const Node& node, std::vector<std::tuple<string, int, double>> &columnNameIndexTransolution)
{
    //sort the solution @wqy0530 lambda expressions might be used here
    for (int i = 0; i < node.LpRelaxSol().varValue.size(); i++){
        if (node.FixedInfo().spIsFixed[node.LpRelaxSol().subprob_no[i]]) {
            columnNameIndexTransolution.push_back(std::make_tuple(node.LpRelaxSol().varName[i], i, -1));//if fixed
        }else{
            columnNameIndexTransolution.push_back(std::make_tuple(node.LpRelaxSol().varName[i], i, node.LpRelaxSol().varValue[i]));
        }
    }

    sort(columnNameIndexTransolution.begin(), columnNameIndexTransolution.end(),
        [&](const std::tuple<string, int, double>& A, const std::tuple<string, int, double>& B) {return (std::get<2>(A) > std::get<2>(B)); });

}


void UpdateNodeNew(Node& node, const FIX& fixedInfos) //fix a columnpool recorded instead of a column
{
    //update fixed info
    node.SetFixedInfo(fixedInfos);

    //the newest fix column is fix a column, other fixed info fix columnpool
    //@wqy all of the columnpools are assigned the fix columnpool, because it is dealed before.
    for (int i = 0; i < node.SubprobNum(); i++){
        //judge whether the current subproblem already has a fixed column @wqy lambda may be used here
        node.SetMasterColumnool(i, fixedInfos.recordColumnpool[i]);
    }
}


void CalculateFixedInfo(std::stack<FIX>& fixedInformationList, const Node& node, int maxDiscrepancy)
{
    //Sort the columns in the node according to the relaxation solution of the master problem(if the relax solution equal to 0, we will not fix this column)
    //std::vector<std::pair<string, int>> columnSortNameIndex;
    //std::vector<std::pair<string,double>> transRMPSloution; //store RMP relax solution after conversion, -1 represent fixed already, other is assignment,pair<name,val>
    std::vector<std::tuple<string, int, double>>columnNameIndexTransolution;
    SortNodeColumn(node, columnNameIndexTransolution);

    int addFixListNum = 0;
    for (int colSortNo = 0; colSortNo < columnNameIndexTransolution.size();  colSortNo++) {
        if (std::get<2>(columnNameIndexTransolution[colSortNo]) <= 0) {
            break;//the column of subprob has fixed || relax solution of masterprob is 0
        }
        FIX fixedInfoCurrent; //fix information of current node
        //feasible after fixed return 1 infeasible after fixed return 0
        if (FindFixedInformation(node, std::get<0>(columnNameIndexTransolution[colSortNo]), fixedInfoCurrent)){
            //@wqy0605test Set fixedinformationlist.size() limit to solve memory out.
            if (fixedInformationList.size() <= 50){
                fixedInformationList.push(fixedInfoCurrent);//if feasible, push
            }
            addFixListNum++;
            if (addFixListNum == maxDiscrepancy) {
                break;//arrive maxDiscrepancy, break
            }
        }
    }
}


int FindFeasSol(const Node& node, SCIP* scipmodel, std::map<int, General_Block> block, double& obj)
{
    // here need to copy an another scipmodel, and kill it at the end of this function
    SCIP* scipModelTemp = nullptr;
    SCIPcreate(&scipModelTemp);
    SCIPcopyOrig(scipmodel, scipModelTemp, nullptr, nullptr, "c", 1, 0, 1, NULL);

    SCIPsetIntParam(scipModelTemp, "limits/solutions", -1);
    SCIP_VAR** varsOri = SCIPgetVars(scipModelTemp);
    int varsOriNum = SCIPgetNVars(scipModelTemp);
    for (int i = 0; i < node.SubprobNum(); i++){
        std::vector<std::pair<int,double>>columnCurrent = node.FixedInfo().fixedColumnValueList[i];
        for (int j = 0; j < columnCurrent.size(); j++) {
            if (!node.VarInt().find(columnCurrent[j].first)->second){
                continue;
            } 
            int fixedValue = int(columnCurrent[j].second);
            string varName;
            double lbNow = 0;
            double ubNow = 0;
            for (const auto& bvars : block[node.FixedInfo().spNoList[i]].bVars){
                if (bvars.Id == columnCurrent[j].first){
                    varName = bvars.Name;
                    lbNow = bvars.Lb;
                    ubNow = bvars.Ub;
                }
            }
            //@wqy just check, delete later
            if (fixedValue - lbNow < -0.000001 || fixedValue - ubNow > 0.0000001){
                cout << "Fix primal problem variable error!" << endl;
                continue;
            }
            for (int i = 0; i < varsOriNum; i++){
                if (SCIPvarGetName(varsOri[i]) == "t_" + varName || SCIPvarGetName(varsOri[i]) == varName){
                    SCIPchgVarLb(scipModelTemp, varsOri[i], fixedValue);
                    SCIPchgVarUb(scipModelTemp, varsOri[i], fixedValue);
                    break;
                }
            }
        }  
    }

    //SCIPwriteOrigProblem(scipModelTemp, "D://finalModel.lp", nullptr, 0);
    SCIPsolve(scipModelTemp);
    if (SCIPgetStatus(scipModelTemp) == SCIP_STATUS_INFEASIBLE){
        return 0;
    }
    obj = SCIPgetSolOrigObj(scipModelTemp, SCIPgetBestSol(scipModelTemp));
    SCIPfree(&scipModelTemp);
    return 1;

}