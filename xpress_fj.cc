/*
 * XPress integration for the Feasibility Jump heuristic.
 */

extern "C"
{

}

#include <cstdio>
#include <string>
#include <chrono>
#include <vector>
#include <cassert>
#include <mutex>
#include <thread>
#include <functional>
#include <atomic>
#include <cmath>
#include <climits>
#include <chrono>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <map>
#include <unordered_map>
#include <memory>
#include "feasibilityjump.hh"

const int NUM_THREADS = 1;

using namespace std;

std::atomic_size_t totalNumSolutionsFound(0);
std::atomic_size_t totalNumSolutionsAdded(0);
std::atomic_bool presolveFinished(false);
std::atomic_bool heuristicFinished(false);
std::chrono::steady_clock::time_point startTime;
const char* mps_filename;
unordered_map<string, int> Var_Ni;

struct _Solution
{
    std::vector<double> assignment;
    bool includesContinuous;
};

std::vector<_Solution> heuristicSolutions;
std::mutex heuristicSolutions_mutex;

std::mutex presolvedProblem_mutex;
std::mutex nonPresolvedProblem_mutex;

struct ProblemInstance
{
    int numCols;
    std::vector<char> varTypes;
    std::vector<double> lb;
    std::vector<double> ub;
    std::vector<double> objCoeffs;

    int numRows;
    int numNonZeros;
    std::vector<char> rowtypes;
    std::vector<double> rhs;
    std::vector<double> rhsrange;
    std::vector<double> lhs;
    std::vector<int> rowStart;
    std::vector<int> colIdxs;
    std::vector<double> colCoeffs;
};

struct FJData
{
    std::vector<int> originalIntegerCols;
    ProblemInstance originalData;
    ProblemInstance presolvedData;
};



std::string inputFilename;
std::string outDir;



ProblemInstance getXPRSProblemData(SCIP*& scip, map<int, SCIP_VAR*>& vars_order)
{
    ProblemInstance data;

    SCIPcreate(&scip);
    SCIPincludeDefaultPlugins(scip);
    SCIPcreateProbBasic(scip, "Prob");
    //D:\C++ Work\FJ_for_improve\FJ_for_improve\miplib2017benchmark
    SCIPreadProb(scip, mps_filename, NULL);

    //获得变量个数
    int num_vars = SCIPgetNVars(scip);
    data.numCols = num_vars;
    //获得变量类型
    data.varTypes = std::vector<char>(data.numCols);
    //获得变量下界
    data.lb = std::vector<double>(data.numCols);
    //获得变量上界
    data.ub = std::vector<double>(data.numCols);
    //获得变量目标系数
    data.objCoeffs = std::vector<double>(data.numCols);
    //注意：xpress根据变量id由小到大加入变量
    SCIP_VAR** vars = SCIPgetVars(scip);
    for (int i = 0; i < num_vars; i++) 
    {
        string name = SCIPvarGetName(vars[i]);
        int id = Var_Ni[name];
        vars_order.insert(pair<int, SCIP_VAR*>(id, vars[i]));
    }
    if (vars_order.size() == 1)
        int tong = 0;
    for (int i = 0; i < num_vars; i++) {
        SCIP_VARTYPE var_type = SCIPvarGetType(vars_order[i]);
        if (var_type == SCIP_VARTYPE_BINARY) {
            data.varTypes[i] = 'B';
        }
        else if (var_type == SCIP_VARTYPE_INTEGER) {
            data.varTypes[i] = 'I';
        }
        else if (var_type == SCIP_VARTYPE_CONTINUOUS) {
            data.varTypes[i] = 'C';
        }
        SCIP_Real lb = SCIPvarGetLbGlobal(vars_order[i]);
        data.lb[i] = lb;
        SCIP_Real ub = SCIPvarGetUbGlobal(vars_order[i]);
        data.ub[i] = ub;
        SCIP_Real coef = SCIPvarGetObj(vars_order[i]);
        data.objCoeffs[i] = coef;
    }

    //获得约束个数
    int num_conss = SCIPgetNConss(scip);
    data.numRows = num_conss;
    SCIP_CONS** conss = SCIPgetConss(scip);
    unsigned success;
    //获得约束类型
    data.rowtypes = std::vector<char>(data.numRows);
    //获得右端项
    data.rhs = std::vector<double>(data.numRows);
    //获得左端项
    data.lhs = std::vector<double>(data.numRows);
    //获得每条约束系数非零变量的个数
    data.rowStart = std::vector<int>(data.numRows + 1);
    data.rowStart[0] = 0;
    //获得每条约束系数非零变量的ID
    data.colIdxs;
    //获得每条约束系数非零变量的系数
    data.colCoeffs;
    for (int i = 0; i < num_conss; i++) {
        SCIP_CONS* cons = conss[i];
        auto Lb = SCIPconsGetLhs(scip, cons, &success);
        auto Ub = SCIPconsGetRhs(scip, cons, &success);
        // 判断约束类型
        if (Lb == Ub)
        {
            data.rowtypes[i] = 'E';
            data.rhs[i] = Lb;
            data.lhs[i] = Ub;
        }
        if (Lb > -SCIPinfinity(scip) && Ub == SCIPinfinity(scip))
        {
            data.rowtypes[i] = 'G';
            data.rhs[i] = Lb;
            data.lhs[i] = Ub;
        }
        if (Lb == -SCIPinfinity(scip) && Ub < SCIPinfinity(scip))
        {
            data.rowtypes[i] = 'L';
            data.rhs[i] = Lb;
            data.lhs[i] = Ub;
        }
        if (Lb > -SCIPinfinity(scip) && Ub < SCIPinfinity(scip) && Lb != Ub)
        {
            data.rowtypes[i] = 'R';
            data.rhs[i] = Lb;
            data.lhs[i] = Ub;
        }
        int nnonz = SCIPgetNVarsLinear(scip, cons);
        data.rowStart[i + 1] = data.rowStart[i] + nnonz;
        SCIP_VAR** vars = SCIPgetVarsLinear(scip, cons);
        SCIP_Real* coefs = SCIPgetValsLinear(scip, cons);
        for (int k = 0; k < nnonz; k++)
        {
            string name = SCIPvarGetName(vars[k]);
            int id = Var_Ni[name];
            data.colIdxs.push_back(id);
            data.colCoeffs.push_back(coefs[k]);
        }
    }
    return data;
}

bool copyDataToHeuristicSolver(FeasibilityJumpSolver &solver, ProblemInstance &data, int relaxContinuous)
{
    printf("initializing FJ with %d vars %d constraints\n", data.numCols, data.numRows);
    for (int colIdx = 0; colIdx < data.numCols; colIdx += 1)
    {
        varType vartype = varType::Continuous;
        if (data.varTypes[colIdx] == 'C')
        {
            vartype = varType::Continuous;
        }
        else if (data.varTypes[colIdx] == 'I')
        {
            vartype = varType::Integer;
        }
        else if (data.varTypes[colIdx] == 'B')
        {
            vartype = varType::Integer;
        }
        else
        {
            printf(FJ_LOG_PREFIX "unsupported variable type '%c' (%d).\n",
                   data.varTypes[colIdx], data.varTypes[colIdx]);
            return false;
        }

        solver.addVar(vartype, data.lb[colIdx], data.ub[colIdx], data.objCoeffs[colIdx]);
    }

    for (int rowIdx = 0; rowIdx < data.numRows; rowIdx += 1)
    {
        rowType rowtype;
        if (data.rowtypes[rowIdx] == 'N')
        {
            continue;
        }
        else if (data.rowtypes[rowIdx] == 'L')
        {
            rowtype = rowType::lte;
        }
        else if (data.rowtypes[rowIdx] == 'G')
        {
            rowtype = rowType::gte;
        }
        else if (data.rowtypes[rowIdx] == 'E')
        {
            rowtype = rowType::equal;
        }
        else if (data.rowtypes[rowIdx] == 'R')
        {
            // For the range constraint, we need two linear inequalities:
            // rhs - range <= lhs <= rhs

            solver.addConstraint(rowType::gte,
                                 data.lhs[rowIdx],
                                 data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
                                 &data.colIdxs[data.rowStart[rowIdx]],
                                 &data.colCoeffs[data.rowStart[rowIdx]],
                                 relaxContinuous);
            solver.addConstraint(rowType::lte,
                                 data.rhs[rowIdx],
                                 data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
                                 &data.colIdxs[data.rowStart[rowIdx]],
                                 &data.colCoeffs[data.rowStart[rowIdx]],
                                 relaxContinuous);
            continue;
        }
        else
        {
            printf(FJ_LOG_PREFIX "unsupported constraint type '%c'. Ignoring constraint.\n", data.rowtypes[rowIdx]);
            return false;
        }

        solver.addConstraint(rowtype,
                             data.rhs[rowIdx],
                             data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
                             &data.colIdxs[data.rowStart[rowIdx]],
                             &data.colCoeffs[data.rowStart[rowIdx]],
                             relaxContinuous);
    }

    return true;
}


// An object containing a function to be executed when the object is destructed.
struct Defer
{
    std::function<void(void)> func;
    Defer(std::function<void(void)> pFunc) : func(pFunc) {};
    ~Defer() { func(); }
};


void mapHeuristicSolution(FJStatus &status, bool usePresolved, SCIP* scip, map<int, SCIP_VAR*>& vars_order, FJData& gFJData)
{

    _Solution s;
    bool conversionOk = false;
    if (usePresolved)
    {
        int tong = 0;
    }
    else
    {
        printf(FJ_LOG_PREFIX "received a solution from non-presolved instance.\n");
        s.assignment = std::vector<double>(status.solution, status.solution + status.numVars);
        s.includesContinuous = true;
        conversionOk = true;
    }

    if (conversionOk)
    {
        {
            std::lock_guard<std::mutex> guard(heuristicSolutions_mutex);
            heuristicSolutions.push_back(s);
            totalNumSolutionsFound += 1;
        }
        //用scip固定整数解再求一下
        vector<SCIP_VAR*> vars;
        for (int i = 0; i < gFJData.originalData.numCols; i += 1)
        {
            vars.push_back(vars_order[i]);
            if (gFJData.originalData.varTypes[i] != 'C')
            {
                int fix_value = std::round(s.assignment[i]);
                SCIP_VAR* fix_var = vars_order[i];
                
                SCIPchgVarLbGlobal(scip, fix_var, fix_value);
                SCIPchgVarUbGlobal(scip, fix_var, fix_value);
            }
        }
        SCIPsolve(scip);
        SCIP_STATUS status = SCIPgetStatus(scip);
        if (status == SCIP_STATUS_OPTIMAL)
        {
            SCIPfreeTransform(scip);
            for (int j = 0; j < gFJData.originalData.numCols; j++)
            {
                string name_real = SCIPvarGetName(vars_order[j]);
                name_for_vars_values[name_real] = SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars_order[j]);
            }
            SCIP_SOL* sol;
            SCIPcreateSol(scip, &sol, NULL);
            for (int i = 0; i < vars.size(); i++)
            {
                string name = SCIPvarGetName(vars[i]);
                SCIPsetSolVal(scip, sol, vars[i], name_for_vars_values[name]);
            }
            unsigned int isFeas = 0;
            SCIP_Bool b = 0;
            SCIP_Bool c = 0;
            SCIPcheckSolOrig(scip, sol, &isFeas, b, c);

            best_sol_now = SCIPgetSolOrigObj(scip, sol);
            std::cout << "FJ find a obj:" << best_sol_now << std::endl;
            FJ_feas_flag = 1;
        }
        else
        {
            SCIPfreeTransform(scip);
            return;
        }
    }
}

const int maxEffort = 100000000;

// Starts background threads running the Feasibility Jump heuristic.
// Also installs the check-time callback to report any feasible solutions
// back to the MIP solver.
void start_feasibility_jump_heuristic(double* init_values, size_t maxTotalSolutions, bool heuristicOnly, SCIP* scip, map<int, SCIP_VAR*>& vars_order, FJData& gFJData, bool relaxContinuous = false, bool exponentialDecay = false, int verbose = 0)
{
    {
        auto allThreadsTerminated = std::make_shared<Defer>([]()
                                                            { heuristicFinished = true; });

        for (int thread_idx = 0; thread_idx < NUM_THREADS; thread_idx += 1)
        {
            auto seed = thread_idx;
            bool usePresolved = false;
            double decayFactor = (!exponentialDecay) ? 1.0 : 0.9999;
            if (scip == nullptr) {
                std::cout << "SCIP is empty!" << std::endl;
            }
            std::thread(
                [init_values, verbose, maxTotalSolutions, usePresolved, seed,
                 relaxContinuous, decayFactor, allThreadsTerminated, &scip, &vars_order, &gFJData]()
                {
                    // Prepare data for the non-presolved version.
                    {
                        std::lock_guard<std::mutex> guard(nonPresolvedProblem_mutex);
                        if (gFJData.originalData.numCols == 0)
                        {
                            gFJData.originalData = getXPRSProblemData(scip, vars_order);
                        }
                    }

                    if (scip == nullptr) {
                        std::cout << "data copied!" << std::endl;
                        std::cout << "SCIP is empty!" << std::endl;
                    }

                    // Produce the presolved solution
                    if (usePresolved)
                    {
                        std::lock_guard<std::mutex> guard(presolvedProblem_mutex);
                    }

                    ProblemInstance &data = usePresolved ? gFJData.presolvedData : gFJData.originalData;
                    FeasibilityJumpSolver solver(seed, verbose, decayFactor);
                    bool copyOk = copyDataToHeuristicSolver(solver, data, relaxContinuous);
                    if (!copyOk)
                        return;
                    std::cout << "目前vars_order.size()为" << vars_order.size() << std::endl;
                    solver.solve(
                        init_values, [maxTotalSolutions, usePresolved, &scip, &vars_order, &gFJData](FJStatus status) -> CallbackControlFlow
                        {
                            double time = std::chrono::duration_cast<std::chrono::milliseconds>(
                                std::chrono::steady_clock::now() - startTime).count() /1000.0;
                            if (vars_order.size() == 1)
                                int tong = 0;
                            // If we received a solution, put it on the queue.
                            if (status.solution != nullptr)
                            {
                                printf("FJSOL %g %g\n", time, status.solutionObjectiveValue);
    
                                mapHeuristicSolution(status, usePresolved, scip, vars_order, gFJData);
                            }
    
                            // If we have enough solutions or spent enough time, quit.
                            auto quitNumSol = totalNumSolutionsFound >= maxTotalSolutions;
                            if(quitNumSol) printf(FJ_LOG_PREFIX "quitting because number of solutions %zd >= %zd.\n", totalNumSolutionsFound.load(), maxTotalSolutions);
                            auto quitEffort = status.effortSinceLastImprovement > maxEffort;
                            if(quitEffort) printf(FJ_LOG_PREFIX "quitting because effort %d > %d.\n", status.effortSinceLastImprovement , maxEffort);
                            
                            auto quit = quitNumSol || quitEffort || heuristicFinished;
                            if (quit)
                                printf(FJ_LOG_PREFIX "effort rate: %g Mops/sec\n", status.totalEffort / time / 1.0e6);
                            return quit ? CallbackControlFlow::Terminate : CallbackControlFlow::Continue; });
                })
                .detach();
        }
    }

    if (heuristicOnly)
    {
        while (!heuristicFinished)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        printf(FJ_LOG_PREFIX "all threads exited.\n");
    }
}

int Use_FJ(unordered_map<string, int> In_Var_Ni, unordered_map<int, double> in_values, const char* filename)
{
    int verbose = 0;
    bool heuristicOnly = false;
    bool relaxContinuous = false;
    bool exponentialDecay = false;
    int timeout = INT32_MAX/2;

    std::string inputPath;
    heuristicOnly = true;

    Var_Ni = In_Var_Ni;
    double* init_values = new double[in_values.size()];
    for (int i = 0; i < in_values.size(); i++)
    {
        init_values[i] = in_values[i];
    }
   
    mps_filename = filename;
    startTime = std::chrono::steady_clock::now();

    SCIP* scip = nullptr;
    map<int, SCIP_VAR*> vars_order;
    FJData gFJData;
    gFJData.originalData.numCols = 0;

    start_feasibility_jump_heuristic(init_values, 1, heuristicOnly, scip, vars_order, gFJData, relaxContinuous, exponentialDecay, verbose);

    std::unordered_map<string, int>().swap(Var_Ni);
    std::vector<_Solution>().swap(heuristicSolutions);
    heuristicFinished = false;
    totalNumSolutionsFound = 0;

    if (FJ_feas_flag == 1)
    {
        SCIPfree(&scip);
    }

    return FJ_feas_flag;
}
