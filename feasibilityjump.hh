#include <algorithm>
#include <functional>
#include <vector>
#include <numeric>
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>
//#include <algorithm>

using namespace std;

#define FJ_LOG_PREFIX "Feasibility Jump: "

static unordered_map<string, double> name_for_vars_values;
static double best_sol_now = 0;
static double FJ_feas_flag = 0;

// 约束类型
enum rowType
{
    equal,
    lte,
    gte,
};
// 变量类型
enum varType
{
    Continuous,
    Integer
};

enum CallbackControlFlow
{
    Terminate,
    Continue,
};

struct FJStatus
{
    int totalEffort;
    int effortSinceLastImprovement;
    int numVars;
    double solutionObjectiveValue;
    double *solution;
};

const double violationTolerance = 1.0e-5;
const double equalityTolerance = 1.0e-5;

// Measures if two doubles are equal within a tolerance of 1.0e-5.
inline bool eq(double a, double b)
{
    return fabs(a - b) < equalityTolerance;
}

struct IdxCoeff
{
    uint32_t idx;
    double coeff;

    IdxCoeff(uint32_t idx, double coeff) : idx(idx), coeff(coeff) {}
};

struct Var
{
    varType vartype; // 变量类型（整数/连续）
    double lb; // 下界
    double ub; // 上界
    double objectiveCoeff; // 该变量在原目标函数的系数 
    std::vector<IdxCoeff> coeffs; // 该变量在约束中的系数
};

struct Constraint // score是约束结构的函数（约束违反量）
{
    rowType sense; // 约束类型，大于等于/小于等于/等于
    double rhs; // 该约束的b
    std::vector<IdxCoeff> coeffs; //约束中各变量的系数
    double weight; // 约束权重
    double incumbentLhs; // incumbentAssignment在该约束条件系数下的加和值( E:aij*xj )
    int32_t violatedIdx; // 所有违反约束列表中，此约束的索引

    // Computes the constraint's contribution to the feasibility score: 
    // If the constraint is satisfied by the given LHS value, returns 0.
    // If the constraint is violated by the given LHS value, returns -|lhs-rhs|.
    double score(double lhs) // Gj
    {
		//if (rhs < -1e19) {
			//return 0.;
		//}
        if (sense == rowType::equal)
            return -fabs(lhs - rhs);
        else if (sense == rowType::lte)
            return -(std::max(0., lhs - rhs));
        else
            return -(std::max(0., rhs - lhs));
    }
};

// A potential new value for a varaiable, including its score.
// 变量的潜在新值，包括其分数
struct Move
{
    double value;
    double score;

    static Move undef()
    {
        Move move;
        move.value = NAN;
        move.score = -std::numeric_limits<double>::infinity();
        return move;
    }
};

// Represents a modification of the LHS in a constraint, for a specific
// variable/constraint combination.The `modifyMove` function below is used to
// update the score of a `Move` to reflect the LHS modification.
struct LhsModification // 表示对特定变量/约束组合的约束中LHS的修改。下面的“modifyMove”函数用于更新“Move”的分数，以反映LHS修改
{
    uint32_t varIdx; // 变量索引
    uint32_t constraintIdx; // 约束索引
    double coeff; // 变量系数
    double oldLhs; // 原LHS
    double newLhs; // 修改后的LHS
};

// Stores the MIP problem, an incumbent assignment, and the set of constraints
// that are violated in the current incumbent assignment. This set is maintained
// when changes are given to the incumbent assignment using `setValue`.
// 存储MIP问题、当前解以及当前解中违反的约束集。当使用“setValue”对现有解进行更改时，将保留此集合。
struct Problem
{
    std::vector<Var> vars; // 变量
    std::vector<Constraint> constraints; // 约束
    std::vector<double> incumbentAssignment; // 当前解
    std::vector<uint32_t> violatedConstraints; //  违反的约束索引
    bool usedRelaxContinuous = false; 

    size_t nNonzeros; // 所有约束中系数非零变量的个数总和，即约束中包含的所有变量的总个数
    double incumbentObjective = NAN; // 当前原问题目标函数
	double incumbentFJObj = std::numeric_limits<double>::infinity(); // 所有约束违反量和
	// 添加变量
    int addVar(varType vartype, double lb, double ub, double objCoeff)
    {
        auto idx = vars.size();
        Var var;
        var.vartype = vartype;
        var.lb = lb;
        var.ub = ub;
        var.objectiveCoeff = objCoeff;
        vars.push_back(var);
        incumbentAssignment.push_back(lb);
        return idx;
    }

	// 添加约束
    int addConstraint(rowType sense, double rhs, int numCoeffs, int *rowVarIdxs, double *rowCoeffs, int relax_continuous)
    {
        if (relax_continuous)
            usedRelaxContinuous = true;

        // If we are relaxing continuous variables, an equality needs to be split into gte and lte.
        if (relax_continuous > 0 && sense == rowType::equal)
            if (std::any_of(rowVarIdxs, rowVarIdxs + numCoeffs, [&](double varIdx)
                            { return vars[varIdx].vartype == varType::Continuous; }))
            {
                addConstraint(rowType::gte, rhs, numCoeffs, rowVarIdxs, rowCoeffs, relax_continuous);
                addConstraint(rowType::lte, rhs, numCoeffs, rowVarIdxs, rowCoeffs, relax_continuous);
                return INT_MAX;
            }

        std::vector<IdxCoeff> coeffs;
        for (int i = 0; i < numCoeffs; i += 1)
        {
            if (relax_continuous > 0 && vars[rowVarIdxs[i]].vartype == varType::Continuous)
            {
                if (sense == rowType::lte)
                {
                    if (rowCoeffs[i] >= 0.)
                        rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
                    else
                        rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
                }
                else if (sense == rowType::gte)
                {
                    if (rowCoeffs[i] >= 0.)
                        rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
                    else
                        rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
                }
                else
                    return INT_MIN;
            }
            else
                coeffs.emplace_back(rowVarIdxs[i], rowCoeffs[i]);
        }

        if (coeffs.empty())
        {
            bool ok;
            if (sense == rowType::lte)
                ok = 0 <= rhs + equalityTolerance;
            else if (sense == rowType::gte)
                ok = 0 + equalityTolerance >= rhs;
            else
                ok = eq(0, rhs);

            return ok ? INT_MAX : INT_MIN;
        }

        int newConstraintIdx = constraints.size();
        for (auto &c : coeffs)
        {
            vars[c.idx].coeffs.emplace_back(newConstraintIdx, c.coeff);
        }

        nNonzeros += coeffs.size();
        Constraint newConstraint;
        newConstraint.coeffs = coeffs;
        newConstraint.incumbentLhs = NAN;
        newConstraint.violatedIdx = -1;
        newConstraint.rhs = rhs;
        newConstraint.sense = sense;
        newConstraint.weight = 1.0;

        constraints.push_back(newConstraint);
        return newConstraintIdx;
    }
	//<==============================================
	// 计算所有约束违反量和
	double calIncumbentObjective(vector<double> current_value, vector<Constraint> tempConstriant) {
		double currentFJObj = 0;
		int no_vc = 0;
		for (size_t cIdx = 0; cIdx < tempConstriant.size(); cIdx += 1) // 遍历约束 cIdx:约束索引
		{
			Constraint& cstr = tempConstriant[cIdx];// cstr:当前约束

			cstr.incumbentLhs = 0.0;
			// 计算当前解的LHS
			for (auto& vc : cstr.coeffs)
				cstr.incumbentLhs += vc.coeff * current_value[vc.idx];
			// this->incumbentFJObj += cstr.weight*(cstr.score(cstr.incumbentLhs));// 约束违反量加和	
			//if (-cstr.score(cstr.incumbentLhs) < 1.0e+19 && -cstr.score(cstr.incumbentLhs) > 0)
			if (cstr.score(cstr.incumbentLhs) < -violationTolerance)
			{
				//this->incumbentFJObj += -cstr.score(cstr.incumbentLhs);// 约束违反量加和
				currentFJObj += -cstr.weight * (cstr.score(cstr.incumbentLhs));// 约束违反量加和
				no_vc++;
			}
			/***else if(-cstr.score(cstr.incumbentLhs) > 1.0e+19)
			{
				int a = -cstr.score(cstr.incumbentLhs);
				int b = 100;
			}***/
			//else
			//{
			//	this->incumbentFJObj += -cstr.score(cstr.incumbentLhs);// 约束违反量加和
			//}
			//cout << cIdx << " -- " << this->incumbentFJObj << endl; //test
		}
		currentFJObj = currentFJObj * no_vc;
		return currentFJObj;
	}
	//==============================================>

	// 初始化问题（变量值、当前原问题目标函数、违反约束索引等）
    void resetIncumbent(double *initialValues)
    {
        // Set the initial values, if given.
        if (initialValues)
            for (size_t i = 0; i < vars.size(); i += 1)
                incumbentAssignment[i] = initialValues[i];
        // std::copy(initialValues, initialValues + vars.size(), incumbentAssignment);

        // Reset the incumbent objective.
        incumbentObjective = 0;
        for (size_t i = 0; i < vars.size(); i += 1)
            incumbentObjective += vars[i].objectiveCoeff * incumbentAssignment[i];

        // Reset the constraint LHSs and the violatedConstraints list.
        violatedConstraints.clear();
        for (size_t cIdx = 0; cIdx < constraints.size(); cIdx += 1) // 遍历约束 cIdx:约束索引
        {
            Constraint &cstr = constraints[cIdx];// cstr:当前约束

            cstr.incumbentLhs = 0.0;
			// 计算当前解的LHS
            for (auto &vc : cstr.coeffs)
                cstr.incumbentLhs += vc.coeff * incumbentAssignment[vc.idx];
			// 判断是否违反约束
            if (cstr.score(cstr.incumbentLhs) < -violationTolerance)
            {	
                cstr.violatedIdx = violatedConstraints.size(); 
                violatedConstraints.push_back(cIdx); // 将违反的约束索引送入violatedConstraints
            }
            else
                cstr.violatedIdx = -1;
        }
		int test = 1;
    }

    // Updates a variable assignment for `varIdx` to `newValue`.
    // Takes a function parameter f that receives a LhsModification
    // for every variable/constraint combination (except for `varIdx` itself)
    // where the LHS of the constraint has changed.
	// 将“varIdx”的变量赋值更新为“newValue”。取一个函数参数f，该参数接收约束的LHS发生变化的每个变量/约束组合（“varIdx”本身除外）的LhsModification。
    template <typename F>
    size_t setValue(uint32_t varIdx, double newValue, F f)
    {
        size_t dt = 0; // 该变量涉及到的全部约束中涉及到的全部变量的个数
        double oldValue = incumbentAssignment[varIdx];
        double delta = (newValue - oldValue); // varIdx变量的变化
        incumbentAssignment[varIdx] = newValue; // 更新当前解
        incumbentObjective += vars[varIdx].objectiveCoeff * delta; // 计算新的原问题目标函数
        // printf("Setting v%d to from %g to value %g\n", varIdx, oldValue, newValue);

        // Update the LHSs of all involved constraints.
		// 对涉及到该变量的所有约束更新LHS
        for (auto &cstrCoeff : vars[varIdx].coeffs) // 遍历该变量涉及到的全部约束
        {
            double oldLhs = constraints[cstrCoeff.idx].incumbentLhs;
            double newLhs = oldLhs + cstrCoeff.coeff * delta;
            constraints[cstrCoeff.idx].incumbentLhs = newLhs;
            double newCost = constraints[cstrCoeff.idx].score(newLhs);// 正为违反，负为满足

            // Add/remove from the violatedConstraints list.
            if (newCost < -violationTolerance && constraints[cstrCoeff.idx].violatedIdx == -1)
            {
                // Became violated.
                constraints[cstrCoeff.idx].violatedIdx = violatedConstraints.size();
                violatedConstraints.push_back(cstrCoeff.idx);
            }
            if (newCost >= -violationTolerance && constraints[cstrCoeff.idx].violatedIdx != -1)
            {
                // Became satisfied.
                auto lastViolatedIdx = violatedConstraints.size() - 1;
                auto lastConstraintIdx = violatedConstraints[lastViolatedIdx];
                auto thisViolatedIdx = constraints[cstrCoeff.idx].violatedIdx;
                std::swap(violatedConstraints[thisViolatedIdx], violatedConstraints[lastViolatedIdx]);
                constraints[lastConstraintIdx].violatedIdx = thisViolatedIdx;
                constraints[cstrCoeff.idx].violatedIdx = -1;
                violatedConstraints.pop_back();
            }

            // Now, report the changes in LHS for other variables.
			// 将该约束新的LHS传递给该约束涉及到的其他变量
            dt += constraints[cstrCoeff.idx].coeffs.size();
            for (auto &varCoeff : constraints[cstrCoeff.idx].coeffs) //  遍历该变量设计到的全部约束中的全部变量
            {
                if (varCoeff.idx != varIdx) // 如果遍历到的变量不是当前操作的变量
                {
                    LhsModification m;
                    m.varIdx = varCoeff.idx; // 变量索引
                    m.constraintIdx = cstrCoeff.idx; // 约束索引
                    m.coeff = varCoeff.coeff; // 变量在该约束下的系数
                    m.oldLhs = oldLhs;
                    m.newLhs = newLhs;
                    f(m);
                }
            }
        }

        return dt;
    }
};

// 更新“Move”的score，以反映LHS修改
// 当前操作变量jump后，涉及到的约束中其他变量move到潜在value后score进行更新
inline void modifyMove(LhsModification mod, Problem &problem, Move &move)
{
    Constraint &c = problem.constraints[mod.constraintIdx]; // 确定约束
    auto incumbent = problem.incumbentAssignment[mod.varIdx]; // 确定varIdx变量当前值
    double oldModifiedLhs = mod.oldLhs + mod.coeff * (move.value - incumbent); // 在当前正在操作的变量进行jump之前，该变量jump到潜在的move_value产生的Lhs
    double oldScoreTerm = c.weight * (c.score(oldModifiedLhs) - c.score(mod.oldLhs)); // 在当前正在操作的变量进行jump之前，该约束乘上权重的move后的分数
    double newModifiedLhs = mod.newLhs + mod.coeff * (move.value - incumbent); // 在当前正在操作的变量进行jump之后，该变量jump到潜在的move_value产生的Lhs
    double newScoreTerm = c.weight * (c.score(newModifiedLhs) - c.score(mod.newLhs)); // 在当前正在操作的变量进行jump之后，该约束乘上权重的move后的分数
    move.score += newScoreTerm - oldScoreTerm; // 在当前操作变量jump后，该变量在该约束上move的的sj，因为其他约束的sj不变，所以score为旧score+当前约束的sj变化量
}

// Stores current moves and computes updated jump values for
// the "Jump" move type.
// 存储当前移动并计算“跳跃”移动类型的更新跳跃值。
class JumpMove
{
    std::vector<Move> moves; // 每个变量潜在value值和score
    std::vector<std::pair<double, double>> bestShiftBuffer; // 论文中的算法1的delta

public:
    void init(Problem &problem) // 初始化move大小
    {
        moves.resize(problem.vars.size());
    }

	// 对varIdx变量，对其move进行操作
    template <typename F>
    void forEachVarMove(int32_t varIdx, F f)
    {
        f(moves[varIdx]);
    }

	// 更新 varIdx变量的jump_value （论文中的算法1）
    void updateValue(Problem &problem, uint32_t varIdx)
    {
        bestShiftBuffer.clear();
        auto varIncumbentValue = problem.incumbentAssignment[varIdx];
        double currentValue = problem.vars[varIdx].lb;
        double currentScore = 0.0;
        double currentSlope = 0.0;

        // printf(" updatevalue lb %g ub %g numcells %d\n",
        //        problem.vars[varIdx].lb,
        //        problem.vars[varIdx].ub, problem.vars[varIdx].coeffs.size());

        for (auto &cell : problem.vars[varIdx].coeffs)
        {
            auto &constraint = problem.constraints[cell.idx];

            std::vector<std::pair<double, double>> constraintBounds;
            if (constraint.sense == rowType::lte)
                constraintBounds.emplace_back(-std::numeric_limits<double>::infinity(), constraint.rhs);
            else if (constraint.sense == rowType::gte)
                constraintBounds.emplace_back(constraint.rhs, std::numeric_limits<double>::infinity());
            else
            {
                constraintBounds.emplace_back(-std::numeric_limits<double>::infinity(), constraint.rhs);
                constraintBounds.emplace_back(constraint.rhs, constraint.rhs);
                constraintBounds.emplace_back(constraint.rhs, std::numeric_limits<double>::infinity());
            }

            for (auto &bound : constraintBounds)
            {
                double residualIncumbent = constraint.incumbentLhs - cell.coeff * varIncumbentValue;

                std::pair<double, double> validRange = {
                    ((1.0 / cell.coeff) * (bound.first - residualIncumbent)),
                    ((1.0 / cell.coeff) * (bound.second - residualIncumbent)),
                };

                if (problem.vars[varIdx].vartype == varType::Integer)
                    validRange = {
                        std::ceil(validRange.first - equalityTolerance),
                        std::floor(validRange.second + equalityTolerance),
                    };

                if (validRange.first > validRange.second)
                    continue;

                if (validRange.first > currentValue)
                {
                    currentScore += constraint.weight * (validRange.first - currentValue);
                    currentSlope -= constraint.weight;
                    if (validRange.first < problem.vars[varIdx].ub)
                        bestShiftBuffer.emplace_back(validRange.first, constraint.weight);
                }

                if (validRange.second <= currentValue)
                {
                    currentScore += constraint.weight * (validRange.second - currentValue);
                    currentSlope += constraint.weight;
                }
                else if (validRange.second < problem.vars[varIdx].ub)
                    bestShiftBuffer.emplace_back(validRange.second, constraint.weight);
            }
        }

        bestShiftBuffer.emplace_back(problem.vars[varIdx].lb, 0);
        bestShiftBuffer.emplace_back(problem.vars[varIdx].ub, 0);
        std::sort(bestShiftBuffer.begin(), bestShiftBuffer.end());

        double bestScore = currentScore;
        double bestValue = currentValue;
        // printf("evaluating best shift buffer size %d \n", bestShiftBuffer.size());
        for (auto &item : bestShiftBuffer)
        {
            currentScore += (item.first - currentValue) * currentSlope;
            currentSlope += item.second;
            currentValue = item.first;

            // printf("bestshift cscore %g cslope %g cval %g bestval %g bestscore %g\n",
            // currentScore,currentSlope, currentValue, bestScore, bestValue
            // );

            if (eq(bestValue, problem.incumbentAssignment[varIdx]) ||
                (!eq(currentValue, problem.incumbentAssignment[varIdx]) && currentScore < bestScore))
            {
                bestScore = currentScore;
                bestValue = currentValue;
            }

            // Slope is always increasing, so if we have a valid value, we can quit
            // as soon as the slope turns nonnegative, since we must already have
            // visited the minimum.
            if (!eq(bestValue, problem.incumbentAssignment[varIdx]) && currentSlope >= 0.)
                break;
        }

        // printf("Setting jump for %d to from %g to %g\n", varIdx, problem.incumbentAssignment[varIdx], moves[varIdx].value);
        moves[varIdx].value = bestValue;
    }
};

class FeasibilityJumpSolver
{
    int verbosity; // 输出等级
    Problem problem; // 问题
    JumpMove jumpMove; // jumpvalue

    std::vector<uint32_t> goodVarsSet; // 好的变量索引（P，算法2）
    std::vector<int32_t> goodVarsSetIdx;  // 好的变量索引集中变量索引的索引

    std::mt19937 rng;

    double bestObjective = std::numeric_limits<double>::infinity(); // 最好的目标函数
	double objectiveWeight = 0.0; // 原问题目标函数的权重 = 0
    size_t bestViolationScore = SIZE_MAX; 
    size_t effortAtLastCallback = 0; 
    size_t effortAtLastImprovement = 0;
    size_t totalEffort = 0;

    double weightUpdateDecay;
    double weightUpdateIncrement = 1.1;

    size_t nBumps;

    // The probability of choosing a random positive-score variable.
    const double randomVarProbability = 0.001;

    // The probability of choosing a variable using a random constraint's
    // non-zero coefficient after updating weights.
	const double randomCellProbability = 0.01;

    // The number of moves to evaluate, if there are many positive-score
    // variables available.
    const size_t maxMovesToEvaluate = 25;


	//<==============================================
	int Gmax = 50;
	int NP = 10;		//种群规模不是越大越好（10、30、50、100、1000）
	int D;
	int indexBest = 0;
	double F = 0.3;	//偏向局部搜索（0.1、0.3、0.5）
	double CR = 0.3;	//偏向局部搜索（0.1、0.3、0.5）
	int rand1, rand2, rand3;
	vector<vector<double>> X, V, U;
	vector<double> Fx, Fu;
	double fitnessAve, fitnessMini;		//先把太大的违反量（1E+10~1E+20）屏蔽了，？
	//==============================================>



public:
	// 构造函数
    FeasibilityJumpSolver(int seed = 0, int _verbosity = 0, double _weightUpdateDecay = 1.0)
    {
        verbosity = _verbosity;
        weightUpdateDecay = _weightUpdateDecay;
        rng = std::mt19937(seed);
    }

    int addVar(varType vartype, double lb, double ub, double objCoeff)
    {
        goodVarsSetIdx.push_back(-1);
        return problem.addVar(vartype, lb, ub, objCoeff);
    }

    int addConstraint(rowType sense, double rhs, int numCoeffs, int *rowVarIdxs, double *rowCoeffs, int relax_continuous)
    {
        return problem.addConstraint(sense, rhs, numCoeffs, rowVarIdxs, rowCoeffs, relax_continuous);
    }

	//<==============================================
	/***double calculateFitness(vector<double> individual, vector<Constraint> currentConstraint)
	{
		for (int j = 0; j < D; j++)
		{
			problem.incumbentAssignment[j] = individual[j];
		}
		int no_vc = problem.calIncumbentObjective(); // 测试计算目标函数
		//printf("current FJObj:  %f\n", problem.incumbentFJObj); // 打印	
		return problem.incumbentFJObj;
	}***/

	//==============================================>



	// 求解
	int solve(double* initialValues, std::function<CallbackControlFlow(FJStatus)> callback)
	{
		assert(callback);
		if (verbosity >= 1)
			printf(FJ_LOG_PREFIX "starting solve. weightUpdateDecay=%g, relaxContinuous=%d  \n", weightUpdateDecay, problem.usedRelaxContinuous);
		init(initialValues);

		double* currentValue = new double[100000];
		std::vector<double> current_weight;
		D = problem.vars.size();

		for (int i = 0; i < D; i += 1) 
		{
			currentValue[i] = problem.incumbentAssignment[i];
		}
		
		bool isterminate = false;

		for (int k = 0; k < 27; k++)
		{
			if (isterminate)
				break;

			for (int i = 0; i < D; i++) {
				goodVarsSetIdx[i] = -1;
			}
			for (size_t i = 0; i < problem.vars.size(); i += 1)
			{
				//if (i == 203)
					//int a = 100;
				//cout << i << endl;
				jumpMove.updateValue(problem, i);
			}
			/***
			int ii = 0;
			for (auto& cIdx : problem.violatedConstraints)// 遍历所有违反约束的索引
			{
				problem.constraints[cIdx].weight = current_weight[ii++];
			}
			
			current_weight.clear();***/
			//<==============================================
			vector<double> left, right, mut;//决策空间上下界		
			
			
			// 确定搜索空间
			left = problem.incumbentAssignment;
			for (int i = 0; i < D; i++)
			{
				auto m = bestMove(i);
				auto Value = m.value;
				if (Value > left[i])
				{
					right.push_back(Value);
				}
				else
				{
					double temp;
					temp = Value;
					right.push_back(left[i]);
					left[i] = temp;
				}
				mut.push_back(Value);
			}
			// 初始化种群
			for (int i = 0; i < NP; i++)
			{
				vector<double> individual;
				if (i == 0)
				{
					individual = problem.incumbentAssignment;
					X.push_back(individual);
					V.push_back(individual);
					U.push_back(individual);
					continue;
				}
				for (int j = 0; j < D; j++)
				{
					individual.push_back(left[j] + (rand() / (RAND_MAX + 1.0)) * (right[j] - left[j]));
					if (problem.vars[j].vartype == Integer)
					{
						individual[j] = round(individual[j]);
					}					
				}
				X.push_back(individual);
				V.push_back(individual);
				U.push_back(individual);
			}

			for (int i = 0; i < NP; i++)
			{
				Fx.push_back(problem.calIncumbentObjective(X[i], problem.constraints));
				Fu.push_back(Fx[i]);
			}
			
			for (int n = 0; n < Gmax; n++)
			{
				// mutation
				for (int i = 0; i < NP; i++)
				{
					do
					{
						rand1 = rand() % NP;
						rand2 = rand() % NP;
						rand3 = rand() % NP;
					} while (rand1 == rand2 || rand3 == rand2 || rand1 == rand3 || rand1 == i || rand2 == i || rand3 == i);

					if ((rand() / (RAND_MAX + 1.0)) < 0.2)
					{
						for (int j = 0; j < D; j++)
						{
							V[i][j] = X[indexBest][j] + F * (X[rand2][j] - X[rand3][j]);
						}
					}
					else
					{
						if ((rand() / (RAND_MAX + 1.0)) < 0.3)
						{
							for (int j = 0; j < D; j++)
							{
								V[i][j] = X[indexBest][j] + (rand() / (RAND_MAX + 1.0)) * 10;
							}
						}
						else
						{
							for (int j = 0; j < D; j++)
							{
								if ((rand() / (RAND_MAX + 1.0)) < 0.1)
								{
									V[i][j] = X[indexBest][j] + (rand() / (RAND_MAX + 1.0)) * 10;
								}
							}
						}
					}
				}
				// check
				for (int i = 0; i < NP; i++)
				{
					for (int j = 0; j < D; j++)
					{
						if (V[i][j]< left[j] || V[i][j] > right[j])
						{
							//V[i][j] = problem.vars[j].lb + (rand() / (RAND_MAX + 1.0)) * (problem.vars[j].ub - problem.vars[j].lb);
							V[i][j] = left[j] + (rand() / (RAND_MAX + 1.0)) * (right[j] - left[j]);
						}
						if (problem.vars[j].vartype == Integer)
						{
							V[i][j] = round(V[i][j]);
						}
					}
				}
				// crossover
				for (int i = 0; i < NP; i++)
				{
					for (int j = 0; j < D; j++)
					{
						if ((rand() / (RAND_MAX + 1.0)) < CR)
						{
							U[i][j] = V[i][j];
						}
						else
						{
							U[i][j] = X[i][j];
						}
					}
				}
				// selelction
				for (int i = 0; i < NP; i++)
				{
					Fu[i] = problem.calIncumbentObjective(U[i], problem.constraints);
					if (Fu[i] <= Fx[i])
					{
						Fx[i] = Fu[i];
						for (int j = 0; j < D; j++)
						{
							X[i][j] = U[i][j];
						}
					}
				}

				fitnessAve = Fx[0];
				fitnessMini = Fx[0];
				indexBest = 0;
				for (int i = 1; i < NP; i++)
				{
					fitnessAve += Fx[i];
					if (fitnessMini > Fx[i])
					{
						fitnessMini = Fx[i];
						indexBest = i;
					}
				}
				fitnessAve = fitnessAve / NP;
				cout << k << "  " << n << "	M: " << fitnessMini << "	A: " << fitnessAve << endl;
				//printf("current FJObj:  %f\n", fitnessAve); // 打印
			}
			
			left.clear(); right.clear(); mut.clear(); 
			//system("pause");
			//==============================================>

			for (int i = 0; i < D; i += 1)
			{
				currentValue[i] = X[indexBest][i];
			}

			for (int i = 0; i < D; i++) {
				goodVarsSetIdx[i] = -1;
			}
			X.clear(); V.clear(); U.clear(); Fx.clear(); Fu.clear();

			init(currentValue);


			for (int step = 0; step < 100001; step += 1)
			{
				if (user_terminate(callback, nullptr)) 
				{
					isterminate = true;
					break;
				}
				if (step % 10000 == 0)
				{
					if (verbosity >= 1)
						printf(FJ_LOG_PREFIX "step %d viol %zd good %zd bumps %zd\n", step+k*100000, problem.violatedConstraints.size(), goodVarsSet.size(), nBumps);
				}

				if (problem.violatedConstraints.size() < bestViolationScore)
				{
					effortAtLastImprovement = totalEffort;
					bestViolationScore = problem.violatedConstraints.size();
				}
				
				if (problem.violatedConstraints.empty() && problem.incumbentObjective < bestObjective)
				{
					effortAtLastImprovement = totalEffort;
					bestObjective = problem.incumbentObjective;
					if (user_terminate(callback, problem.incumbentAssignment.data()))
						isterminate = true;
						break;
				}
				

				if (problem.vars.size() == 0)
				{
					isterminate = true;
					break;
				}

				uint32_t var = selectVariable();
				doVariableMove(var);

				//problem.calIncumbentObjective(initialValues); // 测试计算目标函数
				//printf("current FJObj:  %f\n", problem.incumbentFJObj); // 打印
			}
			for (int i = 0; i < D; i += 1)
				currentValue[i] = problem.incumbentAssignment[i];
			for (auto& cIdx : problem.violatedConstraints)// 遍历所有违反约束的索引
			{
				auto& constraint = problem.constraints[cIdx];
				current_weight.push_back(constraint.weight);
				//cout << cIdx << endl;
			}
		}

        return 0;
    }

private:
	// 初始化问题，对每个变量都寻找一个jumpvalue并计算score，并根据score正负判断是否放入P
    void init(double *initialValues)
    {
        problem.resetIncumbent(initialValues); // 初始化问题
        jumpMove.init(problem); // 初始化move大小
        totalEffort += problem.nNonzeros; 

        // Reset the variable scores.
        goodVarsSet.clear();
		for (size_t i = 0; i < problem.vars.size(); i += 1) 
		{
			//if (i == 203)
				//int a = 100;
			//cout << i << endl;
			resetMoves(i);
		}
           
    }

    uint32_t selectVariable()
    {
        if (!goodVarsSet.empty())
        {
            if (std::uniform_real_distribution<double>(0., 1.)(rng) < randomVarProbability)
                return goodVarsSet[rng() % goodVarsSet.size()];

            auto sampleSize = std::min(maxMovesToEvaluate, goodVarsSet.size());
            totalEffort += sampleSize;

            double bestScore = -std::numeric_limits<double>::infinity();
            uint32_t bestVar = UINT_MAX;
            for (size_t i = 0; i < sampleSize; i++)
            {
                auto setidx = rng() % goodVarsSet.size();
                auto varIdx = goodVarsSet[setidx];
                // assert(goodVarsSetIdx[varIdx] >= 0 && goodVarsSetIdx[varIdx] == setidx);
                Move move = bestMove(varIdx);
                // assert(move.score > equalityTolerance);
                if (move.score > bestScore)
                {
                    bestScore = move.score;
                    bestVar = varIdx;
                }
            }
            assert(bestVar != UINT_MAX);
            return bestVar;
        }

        // Local minimum, update weights.
        updateWeights();

        if (!problem.violatedConstraints.empty())
        {
            size_t cstrIdx = problem.violatedConstraints[rng() % problem.violatedConstraints.size()];
            auto &constraint = problem.constraints[cstrIdx];

            if (std::uniform_real_distribution<double>(0., 1.)(rng) < randomCellProbability)
                return constraint.coeffs[rng() % constraint.coeffs.size()].idx;

            double bestScore = -std::numeric_limits<double>::infinity();
            uint32_t bestVarIdx = UINT_MAX;

            for (auto &cell : constraint.coeffs)
            {
                Move move = bestMove(cell.idx);
                if (move.score > bestScore)
                {
                    bestScore = move.score;
                    bestVarIdx = cell.idx;
                }
            }
            return bestVarIdx;
        }

        // Fallback to random choice.
        return rng() % problem.vars.size();
    }

	// 更新权重并重新计算分数，更新P
    void updateWeights()
    {
        if (verbosity >= 2)
            printf(FJ_LOG_PREFIX "Reached a local minimum.\n");
        nBumps += 1;
        bool rescaleAllWeights = false;
        size_t dt = 0;

		// 已经找到可行解
        if (problem.violatedConstraints.empty())
        {
            objectiveWeight += weightUpdateIncrement;// 原问题目标函数权重+1
            if (objectiveWeight > 1.0e20)
                rescaleAllWeights = true;

            dt += problem.vars.size();
			// 重新计算score
            for (size_t varIdx = 0; varIdx < problem.vars.size(); varIdx += 1)
                forEachMove(
                    varIdx, [&](Move &move)
                    { move.score += weightUpdateIncrement *
                                    problem.vars[varIdx].objectiveCoeff *
                                    (move.value - problem.incumbentAssignment[varIdx]); });
        }
		// 未找到可行解
        else
        {
            for (auto &cIdx : problem.violatedConstraints)// 遍历所有违反约束的索引
            {
                auto &constraint = problem.constraints[cIdx];
                constraint.weight *= weightUpdateIncrement; // 违反的约束权重+1
                if (constraint.weight > 1.0e20)
                    rescaleAllWeights = true;

                dt += constraint.coeffs.size();
                for (auto &cell : constraint.coeffs) // 遍历约束中所有的变量，更新score，更新P
                {
                    forEachMove(
                        cell.idx, [&](Move &move)
                        {
                            double candidateLhs = constraint.incumbentLhs + cell.coeff * (move.value - problem.incumbentAssignment[cell.idx]);
                            double diff = weightUpdateIncrement * (constraint.score(candidateLhs) -
                                                                   constraint.score(constraint.incumbentLhs));
                            move.score += diff; });

                    updateGoodMoves(cell.idx);
                }
            }
        }

        weightUpdateIncrement /= weightUpdateDecay;
        if (rescaleAllWeights)
        {
            weightUpdateIncrement *= 1.0e-20;
            objectiveWeight *= 1.0e-20;

            for (auto &c : problem.constraints)
                c.weight *= 1.0e-20;
            dt += problem.constraints.size();

            for (size_t i = 0; i < problem.vars.size(); i += 1)
                resetMoves(i);
        }

        totalEffort += dt;
    }

	// 返回varIdx变量的move
    Move bestMove(uint32_t varIdx)
    {
        Move best = Move::undef();
        forEachMove(varIdx, [&](Move &move)
                    { if (move.score > best.score)
                        best = move; });

        return best;
    }

	// 跳跃，更新varIdx的jumpvalue，更新其他所有的move.score，更新P
    void doVariableMove(uint32_t varIdx)
    {
        // First, we get the best move for the variable;
        auto m = bestMove(varIdx);
        auto newValue = m.value;
        // assert(!isnan(newValue));

        // Update the incumbent solution.
        // printf("Setting var %d from %g to %g for a score of %g\n", varIdx, oldValue, newValue, m.score);

        totalEffort += problem.setValue(
            varIdx, newValue, [&](LhsModification mod)
            {
                forEachMove(mod.varIdx, [&](Move &m)
                            { modifyMove(mod, problem, m); });
                updateGoodMoves(mod.varIdx); });

        resetMoves(varIdx);
    }

	// 将根据score判断是否放入或移出P（算法2）
    void updateGoodMoves(int32_t varIdx)
    {
        bool anyGoodMoves = bestMove(varIdx).score > 0.; // 判断该变量的score是否为正，正说明跳跃后有改进
        if (anyGoodMoves && goodVarsSetIdx[varIdx] == -1)
        {
            // Became good, add to good set.
            goodVarsSetIdx[varIdx] = goodVarsSet.size();
            goodVarsSet.push_back(varIdx);
        }
        else if (!anyGoodMoves && goodVarsSetIdx[varIdx] != -1)
        {
            // Became bad, remove from good set.
            auto lastSetIdx = goodVarsSet.size() - 1;
            auto lastVarIdx = goodVarsSet[lastSetIdx];
            auto thisSetIdx = goodVarsSetIdx[varIdx];
			if (goodVarsSet.size() < thisSetIdx+1) 
			{
				thisSetIdx = goodVarsSet.size() - 1;
			}
            std::swap(goodVarsSet[thisSetIdx], goodVarsSet[lastSetIdx]);
            goodVarsSetIdx[lastVarIdx] = thisSetIdx;
            goodVarsSetIdx[varIdx] = -1;
            goodVarsSet.pop_back();
        }
    }

	// 对varIdx变量，对其move进行操作
    template <typename F>
    void forEachMove(int32_t varIdx, F f)
    {
        jumpMove.forEachVarMove(varIdx, f);

        // TODO: here, we can add more move types.
        // upDownMove.forEachVarMove(varIdx, f);
    }

	// 寻找新的jumpvalue并计算move.score,并根据score判断是否将变量加入good
    void resetMoves(uint32_t varIdx)
    {
        totalEffort += problem.vars[varIdx].coeffs.size();
        jumpMove.updateValue(problem, varIdx);  // 寻找该索引变量的jumpvalue

        forEachMove(
            varIdx, [&](Move &move)
            {
                move.score = 0.0;
                move.score += objectiveWeight *
                 problem.vars[varIdx].objectiveCoeff * 
                 (move.value - problem.incumbentAssignment[varIdx]); // 分数为（jump值-原值）* 原问题目标函数该变量的系数* 原问题系数 = 0

                for (auto &cell : problem.vars[varIdx].coeffs) // 遍历当前变量涉及到的全部约束
                {
                    auto &constraint = problem.constraints[cell.idx]; // 当前变量涉及到的第idx个约束
                    auto candidateLhs = constraint.incumbentLhs +
                                        cell.coeff *
                                        (move.value - problem.incumbentAssignment[varIdx]); // 跳跃后的该约束Lhs值
                    move.score += constraint.weight *
                     (constraint.score(candidateLhs) - constraint.score(constraint.incumbentLhs)); // 各个约束乘权重的分数（sj  公式8）加和
                } });

        updateGoodMoves(varIdx);
    }

    bool user_terminate(std::function<CallbackControlFlow(FJStatus)> callback, double *solution)
    {
        const int CALLBACK_EFFORT = 500000;
        if (solution != nullptr || totalEffort - effortAtLastCallback > CALLBACK_EFFORT)
        {
            if (verbosity >= 2)
                printf(FJ_LOG_PREFIX "calling user termination.\n");
            effortAtLastCallback = totalEffort;

            FJStatus status;
            status.totalEffort = totalEffort;
            status.effortSinceLastImprovement = totalEffort - effortAtLastImprovement;

            status.solution = solution;
            status.numVars = problem.vars.size();
            status.solutionObjectiveValue = problem.incumbentObjective;

            auto result = callback(status);
            if (result == CallbackControlFlow::Terminate)
            {
                if (verbosity >= 2)
                    printf(FJ_LOG_PREFIX "quitting.\n");
                return true;
            }
        }
        return false;
    }
};

int Use_FJ(unordered_map<string, int> In_Var_Ni, unordered_map<int, double> in_values, const char* filename);
