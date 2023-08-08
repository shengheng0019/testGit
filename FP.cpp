#include "FP.h"

template<class _T>
class equalkey_for_pairs {
private:
    typedef std::pair<std::vector<_T>, _T> Key;
public:
    bool operator()(const Key& k1, const Key& k2) const
    {
        return (k1.first == k2.first) && (fabs(k1.second - k2.second) <= FPPARA::alphadist);
    };
};

template<class _T>
class hash_compare_for_pairs {
private:
    typedef std::pair<std::vector<_T>, _T> Key;
public:
    static const size_t bucket_size = 4;
    static const size_t min_buckets = 20000;
    size_t operator()(const Key& k) const {
        size_t h = 0;
        for (int i = 0; i < k.first.size(); i++)
            h = ((h << 1) | ((h >> 31) & 1)) ^ (long)k.first[i];
        return h;
    };
};


typedef std::pair<double, int> dipair;	// (sigma, index) pair and comparison function
class dipairCmp {
public:	bool operator() (const dipair& p1, const dipair& p2) const { return p1.first > p2.first; }

};


FP::FP(std::string fileName)
        :alpha{ FPPARA::alpha }, maxIter1{FPPARA::maxIter1}, epInt{FPPARA::epInt}, maxit1wi{FPPARA::maxit1wi},
         alpha_quot{ FPPARA::alpha_quot }, trsrnd{ FPPARA::trsrnd}, minChangeNum{ FPPARA::minChangeNum}, trsmin{ FPPARA::trsmin}
{
    //init scip
    scip = NULL;
    if (SCIPcreate(&scip) != SCIP_RETCODE::SCIP_OKAY)
        std::cout << "create model error" << std::endl;
    SCIPincludeDefaultPlugins(scip);
    SCIPreadProb(scip, fileName.c_str(), NULL);
    SCIPwriteOrigProblem(scip, "aflow30a.lp", NULL, false);
    SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), 1);
    offSet =  SCIPgetOrigObjoffset(scip);

    //load vars, conss and var bound
    SCIP_VAR** vars_ori = SCIPgetVars(scip);
    SCIP_CONS** cons_ori = SCIPgetConss(scip);
    int nvars = SCIPgetNVars(scip);
    int nconss = SCIPgetNConss(scip);
    for (int i = 0; i < nvars; i++)
    {
        vars.push_back(vars_ori[i]);
        if (SCIPvarGetType(vars_ori[i]) == SCIP_VARTYPE::SCIP_VARTYPE_CONTINUOUS)
        {
            ;
        }
        else if (SCIPvarGetType(vars_ori[i]) == SCIP_VARTYPE::SCIP_VARTYPE_BINARY)
        {
            intVars.push_back(i);
            isBin.push_back(true);
            double l = floor(SCIPvarGetLbOriginal(vars[i]) + epInt);
            double u = floor(SCIPvarGetUbOriginal( vars[i]) + epInt);
            intVarsLB.push_back(l);
            intVarsUB.push_back(u);
        }
        else if (SCIPvarGetType(vars_ori[i]) == SCIP_VARTYPE::SCIP_VARTYPE_INTEGER)
        {
            intVars.push_back(i);
            isBin.push_back(false);
            double l = floor(SCIPvarGetLbOriginal(vars[i]) + epInt);
            double u = floor(SCIPvarGetUbOriginal(vars[i]) + epInt);
            intVarsLB.push_back(l);
            intVarsUB.push_back(u);
        }
    }

    for (int i = 0; i < nconss; i++)
        conss.push_back(cons_ori[i]);

    //relax the model
    for (int ii = 0; ii < nvars; ii++)
    {
        auto var = vars[ii];
        SCIP_VARTYPE const varType_ = SCIPvarGetType(var);
        if (varType_ == SCIP_VARTYPE_CONTINUOUS)
            continue;

        SCIP_Bool infeasible;
        SCIPchgVarType(scip, var, SCIP_VARTYPE_CONTINUOUS, &infeasible);
    }

    //
    oriObj.clear();
    for (int i = 0; i < vars.size(); i++)
    {
        oriObj.push_back(SCIPvarGetObj(vars[i]));
    }

    roundedIntVars.resize(intVars.size());

    //init auxilary variables and constraints
    auxVars.resize(intVars.size());
    auxLBConss.resize(intVars.size());
    auxUBConss.resize(intVars.size());
}

void FP::IntroduceAux(int intVarIdx)
{
    //introduce auxilary variables and constraints
    SCIP_VAR** varsTemp = new SCIP_VAR * [2];
    varsTemp[0] = vars[intVarIdx];
    SCIP_VAR* tempVar;
    SCIPcreateVarBasic(scip, &tempVar, ("auxVar_" + std::to_string(intVarIdx)).c_str(), 0, SCIPinfinity(scip), 1, SCIP_VARTYPE_CONTINUOUS);
    SCIPaddVar(scip, tempVar);
    varsTemp[1] = tempVar;
    auxVars[intVarIdx] = tempVar;

    double* val_ = new double[2];

    val_[0] = 1, val_[1] = -1;
    SCIP_CONS* tempUBCons;
    SCIPcreateConsBasicLinear(scip, &tempUBCons, ("auxUbCons_" + std::to_string(intVarIdx)).c_str(), 2, varsTemp, val_, -SCIPinfinity(scip), SCIPinfinity(scip));
    SCIPaddCons(scip, tempUBCons);
    auxUBConss[intVarIdx] = tempUBCons;

    val_[0] = 1, val_[1] = 1;
    SCIP_CONS* tempLBCons;
    SCIPcreateConsBasicLinear(scip, &tempLBCons, ("auxLbCons_" + std::to_string(intVarIdx)).c_str(), 2, varsTemp, val_, -SCIPinfinity(scip), SCIPinfinity(scip));
    SCIPaddCons(scip, tempLBCons);
    auxLBConss[intVarIdx] = tempLBCons;

    delete[] val_;
}


FP::~FP(){}


void FP::setLpObj()
{
    consTerm = 0;
    for (int i = 0; i < intVars.size(); i++)
    {
        if (isBin[i])
        {
            SCIP_VAR* scipVar = vars[intVars[i]];

            if (fabs(intVarsLB[i] - roundedIntVars[i]) < epInt)
            {
                SCIPchgVarObj(scip, scipVar, 1);
                consTerm -= intVarsLB[i];
            }
            else if (fabs(intVarsUB[i] - roundedIntVars[i]) < epInt)
            {
                SCIPchgVarObj(scip, scipVar, -1);
                consTerm += intVarsUB[i];
            }
            else
            {
                throw std::string("binary variables rounding failure");
            }
        }
    }

}


void FP::Round()
{
    changed = 0;
    double thrs;
    int etrsrnd = trsrnd; //effective randomization used
    //if (etrsrnd == 4)
    //    etrsrnd = (restarts - s1restarts < 5) ? 0 : rng.getInt(4); //if auto, use standard initially, then choose randomly
    if (etrsrnd)
    {
        thrs = rand() / (((double)RAND_MAX) + 1);
        if (etrsrnd > 2)
        {
            /* random quadratic threshold */
            if (thrs <= 0.5)
                thrs = 2 * thrs * (1 - thrs);
            else
                thrs = 2 * thrs * (thrs - 1) + 1;
        }
        else {/* random threshold */
            if (etrsrnd == 2)
                thrs = thrs / 2 + .25;
        }
    }
    else
        thrs = 0.5; /* standard */

    for (int i = 0; i < intVars.size(); i++)
    {
        if (stage > 1 || isBin[i])  // (stage > 1 || isBin[i])
        {
            double r = floor(relaxedIntVars[i] + thrs);
            if (r < intVarsLB[i] || r > intVarsUB[i])
            {
                if (fabs(intVarsUB[i] - intVarsLB[i]) < 0.00001)
                {
                    r = intVarsUB[i];
                }
                else
                {
                    /*std::cout << "i = " << i << std::endl;
                    std::cout << "relaxedIntVars[i]=" << relaxedIntVars[i] << std::endl;
                    std::cout << "thrs=" << thrs << std::endl;
                    std::cout << "intVars[i].getLB()" << intVarsLB[i] << "intVars[i].getUB()" << intVarsUB[i] << std::endl;*/
                    throw std::string("relax out of bound");
                }
            }

            if (roundedIntVars[i] != r)
            {
                roundedIntVars[i] = r;
                changed++;
            }
        }
        else
            roundedIntVars[i] = floor(relaxedIntVars[i] + 0.5);
    }

}


void FP::Solve(int& found, int stage)
{
    if (!SCIPsolve(scip) == SCIP_RETCODE::SCIP_OKAY)
    {
        throw std::string( "scip solve error!");
    }
    if ((SCIPgetStatus(scip) == SCIP_STATUS_SOLLIMIT || SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL) && stage == 3)
    {
        found = 1;
        std::cout << "SCIP IN FP FOUND FEAS SOL!!!" << std::endl;
        SCIP_SOL* sol;
        sol = SCIPgetBestSol(scip);
        SCIP_Real objval = SCIPgetSolOrigObj(scip, sol);
        std::cout << "obj value:" << objval << std::endl;
    }
    relaxedIntVars.clear();
    for (int i = 0; i < intVars.size(); i++)
    {
        relaxedIntVars.push_back(SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars[intVars[i]]));
    }
    SCIPfreeTransform(scip);
}

void FP::Flip()
{
    //temporarily
    if (!changed)
    {
        std::priority_queue<dipair, std::vector<dipair>, dipairCmp> q;
        unsigned int toBeChanged = minChangeNum * (0.5 + rand() / (((double)RAND_MAX) + 1));
        // populate queue with top <toBeChanged> biggest sigma (and sigma>trsld)
        double sigma, min_sigma = trsmin;
        for (int i = 0; i < intVars.size(); i++)
        {
            if (isBin[i])
            {
                sigma = fabs(roundedIntVars[i] - relaxedIntVars[i]);
                if (sigma > min_sigma || min_sigma == 0)
                {
                    q.push(dipair(sigma, i));
                    if (q.size() > toBeChanged)
                    {
                        q.pop();
                        min_sigma = q.top().first;
                    }
                }
            }
        }
        // change the vars in the queue
        while (!q.empty())
        {
            int i = q.top().second;
            if (roundedIntVars[i] > intVarsLB[i] + 0.5)
                roundedIntVars[i]--;
            else if (roundedIntVars[i] < intVarsUB[i] - 0.5)
                roundedIntVars[i]++;
            changed++;
            q.pop();
        }
    }
}


void FP::SetIntegerValue(std::unordered_map<std::string, double> patialSol)
{
    for (int i = 0; i < intVars.size(); i++)
    {
        std::string varName = SCIPvarGetName(vars[intVars[i]]);
        auto iter = patialSol.find(varName);
        if (iter != patialSol.end())
        {
            roundedIntVars[i] = iter->second;
        }
        else
        {
            std::cout << "variable not found" << std::endl;
        }
    }

    for (auto iter = patialSol.begin(); iter != patialSol.end(); iter++)
    {
        for (int i = 0; i < intVars.size(); i++)
        {
            if (iter->first == SCIPvarGetName(vars[intVars[i]]))
            {
                roundedIntVars[i] = iter->second;
            }
        }
    }
}

double FP::getObjVal()
{
    return SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip));
}

int FP::RUNFP()
{
    int stage = 1;
    std::cout << "\nStage " << stage << "..." << std::endl;
    int missedDecr = 0;
    int maxMissedDecr = maxit1wi;
    double minDelta = std::numeric_limits<double>::infinity();
    int minDeltaIter = 0;

    int iter = 0;
    int found = false;

    std::unordered_set<std::pair<std::vector<double>, double>, hash_compare_for_pairs<double>, equalkey_for_pairs<double> > history;

    //stage 1
    while (iter < maxIter1)
    {
        iter++;
        if (!(iter % 20))
            //std::cout << iter << " iterations" << std::endl;

            alpha *= alpha_quot;

        setLpObj();
        time_t cplex_time_start = time(NULL);

        Solve(found, stage);
        time_t cplex_time_end = time(NULL);


        double fodelta = consTerm + getObjVal();

        prevRoundedIntVars = roundedIntVars;

        Round();

        if (fodelta < minDelta)
        {
            //std::cout << "improved distance by " << 100 * fodelta / minDelta << "percent" << std::endl;

            if (fodelta / minDelta < 0.9)
                missedDecr = 0;
            minDelta = fodelta;
            minDeltaIter = iter;

            bestpoint = roundedIntVars;
        }

        else
            missedDecr++;

        if (fodelta < epInt || isCoincident())
        {
            //std::cout << std::endl << "Binary variables didn't change" << std::endl;
            break;
        }

        if (missedDecr > maxMissedDecr)
        {
            //std::cout << std::endl << "Too many iteration without 10% improvement" << std::endl;
            break;
        }

        Flip();

        while (!history.insert(std::pair<std::vector<double>, double>(roundedIntVars, alpha)).second && restarts < 10 * maxIter1)
        {
            Restart();
        }
        while (history.size() > 2000)
            history.erase(history.begin());

        /*if (restarts > AlgPara::maxRestarts)
        {
            std::cout << std::endl << "Too many restarts" << std::endl;
            break;
        }*/


    }

    found = isCoincidentS2();
    if (found)
    {
        if (CheckFeas())
        {
            std::cout << "integer solution found!" << std::endl;
            std::cout << "obj value: " << GetTempObjVal() << std::endl;
        }
    }
    history.clear();

    //stage 2


    //stage3

    if( !found )
    {
        // enumerate
        stage = 3;
        //std::cout << "\nStage " << stage << "..." << std::endl;
        roundedIntVars = bestpoint;
        //std::cout << "Using best point from iter: " << minDeltaIter << std::endl;
        // make problem MIP again
        RecoverIntVars();
        /*cplex.setParam(IloCplex::IntSolLim, 1);
        cplex.setParam(IloCplex::MIPEmphasis, s3me);
        cplex.setParam(IloCplex::EpRHS, cplex.getDefault(IloCplex::EpRHS));*/
        //cplex.setOut(env.out());
        SetNewObjS2();

        SCIPsetIntParam(scip, "limits/solutions", 1);
        Solve(found, stage);

    }

    if (found)
        return 1;
    else
        return 0;
}

bool FP::CheckFeas()
{
    SCIP_SOL* sol;
    SCIPcreateSol(scip, &sol, NULL);

    for (int i = 0; i < vars.size(); i++)
    {
        SCIPsetSolVal(scip, sol, vars[i], SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars[i]));
    }
    unsigned int isFeas = 0;
    SCIP_Bool b = 0;
    SCIP_Bool c = 0;
    SCIPcheckSolOrig(scip, sol, &isFeas, b, c);
    return isFeas == 1 ? true : false;
}

std::map<std::string, double> FP::OutSol()
{
    std::map<std::string, double> solMap;
    for (int i = 0; i < vars.size(); i++)
    {
        solMap.insert(std::pair<std::string, double>(SCIPvarGetName(vars[i]), SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars[i]) ) );
    }
    return solMap;
}

double FP::GetTempObjVal()
{
    double val = 0;
    for (int i = 0; i < vars.size(); i++)
    {
        val += SCIPgetSolVal(scip, SCIPgetBestSol(scip), vars[i]) * oriObj[i] + offSet;
    }
    return val;
}

void FP::Restart()
{
    restarts++;
    for (int i = 0; i < intVars.size(); i++)
    {
        if (isBin[i])
        {
            double r = rand() / (((double)RAND_MAX) + 1) - .47;
            if (r > 0 && prevRoundedIntVars[i] == roundedIntVars[i])
            {
                double sigma = fabs(roundedIntVars[i] - relaxedIntVars[i]);
                if (sigma + r > 0.5)
                {
                    if (fabs(roundedIntVars[i] - intVarsUB[i]) < 0.00001 && roundedIntVars[i] - 1 >= intVarsLB[i] - 0.00001)
                    {
                        roundedIntVars[i]--;
                        changed++;
                    }
                    else if (fabs(roundedIntVars[i] - intVarsLB[i]) < 0.00001 && roundedIntVars[i] + 1 <= intVarsUB[i] + 0.00001)
                    {
                        roundedIntVars[i]++;
                        changed++;
                    }
                }
            }
        }
    }
}

bool FP::isCoincident()
{
    for (int i = 0; i < intVars.size(); i++)
        if (isBin[i] && fabs(relaxedIntVars[i] - roundedIntVars[i]) > epInt)
            return false;
    return true;
}


bool FP::isCoincidentS2()
{
    for (int i = 0; i < intVars.size(); i++)
        if (fabs(relaxedIntVars[i] - roundedIntVars[i]) > epInt)
            return false;
    return true;
}


void FP::RecoverIntVars()
{
    for (int i: intVars)
    {
        auto var = vars[i];
        SCIP_VARTYPE const varType_ = SCIPvarGetType(var);
        SCIP_Bool infeasible;
        if (varType_ == SCIP_VARTYPE_CONTINUOUS)
            SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible);
        else
            std::cout << "variable type error! it should be a continuous variable" << std::endl;
    }
}

void FP::SetNewObjS2()
{
    consTerm = 0;
    double expr_norm = 0.0;
    for (int i = 0; i < intVars.size(); i++)
    {
        if (!auxVars[i])
        {
            if (roundedIntVars[i] == intVarsLB[i])
            {
                SCIPchgVarObj(scip, vars[intVars[i]], 1);
                consTerm -= intVarsLB[i];
            }
            else if (roundedIntVars[i] == intVarsUB[i])
            {
                SCIPchgVarObj(scip, vars[intVars[i]], 1);
                consTerm += intVarsUB[i];
            }
            else
                IntroduceAux(i);
        }
        if (auxVars[i])
        {
            SCIPchgRhsLinear(scip, auxUBConss[i], roundedIntVars[i]);
            SCIPchgLhsLinear(scip, auxLBConss[i], roundedIntVars[i]);
        }
        expr_norm += 1.0;

    }
    //dist = expr;
    if (stage != 3)
    {
        //expr = (1 - alpha) * expr + alpha * sqrt(expr_norm) * origObjExpr / origobj_norm;
        /*expr = (1 - alpha) * expr / sqrt(expr_norm) + alpha * origObjExpr / origobj_norm;
        expr.normalize();*/
    }
}



