#include "SetVarsBounds.h"

int ToIntUB(double val);

int ToIntLB(double val);


std::pair<std::unordered_map<int, Vars>, std::vector<Constraints>> getAllVarsandCons(std::map<int, General_Block>& blocks)
{
	std::unordered_map<int, Vars> allVars;
	std::vector<Constraints> allCons;

	auto num_of_blocks = blocks.size();
	for (int i = 0; i < num_of_blocks; ++i)
	{
		auto tempBlock = blocks[i];
		allCons.insert(allCons.end(), tempBlock.bCons.begin(), tempBlock.bCons.end());
		for (const auto& tempvar : tempBlock.bVars)
		{
			allVars[tempvar.Id] = tempvar;
		}
	}

	return { allVars,allCons };
}


void SetVarsBounds(std::unordered_map<int, Vars>& vars_, const std::vector<Constraints>& conss_)
{
	std::cout << "setting variables bounds..." << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	
	std::unordered_map<int, Vars> vars;
	vars.insert(vars_.begin(), vars_.end());
	//std::vector<Constraints> conss(conss_);
	std::vector<Constraints> conss;
	for (auto consIter = conss_.begin(); consIter != conss_.end(); consIter++)
	{
		for (auto varIter = consIter->exprDic.begin(); varIter != consIter->exprDic.end(); varIter++)
		{
			if (varIter->first.Lb < -INF || varIter->first.Ub > INF)
			{
				conss.push_back(*consIter);
				break;
			}
		}
	}
	//convert to <= or ==
	for (int i = 0; i < conss.size(); i++)
	{
		if (conss[i].Prop != PROP::geq)
		{
			continue;
		}
		
		for (auto iter = conss[i].exprDic.begin(); iter != conss[i].exprDic.end(); iter++)
		{
			iter->second = -iter->second;
		}
		conss[i].rhs = -conss[i].rhs;
		conss[i].Prop = PROP::leq;
	}

	queue<int> varsToCheck;
	unordered_set<int> ubVars;
	unordered_set<int> lbVars;
	unordered_set<int> varsInQueue;

	int oriUbCounter = 0;
	int oriLbCounter = 0;

	//Init: insert vars with bounds in the set and into the queue
	for (auto iter = vars.begin(); iter != vars.end(); iter++)
	{
		if (iter->second.Ub < INF)		//check if ub is set by the original model
		{
			ubVars.insert(iter->first);
			oriUbCounter++;
		}
		if (iter->second.Lb > -INF)		//check if lb is set by the original model
		{
			lbVars.insert(iter->first);
			oriLbCounter++;
		}
		if (iter->second.Ub < INF || iter->second.Lb > -INF)
		{
			varsToCheck.push(iter->first);
			varsInQueue.insert(iter->first);
		}
	}
	//main alg: 
	while (varsToCheck.size() > 0)
	{
		int varId = varsToCheck.front();

		//to find conss with the specific var and check
		for (auto iter = conss.begin(); iter != conss.end(); )
		{
			if (iter->exprDic.find(vars[varId]) == iter->exprDic.end())
			{
				iter++;
				continue;
			}
			
			// for a given cons...
			double rhs = iter->rhs;
			if (iter->Prop == PROP::leq)	//the type of cons is <=
			{
				/*std::unordered_map<Vars, double>::iterator varIter;
				for (varIter = iter->exprDic.begin(); varIter != iter->exprDic.end(); varIter++)
				{
					if (varIter->second < 0 || ubVars.find(varIter->first.Id) != ubVars.end())
					{
						;
					}
					else if (varIter->second > 0 || lbVars.find(varIter->first.Id) != lbVars.end())
					{
						;
					}
					else
						break;
				}*/
				if (true)
				{
					// \sigma{a_i * x_i} + \sigma{a_j * x_j} <= b	a_i>0, a_j<0
					// each x_i has lb and each x_j has ub
					// then ub of each x_i and lb of each x_j can be inferred
					for (auto varIter = iter->exprDic.begin(); varIter != iter->exprDic.end(); varIter++)
					{
						//not necessary: 
						if (varIter->second > 0 && ubVars.find(varIter->first.Id) != ubVars.end())
							continue;
						if (varIter->second < 0 && lbVars.find(varIter->first.Id) != lbVars.end())	//already with lb, then not necessary to cal lb again
							continue;

						//rules: 
						//x_i <= (b - \sigma{a_i * lb(x_i)} - \sigma{a_j * ub(x_j)}) / a_i
						//x_j >= (b - \sigma{a_i * ub(x_i)} - \sigma{a_i * ub(x_i)}) / a_j
						
						double rhsSum = rhs;
						unordered_map<Vars, double>::iterator varIter_;
						for (varIter_ = iter->exprDic.begin(); varIter_ != iter->exprDic.end(); varIter_++)
						{
							if (varIter_ == varIter)
								continue;
							if (varIter_->second > 0)
							{
								if (lbVars.find(varIter_->first.Id) == lbVars.end())
									break;
								rhsSum -= varIter_->second * vars[varIter_->first.Id].Lb;
							}
							else
							{
								if (ubVars.find(varIter_->first.Id) == ubVars.end())
									break;
								rhsSum -= varIter_->second * vars[varIter_->first.Id].Ub;
							}
						}
						if (varIter_ != iter->exprDic.end())
							break;
						if (varIter->second > 0)
						{
							vars[varIter->first.Id].Ub = ToIntUB(rhsSum / varIter->second);	//add the var ub
							ubVars.insert(varIter->first.Id);				//to record this var has a ub
							/*std::cout << "varId: " << varIter->first.Name << " upper Bound: " << vars[varIter->first.Id].Ub
								<< " from cons: " << iter->Name << endl;*/
							if (varsInQueue.find(varIter->first.Id) == varsInQueue.end())
							{
								varsToCheck.push(varIter->first.Id);			//waiting to be checked to infer the corresponding vars' bounds
								varsInQueue.insert(varIter->first.Id);
							}
						}
						else
						{
							vars[varIter->first.Id].Lb = ToIntLB(rhsSum / varIter->second);	//add the var ub
							lbVars.insert(varIter->first.Id);				//to record this var has a ub
							/*std::cout << "varId: " << varIter->first.Name << " lower Bound: " << vars[varIter->first.Id].Lb
								<< " from cons: " << iter->Name << endl;*/
							if (varsInQueue.find(varIter->first.Id) == varsInQueue.end())
							{
								varsToCheck.push(varIter->first.Id);			//waiting to be checked to infer the corresponding vars' bounds
								varsInQueue.insert(varIter->first.Id);
							}
						}
					}
				}
			}
			else		// ==
			{
				if (true)
				{
					// \sigma{a_i * x_i} + \sigma{a_j * x_j} == b	a_i>0, a_j<0
					// each x_i has lb and each x_j has ub
					// then ub of x_i and lb of x_j can be inferred

					//for each var in a cons
					for (auto varIter = iter->exprDic.begin(); varIter != iter->exprDic.end(); varIter++)
					{
						//not necessary: 
						if (varIter->second > 0 && ubVars.find(varIter->first.Id) != ubVars.end())
							continue;
						if (varIter->second < 0 && lbVars.find(varIter->first.Id) != lbVars.end())	//already with lb, then not necessary to cal lb again
							continue;

						//rules: 
						//x_i <= (b - \sigma{a_i * lb(x_i)} - \sigma{a_j * ub(x_j)}) / a_i
						//x_j >= (b - \sigma{a_i * ub(x_i)} - \sigma{a_i * ub(x_i)}) / a_j

						double rhsSum = rhs;
						unordered_map<Vars, double>::iterator varIter_;
						for (varIter_ = iter->exprDic.begin(); varIter_ != iter->exprDic.end(); varIter_++)
						{
							if (varIter_ == varIter)
								continue;
							if (varIter_->second > 0)
							{
								if (lbVars.find(varIter_->first.Id) == lbVars.end())
									break;
								rhsSum -= varIter_->second * vars[varIter_->first.Id].Lb;
							}
							else
							{
								if (ubVars.find(varIter_->first.Id) == ubVars.end())
									break;
								rhsSum -= varIter_->second * vars[varIter_->first.Id].Ub;
							}
						}
						if (varIter_ != iter->exprDic.end())
							break;
						if (varIter->second > 0)
						{
							vars[varIter->first.Id].Ub = ToIntUB(rhsSum / varIter->second);	//add the var ub
							ubVars.insert(varIter->first.Id);				//to record this var has a ub
							/*std::cout << "varId: " << varIter->first.Name << " upper Bound: " << vars[varIter->first.Id].Ub 
								<< " from cons: " << iter->Name << endl;*/
							if (varsInQueue.find(varIter->first.Id) == varsInQueue.end())
							{
								varsToCheck.push(varIter->first.Id);			//waiting to be checked to infer the corresponding vars' bounds
								varsInQueue.insert(varIter->first.Id);
							}
						}
						else
						{
							vars[varIter->first.Id].Lb = ToIntLB(rhsSum / varIter->second);	//add the var ub
							lbVars.insert(varIter->first.Id);				//to record this var has a ub
							/*std::cout << "varId: " << varIter->first.Name << " lower Bound: " << vars[varIter->first.Id].Lb
								<< " from cons: " << iter->Name << endl;*/
							if (varsInQueue.find(varIter->first.Id) == varsInQueue.end())
							{
								varsToCheck.push(varIter->first.Id);			//waiting to be checked to infer the corresponding vars' bounds
								varsInQueue.insert(varIter->first.Id);
							}
						}
					}
				}
				
				if (true)
				{
					// \sigma{a_i * x_i} + \sigma{a_j * x_j} == b	a_i>0, a_j<0
					// each x_j has lb and each x_i has ub
					// then ub of x_j and lb of x_i can be inferred
					for (auto varIter = iter->exprDic.begin(); varIter != iter->exprDic.end(); varIter++)
					{
						if (varIter->second > 0 && lbVars.find(varIter->first.Id) != lbVars.end())
							continue;
						if (varIter->second < 0 && ubVars.find(varIter->first.Id) != ubVars.end())
							continue;
						//x_i >= (b - \sigma{a_j * lb(x_j)} - \sigma{a_i * ub(x_i)}) / a_i
						double rhsSum = rhs;
						unordered_map<Vars, double>::iterator varIter_;
						for (varIter_ = iter->exprDic.begin(); varIter_ != iter->exprDic.end(); varIter_++)
						{
							if (varIter_ == varIter)
								continue;
							if (varIter_->second > 0)
							{
								if (ubVars.find(varIter_->first.Id) == ubVars.end())
									break;
								rhsSum -= varIter_->second * vars[varIter_->first.Id].Ub;
							}
							else
							{
								if (lbVars.find(varIter_->first.Id) == lbVars.end())
									break;
								rhsSum -= varIter_->second * vars[varIter_->first.Id].Lb;
							}
						}
						if (varIter_ != iter->exprDic.end())
							break;
						if (varIter->second < 0)
						{
							vars[varIter->first.Id].Ub = ToIntUB(rhsSum / varIter->second);	//add the var ub
							ubVars.insert(varIter->first.Id);				//to record this var has a ub
							/*std::cout << "varId: " << varIter->first.Name << " upper Bound: " << vars[varIter->first.Id].Ub
								<< " from cons: " << iter->Name << endl;*/
							if (varsInQueue.find(varIter->first.Id) == varsInQueue.end())
							{
								varsToCheck.push(varIter->first.Id);			//waiting to be checked to infer the corresponding vars' bounds
								varsInQueue.insert(varIter->first.Id);
							}
						}
						else
						{
							vars[varIter->first.Id].Lb = ToIntLB(rhsSum / varIter->second);	//add the var ub
							lbVars.insert(varIter->first.Id);				//to record this var has a ub
							/*std::cout << "varId: " << varIter->first.Name << " lower Bound: " << vars[varIter->first.Id].Lb
								<< " from cons: " << iter->Name << endl;*/
							if (varsInQueue.find(varIter->first.Id) == varsInQueue.end())
							{
								varsToCheck.push(varIter->first.Id);			//waiting to be checked to infer the corresponding vars' bounds
								varsInQueue.insert(varIter->first.Id);
							}
						}
					}
				}

			}
			//eliminate the constraints with all variables having boundings
			auto vIter = iter->exprDic.begin();
			for (; vIter != iter->exprDic.end(); vIter++)
			{
				if (vIter->first.Lb < -INF || vIter->first.Ub > INF)
				{
					iter++;
					break;
				}
			}
			if (vIter == iter->exprDic.end())
				iter = conss.erase(iter);
		}



		varsToCheck.pop();
		varsInQueue.erase(varId);
	}
	
	//count the result of the alg
	
	for (auto iter = vars_.begin(); iter != vars_.end(); iter++)
	{
		if (vars[iter->first].Lb > -INF)
			iter->second.Lb = vars[iter->first].Lb;
		else
			cout << "var: " << iter->second.Name << " has no LB";

		if (vars[iter->first].Ub < INF)
			iter->second.Ub = vars[iter->first].Ub;
		else
			cout << "var: " << iter->second.Name << " has no UB";
	}
	
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "_______________________________________________" << endl;
	std::cout << "finish setting variables bounds with " << duration << " ms" << std::endl;
	std::cout << "original upper bounds: " << oriUbCounter << " inferred bounds:" << ubVars.size() - oriUbCounter << std::endl;
	std::cout << "original lower bounds: " << oriLbCounter << " inferred bounds:" << lbVars.size() - oriLbCounter << std::endl;
	std::cout << "-----------------------------------------------" << endl;
}

int ToIntUB(double val)
{
	return (int)floor(val) + 1;
}

int ToIntLB(double val)
{
	return (int)floor(val);
}

