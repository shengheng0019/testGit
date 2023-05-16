#include "Solution.h"

Solution Solution::operator=(const Solution& sol)
{
	this->objValue = sol.objValue;
	this->varValue.insert(sol.varValue.begin(), sol.varValue.end());
	return *this;
}

Solution::Solution(const Solution& sol)
	: objValue{ sol.objValue }
{
	varValue.insert(sol.varValue.begin(), sol.varValue.end());
}

Solution::Solution()
{

}
//Solution::Solution(int solSize)
//{
//	this->solSize = solSize;
//	varValue = new double[solSize];
//	varId = new int[solSize];
//}

