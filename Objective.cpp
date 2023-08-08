#include "Objective.h"

Objective::Objective(const Objective& obj_)
	:direction{ obj_.direction }
{
	coef.insert(obj_.coef.begin(), obj_.coef.end());
}

//Objective::~Objective()
//{
//	delete[] containedVars;
//}

Objective::Objective(){}

