#include"heuristicTools.h"

double chgInt(double inaccuracy)
{
	double ret = 0;
	ret = std::ceil(inaccuracy) - inaccuracy < EPSFORINT ? std::ceil(inaccuracy) : inaccuracy;
	ret = inaccuracy - std::floor(inaccuracy) < EPSFORINT ? std::floor(inaccuracy) : ret;

	return ret;
}