#include "Constraints.h"

Constraints::Constraints(const Constraints& cons_)
	:Prop{ cons_.Prop }, rhs{cons_.rhs}, Name{cons_.Name}
{
	exprDic.insert(cons_.exprDic.begin(), cons_.exprDic.end());
	idDic.insert(cons_.idDic.begin(), cons_.idDic.end());
}
Constraints& Constraints::operator = (const Constraints& cons)
{
	if (this == &cons)
		return *this;

	//this->Coefficient = cons.Coefficient;
	//this->ContainedVars = cons.ContainedVars;
	this->exprDic = cons.exprDic;
	this->idDic = cons.idDic;
	//this->Id = cons.Id;
	//this->Lb = cons.Lb;
	//this->Ub = cons.Ub;
	this->Prop = cons.Prop;
	this->rhs = cons.rhs;
	this->Name = cons.Name;

	return *this;

}

size_t std::hash<Vars>::operator()(const Vars& vars) const noexcept
{
	return hash<int>()(vars.Id);
}

//Constraints::~Constraints()
//{
//	delete[] ContainedVars;
//}

Constraints::Constraints(){}

