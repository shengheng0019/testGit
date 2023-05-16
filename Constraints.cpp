#include"Constraints.h"

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