#include "Vars.h"


Vars::Vars(const Vars& var_)
	:Id{ var_.Id }, Lb{ var_.Lb }, Ub{ var_.Ub }, Type{ var_.Type }, Name{var_.Name}{}

bool operator==(const Vars& var1, const Vars& var2)
{
	if (var1.Id == var2.Id)
		return true;
	return false;
}

Vars::Vars(){}

Vars Vars::operator=(const Vars& var1)
{
	this->Id = var1.Id;
	this->Lb = var1.Lb;
	this->Ub = var1.Ub;
	this->Type = var1.Type;
	this->Name = var1.Name;
	return *this;
}

