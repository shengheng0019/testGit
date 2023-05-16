#include "Vars.h"

bool operator==(const Vars& var1, const Vars& var2)
{
	if (var1.Id == var2.Id)
		return true;
	return false;
}

Vars& Vars::operator = (const Vars& var)
{
	if (this == &var)
		return *this;
	this->Id = var.Id;
	this->Lb = var.Lb;
	this->Name = var.Name;
	this->Type = var.Type;
	this->Ub = var.Ub;

	return *this;
}