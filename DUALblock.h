#pragma once
#include <vector>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

class DUALBlock
{
public:
    std::vector<SCIP_VAR*> Independent_vars;
    std::vector<SCIP_VAR*> Coupling_vars;
    std::vector<SCIP_CONS*> Cons;
};
