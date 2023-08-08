#pragma once

#include <vector>
#include <string>
#include "json.hpp"

namespace writeJson
{
    class Soljson
    {
    public:
        double obj;
        std::map<std::string,double> solution;    //var:name-value
        std::string source; //source == "L" means from lagrangean ; == "F" means from feasibilityjump
        double duration;
        void to_json(nlohmann::json& j)
        {
            j = nlohmann::json{{"source", this->source}, {"obj", this->obj}, {"solution", this->solution},{"duration", this->duration}};
        }
    };
}
