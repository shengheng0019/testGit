#include "HYPEHypergraph.hpp"
#include "Parsing.hpp"
#include <vector>



auto part::MatrixIntoHypergraph(const std::unordered_map<int, std::vector<uint64_t>> matrix)
    -> part::Hypergraph
{

    part::Hypergraph ret_graph{};
    for (auto element :matrix )
    {
        std::vector<int64_t> temp_vec(element.second.begin(),element.second.end());
        ret_graph.addNodeList(element.first, temp_vec);
    }
    return ret_graph;
}

