#include "HYPEHypergraph.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace part {
auto MatrixIntoHypergraph(const std::unordered_map<int, std::vector<uint64_t>> matrix)
    ->part::Hypergraph;
} // namespace part
