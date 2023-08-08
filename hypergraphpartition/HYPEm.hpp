#pragma once

#include "HYPEHypergraph.hpp"
#include "Parsing.hpp"
#include "Partition.hpp"
#include "Partitioning.hpp"
#include "SSet.hpp"
#include <iostream>
#include <string>
#include <set>




std::vector<part::Partition> RunHYPE(
            const std::unordered_map<int, std::vector<uint64_t>>& matrix,
			std::size_t partitions = 6)
{
    std::size_t sset_size = 50;
    std::size_t nh_expand_candidates = 2;
    bool raw = false;
    bool optput = false;
    double percent_of_edges_ignored = 0;
    part::NodeSelectionMode node_select_mode = part::NodeSelectionMode::NextBest;
    std::uint32_t seed = 0;
    part::NodeHeuristicMode heuristic_calc_method = part::NodeHeuristicMode::Exact;
    //("raw,r","output raw numbers to make it easier to redirect output into files")
    //("output,o","write the final partitions into files")
    //("input,i","input hypergraph file")
    //("format,f","specify the input format of the hypergraph file")
    //("partitions,p","number of partitions")
    //("sset-size,s","maximum size of the secondary set")
    //("nh-expand-candidates,n","number of candidates explored during neighbouhood expantion")
    //("percent-of-edges-ignored,e","how many percent of the biggest edges will be removed")
    //("node-select-mode,m","specifies how the a node will be choosen to when sset is empty")
    //("seed,x","seed used to initialze random number generators if used")
    //("heuristic-calc-method,c","Switch to choose between exact and cached calculation for the node heuristic");
    if (matrix.size()>500000)
    {
        heuristic_calc_method = part::NodeHeuristicMode::Cached;
    }
    if (!raw) 
    {
        std::cout << "----------------------------------------------------------------------------\n"
            << "into "
            << partitions
            << " partitions\n"
            << "max secondary set size: "
            << sset_size
            << "\n"
            << "while the secondary set expantion, the biggest "
            << percent_of_edges_ignored
            << "% of edges will be ignored\n"
            << "----------------------------------------------------------------------------\n";

        std::cout << "parsing graph ...\n";
    }

    part::Hypergraph::setSeed(seed);
    auto begin = std::chrono::steady_clock::now();
    part::Hypergraph graph;
	graph = part::MatrixIntoHypergraph(matrix);

    auto number_of_nodes = graph.getVertices().size();
    auto number_of_edges = graph.getEdges().size();

    auto end = std::chrono::steady_clock::now();
    auto parsing_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
        .count();

    if (!raw) 
    {
        std::cout << "graph parsed in "
            << parsing_time
            << " milliseconds\n"
            << "#Nodes:\t"
            << number_of_nodes
            << "\n"
            << "#HyperEdges:\t"
            << number_of_edges
            << "\n"
            << "----------------------------------------------------------------------------\n";

        std::cout << "partitioning graph\n";
    }


    begin = std::chrono::steady_clock::now();
    auto parts = part::partitionGraph(std::move(graph),
											        partitions,
											        sset_size,
											        nh_expand_candidates,
											        percent_of_edges_ignored,
											        heuristic_calc_method,
											        node_select_mode);
    end = std::chrono::steady_clock::now();

    auto partitioning_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    return parts;
}
