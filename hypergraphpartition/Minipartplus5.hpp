// Copyright (C) 2019 Gabriel Gouvine - All Rights Reserved

#include "hypergraph.hpp"
#include "partitioning_params.hpp"
#include "blackbox_optimizer.hpp"
#include "objective.hpp"
#include "Partition.hpp"

#include <iostream>
#include <iomanip>
#include <unordered_map>





class VM
{
public:
    Index blocks;
    double imbalance = 30.0;
    ObjectiveType objective = ObjectiveType::Cut;
    Index verbosity = 1;
    size_t seed = 0;
    Index pool_size = 2;  //this param can be change , That means there are several solutions
    Index v_cycles = 0;
    double min_c_factor = 1.2;
    double max_c_factor = 3.0;
    Index min_c_nodes = 8;
    double move_ratio = 8.0;
    std::unordered_map<int, std::vector<uint64_t>> matrix;
    std::unordered_map<int, std::vector<uint64_t>> dual_matrix;
    //std::set<uint64_t> num_of_vars;
};


VM parseArguments(Index blocks, std::unordered_map<int, std::vector<uint64_t>> matrix, std::unordered_map<int, std::vector<uint64_t>> matrix_dual)
{
    VM vm;
    vm.blocks = blocks;
    vm.matrix = matrix;
    vm.dual_matrix = matrix_dual;
    return vm;
}


MiniPartHypergraph::Hypergraph readHypergraph(const VM& vm)
{
    auto matrix = vm.matrix;

    Index nNodes = vm.dual_matrix.size();
    Index nHedges = matrix.size();
    MiniPartHypergraph::Hypergraph hg = MiniPartHypergraph::Hypergraph::readFile(matrix, nNodes, nHedges);
    //hg.checkConsistency();
    hg.mergeParallelHedges();
    hg.setupBlocks(
        vm.blocks,
        vm.imbalance / 100.0
    );
    return hg;
}

PartitioningParams readParams(const VM& vm, const MiniPartHypergraph::Hypergraph& hg) {
    PartitioningParams ParParms;
    ParParms.verbosity = vm.verbosity;
    ParParms.seed = vm.seed;
    ParParms.objective = vm.objective;
    ParParms.nSolutions = vm.pool_size;
    ParParms.nCycles = vm.v_cycles;
    ParParms.minCoarseningFactor = vm.min_c_factor;
    ParParms.maxCoarseningFactor = vm.max_c_factor;
    ParParms.minCoarseningNodes = vm.min_c_nodes;
    ParParms.movesPerElement = vm.move_ratio;
    ParParms.nNodes = hg.nNodes();
    ParParms.nHedges = hg.nHedges();
    ParParms.nPins = hg.nPins();
    ParParms.nParts = hg.nParts();
    return ParParms;
}



vector<hgsol::Solution> readInitialSolutions(const VM& vm, const MiniPartHypergraph::Hypergraph& hg) {
    vector<hgsol::Solution> solutions;
    for (hgsol::Solution& sol : solutions) {
        if (sol.nParts() > hg.nParts()) {
	        std::cerr << "The initial solution has more blocks than specified" << std::endl;
            exit(1);
        }
        sol.resizeParts(hg.nParts());
    }
    return solutions;
}


void report(const PartitioningParams&, const MiniPartHypergraph::Hypergraph& hg) {
	std::cout << "Nodes: " << hg.nNodes() << std::endl;
	std::cout << "Edges: " << hg.nHedges() << std::endl;
	std::cout << "Pins: " << hg.nPins() << std::endl;
	std::cout << "Parts: " << hg.nParts() << std::endl;
	std::cout << std::endl;
}

void reportMainMetrics(const PartitioningParams& params, const MiniPartHypergraph::Hypergraph& hg, const hgsol::Solution& sol) {
	std::cout << "Cut: " << hg.metricsCut(sol) << std::endl;
    if (hg.nParts() > 2) {
	    std::cout << "Connectivity: " << hg.metricsConnectivity(sol) << std::endl;
	    std::cout << "Maximum degree: " << hg.metricsMaxDegree(sol) << std::endl;
    }
    if (params.isDaisyChainObj()) {
        if (hg.nParts() > 2) {
	        std::cout << "Daisy-chain distance: " << hg.metricsDaisyChainDistance(sol) << std::endl;
	        std::cout << "Daisy-chain maximum degree: " << hg.metricsDaisyChainMaxDegree(sol) << std::endl;
        }
    }
    if (params.isRatioObj()) {
	    std::cout << "Ratio cut: " << hg.metricsRatioCut(sol) << std::endl;
        if (hg.nParts() > 2) {
	        std::cout << "Ratio connectivity: " << hg.metricsRatioConnectivity(sol) << std::endl;
	        std::cout << "Ratio maximum degree: " << hg.metricsRatioMaxDegree(sol) << std::endl;
        }
	    std::cout << "Ratio penalty: " << 100.0 * (hg.metricsRatioPenalty(sol) - 1.0) << "%" << std::endl;
    }
	std::cout << std::endl;
}

void reportPartitionUsage(const PartitioningParams& params, const MiniPartHypergraph::Hypergraph& hg, const hgsol::Solution& sol) {
    std::vector<Index> usage = hg.metricsPartitionUsage(sol);
    if (params.isRatioObj()) {
        Index totNodeWeight = hg.totalNodeWeight();
        std::cout << "Partition usage:" << std::endl;
        for (Index p = 0; p < hg.nParts(); ++p) {
	        std::cout << "\tPart#" << p << "  \t";
	        std::cout << usage[p] << "\t";
	        std::cout << "(" << 100.0 * usage[p] / totNodeWeight << "%)";
	        std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    else {
	    std::cout << "Partition usage:" << std::endl;
        for (Index p = 0; p < hg.nParts(); ++p) {
	        std::cout << "\tPart#" << p << "  \t";
	        std::cout << usage[p] << "\t/ " << hg.partWeight(p) << "\t";
	        std::cout << "(" << 100.0 * usage[p] / hg.partWeight(p) << "%)\t";
            if (usage[p] > hg.partWeight(p)) std::cout << "(overflow)";
	        std::cout << std::endl;
        }
	    std::cout << std::endl;
    }
}

void reportPartitionDegree(const PartitioningParams& params, const MiniPartHypergraph::Hypergraph& hg, const hgsol::Solution& sol) {
    if (hg.nParts() <= 2) return;

    std::vector<Index> degree = hg.metricsPartitionDegree(sol);
    std::cout << "Partition degrees:" << std::endl;
    for (Index p = 0; p < hg.nParts(); ++p) {
	    std::cout << "\tPart#" << p << "  \t";
	    std::cout << degree[p] << std::endl;
    }
    if (params.isDaisyChainObj()) {
	    std::cout << std::endl;
	    std::cout << "Daisy-chain partition degrees: " << std::endl;
        degree = hg.metricsPartitionDaisyChainDegree(sol);
        for (Index p = 0; p < hg.nParts(); ++p) {
	        std::cout << "\tPart#" << p << "  \t";
	        std::cout << degree[p] << std::endl;
        }
    }
}

void report(const PartitioningParams& params, const MiniPartHypergraph::Hypergraph& hg, const hgsol::Solution& sol) {
    reportMainMetrics(params, hg, sol);
    reportPartitionUsage(params, hg, sol);
    reportPartitionDegree(params, hg, sol);
}

std::unique_ptr<MiniPartObj::Objective> readObjective(const VM& vm) {
    switch (vm.objective) {
    case ObjectiveType::Cut:
        return make_unique<MiniPartObj::CutObjective>();
    case ObjectiveType::Soed:
        return make_unique<MiniPartObj::SoedObjective>();
    case ObjectiveType::MaxDegree:
        return make_unique<MiniPartObj::MaxDegreeObjective>();
    case ObjectiveType::DaisyChainDistance:
        return make_unique<MiniPartObj::DaisyChainDistanceObjective>();
    case ObjectiveType::DaisyChainMaxDegree:
        return make_unique<MiniPartObj::DaisyChainMaxDegreeObjective>();
    case ObjectiveType::RatioCut:
        return make_unique<MiniPartObj::RatioCutObjective>();
    case ObjectiveType::RatioSoed:
        return make_unique<MiniPartObj::RatioSoedObjective>();
    case ObjectiveType::RatioMaxDegree:
        return make_unique<MiniPartObj::RatioMaxDegreeObjective>();
    default:
	    std::cout << "Objective type is not supported." << std::endl;
        exit(1);
    }
}

void initialReport(const MiniPartHypergraph::Hypergraph& hg, const PartitioningParams& params, const vector<hgsol::Solution>& initialSolutions) {
    if (params.verbosity >= 1) {
        report(params, hg);
        if (initialSolutions.size() > 0) {
	        std::cout << "Initial solution:" << std::endl;
        }
        for (const hgsol::Solution& sol : initialSolutions) {
            report(params, hg, sol);
        }
    }
}

void finalReport(const MiniPartHypergraph::Hypergraph& hg, const PartitioningParams& params, const vector<hgsol::Solution>& finalSolutions) {
    if (params.verbosity >= 1) {
        for (const hgsol::Solution& sol : finalSolutions) {
            report(params, hg, sol);
        }
    }
}


std::tuple<VM,MiniPartHypergraph::Hypergraph> MiniPartStage1(std::unordered_map<int, std::vector<uint64_t>> matrix, std::unordered_map<int, std::vector<uint64_t>> matrix_dual, Index k)
{

    auto vm = parseArguments(k, matrix, matrix_dual);

    MiniPartHypergraph::Hypergraph hg = readHypergraph(vm);


    return { vm,hg};
}

std::unordered_map<int, std::set<uint64_t>> MiniPartStage2(const VM& vm, MiniPartHypergraph::Hypergraph& hg,Index k,const std::vector<part::Partition>& initial_parts)
{
    PartitioningParams params = readParams(vm, hg);
    std::unique_ptr<MiniPartObj::Objective> objectivePtr = readObjective(vm);
    vector<hgsol::Solution> initialSolutions = readInitialSolutions(vm, hg);
    hgsol::Solution solution = BlackboxOptimizer::run(hg, params, *objectivePtr, initialSolutions,initial_parts);

    auto parts = solution.ReParts(k);
    std::cout << "partition success ! " << std::endl;
    return parts;
}