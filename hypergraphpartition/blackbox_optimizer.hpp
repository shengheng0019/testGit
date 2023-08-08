#pragma once
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include "common.hpp"
#include "local_search_optimizer.hpp"
#include "objective.hpp"
#include "Partition.hpp"

namespace {

	class CoarseningComparer {
	public:
		CoarseningComparer(const PartitioningParams& params)
			: params_(params) {
		}

		bool operator()(const hgsol::Solution& c1, const hgsol::Solution& c2) const {
			assert(c1.nNodes() == c2.nNodes());
			Index nNodes = c1.nNodes();
			double fac1 = nNodes / (double)c1.nParts();
			double fac2 = nNodes / (double)c2.nParts();
			if (fac1 < params_.minCoarseningFactor || fac2 < params_.minCoarseningFactor)
				return fac1 > fac2;
			if (fac1 > params_.maxCoarseningFactor || fac2 < params_.maxCoarseningFactor)
				return fac1 < fac2;
			double tgt = 0.5 * (params_.maxCoarseningFactor + params_.minCoarseningFactor);
			return abs(fac1 - tgt) < abs(fac2 - tgt);
		}

	private:
		const PartitioningParams& params_;
	};
} // End anonymous namespace


namespace {
	class SolutionHasher {
	public:
		SolutionHasher(const vector<hgsol::Solution>& solutions)
			: solutions_(solutions) {}

		uint64_t operator()(const Index& node) const {
			// FNV hash
			uint64_t magic = 1099511628211llu;
			uint64_t ret = 0;
			for (const hgsol::Solution& solution : solutions_) {
				ret = (ret ^ (uint64_t)solution[node]) * magic;
			}
			return ret;
		}

	private:
		const vector<hgsol::Solution>& solutions_;
	};

	class SolutionComparer {
	public:
		SolutionComparer(const vector<hgsol::Solution>& solutions)
			: solutions_(solutions) {}

		bool operator()(const Index& n1, const Index& n2) const {
			for (const hgsol::Solution& solution : solutions_) {
				if (solution[n1] != solution[n2]) return false;
			}
			return true;
		}

	private:
		const vector<hgsol::Solution>& solutions_;
	};
} // End anonymous namespace


class BlackboxOptimizer {
public:
	static hgsol::Solution run(const MiniPartHypergraph::Hypergraph& hypergraph, const PartitioningParams& params, const MiniPartObj::Objective& objective, const std::vector<hgsol::Solution>& solutions, const std::vector<part::Partition>& initial_parts)
	{
		std::mt19937 rgen(params.seed);
		// Copy because modified in-place
		std::vector<hgsol::Solution> sols = solutions;
		BlackboxOptimizer opt(hypergraph, params, objective, rgen, sols, 0);
		return opt.run(initial_parts);
	}

private:
	BlackboxOptimizer(const MiniPartHypergraph::Hypergraph& hypergraph, const PartitioningParams& params, const MiniPartObj::Objective& objective, std::mt19937& rgen, std::vector<hgsol::Solution>& solutions, Index level) : hypergraph_(hypergraph)
	                                                                                                                                                                                                                           , params_(params)
	                                                                                                                                                                                                                           , objective_(objective)
	                                                                                                                                                                                                                           , rgen_(rgen)
	                                                                                                                                                                                                                           , solutions_(solutions)
	                                                                                                                                                                                                                           , level_(level) {
	}

	hgsol::Solution run(const std::vector<part::Partition>& initial_parts)
	{
		reportStartSearch();
		runInitialPlacement(initial_parts);
		runLocalSearch();
		for (cycle_ = 0; cycle_ < params_.nCycles; ++cycle_) {
			reportStartCycle();
			runVCycle();
			reportEndCycle();
		}
		reportEndSearch();

		return bestSolution();
	}

	hgsol::Solution bestSolution() const
	{
		assert(!solutions_.empty());
		size_t best = 0;
		for (size_t i = 1; i < solutions_.size(); ++i) {
			if (objective_.eval(hypergraph_, solutions_[i]) < objective_.eval(hypergraph_, solutions_[best])) {
				best = i;
			}
		}
		return solutions_[best];
	}

	// this part could be improve , use hype's result to initialize the begin solution
	void runInitialPlacement(const std::vector<part::Partition>& initial_parts)
	{
		for (int i = 0; i < initial_parts.size(); ++i)
		{
			hgsol::Solution solution(hypergraph_.nNodes(), hypergraph_.nParts());
			auto nodes = initial_parts[i].getNodes();
			for (auto node : nodes)
			{
				solution[node] = i;
			}
			solutions_.push_back(solution);
		}
		while ((Index)solutions_.size() < params_.nSolutions) {
			std::uniform_int_distribution<int> partDist(0, hypergraph_.nParts() - 1);
			hgsol::Solution solution(hypergraph_.nNodes(), hypergraph_.nParts());
			for (Index i = 0; i < hypergraph_.nNodes(); ++i) {
				solution[i] = partDist(rgen_);
			}
			solutions_.push_back(solution);
		}
	}
	void runLocalSearch()
	{
		report("Local search");
		for (hgsol::Solution& solution : solutions_) {
			std::unique_ptr<IncrementalObjective> inc = objective_.incremental(hypergraph_, solution);
			LocalSearchOptimizer(*inc, params_, rgen_).run();
			//inc->checkConsistency();
		}
	}
	void runVCycle()
	{
		//checkConsistency();

		if (hypergraph_.nNodes() < params_.minCoarseningNodes * hypergraph_.nParts()) return;
		report("V-cycle step");

		// Pick the best number of solutions for the coarsening
		// If the coarsening is still too large, stop the recursion
		shuffle(solutions_.begin(), solutions_.end(), rgen_);

		vector<hgsol::Solution> coarsenings;
		for (size_t nSols = 1; nSols <= solutions_.size(); ++nSols) {
			coarsenings.emplace_back(computeCoarsening(vector<hgsol::Solution>(solutions_.begin(), solutions_.begin() + nSols)));
		}
		size_t coarseningIndex = min_element(coarsenings.begin(), coarsenings.end(), CoarseningComparer(params_)) - coarsenings.begin();
		hgsol::Solution coarsening = coarsenings[coarseningIndex];
		if (coarsening.nNodes() / (double)coarsening.nParts() < params_.minCoarseningFactor) return;

		MiniPartHypergraph::Hypergraph cHypergraph = hypergraph_.coarsen(coarsening);
		vector<hgsol::Solution> cSolutions;
		for (size_t i = 0; i <= coarseningIndex; ++i) {
			cSolutions.emplace_back(solutions_[i].coarsen(coarsening));
		}
		BlackboxOptimizer nextLevel(cHypergraph, params_, objective_, rgen_, cSolutions, level_ + 1);
		nextLevel.runLocalSearch();
		nextLevel.runVCycle();
		report("Refinement", coarseningIndex + 1);
		for (size_t i = 0; i <= coarseningIndex; ++i) {
			solutions_[i] = cSolutions[i].uncoarsen(coarsening);
			std::unique_ptr<IncrementalObjective> inc = objective_.incremental(hypergraph_, solutions_[i]);
			LocalSearchOptimizer(*inc, params_, rgen_).run();
			//inc->checkConsistency();
		}
		//checkConsistency();
	}

	void report(const std::string& step) const
	{
		report(step, solutions_.size());
	}
	void report(const std::string& step, Index nSols) const
	{
		if (params_.verbosity >= 3) {
			for (int i = 0; i < level_; ++i) std::cout << "  ";
			std::cout << step << ": ";
			std::cout << hypergraph_.nNodes() << " nodes, ";
			std::cout << hypergraph_.nHedges() << " edges, ";
			std::cout << hypergraph_.nPins() << " pins ";
			std::cout << "on " << nSols << " solutions";
			std::cout << std::endl;
		}
	}
	void reportStartCycle() const
	{
		if (params_.verbosity >= 2) {
			std::cout << "Starting V-cycle #" << cycle_ + 1 << std::endl;
		}
	}
	void reportEndCycle() const
	{
		if (params_.verbosity >= 2) {
			hgsol::Solution solution = bestSolution();
			vector<int64_t> obj = objective_.eval(hypergraph_, solution);
			std::cout << "Objectives: ";
			for (size_t i = 0; i < obj.size(); ++i) {
				if (i > 0) std::cout << ", ";
				std::cout << obj[i];
			}
			std::cout << std::endl;
		}
	}
	void reportStartSearch() const{}
	void reportEndSearch() const
	{
		if (params_.verbosity >= 2) {
			std::cout << std::endl;
		}
	}

	void checkConsistency() const
	{
		hypergraph_.checkConsistency();
		for (const hgsol::Solution& solution : solutions_) {
			if (hypergraph_.nNodes() != solution.nNodes())
				throw std::runtime_error("Hypergraph and solutions must have the same number of nodes");
			if (hypergraph_.nParts() != solution.nParts())
				throw std::runtime_error("Hypergraph and solutions must have the same number of partitions");
			solution.checkConsistency();
		}
	}

	static hgsol::Solution computeCoarsening(const std::vector<hgsol::Solution>& solutions)
	{
		assert(solutions.size() >= 1);
		Index nNodes = solutions.front().nNodes();
		std::unordered_map<Index, Index, SolutionHasher, SolutionComparer> coarseningMap(nNodes, SolutionHasher(solutions), SolutionComparer(solutions));
		coarseningMap.reserve(nNodes);

		Index nCoarsenedNodes = 0;
		vector<Index> coarsening(nNodes);
		for (Index node = 0; node < nNodes; ++node) {
			auto p = coarseningMap.emplace(node, nCoarsenedNodes);
			if (p.second) {
				coarsening[node] = nCoarsenedNodes++;
			}
			else {
				coarsening[node] = p.first->second;
			}
		}

		return hgsol::Solution(coarsening);
	}

private:
	const MiniPartHypergraph::Hypergraph& hypergraph_;
	const PartitioningParams& params_;
	const MiniPartObj::Objective& objective_;
	std::mt19937& rgen_;
	std::vector<hgsol::Solution>& solutions_;
	Index level_;
	Index cycle_;
};
