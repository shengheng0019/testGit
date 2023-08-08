
#pragma once

#include "hypergraph.hpp"
#include "incremental_objective.hpp"
#include "partitioning_params.hpp"
#include "move.hpp"

#include <random>
#include <iosfwd>
#include <cassert>
#include <memory>



class LocalSearchOptimizer {
public:
	LocalSearchOptimizer(IncrementalObjective& inc, const PartitioningParams& params, std::mt19937& rgen) : inc_(inc)
		, params_(params)
		, rgen_(rgen) {
	}
	void run()
	{
		assert(inc_.nNodes() > 0);
		init();
		while (totalBudget() > 0) {
			doMove();
		}
	}

private:
	void init()
	{
		moves_.clear();
		double targetCount = params_.movesPerElement * params_.nNodes * (params_.nParts - 1);
		moves_.emplace_back(std::make_unique<hyperMove::VertexMoveRandomBlock>(0.1 * targetCount));
		moves_.emplace_back(std::make_unique<hyperMove::VertexMoveBestBlock>(0.4 * targetCount));
		moves_.emplace_back(std::make_unique<hyperMove::VertexPassRandomBlock>(0.0 * targetCount));
		moves_.emplace_back(std::make_unique<hyperMove::VertexPassBestBlock>(0.0 * targetCount));
		moves_.emplace_back(std::make_unique<hyperMove::VertexSwap>(0.1 * targetCount));
		moves_.emplace_back(std::make_unique<hyperMove::EdgeMoveRandomBlock>(0.1 * targetCount));
		moves_.emplace_back(std::make_unique<hyperMove::VertexAbsorptionPass>(0.3 * targetCount));
	}

	void doMove()
	{
		// Pick a move based on the remaining budget
		std::uniform_int_distribution<size_t> dist(0, totalBudget() - 1);
		size_t roll = dist(rgen_);
		size_t tot = 0;
		for (std::unique_ptr<hyperMove::Move>& mv : moves_) {
			if (mv->budget_ > 0) tot += mv->budget_;
			if (tot > roll) {
				mv->run(inc_, rgen_);
				return;
			}
		}
		assert(false);
	}
	std::int64_t totalBudget() const
	{
		int64_t ret = 0;
		for (const std::unique_ptr<hyperMove::Move>& mv : moves_) {
			if (mv->budget_ > 0) ret += mv->budget_;
		}
		return ret;
	}

private:
	IncrementalObjective& inc_;
	const PartitioningParams& params_;
	std::mt19937& rgen_;
	std::vector<std::unique_ptr<hyperMove::Move> > moves_;
};



