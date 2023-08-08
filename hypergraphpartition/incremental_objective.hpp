// Copyright (C) 2019 Gabriel Gouvine - All Rights Reserved

#pragma once
#include "hypergraph.hpp"

using std::vector;

namespace {
	vector<Index> computePartitionDemands(const MiniPartHypergraph::Hypergraph& hypergraph, const hgsol::Solution& solution) {
		vector<Index> ret(hypergraph.nParts(), 0);
		for (Index node = 0; node < hypergraph.nNodes(); ++node) {
			ret[solution[node]] += hypergraph.nodeWeight(node);
		}
		return ret;
	}

	vector<vector<Index> > computeHedgeNbPinsPerPartition(const MiniPartHypergraph::Hypergraph& hypergraph, const hgsol::Solution& solution) {
		vector<vector<Index> > ret(hypergraph.nHedges());
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			vector<Index> cnt(hypergraph.nParts());
			for (Index node : hypergraph.hedgeNodes(hedge)) {
				++cnt[solution[node]];
			}
			ret[hedge] = cnt;
		}
		return ret;
	}

	vector<Index> computeHedgeDegrees(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<vector<Index> >& hedgeNbPinsPerPartition) {
		vector<Index> ret(hypergraph.nHedges());
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			Index degree = 0;
			for (Index cnt : hedgeNbPinsPerPartition[hedge]) {
				if (cnt != 0) ++degree;
			}
			ret[hedge] = degree;
		}
		return ret;
	}

	vector<Index> computePartitionDegrees(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<Index>& hedgeDegrees, const vector<vector<Index> >& hedgeNbPinsPerPartition) {
		vector<Index> ret(hypergraph.nParts(), 0);
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			if (hedgeDegrees[hedge] > 1) {
				for (Index p = 0; p < hypergraph.nParts(); ++p) {
					if (hedgeNbPinsPerPartition[hedge][p] != 0) {
						ret[p] += hypergraph.hedgeWeight(hedge);
					}
				}
			}
		}
		return ret;
	}

	vector<std::pair<Index, Index> > computeDaisyChainMinMax(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<vector<Index> >& hedgeNbPinsPerPartition) {
		vector<std::pair<Index, Index> > ret(hypergraph.nHedges());
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			Index minPart = hypergraph.nParts() - 1;
			Index maxPart = 0;
			for (Index p = 0; p < hypergraph.nParts(); ++p) {
				if (hedgeNbPinsPerPartition[hedge][p] != 0) {
					minPart = std::min(minPart, p);
					maxPart = std::max(maxPart, p);
				}
			}
			ret[hedge] = std::make_pair(minPart, maxPart);
		}
		return ret;
	}

	Index computeDaisyChainDistance(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<std::pair<Index, Index> >& hedgeMinMax) {
		Index distance = 0;
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			Index minPart = hedgeMinMax[hedge].first;
			Index maxPart = hedgeMinMax[hedge].second;
			distance += hypergraph.hedgeWeight(hedge) * (maxPart - minPart);
		}
		return distance;
	}

	vector<Index> computeDaisyChainPartitionDegrees(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<std::pair<Index, Index> >& hedgeMinMax) {
		vector<Index> ret(hypergraph.nParts(), 0);
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			Index minPart = hedgeMinMax[hedge].first;
			Index maxPart = hedgeMinMax[hedge].second;
			for (Index p = minPart; p < maxPart; ++p) {
				ret[p] += hypergraph.hedgeWeight(hedge);
				ret[p + 1] += hypergraph.hedgeWeight(hedge);
			}
		}
		return ret;
	}

	Index computeCut(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<Index>& hedgeDegrees) {
		Index ret = 0;
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			if (hedgeDegrees[hedge] > 1)
				ret += hypergraph.hedgeWeight(hedge);
		}
		return ret;
	}

	Index computeSoed(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<Index>& hedgeDegrees) {
		Index ret = 0;
		for (Index hedge = 0; hedge < hypergraph.nHedges(); ++hedge) {
			ret += hypergraph.hedgeWeight(hedge) * hedgeDegrees[hedge];
		}
		return ret;
	}

	Index computeSumOverflow(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<Index>& partitionDemands) {
		Index ret = 0;
		for (Index p = 0; p < hypergraph.nParts(); ++p) {
			ret += std::max(partitionDemands[p] - hypergraph.partWeight(p), (Index)0);
		}
		return ret;
	}

	Index countEmptyPartitions(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<Index>& partitionDemands) {
		Index count = 0;
		for (Index d : partitionDemands) {
			if (d == 0) count++;
		}
		return count;
	}

	double computeRatioPenalty(const MiniPartHypergraph::Hypergraph& hypergraph, const vector<Index>& partitionDemands) {
		Index sumDemands = 0;
		for (Index d : partitionDemands)
			sumDemands += d;
		double normalizedDemands = ((double)sumDemands) / partitionDemands.size();
		double productDemands = 1.0;
		for (Index d : partitionDemands) {
			productDemands *= (d / normalizedDemands);
		}
		// Geomean squared
		return 1.0 / pow(productDemands, 2.0 / partitionDemands.size());
	}

	Index computeMaxDegree(const MiniPartHypergraph::Hypergraph&, const vector<Index>& partitionDegrees) {
		return *max_element(partitionDegrees.begin(), partitionDegrees.end());
	}
}

/**
 * Base class for incremental evaluation
 */
class IncrementalObjective {
public:
	IncrementalObjective(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution, Index nObjectives) : hypergraph_(hypergraph)
	                                                                                                                       , solution_(solution)
	                                                                                                                       , objectives_(nObjectives, 0)
	{
		assert(hypergraph_.nNodes() == solution_.nNodes());
		assert(hypergraph_.nParts() == solution_.nParts());
	}

	virtual void move(Index node, Index to) = 0;
	virtual void checkConsistency() const
	{}

	Index nNodes() const { return hypergraph_.nNodes(); }
	Index nHedges() const { return hypergraph_.nHedges(); }
	Index nParts() const { return hypergraph_.nParts(); }
	Index nObjectives() const { return objectives_.size(); }

	const MiniPartHypergraph::Hypergraph& hypergraph() const { return hypergraph_; }
	const hgsol::Solution& solution() const { return solution_; }
	const std::vector<int64_t>& objectives() const { return objectives_; }

	virtual ~IncrementalObjective() {}

protected:
	const MiniPartHypergraph::Hypergraph& hypergraph_;
	hgsol::Solution& solution_;
	std::vector<int64_t> objectives_;
};

class IncrementalCut final : public IncrementalObjective {
public:
	IncrementalCut(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution) : IncrementalObjective(hypergraph, solution, 3)
	{
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		currentCut_ = computeCut(hypergraph, hedgeDegrees_);
		currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
		setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			if (pins[to] == 1 && pins[from] != 0) {
				++hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 2) {
					currentCut_ += hypergraph_.hedgeWeight(hedge);
				}
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (pins[to] != 1 && pins[from] == 0) {
				--hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 1) {
					currentCut_ -= hypergraph_.hedgeWeight(hedge);
				}
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(currentCut_ == computeCut(hypergraph_, hedgeDegrees_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = computeSumOverflow(hypergraph_, partitionDemands_);
		objectives_[1] = currentCut_;
		objectives_[2] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	Index currentCut_;
	Index currentSoed_;
};

class IncrementalSoed final : public IncrementalObjective {
public:
	IncrementalSoed(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution) : IncrementalObjective(hypergraph, solution, 2)
	{
			partitionDemands_ = computePartitionDemands(hypergraph, solution);
			hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
			hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
			currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
			setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			if (pins[to] == 1 && pins[from] != 0) {
				++hedgeDegrees_[hedge];
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (pins[to] != 1 && pins[from] == 0) {
				--hedgeDegrees_[hedge];
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = computeSumOverflow(hypergraph_, partitionDemands_);
		objectives_[1] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	Index currentSoed_;
};

class IncrementalMaxDegree final : public IncrementalObjective {
public:
	IncrementalMaxDegree(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution)
		: IncrementalObjective(hypergraph, solution, 3) {
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		partitionDegrees_ = computePartitionDegrees(hypergraph, hedgeDegrees_, hedgeNbPinsPerPartition_);
		currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
		setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			bool becomesCut = false;
			bool becomesUncut = false;
			if (pins[to] == 1 && pins[from] != 0) {
				++hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 2) {
					becomesCut = true;
				}
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (pins[to] != 1 && pins[from] == 0) {
				--hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 1) {
					becomesUncut = true;
				}
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}

			if (becomesUncut) {
				partitionDegrees_[from] -= hypergraph_.hedgeWeight(hedge);
				partitionDegrees_[to] -= hypergraph_.hedgeWeight(hedge);
			}
			else if (becomesCut) {
				partitionDegrees_[from] += hypergraph_.hedgeWeight(hedge);
				partitionDegrees_[to] += hypergraph_.hedgeWeight(hedge);
			}
			else if (hedgeDegrees_[hedge] >= 2) {
				if (pins[from] == 0) {
					partitionDegrees_[from] -= hypergraph_.hedgeWeight(hedge);
				}
				if (pins[to] == 1) {
					partitionDegrees_[to] += hypergraph_.hedgeWeight(hedge);
				}
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(partitionDegrees_ == computePartitionDegrees(hypergraph_, hedgeDegrees_, hedgeNbPinsPerPartition_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = computeSumOverflow(hypergraph_, partitionDemands_);
		objectives_[1] = computeMaxDegree(hypergraph_, partitionDegrees_);
		objectives_[2] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	std::vector<Index> partitionDegrees_;
	Index currentSoed_;
};

class IncrementalDaisyChainDistance final : public IncrementalObjective {
public:
	IncrementalDaisyChainDistance(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution)
		: IncrementalObjective(hypergraph, solution, 3) {
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		hedgeMinMax_ = computeDaisyChainMinMax(hypergraph, hedgeNbPinsPerPartition_);
		currentDistance_ = computeDaisyChainDistance(hypergraph, hedgeMinMax_);
		currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
		setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			bool reachesPart = pins[to] == 1;
			bool leavesPart = pins[from] == 0;
			if (reachesPart) {
				++hedgeDegrees_[hedge];
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (leavesPart) {
				--hedgeDegrees_[hedge];
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}
			if (reachesPart || leavesPart) {
				Index minAfter = nParts() - 1;
				Index maxAfter = 0;
				for (Index p = 0; p < nParts(); ++p) {
					bool countsAfter = pins[p] != 0;
					if (countsAfter) {
						minAfter = std::min(minAfter, p);
						maxAfter = std::max(maxAfter, p);
					}
				}
				Index minBefore = hedgeMinMax_[hedge].first;
				Index maxBefore = hedgeMinMax_[hedge].second;
				hedgeMinMax_[hedge] = std::make_pair(minAfter, maxAfter);
				currentDistance_ += hypergraph_.hedgeWeight(hedge) * (maxAfter - minAfter - maxBefore + minBefore);
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(hedgeMinMax_ == computeDaisyChainMinMax(hypergraph_, hedgeNbPinsPerPartition_));
		assert(currentDistance_ == computeDaisyChainDistance(hypergraph_, hedgeMinMax_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = computeSumOverflow(hypergraph_, partitionDemands_);
		objectives_[1] = currentDistance_;
		objectives_[2] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	std::vector<std::pair<Index, Index> > hedgeMinMax_;
	Index currentDistance_;
	Index currentSoed_;
};


class IncrementalDaisyChainMaxDegree final : public IncrementalObjective {
public:
	IncrementalDaisyChainMaxDegree(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution)
		: IncrementalObjective(hypergraph, solution, 3) {
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		hedgeMinMax_ = computeDaisyChainMinMax(hypergraph, hedgeNbPinsPerPartition_);
		partitionDegrees_ = computeDaisyChainPartitionDegrees(hypergraph, hedgeMinMax_);
		currentDistance_ = computeDaisyChainDistance(hypergraph, hedgeMinMax_);
		setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			bool reachesPart = pins[to] == 1;
			bool leavesPart = pins[from] == 0;
			if (reachesPart) {
				++hedgeDegrees_[hedge];
			}
			if (leavesPart) {
				--hedgeDegrees_[hedge];
			}
			if (reachesPart || leavesPart) {
				Index minAfter = nParts() - 1;
				Index maxAfter = 0;
				for (Index p = 0; p < nParts(); ++p) {
					bool countsAfter = pins[p] != 0;
					if (countsAfter) {
						minAfter = std::min(minAfter, p);
						maxAfter = std::max(maxAfter, p);
					}
				}
				Index minBefore = hedgeMinMax_[hedge].first;
				Index maxBefore = hedgeMinMax_[hedge].second;
				hedgeMinMax_[hedge] = std::make_pair(minAfter, maxAfter);
				if (minAfter != minBefore || maxAfter != maxBefore) {
					currentDistance_ += hypergraph_.hedgeWeight(hedge) * (maxAfter - minAfter - maxBefore + minBefore);
					// Update degrees
					for (Index p = minBefore; p < maxBefore; ++p) {
						partitionDegrees_[p] -= hypergraph_.hedgeWeight(hedge);
						partitionDegrees_[p + 1] -= hypergraph_.hedgeWeight(hedge);
					}
					for (Index p = minAfter; p < maxAfter; ++p) {
						partitionDegrees_[p] += hypergraph_.hedgeWeight(hedge);
						partitionDegrees_[p + 1] += hypergraph_.hedgeWeight(hedge);
					}
				}
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(hedgeMinMax_ == computeDaisyChainMinMax(hypergraph_, hedgeNbPinsPerPartition_));
		assert(currentDistance_ == computeDaisyChainDistance(hypergraph_, hedgeMinMax_));
		assert(partitionDegrees_ == computeDaisyChainPartitionDegrees(hypergraph_, hedgeMinMax_));
	}

private:
	void setObjective()
	{
		objectives_[0] = computeSumOverflow(hypergraph_, partitionDemands_);
		objectives_[1] = computeMaxDegree(hypergraph_, partitionDegrees_);
		objectives_[2] = currentDistance_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	std::vector<std::pair<Index, Index> > hedgeMinMax_;
	std::vector<Index> partitionDegrees_;
	Index currentDistance_;
};

class IncrementalRatioCut final : public IncrementalObjective {
public:
	IncrementalRatioCut(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution) : IncrementalObjective(hypergraph, solution, 4) {
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		currentCut_ = computeCut(hypergraph, hedgeDegrees_);
		currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
		setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			if (pins[to] == 1 && pins[from] != 0) {
				++hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 2) {
					currentCut_ += hypergraph_.hedgeWeight(hedge);
				}
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (pins[to] != 1 && pins[from] == 0) {
				--hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 1) {
					currentCut_ -= hypergraph_.hedgeWeight(hedge);
				}
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(currentCut_ == computeCut(hypergraph_, hedgeDegrees_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = countEmptyPartitions(hypergraph_, partitionDemands_);
		objectives_[1] = 100.0 * currentCut_ * computeRatioPenalty(hypergraph_, partitionDemands_);
		objectives_[2] = currentCut_;
		objectives_[3] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	Index currentCut_;
	Index currentSoed_;
};

class IncrementalRatioSoed final : public IncrementalObjective {
public:
	IncrementalRatioSoed(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution) : IncrementalObjective(hypergraph, solution, 3) {
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
		setObjective();
	}
	
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			if (pins[to] == 1 && pins[from] != 0) {
				++hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 2) {
					currentCut_ += hypergraph_.hedgeWeight(hedge);
				}
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (pins[to] != 1 && pins[from] == 0) {
				--hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 1) {
					currentCut_ -= hypergraph_.hedgeWeight(hedge);
				}
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = countEmptyPartitions(hypergraph_, partitionDemands_);
		objectives_[1] = 100.0 * currentSoed_ * computeRatioPenalty(hypergraph_, partitionDemands_);
		objectives_[2] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	Index currentCut_;
	Index currentSoed_;
};

class IncrementalRatioMaxDegree final : public IncrementalObjective {
public:
	IncrementalRatioMaxDegree(const MiniPartHypergraph::Hypergraph& hypergraph, hgsol::Solution& solution) : IncrementalObjective(hypergraph, solution, 3) {
		partitionDemands_ = computePartitionDemands(hypergraph, solution);
		hedgeNbPinsPerPartition_ = computeHedgeNbPinsPerPartition(hypergraph, solution);
		hedgeDegrees_ = computeHedgeDegrees(hypergraph, hedgeNbPinsPerPartition_);
		partitionDegrees_ = computePartitionDegrees(hypergraph, hedgeDegrees_, hedgeNbPinsPerPartition_);
		currentSoed_ = computeSoed(hypergraph, hedgeDegrees_);
		setObjective();
	}
	void move(Index node, Index to) override
	{
		assert(to < nParts() && to >= 0);
		Index from = solution_[node];
		if (from == to) return;
		solution_[node] = to;
		partitionDemands_[to] += hypergraph_.nodeWeight(node);
		partitionDemands_[from] -= hypergraph_.nodeWeight(node);

		for (Index hedge : hypergraph_.nodeHedges(node)) {
			vector<Index>& pins = hedgeNbPinsPerPartition_[hedge];
			++pins[to];
			--pins[from];
			bool becomesCut = false;
			bool becomesUncut = false;
			if (pins[to] == 1 && pins[from] != 0) {
				++hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 2) {
					becomesCut = true;
				}
				currentSoed_ += hypergraph_.hedgeWeight(hedge);
			}
			if (pins[to] != 1 && pins[from] == 0) {
				--hedgeDegrees_[hedge];
				if (hedgeDegrees_[hedge] == 1) {
					becomesUncut = true;
				}
				currentSoed_ -= hypergraph_.hedgeWeight(hedge);
			}

			if (becomesUncut) {
				partitionDegrees_[from] -= hypergraph_.hedgeWeight(hedge);
				partitionDegrees_[to] -= hypergraph_.hedgeWeight(hedge);
			}
			else if (becomesCut) {
				partitionDegrees_[from] += hypergraph_.hedgeWeight(hedge);
				partitionDegrees_[to] += hypergraph_.hedgeWeight(hedge);
			}
			else if (hedgeDegrees_[hedge] >= 2) {
				if (pins[from] == 0) {
					partitionDegrees_[from] -= hypergraph_.hedgeWeight(hedge);
				}
				if (pins[to] == 1) {
					partitionDegrees_[to] += hypergraph_.hedgeWeight(hedge);
				}
			}
		}
		setObjective();
	}
	void checkConsistency() const override
	{
		assert(partitionDemands_ == computePartitionDemands(hypergraph_, solution_));
		assert(hedgeNbPinsPerPartition_ == computeHedgeNbPinsPerPartition(hypergraph_, solution_));
		assert(hedgeDegrees_ == computeHedgeDegrees(hypergraph_, hedgeNbPinsPerPartition_));
		assert(partitionDegrees_ == computePartitionDegrees(hypergraph_, hedgeDegrees_, hedgeNbPinsPerPartition_));
		assert(currentSoed_ == computeSoed(hypergraph_, hedgeDegrees_));
	}

private:
	void setObjective()
	{
		objectives_[0] = countEmptyPartitions(hypergraph_, partitionDemands_);
		objectives_[1] = 100.0 * computeMaxDegree(hypergraph_, partitionDegrees_) * computeRatioPenalty(hypergraph_, partitionDemands_);
		objectives_[2] = currentSoed_;
	}

private:
	std::vector<Index> partitionDemands_;
	std::vector<std::vector<Index> > hedgeNbPinsPerPartition_;
	std::vector<Index> hedgeDegrees_;
	std::vector<Index> partitionDegrees_;
	Index currentSoed_;
};





