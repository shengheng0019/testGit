#pragma once


#include <random>

#include "common.hpp"
#include "incremental_objective.hpp"

namespace {

    void tryMoveRandomBlock(IncrementalObjective& inc, std::mt19937& rgen, Index node) {
	    std::uniform_int_distribution<Index> partDist(0, inc.nParts() - 1);
        Index src = inc.solution()[node];
        Index dst = partDist(rgen);

        vector<int64_t> before = inc.objectives();
        inc.move(node, dst);
        vector<int64_t> after = inc.objectives();
        if (before < after) {
            inc.move(node, src);
        }
    }

    void tryMoveBestBlock(IncrementalObjective& inc, std::mt19937& rgen, Index node) {
        Index src = inc.solution()[node];
        Index bestDst = src;
        vector<int64_t> bestObj = inc.objectives();
        for (Index dst = 0; dst < inc.nParts(); ++dst) {
            if (dst == src) continue;
            inc.move(node, dst);
            vector<int64_t> curObj = inc.objectives();
            if (curObj < bestObj) {
                bestObj = curObj;
                bestDst = dst;
            }
        }
        inc.move(node, bestDst);
    }

}

namespace hyperMove
{
    class Move {
    public:
        Move(Index budget) : budget_(budget) {}
        virtual void run(IncrementalObjective& inc, std::mt19937& rgen) = 0;
        virtual ~Move() {}

        std::int64_t budget_;
    };

    class VertexMoveRandomBlock : public Move {
    public:
        VertexMoveRandomBlock(Index budget) : Move(budget) {}
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            std::uniform_int_distribution<Index> nodeDist(0, inc.nNodes() - 1);
            Index node = nodeDist(rgen);
            tryMoveRandomBlock(inc, rgen, node);
            --this->budget_;
        }
    };

    class VertexMoveBestBlock : public Move {
    public:
        VertexMoveBestBlock(Index budget) : Move(budget) {}
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            std::uniform_int_distribution<Index> nodeDist(0, inc.nNodes() - 1);
            Index node = nodeDist(rgen);
            tryMoveBestBlock(inc, rgen, node);
            this->budget_ -= inc.nParts() - 1;
        }
    };

    class VertexPassRandomBlock : public Move {
    public:
        VertexPassRandomBlock(Index budget) : Move(budget) {}
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            vector<Index> nodes;
            nodes.reserve(inc.nNodes());
            for (Index i = 0; i < inc.nNodes(); ++i) {
                nodes.push_back(i);
            }
            shuffle(nodes.begin(), nodes.end(), rgen);
            for (Index node : nodes) {
                if (this->budget_ <= 0) break;
                tryMoveBestBlock(inc, rgen, node);
                --this->budget_;
            }
        }
    };

    class VertexPassBestBlock : public Move {
    public:
        VertexPassBestBlock(Index budget) : Move(budget) {}
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            vector<Index> nodes;
            nodes.reserve(inc.nNodes());
            for (Index i = 0; i < inc.nNodes(); ++i) {
                nodes.push_back(i);
            }
            shuffle(nodes.begin(), nodes.end(), rgen);
            for (Index node : nodes) {
                if (this->budget_ <= 0) break;
                tryMoveBestBlock(inc, rgen, node);
                this->budget_ -= inc.nParts() - 1;
            }
        }
    };

    class EdgeMoveRandomBlock : public Move {
    public:
        EdgeMoveRandomBlock(Index budget) : Move(budget) {
            edgeDegreeCutoff_ = 10;
        }
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            if (inc.nHedges() == 0) {
                --this->budget_;
                return;
            }

            initialStatus_.clear();
            std::uniform_int_distribution<Index> edgeDist(0, inc.nHedges() - 1);
            std::uniform_int_distribution<Index> partDist(0, inc.nParts() - 1);
            Index hedge = edgeDist(rgen);
            Index dst = partDist(rgen);
            if (inc.hypergraph().hedgeNodes(hedge).size() > edgeDegreeCutoff_) {
                --this->budget_;
                return;
            }

            vector<int64_t> before = inc.objectives();
            for (Index node : inc.hypergraph().hedgeNodes(hedge)) {
                Index src = inc.solution()[node];
                inc.move(node, dst);
                initialStatus_.emplace_back(node, src);
            }
            vector<int64_t> after = inc.objectives();
            if (before < after) {
                for (std::pair<Index, Index> status : initialStatus_) {
                    inc.move(status.first, status.second);
                }
            }
            this->budget_ -= inc.hypergraph().hedgeNodes(hedge).size();
        }

    private:
        std::vector<std::pair<Index, Index> > initialStatus_;
        std::size_t edgeDegreeCutoff_;
    };

    class VertexSwap : public Move {
    public:
        VertexSwap(Index budget) : Move(budget) {}
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            --this->budget_;
            std::uniform_int_distribution<Index> nodeDist(0, inc.nNodes() - 1);
            Index n1 = nodeDist(rgen);
            Index n2 = nodeDist(rgen);
            Index p1 = inc.solution()[n1];
            Index p2 = inc.solution()[n2];
            if (p1 == p2) return;

            vector<int64_t> before = inc.objectives();
            inc.move(n1, p2);
            inc.move(n2, p1);
            vector<int64_t> after = inc.objectives();
            if (before < after) {
                inc.move(n1, p1);
                inc.move(n2, p2);
            }
        }
    };

    class VertexAbsorptionPass : public Move {
    public:
        VertexAbsorptionPass(Index budget) : Move(budget) {
            nodeDegreeCutoff_ = 10;
            edgeDegreeCutoff_ = 10;
        }
        void run(IncrementalObjective& inc, std::mt19937& rgen) override
        {
            assert(this->budget_ > 0);
            std::uniform_int_distribution<Index> nodeDist(0, inc.nNodes() - 1);
            std::uniform_int_distribution<Index> partDist(0, inc.nParts() - 1);
            candidates_.clear();
            Index dst = partDist(rgen);
            candidates_.push_back(nodeDist(rgen));

            while (!candidates_.empty() && this->budget_ > 0) {
                Index node = candidates_.back();
                candidates_.pop_back();
                Index src = inc.solution()[node];
                if (src == dst)
                    continue;
                --this->budget_;

                vector<int64_t> before = inc.objectives();
                inc.move(node, dst);
                vector<int64_t> after = inc.objectives();
                if (before < after) {
                    inc.move(node, src);
                }
                else {
                    if (inc.hypergraph().nodeHedges(node).size() <= nodeDegreeCutoff_) {
                        for (Index hEdge : inc.hypergraph().nodeHedges(node)) {
                            if (inc.hypergraph().hedgeNodes(hEdge).size() <= edgeDegreeCutoff_) {
                                for (Index neighbour : inc.hypergraph().hedgeNodes(hEdge)) {
                                    candidates_.push_back(neighbour);
                                }
                            }
                        }
                    }
                }
            }
        }

    private:
        std::vector<Index> candidates_;
        std::size_t nodeDegreeCutoff_;
        std::size_t edgeDegreeCutoff_;
    };
}




