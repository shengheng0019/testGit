
// Copyright (C) 2019 Gabriel Gouvine - All Rights Reserved

#pragma once

#include <iostream>

#include "common.hpp"
#include "solution.hpp"

#include <unordered_map>
#include <unordered_set>

template<typename T>
class Range {
public:
    Range(const T* b, const T* e) : begin_(b), end_(e) {}
    const T* begin() const { return begin_; }
    const T* end() const { return end_; }
    const T* cbegin() const { return begin_; }
    const T* cend() const { return end_; }
    std::size_t size() const { return end_ - begin_; }

    bool operator==(const Range<T>& o) {
        if (size() != o.size()) return false;
        for (std::size_t i = 0, sz = size(); i != sz; ++i) {
            if (begin_[i] != o.begin_[i]) return false;
        }
        return true;
    }
    bool operator!=(const Range<T>& o) { return !operator==(o); }

private:
    const T* begin_;
    const T* end_;
};


namespace MiniPartHypergraph
{
    class Hypergraph {
    public:
        Hypergraph(Index nodeWeights = 1, Index hedgeWeights = 1, Index partWeights = 1)
        {
            nNodes_ = 0;
            nHedges_ = 0;
            nParts_ = 0;
            nNodeWeights_ = nodeWeights;
            nHedgeWeights_ = hedgeWeights;
            nPartWeights_ = partWeights;
            totalNodeWeights_.assign(nodeWeights, 0);
            totalHedgeWeights_.assign(hedgeWeights, 0);
            totalPartWeights_.assign(partWeights, 0);
            hedgeBegin_.push_back(0);
            nodeBegin_.push_back(0);
        }

        Index nNodes() const { return nNodes_; }
        Index nHedges() const { return nHedges_; }
        Index nParts() const { return nParts_; }
        Index nPins() const { return nPins_; }

        Index nNodeWeights() const { return nNodeWeights_; }
        Index nHedgeWeights() const { return nHedgeWeights_; }
        Index nPartWeights() const { return nPartWeights_; }

        Index totalNodeWeight(Index i = 0) const { return totalNodeWeights_[i]; }
        Index totalHedgeWeight(Index i = 0) const { return totalHedgeWeights_[i]; }
        Index totalPartWeight(Index i = 0) const { return totalPartWeights_[i]; }

        Index nodeWeight(Index node, Index i = 0) const { return nodeData_[nodeBegin_[node] + i]; }
        Index hedgeWeight(Index hedge, Index i = 0) const { return hedgeData_[hedgeBegin_[hedge] + i]; }
        Index partWeight(Index part, Index i = 0) const { return partData_[part * nPartWeights_ + i]; }

        Range<Index> hedgeNodes(Index hedge) const {
            const Index* ptr = hedgeData_.data();
            Index b = hedgeBegin_[hedge] + nHedgeWeights_;
            Index e = hedgeBegin_[hedge + 1];
            return Range<Index>(ptr + b, ptr + e);
        }

        Range<Index> nodeHedges(Index node) const {
            const Index* ptr = nodeData_.data();
            Index b = nodeBegin_[node] + nNodeWeights_;
            Index e = nodeBegin_[node + 1];
            return Range<Index>(ptr + b, ptr + e);
        }

        // Metrics
        Index metricsSumOverflow(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nNodeWeights() == 1);
            std::vector<Index> usage = metricsPartitionUsage(solution);
            Index ret = 0;
            for (int i = 0; i < nParts(); ++i) {
                Index ovf = usage[i] - partData_[i];
                if (ovf > 0)
                    ret += ovf;
            }
            return ret;
        }
        Index metricsEmptyPartitions(const hgsol::Solution& solution) const
        {
            std::vector<Index> partitionUsage = metricsPartitionUsage(solution);
            Index count = 0;
            for (Index d : partitionUsage) {
                if (d == 0) count++;
            }
            return count;
        }
        Index metricsCut(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nHedgeWeights() == 1);
            Index ret = 0;
            for (Index hedge = 0; hedge < nHedges(); ++hedge) {
                if (cut(solution, hedge))
                    ret += hedgeWeight(hedge);
            }
            return ret;
        }
        Index metricsSoed(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nHedgeWeights() == 1);
            Index ret = 0;
            for (Index hedge = 0; hedge < nHedges(); ++hedge) {
                ret += hedgeWeight(hedge) * degree(solution, hedge);
            }
            return ret;
        }
        Index metricsConnectivity(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nHedgeWeights() == 1);
            Index ret = 0;
            for (Index hedge = 0; hedge < nHedges(); ++hedge) {
                ret += hedgeWeight(hedge) * (degree(solution, hedge) - 1);
            }
            return ret;
        }
        Index metricsMaxDegree(const hgsol::Solution& solution) const
        {
            std::vector<Index> degree = metricsPartitionDegree(solution);
            return *std::max_element(degree.begin(), degree.end());
        }
        Index metricsDaisyChainDistance(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nHedgeWeights() == 1);
            Index ret = 0;
            for (Index hedge = 0; hedge < nHedges(); ++hedge) {
                Index minPart = nParts() - 1;
                Index maxPart = 0;
                for (Index node : hedgeNodes(hedge)) {
                    minPart = std::min(minPart, solution[node]);
                    maxPart = std::max(maxPart, solution[node]);
                }
                if (minPart >= maxPart) continue;
                ret += hedgeWeight(hedge) * (maxPart - minPart);
            }
            return ret;
        }
        Index metricsDaisyChainMaxDegree(const hgsol::Solution& solution) const
        {
            std::vector<Index> degree = metricsPartitionDaisyChainDegree(solution);
            return *std::max_element(degree.begin(), degree.end());
        }
        double metricsRatioCut(const hgsol::Solution& solution) const
        {
            return metricsCut(solution) * metricsRatioPenalty(solution);
        }
        double metricsRatioSoed(const hgsol::Solution& solution) const
        {
            return metricsSoed(solution) * metricsRatioPenalty(solution);
        }
        double metricsRatioConnectivity(const hgsol::Solution& solution) const
        {
            return metricsConnectivity(solution) * metricsRatioPenalty(solution);
        }
        double metricsRatioMaxDegree(const hgsol::Solution& solution) const
        {
            return metricsMaxDegree(solution) * metricsRatioPenalty(solution);
        }

        double metricsRatioPenalty(const hgsol::Solution& solution) const
        {
            std::vector<Index> partitionUsage = metricsPartitionUsage(solution);
            Index sumUsage = 0;
            for (Index d : partitionUsage)
                sumUsage += d;
            double normalizedUsage = ((double)sumUsage) / partitionUsage.size();
            double productUsage = 1.0;
            for (Index d : partitionUsage) {
                productUsage *= (d / normalizedUsage);
            }
            // Geomean squared
            return 1.0 / pow(productUsage, 2.0 / partitionUsage.size());
        }
        std::vector<Index> metricsPartitionUsage(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            std::vector<Index> usage(nParts(), 0);
            for (int i = 0; i < nNodes(); ++i) {
                assert(solution[i] >= 0 && solution[i] < nParts());
                usage[solution[i]] += nodeWeight(i);
            }
            return usage;
        }
        std::vector<Index> metricsPartitionDegree(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nHedgeWeights() == 1);
            std::vector<Index> degree(nParts(), 0);
            std::unordered_set<Index> parts;
            for (int i = 0; i < nHedges(); ++i) {
                parts.clear();
                for (Index node : hedgeNodes(i)) {
                    parts.insert(solution[node]);
                }
                if (parts.size() > 1) {
                    for (Index p : parts) {
                        degree[p] += hedgeWeight(i);
                    }
                }
            }
            return degree;
        }
        std::vector<Index> metricsPartitionDaisyChainDegree(const hgsol::Solution& solution) const
        {
            assert(solution.nNodes() == nNodes());
            assert(solution.nParts() == nParts());
            assert(nHedgeWeights() == 1);
            std::vector<Index> degree(nParts(), 0);
            std::unordered_set<Index> parts;
            for (int i = 0; i < nHedges(); ++i) {
                Index minPart = nParts() - 1;
                Index maxPart = 0;
                for (Index node : hedgeNodes(i)) {
                    minPart = std::min(minPart, solution[node]);
                    maxPart = std::max(maxPart, solution[node]);
                }
                if (minPart >= maxPart) continue;
                // Count twice for middle partitions
                for (Index p = minPart; p < maxPart; ++p) {
                    degree[p] += hedgeWeight(i);
                    degree[p + 1] += hedgeWeight(i);
                }
            }
            return degree;
        }

        // IO functions
        static Hypergraph readFile(std::unordered_map<int, std::vector<uint64_t>> matrix, Index nNodes, Index nHedges)
        {
            std::cout << "³¬Í¼ÐÅÏ¢: " << nNodes << "," << nHedges << std::endl;
            Hypergraph ret;
            ret.nNodes_ = nNodes;
            ret.nHedges_ = nHedges;
            ret.nParts_ = 0;
            // Read edges
            ret.hedgeBegin_.reserve(nHedges);
            ret.hedgeData_.reserve(3 * nHedges);
            //std::cout << ret.hedgeBegin_.size()<<","<<ret.nodeBegin_.size()<<std::endl;
            //std::cout << nNodes << " " << nHedges << std::endl;

            std::vector<Index> nodes;
            for (Index i = 0; i < nHedges; ++i) {
                Index w = matrix[i].size();
                for (Index j = 0; j < matrix[i].size(); ++j)
                {
                    //std::cout << matrix[i][j]<<",";
                    nodes.push_back(matrix[i][j]);
                }
                //std::cout << i<<" : ";
                //for (auto node : nodes)
                //{
                   // std::cout << node<<",";
                //}
                //std::cout <<std::endl;
                std::sort(nodes.begin(), nodes.end());
                nodes.resize(std::unique(nodes.begin(), nodes.end()) - nodes.begin());
                if (nodes.empty()) throw std::runtime_error("No node on the line");
                ret.hedgeData_.push_back(w);
                ret.hedgeData_.insert(ret.hedgeData_.end(), nodes.begin(), nodes.end());
                ret.hedgeBegin_.push_back(ret.hedgeData_.size());
                nodes.clear();
            }
            // Read node weights
            ret.nodeData_.reserve(nNodes + ret.hedgeData_.size() - nHedges);
            ret.nodeData_.assign(nNodes, 1);
            for (Index i = 1; i <= nNodes; ++i)
            {
                ret.nodeBegin_.push_back(i);
            }

            // Finalize nodes 
            ret.finalize();
            return ret;
        }


        // Coarsening
        Hypergraph coarsen(const hgsol::Solution& coarsening) const
        {
            assert(nNodes() == coarsening.nNodes());
            assert(coarsening.nParts() <= nNodes());
            assert(coarsening.nParts() > 0);
            if (nNodes() == 0) return *this;

            Hypergraph ret(nNodeWeights_, nHedgeWeights_, nPartWeights_);
            ret.nNodes_ = coarsening.nParts();
            ret.nHedges_ = 0;
            ret.nParts_ = nParts_;

            // Hyperedges
            std::vector<Index> pins;
            for (Index hedge = 0; hedge < nHedges_; ++hedge) {
                for (Index node : hedgeNodes(hedge)) {
                    pins.push_back(coarsening[node]);
                }
                sort(pins.begin(), pins.end());
                pins.resize(std::unique(pins.begin(), pins.end()) - pins.begin());
                if (pins.size() > 1) {
                    for (Index i = 0; i < nHedgeWeights_; ++i) {
                        ret.hedgeData_.push_back(hedgeWeight(hedge, i));
                    }
                    ret.hedgeData_.insert(ret.hedgeData_.end(), pins.begin(), pins.end());
                    ret.hedgeBegin_.push_back(ret.hedgeData_.size());
                    ++ret.nHedges_;
                }
                pins.clear();
            }

            // Node weights
            ret.nodeData_.assign(coarsening.nParts() * nNodeWeights_, 0);
            for (Index node = 0; node < nNodes_; ++node) {
                Index coarsened = coarsening[node];
                for (Index i = 0; i < nNodeWeights_; ++i) {
                    ret.nodeData_[coarsened * nNodeWeights_ + i] += nodeWeight(node, i);
                }
            }
            for (Index i = 1; i <= coarsening.nParts(); ++i) {
                ret.nodeBegin_.push_back(i * nNodeWeights_);
            }

            // Partitions
            ret.partData_ = partData_;

            // Finalize
            ret.mergeParallelHedges();

            assert(ret.totalNodeWeights_ == totalNodeWeights_);
            assert(ret.totalPartWeights_ == totalPartWeights_);

            return ret;
        }

        // Modifications
        void setupBlocks(Index nParts, double imbalanceFactor)
        {
            // Setup partitions with weights proportional to the nodes
            if (nPartWeights_ != nNodeWeights_)
                throw std::runtime_error("Unable to generate partition capacities when the number of weights is not the same for parts and nodes");
            nParts_ = nParts;
            if (nParts > 0) {
                partData_.resize(nPartWeights_ * nParts);
                for (Index i = 0; i < nPartWeights_; ++i) {
                    Index totalCapacity = totalNodeWeight(i) * (1.0 + imbalanceFactor);
                    Index partitionCapacity = totalCapacity / nParts;
                    partData_[i] = totalCapacity - partitionCapacity * (nParts - 1);
                    for (Index p = 1; p < nParts; ++p) {
                        partData_[p * nPartWeights_ + i] = partitionCapacity;
                    }
                }
            }
            else {
                partData_.clear();
            }
            finalizePartWeights();
        }


        void mergeParallelHedges()
        {
            std::vector<Index> newHedgeBegin;
            std::vector<Index> newHedgeData;
            newHedgeBegin.push_back(0);

            class HedgeHasher {
            public:
                HedgeHasher(const Hypergraph& hypergraph)
                    : hypergraph_(hypergraph) {}

                size_t operator()(const Index& hedge) const {
                    // FNV hash
                    uint64_t magic = 1099511628211llu;
                    uint64_t ret = 0;
                    for (Index n : hypergraph_.hedgeNodes(hedge)) {
                        ret = (ret ^ (uint64_t)n) * magic;
                    }
                    return ret;
                }

            private:
                const Hypergraph& hypergraph_;
            };

            class HedgeComparer {
            public:
                HedgeComparer(const Hypergraph& hypergraph)
                    : hypergraph_(hypergraph) {}

                bool operator()(const Index& h1, const Index& h2) const {
                    return hypergraph_.hedgeNodes(h1) == hypergraph_.hedgeNodes(h2);
                }

            private:
                const Hypergraph& hypergraph_;
            };

            std::unordered_map<Index, Index, HedgeHasher, HedgeComparer> hedgeWeightMap(nHedges(), HedgeHasher(*this), HedgeComparer(*this));
            for (Index hedge = 0; hedge < nHedges_; ++hedge) {
                Index ind = newHedgeBegin.size() - 1;
                auto p = hedgeWeightMap.emplace(hedge, ind);
                if (p.second) {
                    // No such hyperedge exists yet
                    for (Index i = 0; i < nHedgeWeights_; ++i) {
                        newHedgeData.push_back(hedgeWeight(hedge, i));
                    }
                    for (Index node : hedgeNodes(hedge)) {
                        newHedgeData.push_back(node);
                    }
                    newHedgeBegin.push_back(newHedgeData.size());
                }
                else {
                    // An equivalent hyperedge exists already
                    Index s = newHedgeBegin[p.first->second];
                    for (Index i = 0; i < nHedgeWeights_; ++i) {
                        newHedgeData[s + i] += hedgeWeight(hedge, i);
                    }
                }
            }

            nHedges_ = newHedgeBegin.size() - 1;
            hedgeBegin_ = newHedgeBegin;
            hedgeData_ = newHedgeData;

            finalize();
        }

        void checkConsistency() const
        {
            if (nNodes_ < 0) throw std::runtime_error("Negative number of nodes");
            if (nHedges_ < 0) throw std::runtime_error("Negative number of hedges");
            if (nParts_ < 0) throw std::runtime_error("Negative number of parts");
            if (nPins_ < 0) throw std::runtime_error("Negative number of pins");
            if (nNodeWeights_ < 0) throw std::runtime_error("Negative number of node weights");
            if (nHedgeWeights_ < 0) throw std::runtime_error("Negative number of hedge weights");
            if (nPartWeights_ < 0) throw std::runtime_error("Negative number of part weights");

            if ((Index)nodeBegin_.size() != nNodes_ + 1)
                throw std::runtime_error("Number of node limits and of nodes do not match");
            if ((Index)hedgeBegin_.size() != nHedges_ + 1)
                throw std::runtime_error("Number of hedge limits and of hedges do not match");

            for (size_t i = 0; i + 1 < nodeBegin_.size(); ++i) {
                if (nodeBegin_[i] + nNodeWeights_ > nodeBegin_[i + 1])
                    throw std::runtime_error("Inconsistent node data");
            }
            if (nodeBegin_.front() != 0) throw std::runtime_error("Inconsistent node data begin");
            if (nodeBegin_.back() != (Index)nodeData_.size()) throw std::runtime_error("Inconsistent node data end");

            for (size_t i = 0; i + 1 < hedgeBegin_.size(); ++i) {
                if (hedgeBegin_[i] + nHedgeWeights_ > hedgeBegin_[i + 1])
                    throw std::runtime_error("Inconsistent hedge data");
            }
            if (hedgeBegin_.front() != 0) throw std::runtime_error("Inconsistent hedge data begin");
            if (hedgeBegin_.back() != (Index)hedgeData_.size()) throw std::runtime_error("Inconsistent hedge data end");

            if (nPins_ + nNodeWeights_ * nNodes_ != (Index)nodeData_.size()) throw std::runtime_error("Inconsistent node data size");
            if (nPins_ + nHedgeWeights_ * nHedges_ != (Index)hedgeData_.size()) throw std::runtime_error("Inconsistent hedge data size");
            if (nPartWeights_ * nParts_ != (Index)partData_.size()) throw std::runtime_error("Inconsistent part data size");

            for (Index n = 0; n != nNodes_; ++n) {
                Range<Index> node = nodeHedges(n);
                for (Index hedge : node) {
                    if (hedge < 0 || hedge >= nHedges())
                        throw std::runtime_error("Invalid hedge value");
                }
                std::unordered_set<Index> uq(node.begin(), node.end());
                if (uq.size() != node.size())
                    throw std::runtime_error("Duplicate hedges in a node");
            }
            for (Index h = 0; h != nHedges_; ++h) {
                Range<Index> hedge = hedgeNodes(h);
                for (Index node : hedge) {
                    if (node < 0 || node >= nNodes())
                        throw std::runtime_error("Invalid node value");
                }
                std::unordered_set<Index> uq(hedge.begin(), hedge.end());
                if (uq.size() != hedge.size())
                    throw std::runtime_error("Duplicate nodes in an hedge");
            }
        }

    private:
        void finalize()
        {
            finalizePins();
            finalizeNodes();
            finalizeNodeWeights();
            finalizeHedgeWeights();
            finalizePartWeights();
            //checkConsistency();
        }
        void finalizePins()
        {
            nPins_ = hedgeData_.size() - nHedges_ * nHedgeWeights_;
        }
        void finalizeNodes()
        {
            std::vector<Index> newData;
            std::vector<Index> newBegin;
            newBegin.assign(nNodes_ + 1, 0);

            // Setup node begins : Calculate the number of pins associated with each node.
            for (Index hedge = 0; hedge < nHedges_; ++hedge) {
                for (Index node : hedgeNodes(hedge)) {
                    //std::cout << node<<",";
                    ++newBegin[node + 1];
                }
            }

            // : Set the start location for each node.
            for (Index node = 1; node <= nNodes_; ++node) {
                newBegin[node] += newBegin[node - 1] + nNodeWeights_;
            }
            newData.resize(newBegin.back());
            //std::cout << newData.size() << std::endl;;
            //if(newData.size() != nPins_ + nNodes_ * nNodeWeights_) throw std::runtime_error("Inconsistent node data size");

            // Assign node weights :Go through each node, assigning weights to each node. This step is to store the node weight into the newData vector.
            for (Index node = 0; node < nNodes_; ++node) {
                Index b = newBegin[node];
                for (Index i = 0; i < nNodeWeights_; ++i) {
                    newData[b + i] = nodeWeight(node, i);
                }
            }
            nodeBegin_ = newBegin;

            // Assign pins :Each hyperedge is traversed and the index of the hyperedge is added to the corresponding position of the newData vector for the node associated with each hyperedge. This step is to store the pin into the newData vector
            for (Index hedge = 0; hedge < nHedges_; ++hedge) {
                for (Index node : hedgeNodes(hedge)) {
                    Index ind = --newBegin[node + 1];
                    newData[ind] = hedge;
                }
            }

            nodeData_ = newData;
        }
        void finalizeNodeWeights()
        {
            totalNodeWeights_.assign(nNodeWeights_, 0);
            for (Index i = 0; i < nNodes_; ++i) {
                for (Index j = 0; j < nNodeWeights_; ++j) {
                    totalNodeWeights_[j] += nodeWeight(i, j);    //return node weight , if a node have many weights,through j to index,j default 0
                }
            }
        }
        void finalizeHedgeWeights()
        {
            totalHedgeWeights_.assign(nHedgeWeights_, 0);
            for (Index i = 0; i < nHedges_; ++i) {
                for (Index j = 0; j < nHedgeWeights_; ++j) {
                    totalHedgeWeights_[j] += hedgeWeight(i, j);
                }
            }
        }
        void finalizePartWeights()
        {
            totalPartWeights_.assign(nPartWeights_, 0);
            for (Index i = 0; i < nParts_; ++i) {
                for (Index j = 0; j < nPartWeights_; ++j) {
                    totalPartWeights_[j] += partWeight(i, j);
                }
            }
        }

        bool cut(const hgsol::Solution& solution, Index hedge) const
        {
            std::unordered_set<Index> parts;
            for (Index node : hedgeNodes(hedge)) {
                parts.insert(solution[node]);
            }
            return parts.size() > 1;
        }
        Index degree(const hgsol::Solution& solution, Index hedge) const
        {
            std::unordered_set<Index> parts;
            for (Index node : hedgeNodes(hedge)) {
                parts.insert(solution[node]);
            }
            return parts.size();
        }


    private:
        // Basic stats
        Index nNodes_;
        Index nHedges_;
        Index nParts_;
        Index nPins_;

        // Number of resources for each of them
        Index nNodeWeights_;
        Index nHedgeWeights_;
        Index nPartWeights_;

        // Begin/end in the compressed representation
        std::vector<Index> nodeBegin_;
        std::vector<Index> hedgeBegin_;

        // Data, with the weights then the pins for each node/edge
        std::vector<Index> nodeData_;
        std::vector<Index> hedgeData_;
        std::vector<Index> partData_;

        // Summary stats
        std::vector<Index> totalNodeWeights_;
        std::vector<Index> totalHedgeWeights_;
        std::vector<Index> totalPartWeights_;
    };
}



