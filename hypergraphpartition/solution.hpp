#pragma once
#include <iosfwd>
#include <string>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "common.hpp"
#include <unordered_map>
#include <set>

namespace hgsol
{
	class Solution {
	public:
		Solution(Index nNodes, Index nParts) : parts_(nNodes, 0)
			, nParts_(nParts) {}

		Solution(std::vector<Index> parts)
		{
			parts_ = move(parts);
			if (!parts_.empty())
				nParts_ = *max_element(parts_.begin(), parts_.end()) + 1;
			else
				nParts_ = 1;
		}
		Index nNodes() const { return parts_.size(); }
		Index nParts() const { return nParts_; }

		Index  operator[](Index node) const { return parts_[node]; }
		Index& operator[](Index node) { return parts_[node]; }

		Solution coarsen(const Solution& coarsening) const
		{
			assert(coarsening.nNodes() == nNodes());
			Solution ret(coarsening.nParts(), nParts());
			for (Index node = 0; node < nNodes(); ++node) {
				ret[coarsening[node]] = (*this)[node];
			}
			return ret;
		}
		Solution uncoarsen(const Solution& coarsening) const
		{
			assert(coarsening.nParts() == nNodes());
			Solution ret(coarsening.nNodes(), nParts());
			for (Index node = 0; node < coarsening.nNodes(); ++node) {
				ret[node] = (*this)[coarsening[node]];
			}
			return ret;
		}

		void resizeParts(Index parts)
		{
			if (parts < nParts_)
				throw std::runtime_error("It is only possible to increase the number of blocks");
			nParts_ = parts;
		}

		void checkConsistency() const
		{
			for (Index p : parts_) {
				if (p < 0)
					throw std::runtime_error("Block numbers must be non-negative");
				if (p >= nParts())
					throw std::runtime_error("Block numbers must be smaller than the number of blocks");
			}
		}

		std::unordered_map<int, std::set<uint64_t>> ReParts(int k)
		{
			std::unordered_map<int, std::set<uint64_t>> Remap;
			for (int i = 0; i < parts_.size(); ++i)
			{
				auto p = parts_[i];
				Remap[p].emplace(i);
			}
			return Remap;
		}

	private:
		std::vector<Index> parts_;
		//nodes: 0 1 2 3 4 5 6 7 8 9 10 11 12 13
		//part : 0 0 0 0 0 1 1 1 1 1  2  2  2  2
		Index nParts_;
	};
}
