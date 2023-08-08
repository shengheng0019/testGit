// Copyright (C) 2019 Gabriel Gouvine - All Rights Reserved

#pragma once

#include "common.hpp"
#include "incremental_objective.hpp"
#include <memory>

using std::make_unique;

namespace MiniPartObj
{
	/**
 * Base class for objectives
 */
	class Objective {
	public:
		virtual std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph&, hgsol::Solution&) const = 0;
		virtual std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph&, hgsol::Solution&) const = 0;
		virtual ~Objective() {}
	};

	class CutObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalCut>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsSumOverflow(s), h.metricsCut(s), h.metricsConnectivity(s) };
		}
	};

	class SoedObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalSoed>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsSumOverflow(s), h.metricsConnectivity(s) };
		}
	};

	class MaxDegreeObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalMaxDegree>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsSumOverflow(s), h.metricsMaxDegree(s), h.metricsConnectivity(s) };
		}
	};

	class DaisyChainDistanceObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalDaisyChainDistance>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsSumOverflow(s), h.metricsDaisyChainDistance(s), h.metricsConnectivity(s) };
		}
	};

	class DaisyChainMaxDegreeObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalDaisyChainMaxDegree>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsSumOverflow(s), h.metricsDaisyChainMaxDegree(s), h.metricsDaisyChainDistance(s) };
		}
	};

	class RatioCutObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalRatioCut>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsEmptyPartitions(s), (int64_t)(100.0 * h.metricsRatioCut(s)), h.metricsCut(s), h.metricsConnectivity(s) };
		}
	};

	class RatioSoedObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalRatioSoed>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsEmptyPartitions(s), (int64_t)(100.0 * h.metricsRatioSoed(s)), h.metricsConnectivity(s) };
		}
	};

	class RatioMaxDegreeObjective final : public Objective {
	public:
		std::unique_ptr<IncrementalObjective> incremental(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return make_unique<IncrementalRatioMaxDegree>(h, s);
		}
		std::vector<int64_t> eval(const MiniPartHypergraph::Hypergraph& h, hgsol::Solution& s) const override
		{
			return { h.metricsEmptyPartitions(s), (int64_t)(100.0 * h.metricsRatioMaxDegree(s)), h.metricsConnectivity(s) };
		}
	};
}





