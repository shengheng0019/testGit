#include "Hypergraph.hpp"
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>


auto part::Hypergraph::addVertex(int64_t id)
    -> std::pair<part::Hypergraph::VertexMapIter, bool>
{
    if(auto iter = _vertices.find(id);
       iter != _vertices.end()) {
        return {iter, false};
    } else {
        return _vertices.insert({id, {}});
    }
}

auto part::Hypergraph::addEdge(int64_t id)
    -> std::pair<part::Hypergraph::EdgeMapIter, bool>
{
    if(auto iter = _edges.find(id);
       iter != _edges.end()) {
        return {iter, false};
    } else {
        return _edges.insert({id, {}});
    }
}

auto part::Hypergraph::addEdgeList(const int64_t& vtx,
                                   const std::vector<int64_t>& edge_list)
    -> void
{
    for(auto&& edge : edge_list) {
        connect(vtx, edge);
    }
}


auto part::Hypergraph::addNodeList(const int64_t& edge,
                                   const std::vector<int64_t>& node_list)
    -> void
{
    for(auto&& node : node_list) {
        connect(node, edge);
    }
}

auto part::Hypergraph::connect(const int64_t& vertex,
                               const int64_t& edge) -> void
{
    auto [vertex_iter, dummy1] = addVertex(vertex);
    auto [edge_iter, dummy2] = addEdge(edge);

    edge_iter->second.insert(vertex);
    vertex_iter->second.insert(edge);
}

auto part::Hypergraph::getNodeHeuristicExactly(const int64_t& vtx) const
    -> double
{
    const auto& edges = getEdgesOf(vtx);

    //we need this to not divide by zero later
    if(edges.empty())
        return 0;

    return std::accumulate(std::begin(edges),
                           std::end(edges),
                           std::size_t{0},
                           [this](auto init, auto edge) {
                               return init + getVerticesOf(edge).size() - 1;
                           })
        / edges.size();
}

auto part::Hypergraph::getNodeHeuristicEstimate(const int64_t& vtx) const
    -> double
{
    if(auto iter = _neigbour_map.find(vtx);
       iter != _neigbour_map.end()) {
        return iter->second;
    }

    auto neigs = getNodeHeuristicExactly(vtx);

    _neigbour_map.insert({vtx, neigs});
    return neigs;
}

namespace {

//adds elements of @param from into @param to
//as long as @param to does not have more than @param upto
//elements. Only elements for which the given predicate @param pred is true
template<class From, class To, class Predicate>
auto addAtMostN(From&& from,
                std::size_t upto,
                To&& to,
                Predicate&& pred)
{
    //set end iterator based on wether size(to) + size(from)
    //is bigger than allowed
    auto end = from.begin();
    if(to.size() + from.size() > upto)
        std::advance(end, upto - to.size());
    else
        end = from.end();

    //add elements if predicate is true
    auto iter = from.begin();
    while(iter != end) {
        if(pred(*iter)) {
            to.insert(*iter);
        }
        iter++;
    }

    //return new container with at most @param upto
    //elements
    return std::move(to);
}

} // namespace

auto part::Hypergraph::getSSetCandidates(const int64_t& vtx,
                                         std::size_t n,
                                         std::size_t max_edge_size) const
    -> std::unordered_set<int64_t>
{
    //while we dont have reached the limit
    //and we dont have enough sset-candidates found
    //keep searching for them in bigger edges
    const auto& edges = getEdgesOf(vtx);
    std::unordered_set<int64_t> neigbors;
    std::size_t current_max{2};

    while(current_max < max_edge_size
          && neigbors.size() <= n) 
    {
        //TODO:maybe refactor to something faster
        for(auto&& edge : edges) {
            const auto& vtxs = getVerticesOf(edge);

            // if edge is bigger then max_edge_size then skip
            if(vtxs.size() > current_max)
                continue;

            neigbors = addAtMostN(vtxs,
                                  n,
                                  std::move(neigbors), //i like to move it move it
                                  [&vtx](auto&& elem) {
                                      return elem != vtx;
                                  });

            if(neigbors.size() >= n)
                return neigbors;
        }

        current_max *= 2;
    }
    return neigbors;
}

auto part::Hypergraph::getEdgesOf(const int64_t& vtx) const
    -> const std::unordered_set<int64_t>&
{
    if(auto iter = _vertices.find(vtx); iter != _vertices.end()) {
        return iter->second;
    } else {
        //needed to return a reference
        static const std::unordered_set<int64_t> empty_set;
        return empty_set;
    }
}

auto part::Hypergraph::getVerticesOf(const int64_t& edge) const
    -> const std::unordered_set<int64_t>&
{
    if(auto iter = _edges.find(edge); iter != _edges.end()) {
        return iter->second;
    } else {
        //needed to return a reference
        static const std::unordered_set<int64_t> empty_set;
        return empty_set;
    }
}

auto part::Hypergraph::getEdges() const
    -> const EdgeMap&
{
    return _edges;
}

auto part::Hypergraph::getEdges()
    -> EdgeMap&
{
    return _edges;
}

auto part::Hypergraph::getVertices() const
    -> const VertexMap&
{
    return _vertices;
}

auto part::Hypergraph::getVertices()
    -> VertexMap&
{
    return _vertices;
}

auto part::Hypergraph::getEdgesizeOfPercentBiggestEdge(double percent) const
    -> std::size_t
{
    const auto factor = 1 - percent / 100;
    std::vector<std::size_t> size_vec;
    for(auto&& [key, value] : _edges) {
        size_vec.push_back(value.size());
    }

    std::nth_element(size_vec.begin(),
                     size_vec.begin() + (size_vec.size() - 1) * factor,
                     size_vec.end());

    return size_vec[(size_vec.size() - 1) * factor];
}

auto part::Hypergraph::getRandomNode() const
    -> int64_t
{
    static std::mt19937 engine{Hypergraph::random_seed};
    std::uniform_int_distribution<std::size_t> dist(0, _vertices.size() - 1);

    auto iter = std::begin(_vertices);
    std::advance(iter, dist(engine));

    return iter->first;
}

auto part::Hypergraph::getANode() const
    -> int64_t
{
    return _vertices.begin()->first;
}

auto part::Hypergraph::deleteVertex(int64_t vertex)
    -> void
{
    if(auto vertex_iter = _vertices.find(vertex);
       vertex_iter != _vertices.end()) {
        for(auto&& edge : vertex_iter->second) {
            if(auto edge_iter = _edges.find(edge);
               edge_iter != _edges.end()) {
                edge_iter->second.erase(vertex);

                // delete edges without nodes
                if(edge_iter->second.empty()) {
                    _edges.erase(edge_iter);
                }
            }
        }

        _neigbour_map.erase(vertex);
        _vertices.erase(vertex_iter);
    }
}


auto part::Hypergraph::setSeed(uint32_t seed)
    -> void
{
    Hypergraph::random_seed = seed;
}



//addVertex adds a new vertex to the hypergraph with the given ID.
//addEdge adds a new edge to the hypergraph with the given ID.
//addEdgeList connects the given vertex to a list of edges by calling connect on each pair of(vertex, edge) IDs.
//addNodeList connects the given edge to a list of vertices by calling connect on each pair of(vertex, edge) IDs.
//connect adds a connection between the given vertexand edge, which are assumed to already exist in the hypergraph.
//getNodeHeuristicExactly returns a measure of how "connected" the given vertex is to other vertices in the hypergraph.The heuristic is computed by summing the number of vertices connected to each edge that the given vertex is connected to(minus one, to exclude the given vertex), and dividing by the number of edges.
//getNodeHeuristicEstimate returns an estimate of the same measure as getNodeHeuristicExactly, but memoizes the result so that repeated calls are faster.
//addAtMostN is a helper function that adds elements from one container to another as long as a given predicate is true, up to a maximum number of elements.
//getSSetCandidates returns a set of candidate vertices for a given vertex that are "similar" to the vertex, where "similarity" is defined as being connected to the same set of edges(with some additional conditions).The function searches for candidates in increasing order of edge size, up to a maximum edge size.
//getEdgesOf returns the set of edges that the given vertex is connected to.
//getVerticesOf returns the set of vertices that the given edge is connected to.