#pragma once

#include <future>
#include <unordered_set>
#include <vector>

namespace part {

class Partition
{
public:
    //make partitions move only
    //Partition(Partition&&) = default;
    //Partition() = delete;
    //auto operator=(Partition &&)
    //    -> Partition& = default;
    //auto operator=(const Partition&)
    //    -> Partition& = delete;

    Partition(std::size_t);

    //add node and its hyperedges to the partition
    auto addNode(int64_t node,
                 const std::vector<int64_t>& edges)
        -> void;

    //same as above
    auto addNode(int64_t node,
                 const std::unordered_set<int64_t>& edges)
        -> void;

    //checks if the partitions holds a vertex
    //which is connected to the given edge
    auto hasEdge(int64_t edge) const
        -> bool;
    //return reference to the node set
    auto getNodes() const
        -> const std::unordered_set<int64_t>&;
    //same as above
    auto getNodes()
        -> std::unordered_set<int64_t>&;

    //return reference to the edge set
    auto getEdges() const
        -> const std::unordered_set<int64_t>&;
    //same as above
    auto getEdges()
        -> std::unordered_set<int64_t>&;

    //clears edge and node set
    auto clear()
        -> void;

    //returns the id of the partition
    auto getId() const
        -> std::size_t;

    //returns the number of hyperedges in the partition
    auto numberOfEdges() const
        -> std::size_t;

    //return the number of nodes in the partition
    auto numberOfNodes() const
        -> std::size_t;

    //returns a future which will hold the number of external degrees
    //of the partition to the other given partition
    auto externalDegree(const std::vector<Partition>& parts) const
        -> std::future<std::size_t>;

    auto toString() const
        -> std::string;

private:
    std::size_t _id; //id of the partition
    std::unordered_set<int64_t> _nodes; //nodes in partition
    std::unordered_set<int64_t> _edges; //edges in partition
};

} // namespace part


//Partition(std::size_t)：构造函数，接受一个std::size_t类型的参数，用于指定该分区的ID。
//void addNode(int64_t node, const std::vector<int64_t>& edges)：将节点及其超边添加到分区中，超边以向量的形式提供。
//void addNode(int64_t node, const std::unordered_set<int64_t>& edges)：同上，超边以无序集合的形式提供。
//bool hasEdge(int64_t edge) const：检查分区是否持有连接到给定超边的顶点。
//const std::unordered_set<int64_t>& getNodes() const：返回节点集合的常量引用。
//std::unordered_set<int64_t>& getNodes()：返回节点集合的非常量引用。
//const std::unordered_set<int64_t>& getEdges() const：返回超边集合的常量引用。
//std::unordered_set<int64_t>& getEdges()：返回超边集合的非常量引用。
//void clear()：清空节点集合和超边集合。
//std::size_t getId() const：返回分区的ID。
//std::size_t numberOfEdges() const：返回分区中超边的数量。
//std::size_t numberOfNodes() const：返回分区中节点的数量。
//std::future<std::size_t> externalDegree(const std::vector<Partition>& parts) const：返回一个std::future对象，该对象将在未来持有该分区与另一个给定分区之间的外部度数。
//std::string toString() const：返回表示该分区对象的字符串。
//该类还具有一个构造函数和一个默认构造函数，但没有实现，因此不能使用。成员变量有一个分区ID _id，以及节点和超边的无序集合。