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


//Partition(std::size_t)�����캯��������һ��std::size_t���͵Ĳ���������ָ���÷�����ID��
//void addNode(int64_t node, const std::vector<int64_t>& edges)�����ڵ㼰�䳬����ӵ������У���������������ʽ�ṩ��
//void addNode(int64_t node, const std::unordered_set<int64_t>& edges)��ͬ�ϣ����������򼯺ϵ���ʽ�ṩ��
//bool hasEdge(int64_t edge) const���������Ƿ�������ӵ��������ߵĶ��㡣
//const std::unordered_set<int64_t>& getNodes() const�����ؽڵ㼯�ϵĳ������á�
//std::unordered_set<int64_t>& getNodes()�����ؽڵ㼯�ϵķǳ������á�
//const std::unordered_set<int64_t>& getEdges() const�����س��߼��ϵĳ������á�
//std::unordered_set<int64_t>& getEdges()�����س��߼��ϵķǳ������á�
//void clear()����սڵ㼯�Ϻͳ��߼��ϡ�
//std::size_t getId() const�����ط�����ID��
//std::size_t numberOfEdges() const�����ط����г��ߵ�������
//std::size_t numberOfNodes() const�����ط����нڵ��������
//std::future<std::size_t> externalDegree(const std::vector<Partition>& parts) const������һ��std::future���󣬸ö�����δ�����и÷�������һ����������֮����ⲿ������
//std::string toString() const�����ر�ʾ�÷���������ַ�����
//���໹����һ�����캯����һ��Ĭ�Ϲ��캯������û��ʵ�֣���˲���ʹ�á���Ա������һ������ID _id���Լ��ڵ�ͳ��ߵ����򼯺ϡ�