#pragma once
#include"Solution.h"
#include<unordered_set>
#include"Vars.h"

using namespace std;

class Column : public Solution
{
public:
	bool isAdded;
	Column(const Solution& sol);
	Column(const Column& col);
	Column();
	Column& operator = (const Column& column);
};

bool operator==(const Column& col1, const Column& col2);

void hashCombine(size_t seed, int const& v);

template<>
struct hash<Column>
{
	size_t operator()(const Column& col) const;
};

class ColumnPool
{
public:
	//Vars* var_contained;//这个子问题所包含的变量
	int var_num;//这个子问题包含的变量个数(每个约束变量的并集)
	vector<int>vars_id;//这个子问题包含的变量的序号的并集（可以无序，用的时候自己排）
	unordered_set<Column> columns;
	void AddCol(const Column& col);
	bool DelCol(const Column& col);
	ColumnPool& operator = (const ColumnPool& columnpool);
	ColumnPool();
	~ColumnPool();
};