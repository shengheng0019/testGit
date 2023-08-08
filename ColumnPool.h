#pragma once
#include"Solution.h"
#include<unordered_set>
#include"Vars.h"

//using namespace std;

class Column : public Solution
{
public:
	bool isAdded;
	Column(const Solution& sol);
	Column(const Column& col);
	Column() = default;
	Column& operator = (const Column& column);
};

bool operator==(const Column& col1, const Column& col2);

void hashCombine(size_t seed, int const& v);

template<>
struct std::hash<Column>
{
	size_t operator()(const Column& col)const noexcept;
};

class ColumnPool
{
public:
	//unordered_set<Column> columns;
	//void AddCol(const Column& col);
	//void AddCol(const Solution& sol);
	//bool DelCol(const Column& col);
    int var_num;
	std::vector<int>vars_id;
	std::unordered_set<Column> columns;
	bool AddCol(const Column& col);
    void AddCol(const Solution& sol);
	void AddColumn(const Column& col);         //@wqy using
	bool DelCol(const Column& col);
	ColumnPool& operator = (const ColumnPool& columnpool);
	ColumnPool();
	~ColumnPool() = default;
	ColumnPool(const ColumnPool& columnpool);

	//function
	bool AddColumn(const std::unordered_map<int, double>& column);
	void PushBackColumn(const Column& column) { columns_.push_back(column); }
	double FindColumnValue(int colNo, int varNo)const { return  columns_[colNo].varValue.find(varNo)->second; }

	//inine
	std::vector<Column> Columns()const { return columns_; }

	std::vector<Column> columns_; //vector has the order and I produce the checkRepeat function @getting value with inline function is too slow

private:
};

