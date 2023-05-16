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
	//Vars* var_contained;//����������������ı���
	int var_num;//�������������ı�������(ÿ��Լ�������Ĳ���)
	vector<int>vars_id;//�������������ı�������ŵĲ��������������õ�ʱ���Լ��ţ�
	unordered_set<Column> columns;
	void AddCol(const Column& col);
	bool DelCol(const Column& col);
	ColumnPool& operator = (const ColumnPool& columnpool);
	ColumnPool();
	~ColumnPool();
};