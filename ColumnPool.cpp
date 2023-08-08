#include "ColumnPool.h"
#include <cmath>

ColumnPool::ColumnPool()
{
    var_num = 0;
}

Column::Column(const Solution& sol)
    :Solution(sol) , isAdded{false}
{
    ;
}

Column::Column(const Column& col)
    :Solution{col}, isAdded{ false }
{
    ;
}

Column& Column::operator = (const Column& column)
{
    if (this == &column)
        return *this;
    this->isAdded = column.isAdded;
    this->objValue = column.objValue;
    this->solSize = column.solSize;
    this->varValue = column.varValue;
    return *this;

}

bool operator==(const Column& col1, const Column& col2)
{
    
    /*if (col1.varValue.size() != col2.varValue.size())
        return false;
    for (auto iter1 = col1.varValue.begin(); iter1 != col1.varValue.end(); iter1++)
    {
        if (fabs(iter1->second - col2.varValue.find(iter1->first)->second) > 0.000001)
        {
            return false;
        }
    }
    return true;*/
    if (col1.varValue.size() != col2.varValue.size())
        return false;
    for (auto iter1 = col1.varValue.begin(); iter1 != col1.varValue.end(); iter1++)
    {
        if (std::fabs(col2.varValue.find(iter1->first)->second - iter1->second) > 0.0001)
        {
            return false;
        }
    }
    return true;
}

void hashCombine(size_t seed, int const& v)
{
	seed ^= std::hash<int>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

size_t hash<Column>::operator()(const Column& col) const noexcept
{
    size_t seed = 0;
    for (auto iter = col.varValue.begin(); iter != col.varValue.end(); iter++)
    {
        hashCombine(seed, iter->second);
    }
    return seed;
}

bool ColumnPool::AddCol(const Column& col)
{
    return columns.insert(col).second;
}

void ColumnPool::AddColumn(const Column& col)
{
    columns_.push_back(col);
}

void ColumnPool::AddCol(const Solution& sol)
{
    columns.insert(sol);
}

bool ColumnPool::DelCol(const Column& col)
{
    if(columns.find(col) != columns.end())
        columns.erase(col);
    return true;
}

ColumnPool& ColumnPool::operator = (const ColumnPool& columnpool)
{
    if (this == &columnpool)
        return *this;
    this->columns = columnpool.columns;
    this->vars_id = columnpool.vars_id;
    this->var_num = columnpool.var_num;
    this->columns_ = columnpool.columns_;

    return *this;

}

ColumnPool::ColumnPool(const ColumnPool& columnpool)
{
    columns = columnpool.columns;
    vars_id = columnpool.vars_id;
    var_num = columnpool.var_num;
    columns_ = columnpool.columns_;
}
