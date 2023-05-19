#include "ColumnPool.h"

ColumnPool::ColumnPool()
{
    
}

ColumnPool::~ColumnPool()
{

}
Column::Column(const Solution& sol)
    :Solution(sol), isAdded{ false }
{
    ;
}

Column::Column(const Column& col)
    :Solution{ col }, isAdded{ false }
{
    ;
}
Column::Column()
{

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

    if (col1.varValue.size() != col2.varValue.size())
        return false;
    //for (auto iter1 = col1.varValue.begin(), iter2 = col2.varValue.begin(); iter1 != col1.varValue.end(); iter1++, iter2++)
    //{
    //    if (fabs(iter1->second - iter2->second)>0.0001)
    //    {
    //        return false;
    //    }
    //}
    for (auto iter1 = col1.varValue.begin(); iter1 != col1.varValue.end(); iter1++)
    {
        if (fabs(col2.varValue.find(iter1->first)->second - iter1->second) > 0.0001)
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

size_t hash<Column>::operator()(const Column& col) const
{
    size_t seed = 0;
    for (auto iter = col.varValue.begin(); iter != col.varValue.end(); iter++)
    {
        hashCombine(seed, iter->second);
    }
    return seed;
}

void ColumnPool::AddCol(const Column& col)
{
    columns.insert(col);
}
bool ColumnPool::DelCol(const Column& col)
{

}

ColumnPool& ColumnPool::operator = (const ColumnPool& columnpool)
{
    if (this == &columnpool)
        return *this;
    this->columns = columnpool.columns;
    //this->var_contained = columnpool.var_contained;
    //this->var_num = columnpool.var_num;
    this->vars_id = columnpool.vars_id;
    this->var_num = columnpool.var_num;

    return *this;

}
