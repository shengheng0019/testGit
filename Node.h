#pragma once
#include "ColumnPool.h"
#include "MasterProb.h"
#include "SubProbManagement.h"



class Node
{
public:

	int index; //��ǰ�ڵ�����

	RMPSOL LPrelax_sol;//lp�ɳڽ�
	double LPrelax_obj;//lp�ɳ�Ŀ��ֵ
	int fixed_sp_no;//�̶����������ĸ�������
	Column fixed_column;//���α��̶�����
	vector<FIX> Fixed_Info;//�̶�����Ϣ�����ɱ�

	int subprob_num;//����������

	MasterProb MasterProblem; //������
	SubProb* SubProblem; //������

	//����һЩԼ���࣬Ŀ���ǿ���ֱ�Ӹ���node����ģ�ͣ��ݹ����ʱֱ���޸�node��Լ����Ϳ���
	Objective obj;
	Constraints* couplingCons;
	int couplingCons_num;
	Constraints* constraints;
	int constraints_num;
	ColumnPool* columnpool;


	Node();
	Node(ColumnPool* columnpool);
	Node(Objective obj, Constraints* couplingCons, int couplingCons_num,
		Constraints* constraints, int constraints_num, ColumnPool* columnpool, int subprob_num);
	Node(const Node& node);
	void free();
	~Node();
	Node& operator = (const Node& node);

private:


};


