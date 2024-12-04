#ifndef BRANCH_AND_PRICE_H
#define BRANCH_AND_PRICE_H

#include <ilcplex/ilocplex.h>

#include <iostream>
#include <utility>
#include <vector>

#include "Data.h"

using namespace std;

class BranchAndPrice {
	Data* data;
	vector<pair<int, int>> forbidden;
	vector<pair<int, int>> together;
	IloEnv env;
	IloCplex rmp;
	IloNumVarArray lambda;
	IloRangeArray partition_constraint;
	IloObjective master_objective;

   public:
	BranchAndPrice(BranchAndPrice* bp, int a, int b, bool forbid = true);
	pair<BranchAndPrice*, BranchAndPrice*> Solve();
	BranchAndPrice(Data* data, vector<pair<int, int>> forbidden = vector<pair<int, int>>(),
	               vector<pair<int, int>> together = vector<pair<int, int>>());
	void masterProblem();
};
#endif  // !BRANCH_AND_PRICE_H