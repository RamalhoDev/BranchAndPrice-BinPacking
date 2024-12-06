#ifndef BRANCH_AND_PRICE_H
#define BRANCH_AND_PRICE_H

#include <iostream>
#include <utility>
#include <vector>

#include "MasterProblem.h"

using namespace std;

class BranchAndPrice {
	Data* data;
	vector<pair<int, int>> forbidden;
	vector<pair<int, int>> together;
	vector<int> lambdasToForbid;
	MasterProblem* mp;
	vector<vector<bool>>* lambdaItems;

   public:
	BranchAndPrice(Data* data, MasterProblem* mp);
	BranchAndPrice(BranchAndPrice* bp, int a, int b, bool forbid = true);
	pair<BranchAndPrice*, BranchAndPrice*> Solve();
};
#endif  // !BRANCH_AND_PRICE_H