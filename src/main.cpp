#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "BranchAndPrice.h"
#include "Data.h"

using namespace std;

int main(int argc, char** argv) {
	const double M = 1e6;
	vector<int> weight = {2, 1, 3, 3, 5};
	int capacity = 7;
	int n = weight.size();

	if (argc != 2) {
		printf("Usage:\n./bin instance\n");
		return 0;
	}

	Data data;
	data.readData(argv[1]);

	vector<BranchAndPrice*> bp = {new BranchAndPrice(&data)};

	while (bp.empty()) {
		auto currBP = bp.back();
		bp.pop_back();
		auto nodes = currBP->Solve();
		if (nodes.first) {
			bp.push_back(nodes.first);
			bp.push_back(nodes.second);
		}

		delete nodes.first;
		delete nodes.second;
	}

	// itens 1 and 2 are together only on columns 5 and 11
	// lambda[best.first].setUB(0.0);
	// lambda[best.second].setUB(0.0);

	// // to allow them again:
	// lambda[best.first].setUB(1);
	// lambda[best.second].setUB(1);

	return 0;
}