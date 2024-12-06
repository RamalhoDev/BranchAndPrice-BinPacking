#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "BranchAndPrice.h"
#include "Data.h"

using namespace std;

int main(int argc, char** argv) {
	Data data;
	data.readData(argv[1]);

	IloEnv env;

	auto mp = MasterProblem(&data, env);
	mp.createModel();
	vector<BranchAndPrice*> bp = {new BranchAndPrice(&data, &mp)};

	while (!bp.empty()) {
		auto currBP = bp.back();
		bp.pop_back();
		auto nodes = currBP->Solve();
		if (nodes.first) {
			bp.push_back(nodes.first);
			bp.push_back(nodes.second);
		}

		delete currBP;
	}
	cout << mp.getBest() << endl;
	env.end();

	// itens 1 and 2 are together only on columns 5 and 11
	// lambda[best.first].setUB(0.0);
	// lambda[best.second].setUB(0.0);

	// // to allow them again:
	// lambda[best.first].setUB(1);
	// lambda[best.second].setUB(1);

	return 0;
}