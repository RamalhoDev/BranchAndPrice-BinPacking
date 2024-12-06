#include "BranchAndPrice.h"

BranchAndPrice::BranchAndPrice(Data* data, MasterProblem* mp) {
	this->data = data;
	this->mp = mp;
	this->forbidden = vector<pair<int, int>>();
	this->together = vector<pair<int, int>>();
	this->lambdaItems = new vector<vector<bool>>();
	for (size_t i = 0; i < data->getNItems(); i++) {
		vector<bool> items = vector<bool>();
		for (size_t j = 0; j < data->getNItems(); j++) {
			if (i == j) items.push_back(true);
			items.push_back(false);
		}
		lambdaItems->push_back(items);
	}

	this->lambdasToForbid = vector<int>();
}

pair<BranchAndPrice*, BranchAndPrice*> BranchAndPrice::Solve() {
	int n = data->getNItems();
	vector<vector<double>> z = vector<vector<double>>(n, vector<double>(n, 0));
	mp->setForbiddenLambdas(&lambdasToForbid);

	pair<BranchAndPrice*, BranchAndPrice*> bestBP;
	double cost = mp->solve();

	if (mp->getCurrent() > mp->getBest() || cost < 0) {
		return bestBP;
	}
	cout << "Initial lower bound: " << cost << endl;

	cout << "Initial solution: ";
	for (size_t j = 0; j < mp->getLambdaCounter(); j++) {
		cout << mp->getLambdaValue(j) << " ";
	}
	cout << endl;

	while (true) {
		// Get the dual variables
		IloNumArray pi = mp->getDuals();

		IloModel pricing_model(mp->env);

		IloNumVarArray x(mp->env, n, 0, 1, IloNumVar::Bool);
		IloExpr sum_obj(mp->env, 1);
		IloRange r;
		IloExpr sumX(mp->env);

		for (int i = 0; i < n; i++) {
			char var_name[50];
			sprintf(var_name, "x_%d", i);

			x[i].setName(var_name);
			sum_obj -= pi[i] * x[i];
			sumX += data->getItemWeight(i) * x[i];
		}

		char name[30];
		sprintf(name, "capacity");
		r = (sumX <= data->getBinCapacity());
		r.setName(name);

		pricing_model.add(r);
		for (size_t i = 0; i < this->forbidden.size(); i++) {
			IloRange forbid_ctr;
			IloExpr forbid(mp->env);
			forbid += x[forbidden[i].first] + x[forbidden[i].second];
			char constr[30];
			sprintf(constr, "forbid_%d_%d", forbidden[i].first, forbidden[i].second);
			forbid_ctr = (forbid <= 1);
			forbid_ctr.setName(name);

			pricing_model.add(forbid_ctr);
		}

		for (size_t i = 0; i < this->together.size(); i++) {
			IloRange tog_ctr;
			IloExpr tog(mp->env);
			tog += x[together[i].first] - x[together[i].second];
			char constr[30];
			sprintf(constr, "together_%d_%d", together[i].first, together[i].second);
			tog_ctr = (tog == 0);
			tog_ctr.setName(name);

			pricing_model.add(tog_ctr);
		}

		IloObjective objective = IloMinimize(mp->env, sum_obj);
		pricing_model.add(objective);

		IloCplex pricing_problem(pricing_model);
		pricing_problem.setParam(IloCplex::Param::Threads, 1);
		pricing_problem.setOut(mp->env.getNullStream());
		pricing_problem.solve();
		if (pricing_problem.getStatus() == IloAlgorithm::Infeasible) {
			pricing_problem.end();
			return bestBP;
		}
		if (pricing_problem.getObjValue() < -1e-5) {
			cout << "Reduced cost is equal to " << pricing_problem.getObjValue() << ", which is less than 0..." << endl;

			IloNumArray entering_col(mp->env, n);

			pricing_problem.getValues(x, entering_col);
			vector<bool> items = vector<bool>();
			for (size_t i = 0; i < n; i++) {
				bool insertItem = (entering_col[i] < 0.5 ? false : true);
				items.push_back(insertItem);
			}
			lambdaItems->push_back(items);

			// Add the column to the master problem
			mp->addLambda(entering_col);
			// cout << "Solving the RMP again..." << endl;

			cost = mp->solve();
			pricing_problem.end();
			if (cost < 0) break;
		} else {
			pricing_problem.end();
			break;
		}
	}

	if (mp->getCurrent() > mp->getBest()) {
		return bestBP;
	}

	double min = 10000;
	pair<int, int> best = make_pair(-1, -1);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			for (size_t k = 0; k < mp->getLambdaCounter(); k++) {
				if (lambdaItems->at(k)[i] && lambdaItems->at(k)[j]) {
					double valueLambda = mp->getLambdaValue(k);
					z[i][j] += valueLambda;
					z[j][i] += valueLambda;
				}
			}

			double frac = abs(z[i][j] - 0.5);
			if (frac < min) {
				min = frac;
				best = make_pair(i, j);
			}
		}
	}

	mp->unsetForbiddenLambdas(&lambdasToForbid);
	cout << abs(min - 0.5) << endl;
	if (abs(min - 0.5) <= 0.0000001) {
		if (mp->getCurrent() < mp->getBest()) {
			mp->setBest(mp->getCurrent());
		}
		return bestBP;
	}

	if (best.first != -1)
		bestBP = make_pair(new BranchAndPrice(this, best.first, best.second, true), new BranchAndPrice(this, best.first, best.second, false));

	return bestBP;
}

BranchAndPrice::BranchAndPrice(BranchAndPrice* bp, int a, int b, bool forbid) {
	this->forbidden = bp->forbidden;
	this->together = bp->together;
	this->mp = bp->mp;
	this->data = bp->data;
	this->lambdaItems = bp->lambdaItems;
	this->lambdasToForbid = bp->lambdasToForbid;

	if (forbid) {
		for (size_t i = 0; i < mp->getLambdaCounter(); i++) {
			if (lambdaItems->at(i)[a] && lambdaItems->at(i)[b]) {
				this->lambdasToForbid.push_back(i);
			}
		}

		this->forbidden.push_back(make_pair(a, b));
	} else {
		for (size_t i = 0; i < mp->getLambdaCounter(); i++) {
			if ((lambdaItems->at(i)[a] && !lambdaItems->at(i)[b]) || (!lambdaItems->at(i)[a] && lambdaItems->at(i)[b])) {
				this->lambdasToForbid.push_back(i);
			}
		}
		this->together.push_back(make_pair(a, b));
	}
}