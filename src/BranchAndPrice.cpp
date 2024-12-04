#include "BranchAndPrice.h"

BranchAndPrice::BranchAndPrice(Data* data, vector<pair<int, int>> forbidden, vector<pair<int, int>> together) {
	this->data = data;
	this->forbidden = forbidden;
	this->together = together;
	masterProblem();
}

void BranchAndPrice::masterProblem() {
	const double M = 1e6;
	int n = data->getNItems();

	IloModel master_model(env);

	lambda = IloNumVarArray(env, n, 0, IloInfinity);

	IloExpr sum_obj(env);
	partition_constraint = IloRangeArray(env);

	for (int i = 0; i < n; i++) {
		char var_name[50];
		sprintf(var_name, "y%d", i);

		lambda[i].setName(var_name);
		sum_obj += M * lambda[i];

		partition_constraint.add(lambda[i] == 1);
	}

	master_model.add(partition_constraint);

	IloObjective master_objective = IloMinimize(env, sum_obj);
	master_model.add(master_objective);

	for (size_t i = 0; i < forbidden.size(); i++) {
		auto nodesForb = forbidden[i];
		lambda[nodesForb.first].setUB(0.0);
		lambda[nodesForb.second].setUB(0.0);
	}

	rmp = IloCplex(master_model);

	rmp.setOut(env.getNullStream());  // disables CPLEX log
}

pair<BranchAndPrice*, BranchAndPrice*> BranchAndPrice::Solve() {
	int n = data->getNItems();
	vector<vector<double>> z = vector<vector<double>>(n, vector<double>(n, 0));
	rmp.solve();

	cout << "Initial lower bound: " << rmp.getObjValue() << endl;

	cout << "Initial solution: ";
	for (size_t j = 0; j < lambda.getSize(); j++) {
		cout << rmp.getValue(lambda[j]) << " ";
	}
	cout << endl;

	int lambda_counter = n;
	auto lambdaItems = vector<vector<bool>>(n, vector<bool>());
	while (true) {
		// Get the dual variables
		IloNumArray pi(env, n);

		IloRangeArray partition_constraint(env);
		rmp.getDuals(pi, partition_constraint);

		for (size_t i = 0; i < n; i++) {
			cout << "Dual variable of constraint " << i << " = " << pi[i] << endl;
		}

		// Build and solve the pricing problem

		IloModel pricing_model(env);

		IloNumVarArray x(env, n, 0, 1, IloNumVar::Bool);
		IloExpr sum_obj(env, 1);
		IloRange r;
		IloExpr sumX(env);

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
			IloExpr forbid(env);
			forbid += x[forbidden[i].first] + x[forbidden[i].second];
			char constr[30];
			sprintf(constr, "forbid_%d_%d", forbidden[i].first, forbidden[i].second);
			forbid_ctr = (forbid <= 1);
			forbid_ctr.setName(name);

			pricing_model.add(forbid_ctr);
		}

		for (size_t i = 0; i < this->together.size(); i++) {
			IloRange tog_ctr;
			IloExpr tog(env);
			tog += x[together[i].first] + x[together[i].second];
			char constr[30];
			sprintf(constr, "together_%d_%d", together[i].first, together[i].second);
			tog_ctr = (tog <= 2);
			tog_ctr.setName(name);

			pricing_model.add(tog_ctr);
		}

		IloObjective objective = IloMinimize(env, sum_obj);
		pricing_model.add(objective);

		IloCplex pricing_problem(pricing_model);

		pricing_problem.solve();
		if (pricing_problem.getObjValue() < -1e-5) {
			cout << "Reduced cost is equal to " << pricing_problem.getObjValue() << ", which is less than 0..." << endl;

			IloNumArray entering_col(env, n);

			pricing_problem.getValues(x, entering_col);
			vector<bool> items = vector<bool>();
			cout << endl << "Entering column:" << endl;
			for (size_t i = 0; i < n; i++) {
				bool insertItem = (entering_col[i] < 0.5 ? false : true);
				items.push_back(insertItem);
				cout << entering_col[i] << " " << insertItem << endl;
			}
			cout << endl;
			lambdaItems.push_back(items);

			// Add the column to the master problem
			char var_name[50];
			sprintf(var_name, "y%d", lambda_counter++);
			IloNumVar new_lambda(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
			new_lambda.setName(var_name);

			lambda.add(new_lambda);

			cout << "Solving the RMP again..." << endl;

			rmp.solve();
		} else {
			cout << "No column with negative reduced costs found. The current basis is optimal" << endl;
			cout << "Final master problem: " << endl;
			system("cat model.lp");
			break;
		}
	}

	double min = 1;
	pair<int, int> best = make_pair(-1, -1);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			if (lambdaItems[lambda_counter][i] && lambdaItems[lambda_counter][j]) {
				z[i][j] += rmp.getValue(lambda[lambda_counter]);
				z[j][i] += rmp.getValue(lambda[lambda_counter]);
			}

			double frac = abs(z[i][j] - 0.5);
			if (frac < min) {
				min = frac;
				best = make_pair(i, j);
			}
		}
	}

	for (size_t i = 0; i < forbidden.size(); i++) {
		auto nodesForb = forbidden[i];
		lambda[nodesForb.first].setUB(1.0);
		lambda[nodesForb.second].setUB(1.0);
	}

	pair<BranchAndPrice*, BranchAndPrice*> bestBP;
	if (best.first)
		bestBP = make_pair(new BranchAndPrice(this, best.first, best.second, true), new BranchAndPrice(this, best.first, best.second, false));

	return bestBP;
}

BranchAndPrice::BranchAndPrice(BranchAndPrice* bp, int a, int b, bool forbid) {
	this->forbidden = bp->forbidden;
	this->together = bp->together;

	this->data = bp->data;

	if (forbid)
		this->forbidden.push_back(make_pair(a, b));
	else
		this->together.push_back(make_pair(a, b));
	masterProblem();
}