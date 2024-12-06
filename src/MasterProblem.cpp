#include "MasterProblem.h"

MasterProblem::MasterProblem(Data* data, IloEnv env) {
	this->data = data;
	this->env = env;
	this->lambdaCounter = data->getNItems();
	createModel();
}

void MasterProblem::setForbiddenLambdas(vector<int>* forbidden) {
	for (size_t i = 0; i < forbidden->size(); i++) {
		auto nodesForb = forbidden->at(i);
		lambda[nodesForb].setUB(0.0);
	}
}

void MasterProblem::unsetForbiddenLambdas(vector<int>* forbidden) {
	for (size_t i = 0; i < forbidden->size(); i++) {
		auto nodesForb = forbidden->at(i);
		lambda[nodesForb].setUB(1.0);
	}
}

IloNumArray MasterProblem::getDuals() {
	IloNumArray pi(env, data->getNItems());

	rmp.getDuals(pi, partition_constraint);
	return pi;
}

void MasterProblem::addLambda(IloNumArray& column) {
	char var_name[50];
	sprintf(var_name, "y%d", lambdaCounter++);
	IloNumVar new_lambda(master_objective(1) + partition_constraint(column), 0, IloInfinity);
	new_lambda.setName(var_name);

	lambda.add(new_lambda);
}

void MasterProblem::createModel() {
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

	master_objective = IloMinimize(env, sum_obj);
	master_model.add(master_objective);

	rmp = IloCplex(master_model);

	rmp.setOut(env.getNullStream());  // disables CPLEX log
}

double MasterProblem::solve() {
	rmp.solve();
	this->current = rmp.getObjValue();
	if (rmp.getCplexStatus() == IloAlgorithm::Infeasible) {
		return -1;
	}

	return current;
}
