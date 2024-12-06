#ifndef MASTER_PROBLEM_H
#define MASTER_PROBLEM_H

#include <ilcplex/ilocplex.h>

#include "Data.h"
using namespace std;

class MasterProblem {
	IloCplex rmp;
	IloNumVarArray lambda;
	IloRangeArray partition_constraint;
	IloObjective master_objective;
	int lambdaCounter;
	Data* data;
	double best = 10000000000;
	double current;

   public:
	IloEnv env;
	MasterProblem(Data* data, IloEnv env);
	void setForbiddenLambdas(vector<int>* forbidden = NULL);
	void unsetForbiddenLambdas(vector<int>* forbidden = NULL);

	void createModel();
	void addLambda(IloNumArray& column);
	IloNumArray getDuals();
	double solve();
	double getLambdaValue(int counter) {
		double value = rmp.getValue(lambda[counter]);

		return value;
	}
	double getBest() { return best; }
	void setBest(double newBest) { best = newBest; }
	int getLambdaCounter() { return lambdaCounter; }
	double getCurrent() { return current; }
};

#endif  // !MASTER_PROBLEM_H