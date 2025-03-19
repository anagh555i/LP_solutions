#include "IndependentSet.h"
#include "graph.h"
#include<iostream>
#include<vector>
#include "gurobi_c++.h"

using namespace std;

void IndependentSet() {
	int n, ne,i;
	vector<vector<int>> graph = makeEdgeSet(n, ne);

	try {

		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		vector<GRBVar> X; // for each vertex x=1 => presence of vertex

		for (i = 0;i < n;i++) X.push_back(model.addVar(0, 1, 0, GRB_BINARY, "X" + to_string(i)));

		for (auto e : graph) { //
			GRBLinExpr constraint = X[e[0]] + X[e[1]];
			model.addConstr(constraint <= 1);
		}

		GRBLinExpr objective = 0;
		for (i = 0;i < n;i++) objective += X[i];
		model.setObjective(objective, GRB_MAXIMIZE);

		model.optimize();

		int fin = 0;
		cout << "Independent Set:";
		for (i = 0;i < n;i++) {
			if (X[i].get(GRB_DoubleAttr_X) == 1) {
				cout << i << " ";
				fin++;
			}
		}
		cout << endl;
		cout << "size of Independent set: " << fin << endl;

	}
	catch (GRBException& e) {
		std::cerr << "Error code = " << e.getErrorCode() << "\n";
		std::cerr << e.getMessage() << "\n";
	}
	catch (...) {
		std::cerr << "Exception during optimization" << "\n";
	}
}