#include "BicliquePartition.h"
#include "graph.h"
#include<iostream>
#include<vector>
#include "gurobi_c++.h"

using namespace std;

void bicliquePartition() {
	int n, ne, i,j;
	vector<vector<int>> graph=makeAdjecencyMatrix(n,ne);
	vector<pair<vector<int>, vector<int>>> bicliques = makeAllBicliques(graph, n);
	int nb = bicliques.size();
	int type = 0;
	/*for (auto it : bicliques) {
		cout << "***********************" << endl;
		for (auto v : it.first) cout << v << " ";
		cout << endl;
		for (auto v : it.second) cout << v << " ";
		cout << endl;
	}*/
	int in = 0;
	while (in < 3) {
		cout << "*********************** BICLIQUE PARTITION ****************************" << endl;
		cout << "0)biclique partition\n1)1,2-biclique partition\n2)1,3-biclique partition\n3)exit\n";
		cin >> in;
		try {

			GRBEnv env = GRBEnv(true);
			env.start();
			GRBModel model = GRBModel(env);

			vector<GRBVar> B; // for each b belongs to B, b represents the presence of a biclique

			for (i = 0;i < nb;i++) B.push_back(model.addVar(0, 1, 0, GRB_BINARY, "b" + to_string(i)));

			vector<vector<GRBLinExpr>> constraints(n, vector<GRBLinExpr>(n, 0));
			GRBLinExpr objective = 0;
			for (i = 0;i < nb;i++) {
				objective += B[i];
				for (auto b1 : bicliques[i].first) {
					for (auto b2 : bicliques[i].second) {
						constraints[b1][b2] += B[i];
						constraints[b2][b1] += B[i];
					}
				}
			}
			for (i = 0;i < n;i++) {
				for (j = i;j < n;j++) {
					if (graph[i][j]) {
						switch (in) {
						case PARTITION_1:
							model.addConstr(constraints[i][j] == 1);
							break;
						case PARTITION_12:
							model.addConstr(constraints[i][j] >= 1);
							model.addConstr(constraints[i][j] <= 2);
							break;
						case PARTITION_13: // contraints[i][j]==2x+1
							GRBVar x = model.addVar(0, 1, 0, GRB_BINARY, "x" + to_string(i) + to_string(j));
							model.addConstr(constraints[i][j] == 2*x + 1);
							break;
						}
					}
				}
			}
			model.setObjective(objective, GRB_MINIMIZE);
			model.optimize();

			cout << "Size of biclique Partition:" << model.get(GRB_DoubleAttr_ObjVal) << endl;
			cout << "Biclique Partition:" << endl;
			for (i = 0;i < nb;i++) {
				if (B[i].get(GRB_DoubleAttr_X) == 1) {
					cout << "**************************************" << endl;
					for (auto v : bicliques[i].first) cout << v << " ";
					cout << endl;
					for (auto v : bicliques[i].second) cout << v << " ";
					cout << endl;
				}
			}

		}
		catch (GRBException& e) {
			std::cerr << "Error code = " << e.getErrorCode() << "\n";
			std::cerr << e.getMessage() << "\n";
		}
		catch (...) {
			std::cerr << "Exception during optimization" << "\n";
		}
	}
}