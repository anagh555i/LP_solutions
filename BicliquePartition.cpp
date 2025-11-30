#include "BicliquePartition.h"
#include "graph.h"
#include<iostream>
#include<vector>
#include "gurobi_c++.h"

//using namespace std;

class MyCallback : public GRBCallback {
public:
	vector<GRBVar> vars;
	vector<pair<vector<int>, vector<int>>> bicliques;

	MyCallback(const vector<GRBVar>& v, const vector<pair<vector<int>, vector<int>>> &b) : vars(v), bicliques(b) {}

protected:
	void callback() override {
		if (where == GRB_CB_MIPSOL) {
			double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
			std::cout << "Found feasible solution, obj = " << obj << "\n";

			// get solution values
			std::vector<double> sol(vars.size());
			// the overload returns a heap-allocated array; or use getSolution one-by-one
			double* x = getSolution(vars.data(), (int)vars.size());
			for (int i = 0; i < vars.size(); i++) {
				double val = x[i];
				if (std::abs(val) > 1e-9) {  // filter out 0
					//std::cout << "var[" << i << "] = " << val << std::endl;
					cout << "**************************************" << endl;
					for (auto v : bicliques[i].first) cout << v << " ";
					cout << endl;
					for (auto v : bicliques[i].second) cout << v << " ";
					cout << endl;
				}
			}
			delete[] x;  // free the array returned by getSolution

			std::cout << "-----\n";
		}
	}
};

class MyCallbackV2 : public GRBCallback {
public:
	vector<GRBVar> xc;
	vector<vector<GRBVar>> ac;
	vector<vector<GRBVar>> bc;
	int nc;

	MyCallbackV2(const vector<GRBVar>& x, const vector<vector<GRBVar>> &a, const vector<vector<GRBVar>>& b,const int n) : ac(a), bc(b), xc(x), nc(n) {}

protected:
	void callback() override {
		if (where == GRB_CB_MIPSOL) {
			double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
			std::cout << "Found feasible solution, obj = " << obj << "\n";

			// get solution values
			//std::vector<double> sol(vars.size());
			// the overload returns a heap-allocated array; or use getSolution one-by-one
			double* x = getSolution(xc.data(), (int)xc.size());
			for (int k = 0; k < xc.size(); k++) {
				double val = x[k];
				if (std::abs(val) > 1e-9) {  // filter out 0
					//std::cout << "var[" << i << "] = " << val << std::endl;
					cout << "**************************************" << endl;
					for (int i = 0;i < nc;i++) if (getSolution(ac[i][k]) > 0.5) cout << i << " ";
					cout << endl;
					for (int i = 0;i < nc;i++) if (getSolution(bc[i][k]) > 0.5) cout << i << " ";
					cout << endl;
				}
			}
			delete[] x;  // free the array returned by getSolution

			std::cout << "-----\n";
		}
	}
};


void bicliquePartitionV2() {
	int n, ne, i,j,k;
	vector<vector<int>> graph=makeAdjecencyMatrix(n,ne);
	int in = 0,m;
	cout << "Bound on bicliques, m:";
	cin >> m;
	while (in < 3) {
		cout << "*********************** BICLIQUE PARTITION ****************************" << endl;
		cout << "0)biclique partition\n1)1,2-biclique partition\n2)1,3-biclique partition\n3)odd-biclique partition\n4)exit\n";
		cin >> in;
		if (in == 4) break;
		try {

			GRBEnv env = GRBEnv(true);
			env.start();
			GRBModel model = GRBModel(env);

			vector<vector<GRBVar>> a(n,vector<GRBVar>(m));
			vector<vector<GRBVar>> b(n, vector<GRBVar>(m));
			vector<vector<vector<GRBVar>>> z(n, vector<vector<GRBVar>>(n, vector<GRBVar>(m)));
			vector<GRBVar> x(m); 
			GRBLinExpr Obj = 0;

			for (k = 0;k < m;k++) {
				x[k] = model.addVar(0, 1, 0, GRB_BINARY, "x" + to_string(k));
				Obj += x[k];
			}
			for (i = 0;i < n;i++) for (k = 0;k < m;k++) {
				a[i][k] = model.addVar(0, 1, 0, GRB_BINARY, "a" + to_string(i) + to_string(k));
				b[i][k] = model.addVar(0, 1, 0, GRB_BINARY, "b" + to_string(i) + to_string(k));
			}
			for (i = 0;i < n;i++) for (j = 0;j < n;j++) for (k = 0;k < m;k++) {
				if(i!=j) z[i][j][k] = model.addVar(0, 1, 0, GRB_BINARY, "z" + to_string(i) + to_string(j) + to_string(k));
			}

			for (i = 0;i < n;i++) for (k = 0;k < m;k++) {
				model.addConstr(a[i][k]+b[i][k]<=1); 

				model.addConstr(x[k] >= a[i][k]);
				model.addConstr(x[k] >= b[i][k]);
			}

			for (i = 0;i < n;i++) for (j = 0;j < n;j++) for (k = 0;k < m;k++) {
				if (i == j) continue;
				model.addConstr(z[i][j][k] <= a[i][k]);
				model.addConstr(z[i][j][k] <= b[j][k]);
				model.addConstr(z[i][j][k]+1 >= a[i][k]+b[j][k]);
			}

			for (i = 0;i < n;i++) for (j = 0;j < n;j++) {
				if (i == j) continue;
				GRBLinExpr constr = 0; 
				for (k = 0;k < m;k++) constr += (z[i][j][k] + z[j][i][k]);
				//GRBVar y = model.addVar(0, 1, 0, GRB_BINARY, "y" + to_string(i) + " " + to_string(j));
				GRBVar y;
				switch (in) {
					case 0: // biclique partition
						model.addConstr(constr == 1);
					break;
					case 1: // 1,2-partition 
						model.addConstr(constr >= 1);
						model.addConstr(constr <= 2);
					break;
					case 2: // 1,3-partition
						y = model.addVar(0, 1, 0, GRB_BINARY, "y" + to_string(i) + " " + to_string(j));
						model.addConstr(constr == 2*y + 1);
					break;
					case 3: // odd-partition
						y = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER, "y" + to_string(i) + " " + to_string(j));
						model.addConstr(constr == 2 * y + 1);
					break;
				}
			}

			model.setObjective(Obj, GRB_MINIMIZE);

			MyCallbackV2 cb(x,a,b,n);
			model.setCallback(&cb);

			model.optimize();

			cout << "Size of biclique Partition:" << model.get(GRB_DoubleAttr_ObjVal) << endl;
			cout << "Biclique Partition:" << endl;
			for (k = 0;k < m;k++) {
				if (x[k].get(GRB_DoubleAttr_X) == 1) {
					cout << "**************************************" << endl;
					for (i = 0;i < n;i++) if (a[i][k].get(GRB_DoubleAttr_X) == 1) cout << i << " ";
					cout << endl;
					for (j = 0;j < n;j++) if (b[j][k].get(GRB_DoubleAttr_X) == 1) cout << j << " ";
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

void bicliquePartition() {
	int n, ne, i, j;
	vector<vector<int>> graph = makeAdjecencyMatrix(n, ne);
	vector<pair<vector<int>, vector<int>>> bicliques = makeAllBicliques(graph, n);
	int nb = bicliques.size();
	cout << "No. of bicliques " << nb << endl;
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
							model.addConstr(constraints[i][j] >= 4);
							//model.addConstr(constraints[i][j] == 1);
							break;
						case PARTITION_12:
							model.addConstr(constraints[i][j] >= 1);
							model.addConstr(constraints[i][j] <= 2);
							break;
						case PARTITION_13: // contraints[i][j]==2x+1
							GRBVar x = model.addVar(0, 1, 0, GRB_BINARY, "x" + to_string(i) + " " + to_string(j));
							model.addConstr(constraints[i][j] == 2 * x + 1);
							break;
						}
					}
				}
			}
			model.setObjective(objective, GRB_MINIMIZE);

			//MyCallback cb(B,bicliques);
			//model.setCallback(&cb);

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
