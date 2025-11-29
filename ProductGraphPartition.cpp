#include "ProductGraphPartition.h"
#include "gurobi_c++.h"

using namespace std;

class MyCallback : public GRBCallback {
public:
	vector<GRBVar> vars;
	vector<g4partite> KnKm;

	MyCallback(const vector<GRBVar>& v, vector<g4partite> & b) : vars(v), KnKm(b) {}

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
					g4partite it = KnKm[i];
					cout << "***********************************" << endl;
					for (int i = 0;it.first && i < 32;i++, it.first = (it.first >> 1)) if (it.first & 1) cout << i << " ";
					cout << endl;
					for (int i = 0;it.second && i < 32;i++, it.second = (it.second >> 1)) if (it.second & 1) cout << i << " ";
					cout << endl;
					for (int i = 0;it.third && i < 32;i++, it.third = (it.third >> 1)) if (it.third & 1) cout << i << " ";
					cout << endl;
					for (int i = 0;it.fourth && i < 32;i++, it.fourth = (it.fourth >> 1)) if (it.fourth & 1) cout << i << " ";
					cout << endl;
				}
			}
			delete[] x;  // free the array returned by getSolution

			std::cout << "-----\n";
		}
	}
};


void partition_4partite4graph() {
	/*
		Partitions product graph Kn X Km into 4-partite-4-graphs
		
		Graph representation:
			=> an integer pair represent a biclique, each integer being the bit mask of vertex set of either left or right partition
			=> 2 integer pairs represent a 4-partite-4-graph

	*/
	int n, m;
	cout << "Partitions product graph Kn X Km into 4-partite-4-graphs" << endl;
	cout << "n = ";
	cin >> n;
	cout << "m = ";
	cin >> m;

	vector<pair<uint, uint>> Kn = generateBicliques(n);
	vector<pair<uint, uint>> Km = generateBicliques(m);
	/*for (auto it : Kn) {
		cout << "*********" << endl;
		for (int i = 0;it.first && i < 32;i++, it.first = (it.first >> 1)) if (it.first & 1) cout << i << " ";
		cout << endl;
		for (int i = 0;it.second && i < 32;i++, it.second = (it.second >> 1)) if (it.second & 1) cout << i << " ";
		cout << endl;
		cout << "*********" << endl;
	}*/
	cout << "size of Kn's: " << Kn.size() << endl;
	cout << "size of Km's: " << Km.size() << endl;

	vector<g4partite> KnKm;
	for (auto kn : Kn) {
		for (auto km : Km) {
			KnKm.push_back(g4partite(kn.first,kn.second,km.first,km.second));
		}
	}
	cout << "size of product partitions: " << KnKm.size() << endl;
	/*for (auto it : KnKm) {
		cout << "*********" << endl;
		for (int i = 0;it.first && i < 32;i++, it.first = (it.first >> 1)) if (it.first & 1) cout << i << " ";
		cout << endl;
		for (int i = 0;it.second && i < 32;i++, it.second = (it.second >> 1)) if (it.second & 1) cout << i << " ";
		cout << endl;
		for (int i = 0;it.third && i < 32;i++, it.third = (it.third >> 1)) if (it.third & 1) cout << i << " ";
		cout << endl;
		for (int i = 0;it.fourth && i < 32;i++, it.fourth = (it.fourth >> 1)) if (it.fourth & 1) cout << i << " ";
		cout << endl;
		cout << "*********" << endl;
	}*/

	try
	{
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		vector<GRBVar> B; // for each b belongs to B, b represents the presence of a 4partite4graph

		int nb = KnKm.size();
		for (int i = 0;i < nb;i++) B.push_back(model.addVar(0, 1, 0, GRB_BINARY, "b" + to_string(i)));

		vector<vector<vector<vector<GRBLinExpr>>>> constraints(n, vector<vector<vector<GRBLinExpr>>>(n, vector<vector<GRBLinExpr>>(m, vector<GRBLinExpr>(m,0))));
		GRBLinExpr objective = 0;

		for (int i = 0;i < nb;i++) {
			objective += B[i];
			g4partite g = KnKm[i];
			// yeah edge is a-b-c-d.
			int a, b, c, d;
			int ga, gb, gc, gd;
			ga = g.first;
			for (a = 0;ga;a++,ga=ga>>1) {
				if ((ga & 1) == 0) continue;
				gb = g.second;
				for (b = 0; gb; b++, gb = gb >> 1) {
					if ((gb & 1) == 0) continue;
					gc = g.third;
					for (c = 0; gc; c++, gc = gc >> 1) {
						if ((gc & 1) == 0) continue;
						gd = g.fourth;
						for (d = 0; gd; d++, gd = gd >> 1) {
							if ((gd & 1) == 0) continue;
							//cout << a << " " << b << " " << c << " " << d << " **** " << i << endl;
							constraints[min(a,b)][max(a,b)][min(c,d)][max(c,d)] += B[i];
						}
					}
				}
			}
		}
		int testCount = 0;
		for (int a = 0;a < n;a++) {
			for (int b = 0;b < n;b++) {
				for (int c = 0;c < m;c++) {
					for (int d = 0;d < m;d++) {
						if (constraints[a][b][c][d].size() == 0 && constraints[a][b][c][d].getConstant() == 0.0) {
							testCount++;
							continue;
						}
						model.addConstr(constraints[a][b][c][d] == 1);
					}
				}
			}
		}
		cout << testCount << endl;
		/*double upperBound, lowerBound;
		cout << "Upper Bound: ";
		cin >> upperBound;
		cout << "Lower Bound: ";
		cin >> lowerBound;*/

		//model.addConstr(objective <= upperBound);
		//model.addConstr(objective >= lowerBound);
		//model.addConstr(objective == upperBound);

		//model.set(GRB_IntParam_SolutionLimit, 1);
		//model.set(GRB_IntParam_MIPFocus, 2);
		//model.set(GRB_IntParam_Presolve, 2);
		//model.set(GRB_IntParam_Symmetry, 2);
		//model.set(GRB_IntParam_Cuts, 2);
		//model.set(GRB_DoubleParam_MIPGap, 0.05);

		model.setObjective(objective, GRB_MINIMIZE);
		//MyCallback cb(B, KnKm); // use if you want to print intermediate solutions
		//model.setCallback(&cb);
		model.optimize();
		cout << "Size of partitions:" << model.get(GRB_DoubleAttr_ObjVal) << endl;
		cout << "Biclique Partition:" << endl;
		for (int i = 0;i < nb;i++) {
			if (B[i].get(GRB_DoubleAttr_X) == 1) {
				g4partite it = KnKm[i];
				cout << "***********************************" << endl;
				for (int i = 0;it.first && i < 32;i++, it.first = (it.first >> 1)) if (it.first & 1) cout << i << " ";
				cout << endl;
				for (int i = 0;it.second && i < 32;i++, it.second = (it.second >> 1)) if (it.second & 1) cout << i << " ";
				cout << endl;
				for (int i = 0;it.third && i < 32;i++, it.third = (it.third >> 1)) if (it.third & 1) cout << i << " ";
				cout << endl;
				for (int i = 0;it.fourth && i < 32;i++, it.fourth = (it.fourth >> 1)) if (it.fourth & 1) cout << i << " ";
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

vector<pair<uint, uint>> generateBicliques(int n) {
	/*
		Generates all bicliques from complete graph Kn

		Graph representation:
			=> an integer pair represent a biclique, each integer being the bit mask of vertex set of either left or right partition
	*/

	vector<pair<uint, uint>> result,gen;

	int i,currSize,j;
	for (i = 0;i < n;i++) {
		currSize = gen.size();
		for (j = 0;j < currSize;j++) {
			gen.push_back({
				gen[j].first | (1 << i),
				gen[j].second
				});
			gen.push_back({
				gen[j].first ,
				gen[j].second | (1 << i)
				});
		}
		gen.push_back({ 1 << i,0 });
	}
	for (pair<uint, uint> it : gen) if (it.second != 0) result.push_back(it);

	return result;
}

vector<rpartite> generateRpartiteRgraphs(int n, int n1, int m, int m1) {
	int r = n1 + m1;
	vector<uint> curr(r, 0);
	vector<rpartite> buff(1, rpartite(r, curr));
	int currSize = 1;
	for (int i = 0;i < n;i++) {
		int kn = currSize;
		for (int k = 0;k < kn;k++) {
			curr = buff[k].partitions;
			for (int j = 0;j < n1;j++) {
				curr[j] = (curr[j] ^ (1 << i));
				buff.push_back(rpartite(r, curr));
				currSize++;
				curr[j] = (curr[j] ^ (1 << i));
			}
		}
	} 
	vector<rpartite> buff2;
	for (rpartite& p : buff) if (p.isValidPartitoning(n1)) buff2.push_back(p); 
	currSize = buff2.size(); 
	//cout << "generated till buff2 of size=" << currSize << endl;
	for (int i = 0;i < m;i++) {
		int kn = currSize;
		for (int k = 0;k < kn;k++) {
			curr = buff2[k].partitions;
			for (int j = n1;j < r;j++) {
				curr[j] = (curr[j] ^ (1 << i));
				buff2.push_back(rpartite(r, curr));
				currSize++;
				curr[j] = (curr[j] ^ (1 << i));
			}
		}
	}
	vector<rpartite> result;  
	for (rpartite& p : buff2) {
		if (!p.isValidPartitoning()) continue;
		p.sortPartitions(n1);
		result.push_back(p);
	}
	sort(result.begin(), result.end());
	result.erase(unique(result.begin(), result.end()), result.end());
	return result;
}


void productPartition() {
	int n, m,n1,m1;
	cout << "Partition product graph Kn X Km into r-partite-r-graphs" << endl << "n1-edges from Kn and m1-edges from Km, where r=n1+m1" << endl;
	cout << "n = ";
	cin >> n;
	cout << "n1 = ";
	cin >> n1;
	cout << "m = ";
	cin >> m;
	cout << "m1 = ";
	cin >> m1;
	int r = n1 + m1;
	vector<rpartite> rpartites = generateRpartiteRgraphs(n, n1, m, m1); 
	cout << "Size of partition = " << rpartites.size() << endl;
	/*for (rpartite& g : rpartites) {
		g.printPartitions();
		cout << "********************************" << endl;
	}*/
	try
	{
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		vector<GRBVar> R; // for each b belongs to B, b represents the presence of a RpartiteRgraph

		int nr = rpartites.size();
		for (int i = 0;i < nr;i++) R.push_back(model.addVar(0, 1, 0, GRB_BINARY, "r" + to_string(i)));

		//unordered_map<pair<int, int>, GRBLinExpr> constraints;
		unordered_map<lluint, GRBLinExpr> constraints;

		//vector<vector<vector<vector<GRBLinExpr>>>> constraints(n, vector<vector<vector<GRBLinExpr>>>(n, vector<vector<GRBLinExpr>>(m, vector<GRBLinExpr>(m, 0))));
		GRBLinExpr objective = 0;

		for (int i = 0;i < nr;i++) {
			objective += R[i];
			rpartite rp = rpartites[i];
			//  edges are of form first-second. where first has n1 vertices and second has m1 vertices 
			// generate all edges 
			queue<uint> edges1;
			edges1.push(0);
			int prevSize = 1;
			for (int j = 0;j < n1;j++) {
				int currSize = 0;
				while (prevSize--) {
					uint edge = edges1.front();
					edges1.pop();
					uint p = rp.partitions[j];
					for (int s = 0;p;p = (p >> 1), s++) if (p & 1) {
						edges1.push(edge | (1 << s));
						currSize++;
					}
				}
				prevSize = currSize;
			}
			queue<uint> edges2;
			edges2.push(0);
			prevSize = 1;
			for (int j = n1;j < r;j++) {
				int currSize = 0;
				while (prevSize--) {
					uint edge = edges2.front();
					edges2.pop();
					uint p = rp.partitions[j];
					for (int s = 0;p;p = (p >> 1), s++) if (p & 1) {
						edges2.push(edge | (1 << s));
						currSize++;
					}
				}
				prevSize = currSize;
			}
			//rp.printPartitions();
			//cout << "edges1.size()=" << edges1.size() << endl;
			//cout << "edges2.size()=" << edges2.size() << endl;
			while (!edges1.empty()) {
				uint e1 = edges1.front(); 
				edges1.pop(); 
				prevSize = edges2.size();
				while (prevSize--) {
					uint e2 = edges2.front();
					//cout << e1 << " " << e2 << endl;
					edges2.pop(); 
					constraints[(lluint(e1) << 32) | e2] += R[i];
					edges2.push(e2);
				}
			}
		}

		// for each constraint in each edge, add constraint to model
		for (auto& kv : constraints) {
			//cout << "##################" << endl;
			//cout << kv.first << endl;
			model.addConstr(kv.second == 1);
		}

		model.setObjective(objective, GRB_MINIMIZE);
		model.optimize();

		cout << "Size of partitions:" << model.get(GRB_DoubleAttr_ObjVal) << endl;
		cout << "Biclique Partition:" << endl;
		for (int i = 0;i < nr;i++) {
			if (R[i].get(GRB_DoubleAttr_X) == 1) {
				rpartite rp = rpartites[i];
				cout << "*****************************************************************" << endl;
				rp.printPartitions();
				cout << "*****************************************************************" << endl;
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