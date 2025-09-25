#include "ProductGraphPartition.h"
#include "gurobi_c++.h"


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
		model.set(GRB_IntParam_MIPFocus, 2);
		model.set(GRB_IntParam_Presolve, 2);
		model.set(GRB_IntParam_Symmetry, 2);
		model.set(GRB_IntParam_Cuts, 2);
		model.set(GRB_DoubleParam_MIPGap, 0.05);

		model.setObjective(objective, GRB_MINIMIZE);
		model.optimize();
		cout << "Size of partitions:" << model.get(GRB_DoubleAttr_ObjVal) << endl;
		cout << "Biclique Partition:" << endl;
		for (int i = 0;i < nb;i++) {
			if (B[i].get(GRB_DoubleAttr_X) == 1) {
				auto it = KnKm[i];
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


void productPartition() {
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
			KnKm.push_back(g4partite(kn.first, kn.second, km.first, km.second));
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
}