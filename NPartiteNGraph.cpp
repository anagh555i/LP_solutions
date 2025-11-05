#include "NPartiteNGraph.h"
#include "gurobi_c++.h"

#include <map>
#include <algorithm>



vector<NPartiteNGraph> createAllNPartiteNGraphs(int numOfParts,int numOfVert) {
	vector<NPartiteNGraph> gen;
	for (int v = 0; v < numOfVert; v++) {
		vector<NPartiteNGraph> temp;
		for (auto& nPartiteNGraph : gen) {
			bool emptyPart = false;
			for (int part = 0; part < numOfParts; part++) {
				NPartiteNGraph newNPartiteNGraph = nPartiteNGraph;
				if (newNPartiteNGraph.isPartEmpty(part)) {
					if (emptyPart) continue;
					emptyPart = true;
				}
				newNPartiteNGraph.addVertex(part, v);
				temp.push_back(newNPartiteNGraph);
			}
		}
		if (!temp.empty()) {
			gen.insert(gen.end(), temp.begin(), temp.end());
		}
		NPartiteNGraph newNPartiteNGraph(numOfParts);
		newNPartiteNGraph.addVertex(0, v);
		gen.push_back(newNPartiteNGraph);
	}
	vector<NPartiteNGraph> res;
	for (auto& nPartiteNGraph : gen) {
		bool emptyPart = false;
		for (int part = 0; part < numOfParts; part++) {
			if (nPartiteNGraph.isPartEmpty(part)) {
				emptyPart = true;
				break;
			}
		}
		if (!emptyPart) {
			res.push_back(nPartiteNGraph);
		}
	}
	return res;
}

int getSetBits(uint x) {
	int count = 0;
	while (x > 0) {
		count += (x & 1);
		x >>= 1;
	}
	return count;
}

void populateConstraints(int part,uint curEdge,NPartiteNGraph& nPartiteNGraph,GRBVar& nPartiteNGraphVar,map<uint,GRBLinExpr>& constraints) {
	if (getSetBits(curEdge) == nPartiteNGraph.n) {
		constraints[curEdge] += nPartiteNGraphVar;
		return;
	}

	if (part == nPartiteNGraph.n) {
		return;
	}

	for (int v = 0; v < MAX_VERTICES; v++) {
		if (nPartiteNGraph.isVertexPresent(part, v)) {
			populateConstraints(part + 1,curEdge | (1 << v), nPartiteNGraph, nPartiteNGraphVar, constraints);
		}
	}
}

void partitionNUniformHyperGraphIntoNPartiteNGraphs() {
	int numOfVert,numOfPart;
	cout << "Enter the numOfParts in the N Partite N Graph" << endl;
	cin >> numOfPart;
	cout << "Enter the numOfVert in the M Uniform Hypergraph" << endl;
	cin >> numOfVert;

	vector<NPartiteNGraph> nPartiteNGraphs = createAllNPartiteNGraphs(numOfPart, numOfVert);
	int numOfNPartiteNGraphs = nPartiteNGraphs.size();

	cout << "Number of N Partite N Graphs: " << numOfNPartiteNGraphs << endl;

	try {
		GRBEnv env = GRBEnv(true);
		env.start();
		GRBModel model = GRBModel(env);

		vector<GRBVar> nPartiteNGraphVars;
		for (int i = 0; i < numOfNPartiteNGraphs; i++) {
			nPartiteNGraphVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "npn" + to_string(i)));
		}
		
		map<uint, GRBLinExpr> constraints;
		GRBLinExpr objective = 0;
		for (int i = 0; i < numOfNPartiteNGraphs; i++) {
			objective += nPartiteNGraphVars[i];
			populateConstraints(0,0, nPartiteNGraphs[i], nPartiteNGraphVars[i],constraints);
		}

		for (auto& constraintMapping : constraints) {
			GRBLinExpr constr = constraintMapping.second;
			model.addConstr(constr == 1);
		}

		model.set(GRB_IntParam_MIPFocus, 2);
		model.set(GRB_IntParam_Presolve, 2);
		model.set(GRB_IntParam_Symmetry, 2);
		model.set(GRB_IntParam_Cuts, 2);
		model.set(GRB_DoubleParam_MIPGap, 0.05);

		model.setObjective(objective, GRB_MINIMIZE);
		model.optimize();
		cout << "Size of partitions:" << model.get(GRB_DoubleAttr_ObjVal) << endl;
		printf("%d Partite %d Partition of %d Uniform HyperGraph in %d Vertices\n", numOfPart, numOfPart, numOfPart, numOfVert);
		for (int i = 0; i < numOfNPartiteNGraphs; i++) {
			if (nPartiteNGraphVars[i].get(GRB_DoubleAttr_X) == 1) {
				nPartiteNGraphs[i].printNPartiteNGraph();
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


