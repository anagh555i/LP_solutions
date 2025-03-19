#include "gurobi_c++.h"
#include <iostream>
#include "test.h"
#include <vector>

void test() {
    try {
        // Initialize environment and model
        GRBEnv env = GRBEnv(true);
        env.start();
        GRBModel model = GRBModel(env);

        // Graph: Bipartite edges (U,V) -> E
        std::vector<std::pair<int, int>> edges = { {0, 3}, {0, 4}, {1, 3}, {2, 4} };
        int U = 3, V = 2;  // Set sizes

        // Variable map for edges
        std::vector<GRBVar> edgeVars;
        for (int i = 0; i < edges.size(); i++) {
            edgeVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "x_" + std::to_string(i)));
        }

        // Objective function: Maximize sum of selected edges
        GRBLinExpr obj = 0;
        for (auto& var : edgeVars) obj += var;
        model.setObjective(obj, GRB_MAXIMIZE);

        // Constraints: Each node can be in at most one edge
        for (int u = 0; u < U; u++) {
            GRBLinExpr sum = 0;
            for (int i = 0; i < edges.size(); i++)
                if (edges[i].first == u) sum += edgeVars[i];
            model.addConstr(sum <= 1);
        }

        for (int v = 0; v < V; v++) {
            GRBLinExpr sum = 0;
            for (int i = 0; i < edges.size(); i++)
                if (edges[i].second == v + U) sum += edgeVars[i];
            model.addConstr(sum <= 1);
        }

        // Optimize model
        model.optimize();

        // Print results
        for (int i = 0; i < edges.size(); i++) {
            if (edgeVars[i].get(GRB_DoubleAttr_X) > 0.5) {
                std::cout << "Edge " << edges[i].first << " - " << edges[i].second << " is selected.\n";
            }
        }

    }
    catch (GRBException& e) {
        std::cerr << "Error: " << e.getMessage() << std::endl;
    }

    try {
        // Create environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "gurobi.log");
        env.start();

        // Create model
        GRBModel model = GRBModel(env);

        // Create variables
        GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
        GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
        GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

        // Set objective: maximize x + y + 2z
        model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

        // Add constraint: x + 2y + 3z <= 4
        model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

        // Add constraint: x + y >= 1
        model.addConstr(x + y >= 1, "c1");

        // Optimize model
        model.optimize();

        // Output results
        std::cout << "Optimal solution:\n";
        std::cout << "x = " << x.get(GRB_DoubleAttr_X) << "\n";
        std::cout << "y = " << y.get(GRB_DoubleAttr_X) << "\n";
        std::cout << "z = " << z.get(GRB_DoubleAttr_X) << "\n";

    }
    catch (GRBException& e) {
        std::cerr << "Error code = " << e.getErrorCode() << "\n";
        std::cerr << e.getMessage() << "\n";
    }
    catch (...) {
        std::cerr << "Exception during optimization" << "\n";
    }

}
