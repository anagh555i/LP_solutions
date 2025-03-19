#include "graph.h"
#include<vector>
#include<iostream>

using namespace std;

vector<vector<int>> makeEdgeSet(int& n, int& ne) {
	vector<vector<int>> graph;
	cout << "No of vertices:";
	cin >> n;
	cout << "No of edges:";
	cin >> ne;
	int i, u, v;
	cout << "edge sets u v:\n";
	for (i = 0;i < ne;i++) {
		cin >> u >> v;
		graph.push_back({ u,v });
	}
	return graph;
}