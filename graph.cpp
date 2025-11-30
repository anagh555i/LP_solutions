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

vector<vector<int>> makeAdjecencyMatrix(int& n, int& ne) {
	cout << "No of vertices:";
	cin >> n;
	vector<vector<int>> graph(n,vector<int>(n,0));
	cout << "No of edges(-1 for complete graphs):";
	cin >> ne;
	//if (ne == ((n * (n - 1)) >> 1)) { // instantly return the complete graph
	if (ne == -1 || ne == ((n * (n - 1)) >> 1)) { // instantly return the complete graph
		for (int i = 0;i < n;i++) {
			for (int j = i + 1;j < n;j++) {
				graph[i][j] = 1;
				graph[j][i] = 1;
			}
		}
		return graph;
	}
	int i, u, v;
	cout << "edge sets u v:\n";
	for (i = 0;i < ne;i++) {
		cin >> u >> v;
		graph[u][v] = 1;
		graph[v][u] = 1;
	}
	return graph;
}

bool isConnectingAll(int i, vector<int>& vertices, vector<vector<int>> &graph, int n) {
	for (auto v : vertices) {
		if (graph[i][v] == 0) return false;
	}
	return true;
}
vector<pair<vector<int>, vector<int>>> makeAllBicliques(vector<vector<int>> graph,int n) {
	// graph is an adjacency matrix for an undirected graph;
	vector<pair<vector<int>, vector<int>>> result,curr;
	int i, j, k, m=0;
	for (i = 0;i < n;i++) {
		pair<vector<int>, vector<int>> p( vector<int>(1,i), vector<int>(0));
		curr = vector<pair<vector<int>, vector<int>>>(1, p);
		for (j = i + 1;j < n;j++) {
			int m = curr.size();
			for (k = 0;k < m;k++) {
				p = curr[k];
				if (isConnectingAll(j, p.first,graph,n)) {
					p.second.push_back(j);
					curr.push_back(p);
					p.second.pop_back();
				}
				if (isConnectingAll(j, p.second,graph,n)) {
					p.first.push_back(j);
					curr.push_back(p);
				}
			}
		}
		for (pair<vector<int>, vector<int>> it : curr) result.push_back(it);
	}
	return result;
}