#pragma once
#include<vector>

using namespace std;

vector<vector<int>> makeEdgeSet(int &n, int &ne);
vector<vector<int>> makeAdjecencyMatrix(int& n, int& ne);
bool isConnectingAll(int i, vector<int>& vertices, vector<vector<int>>& graph, int n);
vector<pair<vector<int>, vector<int>>> makeAllBicliques(vector<vector<int>> graph, int n);