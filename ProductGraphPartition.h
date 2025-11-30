#pragma once
#include<iostream>
#include<vector> 
#include<unordered_map>
#include<queue>
#include<algorithm>

using namespace std;

#define uint unsigned int
#define lluint long long unsigned int

class g4partite {
public:
	uint first;
	uint second;
	uint third;
	uint fourth;
	g4partite(uint a, uint b, uint c, uint d) {
		first = a;
		second = b;
		third = c;
		fourth = d;
	}
};

class rpartite {
public: 
	vector<uint> partitions; 
	rpartite(int n, vector<uint>& g) {
		partitions = vector<uint>(n);
		if (n != g.size()) {
			cout<<"Error: mismatch in size of g and n"<<endl;
			return;
		}
		for (int i = 0;i < n;i++) partitions[i] = g[i]; 
	} 
	bool operator<(const rpartite& other) {
		return partitions < other.partitions;
	} 
	bool operator==(const rpartite& other) {
		return partitions == other.partitions;
	}
	void sortPartitions(int n1) {
		int n = partitions.size(); 
		//cout << "n1=" << n1 << endl;
		//cout << "Before sort" << endl; 
		//printPartitions();
		sort(partitions.begin(), partitions.begin() + n1);
		sort(partitions.begin()+n1, partitions.end());
		//cout << "After sort" << endl;
		//printPartitions();
	}
	bool isValidPartitoning() {
		for (uint p : partitions) if (p == 0) return false; 
		return true;
	}
	bool isValidPartitoning(int n1) {
		for (int i = 0;i < n1;i++) if (partitions[i] == 0) return false;
		return true;
	}
	void printPartitions() {
		if (!isValidPartitoning()) {
			cout << "Invalid Partitioning * ";
		}
		for (uint p : partitions) { 
			for (int i = 0;i < 32;i++) if (p & (1 << i)) cout << i << " ";
			cout << endl;
		}
	}
};


void partition_4partite4graph();
vector< pair<uint, uint> > generateBicliques(int n);

void productPartition();
vector<rpartite> generateRpartiteRgraphs(int n, int n1, int m, int m1);
