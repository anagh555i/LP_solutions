#pragma once
#include <vector>
#include <iostream>

using namespace std;


#define uint unsigned int
#define MAX_VERTICES 32

class NPartiteNGraph {
public: 
	int n;
	vector<uint> parts;

	NPartiteNGraph(int n) {
		this->n = n;
		this->parts.resize(n);
	}

	NPartiteNGraph(int n, vector<uint> parts) {
		this->n = n;
		this->parts = parts;
	}

	void addVertex(int part, int i) {
		parts[part] = parts[part] | (1 << i);
	}

	bool isPartEmpty(int part) {
		return parts[part] == 0;
	}

	bool isVertexPresent(int part, int v) {
		return ((parts[part] >> v) & 1) == 1;
	}

	void printNPartiteNGraph() {
		cout << "============================================================" << endl;
		for(int part = 0;part < n;part++){
			for (int v = 0; v < MAX_VERTICES; v++) {
				if (isVertexPresent(part,v)) {
					cout << (v + 1) << " ";
				}
			}
			cout << endl;
		}
		cout << "============================================================" << endl;
	}

	int getPartSize(int part) {
		int count = 0;
		for (int v = 0; v < MAX_VERTICES; v++) {
			if (isVertexPresent(part, v)) {
				count++;
			}
		}
		return count;
	}

	int getTotalSize() {
		int count = 0;
		for (int part = 0; part < n; part++) {
			count += getPartSize(part);
		}
		return count;
	}
};

vector<NPartiteNGraph> createAllNPartiteNGraphs(int numOfParts, int numOfVert);
void partitionNUniformHyperGraphIntoNPartiteNGraphs();