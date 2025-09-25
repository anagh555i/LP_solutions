#pragma once
#include<iostream>
#include<vector>

using namespace std;

#define uint unsigned int

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

void partition_4partite4graph();
vector< pair<uint, uint> > generateBicliques(int n);
void productPartition();
