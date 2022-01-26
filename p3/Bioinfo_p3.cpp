// Bioinfo_p3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <new>
#include <stdexcept>
#include <string>
#include <cstring>
#include <iostream>
#include <queue>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>

#define MAXDIST 999999

using namespace std;

map<string, string> oghash;
map<string, map<string, int>> hash_clusters;
map<string, string> Newick_hash;


fstream input;

int main()
{
	vector <string> codes;
	vector <vector <int>> dist_matrix;
	int numOTU;
	string file_str;
	vector <string> arrayclusters;

	cout << "Enter the name of an input file.\n";
	cin >> file_str;
	input.open(file_str);
	if (input.fail()) {
		cout << "Failed to open file " << file_str << "\n";
		exit(1);
	}
	//Get number of OTU
	string OTU;
	getline(input, OTU);
	numOTU = stoi(OTU, nullptr, 10);
	//Get keys
	string t_codes;
	getline(input, t_codes);
	char* tc = strdup(t_codes.c_str());
	char* tok;
	char* token;
	token = strtok(tc, "\r");
	tok = strtok(token, "\t");
	while (tok != NULL) {
		string token = tok;
		codes.push_back(token);
		tok = strtok(NULL, "\t");
	}
	//Initialize hash codes and populate distance matrix
	for (int i = 0; i < numOTU; i++) {
		Newick_hash[codes[i]] = codes[i];
		oghash[codes[i]] = codes[i];
		vector <int> vtemp;
		string temp;

		getline(input, temp);
		char* t = strdup(temp.c_str());
		char* stupid = strtok(t, "\r");
		tok = strtok(stupid, "\t");
		while (tok != NULL) {
			string token = tok;
			tok = strtok(NULL, "\t");
			vtemp.push_back(stoi(token, nullptr, 10));
		}
		dist_matrix.push_back(vtemp);
	}
	//Print dist_matrix while initializing hash clusters
	cout << "Codes:\n";
	int cs = codes.size();
	for (int i = 0; i < numOTU; i++) {
		printf("%s ", codes[i].c_str());
	}
	cout << "\n";
	for (int i = 0; i < numOTU; i++) {
		printf("%s ", codes[i].c_str());
		for (int j = 0; j < numOTU; j++) {
			cout << dist_matrix[i][j] << " ";
			hash_clusters[codes[i]][codes[j]] = dist_matrix[i][j];
		}
		cout << "\n";
	}
	//agglomerative clustering
	int num_clusters = numOTU;
	while (num_clusters > 2) {

		for (auto const& element : hash_clusters) {
			arrayclusters.push_back(element.first);
		}
		//calculate R values
		vector <int> arrayRvalues;
		for (int i = 0; i < num_clusters; i++) {
			int temp = 0;
			for (int j = 0; j < num_clusters; j++) {
				temp = temp + hash_clusters[arrayclusters[i]][arrayclusters[j]];
			}
			arrayRvalues.push_back(temp / (num_clusters - 2));
		}

		cout << "R values\n";
		for (int i = 0; i < num_clusters; i++) {
			cout << arrayclusters[i] << " " << arrayRvalues[i] << " ";
		}
		cout << "\n";

		//Determine merging cluster using transformed distances
		int smallest = MAXDIST;
		int smallestI = 0;
		int smallestJ = 0;
		cout << "TD matrix:\n";
		for (int i = 0; i < num_clusters - 1; i++) {
			for (int j = i + 1; j < num_clusters; j++) {
				//debug
				//cout << "cluster " << i << ": " << arrayclusters[i] << " cluster " << j << ": " << arrayclusters[j] << " distance: " << oghash[arrayclusters[i]][arrayclusters[j]] << "\n";
				int tempTD = hash_clusters[arrayclusters[i]][arrayclusters[j]] -
					arrayRvalues[i] - arrayRvalues[j];
				//check if smallest transformed distance
				cout << "TD(" << arrayclusters[i] << ", " << arrayclusters[j] <<
					") = " << tempTD << ", ";
				if (tempTD < smallest) {
					smallest = tempTD;
					smallestI = i;
					smallestJ = j;
				}

			}
			cout << "\n";
		}
		string clusterI = arrayclusters[smallestI];
		string clusterJ = arrayclusters[smallestJ];
		string merge = clusterI;
		merge.append(clusterJ);

		//calculate branch length
		int branch1 = (hash_clusters[clusterJ][clusterI] + arrayRvalues[smallestI]
			- arrayRvalues[smallestJ]) / 2;
		int branch2 = (hash_clusters[clusterJ][clusterI] +
			arrayRvalues[smallestJ] - arrayRvalues[smallestI]) / 2;
		cout << "merging " << clusterI << " and " << clusterJ << "\n";
		cout << "Distance between " << clusterI << " and ancestral node = " << branch1 << "\n";
		cout << "Distance between " << clusterJ << " and ancestral node = " << branch2 << "\n";


		for (int i = 0; i < num_clusters; i++) {
			if ((arrayclusters[i] != clusterI) && (arrayclusters[i] != clusterJ)) {
				int d1 = hash_clusters[arrayclusters[i]][clusterI];
				int d2 = hash_clusters[arrayclusters[i]][clusterJ];
				//calculate new distance value
				hash_clusters[merge][arrayclusters[i]] = (d1 + d2 -
					hash_clusters[clusterI][clusterJ]) / 2;
				hash_clusters[arrayclusters[i]][merge] = hash_clusters[merge][arrayclusters[i]];

			}
		}

		string newick;
		newick.append("(");
		newick.append(Newick_hash[clusterI]);
		newick.append(":");
		newick.append(to_string(branch1));
		newick.append(",");
		newick.append(Newick_hash[clusterJ]);
		newick.append(":");
		newick.append(to_string(branch2));
		newick.append(")");
		Newick_hash[merge] = newick;


		//erase the clusters that were merged from original map
		for (int i = 0; i < num_clusters; i++) {
			hash_clusters[clusterI].erase(arrayclusters[i]);
			hash_clusters[clusterJ].erase(arrayclusters[i]);
		}
		hash_clusters.erase(clusterI);
		hash_clusters.erase(clusterJ);
		Newick_hash.erase(clusterI);
		Newick_hash.erase(clusterJ);

		//create new array of cluster values
		arrayclusters.clear();
		for (auto const& element : hash_clusters) {
			arrayclusters.push_back(element.first);
		}
		num_clusters--;
	}
	//print info
	cout << "Distance between remaining clusters: " <<
		hash_clusters[arrayclusters[0]][arrayclusters[1]] << ", " <<
		hash_clusters[arrayclusters[1]][arrayclusters[0]] << "\n";
		
	cout << "Newick format: \n";
	for (int j = 0; j < num_clusters; j++) {
		cout << Newick_hash[arrayclusters[j]] << " " << "\n";
	}
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
