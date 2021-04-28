/*
Convert the multi-valued ALARM Bayes net into a binary Bayes net.
Input (from stdin): a BIF file describing ALARM.
Output (to stdout): MATLAB code defining the graph structure (parent), and the conditional probabilities (p).
This code only works when the maximum in-degree of the input Bayes net is at most 4.
*/

#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define FOR(i,a,b) for(int i=int(a);i<=int(b);++i)
#define REP(i,n) FOR(i,0,(n)-1)
#define PB push_back

#define maxNode 80
#define maxVal 4
#define maxConf (1<<7)

bool isnumber(const string &s) {
	REP(i, s.size()) {
		if (!isdigit(s[i])) return 0;
	}
	return 1;
}

string erase_char(string s, char c) {
	s.erase(remove(s.begin(), s.end(), c), s.end());
	return s;
}

int nbits(int x) {
	if (x <= 2) return 1;
	return 2;
}

int nVal[maxNode], vis[maxNode], newID[maxNode];
double pr[maxNode][maxVal][maxConf], newPr[maxNode][maxVal][maxConf];

int main() {
	map<string, int> node2id;
	map<string, int> val2id[maxNode];
	vector<int> parent[maxNode];
	string s;
	int n = 0;
	while (cin >> s) {
		if (s == "variable") {
			string node_name;
			cin >> node_name;
			node2id[node_name] = ++n;
			while (cin >> s) {
				if (isnumber(s)) {
					istringstream iss(s);
					iss >> nVal[n];
					break;
				}
			}
			while (cin >> s) {
				if (isalpha(s[0])) {
					int k = 0;
					REP(i, nVal[n]) {
						if (i) cin >> s;
						s = erase_char(s, ',');
						val2id[n][s] = k++;
					}
					break;
				}
			}
		} else if (s == "probability") {
			int v;
			// Read v.
			while (cin >> s) {
				if (isalpha(s[0])) {
					v = node2id[s];
					// Read parents of v.
					while (cin >> s) {
						if (s[0] == ')') break;
						if (isalpha(s[0])) {
							s = erase_char(s, ',');
							parent[v].PB(node2id[s]);
						}
					}
					break;
				}
			}
			// Read conditional probabilities.
			//   (parent config): pr[0], pr[1], pr[2] ...
			int totBits = 0;
			REP(i, parent[v].size()) {
				int u = parent[v][i];
				totBits += nbits(nVal[u]);
			}
			while (cin >> s) {
				if (s[0] == '}') break;
				if (s[0] == '(' || isalpha(s[0])) { 
					int config = 0, mul = totBits;
					REP(i, parent[v].size()) {
						int u = parent[v][i];
						mul -= nbits(nVal[u]);
						if (i) cin >> s;
						s = erase_char(s, '(');
						s = erase_char(s, ',');
						s = erase_char(s, ')');
						int val_parent_i = val2id[u][s];
						config += val_parent_i << mul;
					}
					REP(val_v, nVal[v]) {
						scanf("%lf,", &pr[v][val_v][config]);
					}
				}
			}
		}
	}
	
	// Augment all 3-value nodes to 4-value.
	FOR(u, 1, n) if (nVal[u] == 3) {
		nVal[u] = 4;
		// Given a parent configuration of u, split the probability of drawing '10' to drawing '10' or '11'.
		REP(c, maxConf) pr[u][2][c] *= 0.5;
		memcpy(pr[u][3], pr[u][2], sizeof(pr[u][3]));
		// When u is a parent, treat u=='11' the same as u=='10'.
		FOR(v, 1, n) {
			REP(i, parent[v].size()) if (parent[v][i] == u) {
				int numBits = 0;
				FOR(j, i+1, parent[v].size()-1) numBits += nbits(nVal[parent[v][j]]);
				int mask = (1 << numBits) | (1 << (numBits+1));
				REP(val_v, nVal[v])
				REP(c, maxConf) if ((c & mask) == mask) {
					pr[v][val_v][c] = pr[v][val_v][c ^ (1 << numBits)];
				}
			}
		}
	}
	
	// Split each 4-valued node into 2 nodes.
	memset(vis, 0, sizeof(vis));
	int newN = 0, visited = 0;
	while (visited < n) {
		FOR(v, 1, n) if (!vis[v]) {
			bool ok = 1;
			REP(i, parent[v].size()) if (!vis[parent[v][i]]) ok = 0;
			if (ok) {
				newID[v] = newN + 1;
				newN += nbits(nVal[v]);
				++visited;
				vis[v] = 1;
			}
		}
	}
	
	vector<int> newParent[maxNode];
	FOR(v, 1, n) {
		if (nVal[v] <= 2) {
			// Add v's parents in the new graph.
			REP(i, parent[v].size()) {
				int u = parent[v][i];
				newParent[newID[v]].PB(newID[u]);
				if (nVal[u] > 2) newParent[newID[v]].PB(newID[u] + 1);
			}
			REP(i, 2) memcpy(newPr[newID[v]][i], pr[v][i], sizeof(pr[v][i]));
		} else {
			// newID[u] is the high bit of u, newID[u]+1 is the low bit of u.
			int vH = newID[v];
			// vH depends only on v's parents.
			REP(i, parent[v].size()) {
				int u = parent[v][i];
				newParent[vH].PB(newID[u]);
				if (nVal[u] > 2) newParent[vH].PB(newID[u] + 1);
			}
			// newPr[vH][0][c] = pr[v][00][c] + pr[v][01][c].
			REP(val_vH, 2)
			REP(c, maxConf) {
				newPr[vH][val_vH][c] = pr[v][val_vH*2][c] + pr[v][val_vH*2+1][c];
			}
			// vL depends on v's parents and vH.
			int vL = newID[v] + 1;
			REP(i, parent[v].size()) {
				int u = parent[v][i];
				newParent[vL].PB(newID[u]);
				if (nVal[u] > 2) newParent[vL].PB(newID[u] + 1);
			}
			newParent[vL].PB(vH);
			// newPr[vL][0][c 0] = pr[v][00][c] / (pr[v][00][c] + pr[v][01][c]).
			REP(val_vH, 2)
			REP(val_vL, 2) 
			REP(c, 1 << newParent[vH].size()) {
				if (newPr[vH][val_vH][c] < 1e-8) {
					newPr[vL][val_vL][c*2 + val_vH] = 0.5;
				} else {
					newPr[vL][val_vL][c*2 + val_vH] = pr[v][val_vH*2 + val_vL][c] / newPr[vH][val_vH][c];
				}
			}
		}
	}
	
	// Output for MATLAB.
	printf("parent = {");
	FOR(v, 1, newN) {
		if (v > 1) printf(";\n\t");
		printf("[");
		REP(i, newParent[v].size()) {
			int u = newParent[v][i];
			if (i) printf(", ");
			printf("%d", u);
		}
		printf("]");
	}
	puts("};");
	
	printf("p = [");
	FOR(v, 1, newN) {
		REP(c, 1 << newParent[v].size()) {
			if (c) printf(" ");
			printf("%.6lf", newPr[v][1][c]);
		}
		if (v < newN) printf("...\n\t");
		else puts("];");
	}
}