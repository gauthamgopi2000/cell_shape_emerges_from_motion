#ifndef VERLETLIST_H_
#define VERLETLIST_H_

using std::vector;

void VerletList(int Na, vector<int> &Ns, int Nhalf, int Nall,
				vector<double> &Da, double r_cut,
				vector<int> &idx_start, vector<int> &idx_end,
				vector<double> &pos_a, vector<double> &pos_save,
				vector<vector <int>> &VL_aa, vector<int> &VL_counter,
				int first_call, double box_length);

#endif