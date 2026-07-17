#ifndef ENERGY_DPM_H_
#define ENERGY_DPM_H_

using std::vector;

double Energy_DPM(int Na, vector<int> &Ns, int Nhalf, int Nall, double Rb,
				  vector<double> &L0, vector<double> &D0, vector<double> &A0,
				  vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				  double Kp, double KA, double KAA, double KAA_a, double Kw,
				  vector<int> &idx_start, vector<int> &idx_end,
				  vector<int> &ift, vector<int> &jft,
				  vector<double> &pos, vector<double> &force);

double Force_DPM(int Na, vector<int> &Ns, int Nhalf, int Nall,
				double Rb, vector<double> &D0, double A0,
				vector<double> &Da, double a_cut_c,
				double Kp, double KA, double KAA, double KAA_a, double Kw,
				vector<int> &idx_start, vector<int> &idx_end,
				vector<int> &ift, vector<int> &jft,
				vector<double> &pos, vector<vector <int>> &VL_aa,
				vector<int> &VL_counter, vector<double> &F_a, vector<double> &aforces, 
				 vector<double> &L0, vector<double> &lk, vector<double> &A, double Kb, int nt,
				  int dstep, int bending_added, double box_length, int bending_type, vector<double> &Kl);

double Pressure_DPM(int Na, int Nhalf, double Rb, vector<int> &idx_start, vector<int> &idx_end,
					vector<double> &D0, double Kw, vector<double> &pos);

#endif