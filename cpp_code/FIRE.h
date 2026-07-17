#ifndef FIRE_H_
#define FIRE_H_

using std::vector;

double FIRE_DPM_RP(int Na, vector<int> &Ns, int Nhalf, int Nall,
				   double Rb, vector<double> &L0, vector<double> &D0, vector<double> &A0,
				   vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				   double Kp, double KA, double KAA, double KAA_a, double Kw,
				   vector<int> &idx_start, vector<int> &idx_end,
				   vector<int> &ift, vector<int> &jft, vector<double> &pos,
				   double Fthresh, double dt_fire, int Nt_fire);


double FIRE_DPM_RP_VL(int Na, vector<int> &Ns, int Nhalf, int Nall,
				   	  	double Rb, vector<double> &D0, double A0,
						vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				   		double Kp, double KA, double KAA, double KAA_a, double Kw,
				   		vector<int> &idx_start, vector<int> &idx_end,
						vector<int> &ift, vector<int> &jft,
						vector<double> &pos, double Fthresh, double dt_fire, int Nt_fire,  
				 		vector<double> &L0, double Kb, 
						int dstep, int bending_added, double box_length, int bending_type, double l0);

					

double Verlet_Int(int Na, vector<int> &Ns, int Nhalf, int Nall,
				  double Rb, vector<double> &D0, double A0,
				  vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				  double Kp, double KA, double KAA, double KAA_a, double Kw,
				  vector<int> &idx_start, vector<int> &idx_end,
				  vector<int> &ift, vector<int> &jft,
				  vector<double> &pos, double Fthresh, double dt_fire, int Nt_verlet, vector<double> &Vel, int dstep,
				  int restartstep, double Ti, double B, 
				  vector<int> &inflate_ids, vector<double> &L0, double drag_l, double Alf_set, vector<int> &deflation_ids, 
				  int n_deflate, double drag_a, int flag_inflate, int flag_deflate, double kA0, double Kb,
				  double kal_A, int bending_added, int perimeter_plasticity_flag, int area_plasticity_flag, double l0min, int area_plasticity_all, 
				  int whole_dp_plasticity, double box_length, double delta, int temperature_on, double pressure_crit,
				  int active_yes, double f0, double tau, vector<double> &theta, int activity_onset_time, int global_step,
				  double phi_crit, double KE_thresh,double dA0inf, int checking_time, double phi_slow, int max_relaxation_time,
				  int jamming_yes, double l0max, double dphi, double phi_start, int com, int bending_type, double l0);


#endif