#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

#include "VerletList.h"
#include "Energy_DPM.h"

using std::vector;

double FIRE_DPM_RP(int Na, vector<int> &Ns, int Nhalf, int Nall,
				   double Rb, vector<double> &L0, vector<double> &D0, vector<double> &A0,
				   vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				   double Kp, double KA, double KAA, double KAA_a, double Kw,
				   vector<int> &idx_start, vector<int> &idx_end,
				   vector<int> &ift, vector<int> &jft, vector<double> &pos,
				   double Fthresh, double dt_fire, int Nt_fire)
{
	int 		n, nt = 0;
	vector<double>	Acc(Nall, 0.0);
	vector<double>	Vel(Nall, 0.0);
	double 		P, Vel_norm, Acc_norm, Acc_max, Acc_abs;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	double 		Pw;
	Pw = Energy_DPM(Na, Ns, Nhalf, Nall, Rb, L0, D0, A0, Da, alpha_0, a_cut_c,
					Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft, pos, Acc);
	
	Acc_max = std::abs(Acc[0]);
	for (n = 1; n < Nall; n++){
		Acc_abs = std::abs(Acc[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	if (Acc_max < Fthresh){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return Pw;
	}

	for (nt = 1; nt < Nt_fire; nt++){
		P = 0.0;
		for (n = 0; n < Nall; n++)
			P += Acc[n] * Vel[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max){
				break;
			}
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
				for (n = 0; n < Nall; n++){
					pos[n] -= dt_half * Vel[n];
					Vel[n] = 0.0;
				}
			}
		}

		// MD using Verlet method
		for (n = 0; n < Nall; n++)
			Vel[n] += dt_half * Acc[n];
		Vel_norm = 0.0;
		Acc_norm = 0.0;
		for (n = 0; n < Nall; n++){
			Vel_norm += Vel[n] * Vel[n];
			Acc_norm += Acc[n] * Acc[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Nall; n++){
			Vel[n] = delta_a * Vel[n] + v_rsc * Acc[n];
			pos[n] += dt * Vel[n];
		}
		Pw = Energy_DPM(Na, Ns, Nhalf, Nall, Rb, L0, D0, A0, Da, alpha_0, a_cut_c,
						Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft, pos, Acc);
		for (n = 0; n < Nall; n++)
			Vel[n] += dt_half * Acc[n];

		Acc_max = std::abs(Acc[0]);
		for (n = 1; n < Nall; n++){
			Acc_abs = std::abs(Acc[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		if (Acc_max < Fthresh){
			break;
		}
		if (nt % 10000 == 0) printf("nt: %d  Max Force: %.5e\n", nt, Acc_max);
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
	return Pw;
}



double FIRE_DPM_RP_VL(int Na, vector<int> &Ns, int Nhalf, int Nall,
				   	  	double Rb, vector<double> &L0, vector<double> &D0, vector<double> &A0,
						vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				   		double Kp, double KA, double KAA, double KAA_a, double Kw,
				   		vector<int> &idx_start, vector<int> &idx_end,
						vector<int> &ift, vector<int> &jft,
						vector<double> &pos, double Fthresh, double dt_fire, int Nt_fire)
{
	int 		n, nt = 0;
	vector<double>	Acc(Nall, 0.0);
	vector<double> 	Vel(Nall, 0.0);
	double 		P, Vel_norm, Acc_norm, Acc_max, Acc_abs;
	// FIRE parameters
	int 		N_delay = 20;
	int 		N_pn_max = 2000;
	double 		f_inc = 1.1;
	double 		f_dec = 0.5;
	double 		a_start = 0.15;
	double 		f_a = 0.99;
	double 		dt_max = 10.0 * dt_fire;
	double 		dt_min = 0.05 * dt_fire;
	int 		initialdelay = 1;
	double 		v_rsc;
	// FIRE Initialization
	double 		dt = dt_fire;
	double 		dt_half = dt / 2.0;
	double 		a_fire = a_start;
	double 		delta_a = 1.0 - a_fire;
	int 		N_pp = 0, N_pn = 0;

	// Verlet list parameters
	int 		first_call = 1;
	vector<vector <int>>	VL_aa(20 * Nhalf, vector<int> (4, -1));
	vector<int>		 		VL_counter_all(1, 0);
	vector<double>			pos_save(Nall, 0.0);
	double 		r_cut = D0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > D0[n]) ? r_cut : D0[n];
	r_cut *= a_cut_c;

	double Pw;
	
	VerletList(Na, Ns, Nhalf, Nall, Da, r_cut, idx_start, idx_end,
			   pos, pos_save, VL_aa, VL_counter_all, first_call);
	Pw = Force_DPM(Na, Ns, Nhalf, Nall, Rb, L0, D0, A0, Da, alpha_0, a_cut_c,
					Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
					pos, VL_aa, VL_counter_all, Acc);
	/*
	Eval = Energy_DPM(Na, Ns, Nhalf, Nall, Lx, Ly, L0, D0, A0, Da,
				  		Kp, KA, KC, Kb, KAC, KCC, idx_start, idx_end,
				  		ift, jft, pos_all, Nc, Dc, Acc);
	*/
	Acc_max = std::abs(Acc[0]);
	for (n = 1; n < Nall; n++){
		Acc_abs = std::abs(Acc[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	if (Acc_max < Fthresh){
		//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return Pw;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++){
		P = 0.0;
		for (n = 0; n < Nall; n++)
			P += Acc[n] * Vel[n];
		if (P > 0.0){
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay){
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else {
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay)){
				if (f_dec * dt > dt_min){
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
				for (n = 0; n < Nall; n++){
					pos[n] -= dt_half * Vel[n];
					Vel[n] = 0.0;
				}
			}
		}

		// MD using Verlet method
		for (n = 0; n < Nall; n++)
			Vel[n] += dt_half * Acc[n];
		Vel_norm = 0.0;
		Acc_norm = 0.0;
		for (n = 0; n < Nall; n++){
			Vel_norm += Vel[n] * Vel[n];
			Acc_norm += Acc[n] * Acc[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Nall; n++){
			Vel[n] = delta_a * Vel[n] + v_rsc * Acc[n];
			pos[n] += dt * Vel[n];
		}
		
		VerletList(Na, Ns, Nhalf, Nall, Da, r_cut, idx_start, idx_end,
			   		pos, pos_save, VL_aa, VL_counter_all, first_call);
		Pw = Force_DPM(Na, Ns, Nhalf, Nall, Rb, L0, D0, A0, Da, alpha_0, a_cut_c,
						Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
						pos, VL_aa, VL_counter_all, Acc);
		/*
		Eval = Energy_DPM(Na, Ns, Nhalf, Nall, Lx, Ly, L0, D0, A0, Da,
				  		Kp, KA, KC, Kb, KAC, KCC, idx_start, idx_end,
				  		ift, jft, pos_all, Nc, Dc, Acc);
		*/
		for (n = 0; n < Nall; n++)
			Vel[n] += dt_half * Acc[n];

		Acc_max = std::abs(Acc[0]);
		for (n = 1; n < Nall; n++){
			Acc_abs = std::abs(Acc[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		if (Acc_max < Fthresh)
			break;
		// if (nt % 1 == 0) printf("nt: %d  Max Force: %.5e\n", nt, Acc_max);
	}
	//printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
	return Pw;
}