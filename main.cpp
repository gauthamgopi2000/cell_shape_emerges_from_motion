#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <stack>

#include "FIRE.h"
#include "main.h"
#include "VerletList.h"
#include "Energy_DPM.h"

using std::iota;
using std::vector;

int main(int argc, char *argv[])
{

	if (argc < 8)
	{
		fprintf(stderr, "Usage: %s KA KB Kp B Nt F0 tau \n", argv[0]);
		exit(1);
	}

	FILE *otptfl1;
	char dongdir[300], outfile1[300], restartfile[300];

	// DP PARAMETERS
	int Na = 50;  // Number of Dps
	int Nss = 50; // Number of vertices
	double delta = 1.1;
	int Alf_int = 100;
	double l0 = 1.;
	double d0 = delta*l0;
	double kal_A = double(Alf_int) / 100.;
	double KAA_a_ratio = .01;
	double a_cut_c = 1.;
	// double a_cut_c = 1.;
	int pid = 1;
	int Nsb = Nss;
	int Ntot = (Nss + Nsb) * Na;
	int Nhalf = Ntot / 2;
	int Na_half = Na / 2;
	double Alf_set = (double)Alf_int / 100.0;
	double D0_c = 1.;
	double Rb = 10000.;
	double a_cut_c0 = 1.;
	double P;
	double Alf_b = (double)Nsb / PI * std::tan(PI / (double)Nsb);
	double Alf_s = (double)Nss / PI * std::tan(PI / (double)Nss);

	// SPRING CONSTANTS
	double KA = (double)atof(argv[1]);
	double KAA = 1000.;
	double KAA_a;
	double Kw = KAA / 4.0;
	double Kb = (double)atof(argv[2]);
	int bending_added = 1; // is there bending force
	int bending_type = 2;
	double Kp = (double)atof(argv[3]); // 1   perimeter spring
	vector<double> Kl(Na, Kp);
	KAA_a = KAA_a_ratio * KAA;

	sprintf(dongdir, "./");

	// PLASTICITY

	int perimeter_plasticity_flag = 0; // should the perimeter be plastic or not
	int whole_dp_plasticity = 1;
	int area_plasticity_flag = 0;
	int area_plasticity_all = 1; // is area plasticity on all
	double drag_l = (double)atof(argv[8]);			 // drag on perimeter plasticity
	double l0min = 0.75;		 // lower bound on l0
	double drag_a = 10000.;		 // drag on area plasticity
	double l0max = 1.5;			 // upper bound on l0

	// OTHER SIMULATION PARAMETERS
	int resume_simulation = 1; // whether to resume a simulation from a restart file
	int dstep = 100000;		   // how often to output coords
	int restartstep = 500000;  // how ofter to output restart files
	int global_step = 2000;
	double B = (double)atof(argv[4]); // drag
	double Ti = 0.0;				  // temperature
	int temperature_on = 0;			  // should temperature be on
	double T = Ti;
	int param_override = 1; // whether these input parameters should override those from
	vector<int> inflate_ids{};
	int flag_deflate = 0; // should deflation happen after ablation occurs.
	vector<int> deflation_ids = {};
	int n_deflate = deflation_ids.size();
	double kA0 = -1.0e-3;
	double box_length;
	double phi_crit = 0.62;
	double dt_verlet = 0.01;
	double dt_fire = 0.001;
	int Nt_verlet = (int)atoi(argv[5]);
	int Nt_fire = 100000;
	int com = 1; // correct for center of mass velocity

	// ACTIVE BROWNIAN PARTICLES
	int active_yes = 1;					// is activity on
	double f0 = (double)atof(argv[6]);	// driving force for active brownian particles
	double tau = (double)atof(argv[7]); // persistence time
	int activity_onset_time = -1;		// set to -1 for no activity

	// JAMMING PARAMETERS
	double phi_start = 0.4; // starting packing fraction before growth
	box_length = std::sqrt((double)Na * ((double)Nss * (double)Nss * l0 * l0 / (4. * PI * Alf_set) + (PI / 8.) *d0*d0*((double)Nss + 2)) / phi_start);
	// box_length = 1000.;
	int jamming_yes = 0;  // is this a jamming simulation
	int flag_inflate = 0; // should inflation of particles happen
	// double pressure_crit = 1.0e-3;//current value for jamming algorithm
	double pressure_crit = 2.0e-2;
	double KE_thresh = 1.0e-10;
	// double dA0inf = 1.0e-4; //current value for jamming algorithm
	double dA0inf = 2.0e0;
	double phi_slow = 0.85;
	double checking_time = 40000;
	double max_relaxation_time = 1000;
	double dphi = 0.0001;

	char Alf_str[200];
	std::ofstream ofile("./coords.csv"); // Gautham add
	int t = 0;
	sprintf(Alf_str, "%02A1P%d", Alf_int - 100);
	char cmd[300];
	// sprintf(cmd, "mkdir -p %s%s", dongdir, Alf_str);
	std::system(cmd);
	sprintf(outfile1, "%s%s/Pos_%05d_ac_%.2f_Ka_%.2f_Kp_%.2f_N_%04d_Kb_%.2f.txt", dongdir, Alf_str, pid, a_cut_c, KAA_a_ratio, Kp, Na, Kb);
	std::ifstream pfile(outfile1);
	sprintf(restartfile, "restart_regular_KAA_1000.00_KAA_a_10.00_KA_100000.00_Kp_1.00_Kb_0.004_Na_50_t_3400999_phi_0.740.txt"); // name of restart file to open
	std::ifstream rfile(restartfile);
	vector<double> Alf(Na, 0.0);
	vector<double> Da(Na, 0.0), L0(Na, 0.0), D0(Na, 0.0), alpha_0(Nhalf, 0.0), L0i(Nhalf, 0.0);
	double A0;
	vector<int> Ns(Na, 0);
	vector<int> idx_start(Na, 0), idx_end(Na, 0);
	vector<int> ift(Nhalf, 0), jft(Nhalf, 0);
	int countstart = 0, countend;
	int i, na, ns;
	vector<double> pos(Ntot, 0.0);
	vector<double> Vel(Ntot, 0.0);
	vector<double> theta(Na, 0.0);
	// initializing theta with random angles
	for (na = 0; na < Na; na++)
	{
		theta[na] = 2. * PI * ((double)rand()) / RAND_MAX;
	}

	double Fthresh = std::pow(10.0, -10);
	double Pt = std::pow(10.0, -5); // 10-5

	for (na = 0; na < Na; na++)
	{
		// if (na < Na_half)
		//{
		Ns[na] = Nsb;
		Da[na] = D0_c / std::sin(PI / Nsb);
		// D0[na] = delta * D0_c;
		// L0[na] = D0_c;

		A0 = (double)Nss * (double)Nss * l0 * l0 / (4.0 * PI * Alf_set);
		// A0 = (1./(1.+PI*PI*delta*delta*Alf_set*(double(Nsb)+2.)/(2.*(double)Nsb*(double)Nsb)))*phi_start * box_length * box_length / (double)Na;
		// starting with a constant effective packing fraction
		// L0[na] = std::sqrt(4. * PI * A0 * Alf_set) / (double)Nsb;
		L0[na] = l0;
		D0[na] = delta * L0[na];
		// alpha_0[na] = 2.0 * PI / Nsb;
		if (na > 0)
			countstart += Nsb;
		countend = countstart + Nsb - 1;
		idx_start[na] = countstart;
		idx_end[na] = countend;

		for (ns = countstart; ns <= countend; ns++)
		{
			if (ns == countstart)
			{
				ift[ns] = ns + 1;
				jft[ns] = countend;
			}
			else if (ns == countend)
			{
				ift[ns] = countstart;
				jft[ns] = ns - 1;
			}
			else
			{
				ift[ns] = ns + 1;
				jft[ns] = ns - 1;
			}
		}
		//}
		// else
		// {
		// 	Ns[na] = Nss;
		// 	Da[na] = D0_c / std::sin(PI / Nss);
		// 	//D0[na] = delta * D0_c;
		// 	//L0[na] = D0_c;
		// 	//A0[na] = Nss * Nss * D0_c * D0_c / 4.0 / PI / Alf_set / Alf_s;
		// 	// alpha_0[na] = 2.0 * PI / Nss;
		// 	A0[na] = phi_start*box_length*box_length/(double)Na;
		// 	L0[na] = std::sqrt(4.*PI*A0[na]*Alf_set)/(double)Nsb;
		// 	D0[na] = delta*L0[na];
		// 	if (na == Na_half)
		// 		countstart += Nsb;
		// 	else
		// 		countstart += Nss;
		// 	countend = countstart + Nss - 1;
		// 	idx_start[na] = countstart;
		// 	idx_end[na] = countend;
		// 	for (ns = countstart; ns <= countend; ns++)
		// 	{
		// 		if (ns == countstart)
		// 		{
		// 			ift[ns] = ns + 1;
		// 			jft[ns] = countend;
		// 		}
		// 		else if (ns == countend)
		// 		{
		// 			ift[ns] = countstart;
		// 			jft[ns] = ns - 1;
		// 		}
		// 		else
		// 		{
		// 			ift[ns] = ns + 1;
		// 			jft[ns] = ns - 1;
		// 		}
		// 	}
		// }
	}

	for (na = 0; na < Na; na++)
	{
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			alpha_0[ns] = 2.0 * PI / Ns[na];
		}
	}

	for (int na = 0; na < Na; na++)
	{
		for (int ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			L0i[ns] = L0[na];
		}
	}

	if (resume_simulation == 1)
	{

		if (rfile.good())
		{
			rfile >> Rb;
			for (int ntot = 0; ntot < Ntot; ntot++)
			{
				rfile >> pos[ntot];
			}
			for (int ntot = 0; ntot < Ntot; ntot++)
			{
				rfile >> Vel[ntot];
			}
			for (int na = 0; na < Na; na++)
			{
				rfile >> A0;
			}
			for (int na = 0; na < Na; na++)
			{
				for (ns = idx_start[na]; ns <= idx_end[na]; ns++){

					rfile >> L0i[ns];
					

				}
				L0[na]  = L0i[idx_start[na]];// made this way to allow backwards compatibility for reading restart files. 
				
		
			}
			for (int na = 0; na < Na; na++)
			{
				rfile >> D0[na];
			}

			rfile >> box_length;
		}

		// double F = FIRE_DPM_RP_VL(Na, Ns, Nhalf, Ntot, Rb, D0, A0, Da, alpha_0, a_cut_c, Kp, KA, KAA, KAA_a, Kw, idx_start,
		// 					idx_end, ift, jft, pos, Fthresh, dt_fire, Nt_fire, L0i, Kb, dstep, bending_added,
		// 					box_length, bending_type, l0);

		P = Verlet_Int(Na, Ns, Nhalf, Ntot, Rb, D0, A0, Da, alpha_0, a_cut_c,
					   Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
					   pos, Fthresh, dt_verlet, Nt_verlet, Vel, dstep, restartstep, Ti, B,
					   inflate_ids, L0, drag_l, Alf_set, deflation_ids, n_deflate, drag_a, flag_inflate, flag_deflate,
					   kA0, Kb, kal_A, bending_added, perimeter_plasticity_flag, area_plasticity_flag,
					   l0min, area_plasticity_all, whole_dp_plasticity, box_length, delta, temperature_on, pressure_crit,
					   active_yes, f0, tau, theta, activity_onset_time, global_step, phi_crit, KE_thresh, dA0inf,
					   checking_time, phi_slow, max_relaxation_time, jamming_yes, l0max, dphi, phi_start, com, bending_type, l0);
		return 0;
	}
	else
	{

		// ARRANGING PARTICLES ON A RECTANGULAR LATTICE

		vector<double> cell_centers(2 * Na, 0.0);
		double cushion = 10.;
		double rminimum = std::sqrt(A0 / PI);
		double min_distance2 = 4. * rminimum * rminimum + 4. * rminimum * D0[0] + D0[0] * D0[0] + cushion;
		bool overlap;

		for (int na = 0; na < Na; na++)
		{

			do
			{
				overlap = false;
				cell_centers[na] = box_length * (((double)rand()) / RAND_MAX - 0.5);
				cell_centers[na + Na] = box_length * (((double)rand()) / RAND_MAX - 0.5);
				int nb = 0;
				while (nb < na && !overlap)
				{
					double dx = pbc_separation(cell_centers[na] - cell_centers[nb], box_length);
					double dy = pbc_separation(cell_centers[na + Na] - cell_centers[nb + Na], box_length);
					if (dx * dx + dy * dy < min_distance2)
					{
						overlap = true;
					}
					nb++;
				}
			} while (overlap);
		}
		double dgamma = 2. * PI / (double)Nss;

		for (int na = 0; na < Na; na++)
		{

			double gamma = 2. * PI * ((double)rand()) / RAND_MAX;

			for (int ns = idx_start[na]; ns <= idx_end[na]; ns++)
			{

				pos[ns] = cell_centers[na] + rminimum * cos(gamma);
				pos[ns + Nhalf] = cell_centers[na + Na] + rminimum * sin(gamma);
				gamma += dgamma;
			}
		}
	}

	double F = FIRE_DPM_RP_VL(Na, Ns, Nhalf, Ntot, Rb, D0, A0, Da, alpha_0, a_cut_c, Kp, KA, KAA, KAA_a, Kw, idx_start,
							  idx_end, ift, jft, pos, Fthresh, dt_fire, Nt_fire, L0, Kb, dstep, bending_added,
							  box_length,bending_type, l0);

	P = Verlet_Int(Na, Ns, Nhalf, Ntot, Rb, D0, A0, Da, alpha_0, a_cut_c,
				   Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
				   pos, Fthresh, dt_verlet, Nt_verlet, Vel, dstep, restartstep, Ti, B,
				   inflate_ids, L0, drag_l, Alf_set, deflation_ids, n_deflate, drag_a, flag_inflate, flag_deflate,
				   kA0, Kb, kal_A, bending_added, perimeter_plasticity_flag, area_plasticity_flag,
				   l0min, area_plasticity_all, whole_dp_plasticity, box_length, delta, temperature_on, pressure_crit,
				   active_yes, f0, tau, theta, activity_onset_time, global_step, phi_crit, KE_thresh, dA0inf, checking_time,
				   phi_slow, max_relaxation_time, jamming_yes, l0max, dphi, phi_start, com, bending_type, l0);

	return 0;
}

// FUNCTIONS

double pbc_separation(double dl, double box_length)
{
	dl -= box_length * round(dl / box_length);

	return dl;
}

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
				  double kal_A, int bending_added, int perimeter_plasticity_flag, int area_plasticity_flag, double l0min, int area_plasticity_all, int whole_dp_plasticity, double box_length, double delta, int temperature_on, double pressure_crit,
				  int active_yes, double f0, double tau, vector<double> &theta, int activity_onset_time, int global_step,
				  double phi_crit, double KE_thresh, double dA0inf, int checking_time, double phi_slow, int max_relaxation_time,
				  int jamming_yes, double l0max, double dphi, double phi_start, int com, int bending_type, double l0)
{
	int n, nt = 0;
	vector<double> Acc(Nall, 0.0);
	double dt = dt_fire;
	double dt_half = dt / 2.0;
	double dt_root = std::sqrt(dt);
	double root_3 = std::sqrt(3.);
	// Verlet list parameters
	int first_call = 1;
	vector<vector<int>> VL_aa(100 * Nhalf, vector<int>(4, -1));
	vector<int> VL_counter_all(1, 0);
	vector<double> pos_save(Nall, 0.0);
	vector<double> Vel_half(Nall, 0.0);
	vector<double> lk(Nhalf, 0.0);
	vector<double> A(Na, 0.0);
	vector<double> Kl(Na, Kp);
	double r_cut = D0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > D0[n]) ? r_cut : D0[n];
	r_cut *= a_cut_c;

	double Pw;
	double deltal0;

	vector<double> Acc0(Nall, 0.0);
	vector<double> Acc_mid(Nall, 0.0);
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.);
	double T = Ti;
	double sigma = sqrt(2. * T * B);
	double eta1, eta2, xi1, xi2;
	vector<double> eta(Nall, 0.0);
	vector<double> xi(Nall, 0.0);

	char restart_ablate[200];
	char restart_folder[200];
	sprintf(restart_ablate, "restart_ablate.txt");
	sprintf(restart_folder, "restart_files");
	char cmd[300];
	sprintf(cmd, "mkdir -p %s", restart_folder);
	sprintf(cmd, "rm -r %s", restart_folder);
	sprintf(cmd, "mkdir -p %s", restart_folder);
	std::system(cmd);
	FILE *otptfl1;
	char coords[200];
	sprintf(coords, "coords.csv");
	FILE *coords_file;
	char restart_regular[200];

	FILE *regular;
	char restart_wall[200];
	sprintf(restart_wall, "restart_wall2.txt");
	char log[200]; // to output simulation parameters
	sprintf(log, "log.csv");
	FILE *logfile;

	// creating array to store interparticle forces for pressure
	vector<double> aforces(Nall, 0.0);

	VerletList(Na, Ns, Nhalf, Nall, Da, r_cut, idx_start, idx_end,
			   pos, pos_save, VL_aa, VL_counter_all, first_call, box_length);

	Pw = Force_DPM(Na, Ns, Nhalf, Nall, Rb, D0, A0, Da, a_cut_c,
				   Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
				   pos, VL_aa, VL_counter_all, Acc, aforces, L0, lk, A, Kb, nt, dstep, bending_added,
				   box_length, bending_type, Kl);

	int out_t = 0;

	std::ofstream ofile2("./pressure.csv");
	std::ofstream ofile3("./energies.csv");
	first_call = 0; // it is no longer first call of verlets

	int energy_step = 10000;

	// double time_diff  = double(Nt_fire);

	double dTemp = 0.;

	double dkb = 1.0e-5;
	double kbmax = 10.;
	int bend_increase_flag = 0;
	double press;
	double G;
	double phi; // packing fraction
	int time_elapsed = 0;
	int time_stationary = 0;
	double phi_curr = phi_start;

	// setting initial velocity

	for (int nt = 0; nt < Nt_verlet; nt++)

	{

		// outputting simulation data every dstep timesteps

		if (nt % dstep == 0)
		{
			if (nt == 0)
			{
				coords_file = fopen(coords, "w");

				for (int i = 0; i < Nhalf; i++)
				{
					int type = (int)(i + Ns[0]) / Ns[0];
					fprintf(coords_file, "%d, %d, %d, %f, %f, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f\n", nt, i + 1, type, pos[i], pos[i + Nhalf], 0.,
							D0[type - 1] / 2., Nhalf, A[type - 1], A0, lk[i], L0[type-1], alpha_0[i], Kb, box_length);
				}
				logfile = fopen(log, "w");
				fprintf(logfile, "KAA, KAA_a, a_cut_c, B, KA, Kp, Kb, f0, tau, drag_l, l0min, l0max, dA0inf, phi_slow, checking_time, max_relaxation_time\n");
				fprintf(logfile, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d\n",
						KAA, KAA_a, a_cut_c, B, KA, Kp, Kb, f0, tau, drag_l, l0min, l0max, dA0inf, phi_slow, checking_time, max_relaxation_time);
				fclose(logfile);
			}

			else
			{
				coords_file = fopen(coords, "a");
				out_t += 1;

				for (int i = 0; i < Nhalf; i++)
				{
					int type = (int)(i + Ns[0]) / Ns[0];
					fprintf(coords_file, "%d, %d, %d, %f, %f, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f\n", nt, i + 1, type, pos[i], pos[i + Nhalf], 0.,
							D0[type - 1] / 2., Nhalf, A[type - 1], A0, lk[i], L0[type-1], alpha_0[i], Kb, box_length);
				}
			}
			fclose(coords_file);
		}
		double KE = 0.;
		if (nt % energy_step == 0)
		{

			ofile3 << nt << ", " << KE << ", " << Pw << ", " << KE + Pw << ", " << T << std::endl;
		}

		r_cut = D0[0];

		for (int na = 1; na < Na; na++)
		{

			r_cut = (r_cut > D0[na]) ? r_cut : D0[na];
		}

		r_cut *= a_cut_c;

		KE = 0.;

		sigma = sqrt(2. * T * B);

		// outputting simulation data every few timesteps
		if (nt % restartstep == 0)
		{

			sprintf(restart_regular, "%s/restart_regular_KAA_%.2f_KAA_a_%.2f_KA_%.2f_Kp_%.2f_Kb_%.3f_Na_%d_t_%d.txt", restart_folder, KAA, KAA_a, KA, Kp, Kb, Na, nt);

			restart_write(restart_regular, Rb, pos, Vel, A0, L0, KA, Kp, Na, Nall, Nhalf, Kb, D0, box_length);
		}

		if (jamming_yes)
		{

			flag_inflate = 0;

			KE = 0.;

			for (int n = 0; n < Nall; n++)
			{
				KE += 0.5 * Vel[n] * Vel[n];
			}

			KE /= (double)Nhalf;

			press = 0.;

			for (int n = 0; n < Nall; n++)
			{
				press -= 0.5 * pos[n] * aforces[n] / (box_length * box_length);
			}

			// if ((KE < KE_thresh) && (press < pressure_crit))
			// {
			// 	flag_inflate = 1;
			// 	time_stationary = 0;
			// }
			// else
			// {
			// 	flag_inflate = 0;
			// }

			if (flag_inflate == 0)
			{
				time_stationary += 1;
				if (time_stationary == max_relaxation_time)
				{
					flag_inflate = 1;
					time_stationary = 0;
				}
			}
		}

		if (flag_inflate == 1)
		{
			sprintf(restart_regular, "%s/restart_regular_KAA_%.2f_KAA_a_%.2f_KA_%.2f_Kp_%.2f_Kb_%.3f_Na_%d_t_%d_phi_%.3f.txt", restart_folder, KAA,
					KAA_a, KA, Kp, Kb, Na, nt, phi_curr);
			restart_write(restart_regular, Rb, pos, Vel, A0, L0, KA, Kp, Na, Nall, Nhalf, Kb, D0, box_length);
			phi_curr += dphi;

			double l0 = L0[0]; //only makes sense for jamming simulations where all L0s are the same
			double nvertex = Ns[0]; // assuming all dps have same no. of vertices.
			double d0 = delta*l0;

			// A0 = (1./(1.+PI*PI*delta*delta*Alf_set*(double(Nsb)+2.)/(2.*(double)Nsb*(double)Nsb)))*phi_curr * box_length * box_length / (double)Na;
			box_length = std::sqrt((double)Na * (nvertex * nvertex * l0 * l0 / (4. * PI * Alf_set) + (PI / 8.) * d0*d0*(nvertex + 2)) / phi_curr);
			first_call = 1; // rebuild verlet list if box becomes smaller.

			// A0[na] += dA0inf;
			// double perimeter = std::sqrt(4. * PI * Alf_set * A0);

			// double l0curr = perimeter / nvertex;
			// for (int ns = idx_start[na]; ns <= idx_end[na]; ns++)
			// {
			// 	L0i[ns] = l0curr;
			// }

			// D0[na] = l0curr * delta;

			time_elapsed += 1;
		}

		/////TEMPORARY//////
		// press = 0.;

		// for (int n = 0; n < Nall; n++)
		// {
		// 	press -= 0.5 * pos[n] * aforces[n] / (box_length * box_length);
		// }

		// if (press > pressure_crit)
		// {
		// 	flag_inflate = 0;
		// 	perimeter_plasticity_flag = 1;
		// }

		// calculating pressure from virial

		if (nt % global_step == 0)
		{
			///////////////////

			// kinetic energy calculation

			press = 0.;
			phi = 0.;

			for (int na = 0; na < Na; na++)
			{
				phi += A0 + 0.125 * (Ns[na] + 2.) * PI * D0[na] * D0[na];
			}
			phi /= (box_length * box_length);

			for (int n = 0; n < Nall; n++)
			{
				press -= 0.5 * pos[n] * aforces[n] / (box_length * box_length);
			}
			ofile2 << nt << ", " << press << ", " << phi << std::endl;
		}
		if (nt == activity_onset_time)
		{

			active_yes = 1;
		}

		int particle_id = 0;

		// if (temperature_on)
		// {

		// 	for (int n = 0; n < Nhalf; n++)
		// 	{

		// 		if (n % Ns[0] == 0)
		// 		{ /// asuming all dps have same no. of vertices

		// 			eta1 = distribution(generator);
		// 			eta2 = distribution(generator);
		// 			xi1 = distribution(generator);
		// 			xi2 = distribution(generator);
		// 		}

		// 		xi[n] = xi1;
		// 		xi[n + Nhalf] = xi2;
		// 		eta[n] = eta1;
		// 		eta[n + Nhalf] = eta2;
		// 	}
		// }

		for (int n = 0; n < Nall; n++)
		{

			// cicotti integrator
			Vel_half[n] = Vel[n] + 0.5 * dt * Acc[n] - 0.5 * dt * B * Vel[n] + 0.5 * dt_root * sigma * xi[n] - (1. / 8.) * dt * dt * B * (Acc[n] - B * Vel[n]) - 0.25 * dt * dt_root * B * sigma * (0.5 * xi[n] + (1. / root_3) * eta[n]);

			pos[n] = pos[n] + dt * Vel_half[n] + dt * dt_root * sigma * (1. / (2. * root_3)) * eta[n];
		}
		if (perimeter_plasticity_flag)
		{

			for (int na = 0; na < Na; na++)
			{

				int start = idx_start[na];
				int end = idx_end[na];
				double l0curr;
				for (int lid = start; lid <= end; lid++)
				{

					deltal0 = -(Kp / drag_l) * (L0[na] - lk[lid]) * dt;

					// if(deltal0<0.){
					// 	if(L0i[lid]>l0min){
					// 		L0i[lid] += deltal0;   //always allow the perimeter to get shorter
					// 	}
					// }
					// else{
					// 	if (L0i[lid]<l0max){

					// 		L0i[lid] += deltal0;
					// 	}
					// }
				}
			}
		}
		else
		{
			if (whole_dp_plasticity)
			{

				for (int na = 0; na < Na; na++)
				{

					int start = idx_start[na];
					int end = idx_end[na];
					double l0curr;
					for (int lid = start; lid <= end; lid++)
					{

						deltal0 = (Kp / drag_l) * (lk[lid]-L0[na]) * dt;

						L0[na] += deltal0 / Ns[na]; 

						//D0[na] += deltal0 / Ns[na];
					}
					Kl[na] = L0[na]*L0[na]*Kp;
				}
			}
		}

		Pw = Force_DPM(Na, Ns, Nhalf, Nall, Rb, D0, A0, Da, a_cut_c,
					   Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
					   pos, VL_aa, VL_counter_all, Acc, aforces, L0, lk, A, Kb, nt, dstep, bending_added,
					   box_length, bending_type, Kl);

		if (active_yes)
		{

			for (int na = 0; na < Na; na++)
			{

				for (int ns = idx_start[na]; ns <= idx_end[na]; ns++)
				{

					Acc[ns] += f0 * std::cos(theta[na]);
					Acc[ns + Nhalf] += f0 * std::sin(theta[na]);

				}
				G = distribution(generator);
				theta[na] += std::sqrt(2. * dt / tau) * G;
			}
		}

		VerletList(Na, Ns, Nhalf, Nall, Da, r_cut, idx_start, idx_end,
				   pos, pos_save, VL_aa, VL_counter_all, first_call, box_length);
		first_call = 0;

		for (int n = 0; n < Nall; n++)
		{
			Vel[n] = Vel_half[n] + 0.5 * dt * Acc[n] - 0.5 * dt * B * Vel_half[n] + 0.5 * dt_root * sigma * xi[n] - (1. / 8.) * dt * dt * B * (Acc[n] - B * Vel_half[n]) - 0.25 * dt * dt_root * B * sigma * (0.5 * xi[n] + (1. / root_3) * eta[n]);
		}

		if (com)
		{
			// correct for center of mass velocity
			double com_vel_x = 0.;
			double com_vel_y = 0.;
			for (int n = 0; n < Nhalf; n++)
			{
				com_vel_x += Vel[n];
				com_vel_y += Vel[n + Nhalf];
			}
			com_vel_x /= (double)Nhalf;
			com_vel_y /= (double)Nhalf;
			for (int n = 0; n < Nall; n++)
			{
				Vel[n] -= com_vel_x;
				Vel[n + Nhalf] -= com_vel_y;
			}
		}
	}

	ofile2.close();

	return Pw;
}

double Force_DPM(int Na, vector<int> &Ns, int Nhalf, int Nall,
				 double Rb, vector<double> &D0, double A0,
				 vector<double> &Da, double a_cut_c,
				 double Kp, double KA, double KAA, double KAA_a, double Kw,
				 vector<int> &idx_start, vector<int> &idx_end,
				 vector<int> &ift, vector<int> &jft,
				 vector<double> &pos, vector<vector<int>> &VL_aa,
				 vector<int> &VL_counter, vector<double> &F_a, vector<double> &aforces,
				 vector<double> &L0, vector<double> &lk, vector<double> &A, double Kb, int nt,
				 int dstep, int bending_added, double box_length, int bending_type, vector<double> &Kl)
{
	int ns, ns_y, ms, na, ma;
	double lxi, lyi, lki, dli, dlj, deltal;
	double F;
	double Pwall = 0.0;

	vector<double> xa(Nhalf, 0.0);
	vector<double> ya(Nhalf, 0.0);
	vector<double> lx(Nhalf, 0.0);
	vector<double> ly(Nhalf, 0.0);
	vector<double> lk_sq(Nhalf, 0.0);
	// vector<double> A(Na, 0.0);
	for (int na = 0; na < Na; na++)
	{
		A[na] = 0.0;
	}
	vector<double> dxA(Nhalf, 0.0);
	vector<double> dyA(Nhalf, 0.0);
	// vector<double> lk(Nhalf, 0.0);
	double dAa;

	int ifts, jfts, ifts_next, jfts_next;
	double d0;
	int countstart = 0;
	int countend = 0;
	double Eval = 0.;
	double l0curr;

	for (ns = 0; ns < Nhalf; ns++)
	{
		xa[ns] = pos[ns];
		ya[ns] = pos[ns + Nhalf];
		F_a[ns] = 0.0;
		F_a[ns + Nhalf] = 0.0;
		aforces[ns] = 0.0;
		aforces[ns + Nhalf] = 0.0;
	}

	// Perimeter Force
	double L0_a, L0_sq;
	double Ns_a;
	for (na = 0; na < Na; na++)
	{
		// L0_a = L0[na];
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			ifts = ift[ns];
			lxi = xa[ifts] - xa[ns];
			lyi = ya[ifts] - ya[ns];
			lx[ns] = lxi;
			ly[ns] = lyi;
			lk_sq[ns] = lxi * lxi + lyi * lyi;
			lk[ns] = std::sqrt(lk_sq[ns]);
			// dli = lki / L0_c - 1.0;
			// Eval += Kp * dli * dli / double(Ns[nc]) / 2.0;
		}
	}

	for (na = 0; na < Na; na++)
	{
		Ns_a = double(Ns[na]);
		// L0_a = L0[na];
		double spring_const = Kl[na];
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{

			jfts = jft[ns];
			double L0r = L0[na];
			double L0l = L0[na];
			dlj = (L0r / lk[jfts] - 1.0) / (L0r * L0r);
			dli = (1.0 - L0l / lk[ns]) / (L0l * L0l);
			deltal = (1.0 - lk[ns] / L0l);
			Eval += 0.5 * Kp * deltal * deltal / Ns_a;
			F_a[ns] += spring_const / Ns_a * (dlj * lx[jfts] + dli * lx[ns]);
			F_a[ns + Nhalf] += spring_const / Ns_a * (dlj * ly[jfts] + dli * ly[ns]);
		}
	}

	double A0_a, KA_a, grad_x, grad_y, dxx, dyy, xcurr, xcurri, ycurr, ycurri;

	// Area Force

	for (na = 0; na < Na; na++)
	{
		A0_a = A0;
		KA_a = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * KA;
		// printf("na: %d  Ka: %.5e\n", na, KA_a);
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			ifts = ift[ns];
			jfts = jft[ns];
			dxx = xa[ns] - xa[ifts];
			dyy = ya[ns] - ya[ifts];
			xcurr = xa[ns];
			xcurri = xa[ifts];
			ycurr = ya[ns];
			ycurri = ya[ifts];

			if (dxx > 0.5 * box_length)
			{
				xcurri += box_length;
			}
			else if (dxx < -0.5 * box_length)
			{
				xcurri -= box_length;
			}
			if (dyy > 0.5 * box_length)
			{
				ycurri += box_length;
			}
			else if (dyy < -0.5 * box_length)
			{
				ycurri -= box_length;
			}

			A[na] += (xcurr * ycurri - ycurr * xcurri);
			grad_x = xa[ifts] - xa[jfts];
			grad_y = ya[jfts] - ya[ifts];
			grad_x = pbc_separation(grad_x, box_length);
			grad_y = pbc_separation(grad_y, box_length);
			// dxA[ns] = (xa[ifts] - xa[jfts]) / A0_a;
			// dyA[ns] = (ya[jfts] - ya[ifts]) / A0_a;
			dxA[ns] = grad_x / A0_a;
			dyA[ns] = grad_y / A0_a;
		}
		A[na] /= 2.0;
		dAa = 0.5 * KA_a * (A[na] / A0_a - 1.0);
		Eval += dAa * (A[na] / A0_a - 1.0);
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			F_a[ns] += dAa * dyA[ns];
			F_a[ns + Nhalf] += dAa * dxA[ns];
		}
	}
	if (bending_added)
	{
		if (bending_type == 1)
		{

			for (na = 0; na < Na; na++)
			{

				for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
				{
					l0curr = L0[na];
					ifts = ift[ns];
					ifts_next = ift[ifts];
					jfts = jft[ns];
					jfts_next = jft[jfts];
					F_a[ns] -= Kb * (-4. * (xa[ifts] + xa[jfts]) + 6. * xa[ns] + xa[ifts_next] + xa[jfts_next]) / (l0curr * l0curr * l0curr * l0curr);
					F_a[ns + Nhalf] -= Kb * (-4. * (ya[ifts] + ya[jfts]) + 6. * ya[ns] + ya[ifts_next] + ya[jfts_next]) / (l0curr * l0curr * l0curr * l0curr);
				}
			}
		}
		else
		{
			vector<double> theta(Nhalf, 0.0);
			vector<double> d_alpha(Nhalf, 0.0);
			vector<double> cos_theta(Nhalf, 0.0);
			vector<double> sin_theta(Nhalf, 0.0);
			vector<double> alpha(Nhalf, 0.0);
			double alphai, d_alpha_i, d_alpha_j, Kb_s, alphai0;
			for (ns = 0; ns < Nhalf; ns++)
			{
				theta[ns] = std::atan2(ly[ns], lx[ns]);
				cos_theta[ns] = lx[ns] / lk_sq[ns];
				sin_theta[ns] = ly[ns] / lk_sq[ns];
			}

			for (na = 0; na < Na; na++)
			{

				for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
				{
					alphai = theta[ns] - theta[jft[ns]];
					alpha[ns] = alphai - std::round(alphai / 2.0 / PI) * 2.0 * PI;
					d_alpha[ns] = alphai - std::round(alphai / 2.0 / PI) * 2.0 * PI;
				}
			}

			for (na = 0; na < Na; na++)
			{

				for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
				{
					Kb_s = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * Kb;
					// printf("ns: %d  Kb: %.5e\n", ns, Kb_s);
					ifts = ift[ns];
					jfts = jft[ns];
					d_alpha_i = d_alpha[ns] - d_alpha[ifts];
					d_alpha_j = d_alpha[ns] - d_alpha[jfts];
					F_a[ns] -= Kb_s * (d_alpha_i * sin_theta[ns] + d_alpha_j * sin_theta[jfts]);
					F_a[ns + Nhalf] += Kb_s * (d_alpha_i * cos_theta[ns] + d_alpha_j * cos_theta[jfts]);
				}
			}
		}
	}

	// Cell-cell Force

	double Dnm, Dnm_sq, Dnm_a, Dnm_h, Dnm_a_sq;
	double dFx, dFy;
	double dx13, dy13;
	double dnm, dd;
	int vl_idx;
	double a_cut_h = 0.5 + 0.5 * a_cut_c;
	double shift;

	for (vl_idx = 0; vl_idx < VL_counter[0]; vl_idx++)
	{
		ns = VL_aa[vl_idx][0];
		ms = VL_aa[vl_idx][1];
		na = VL_aa[vl_idx][2];
		ma = VL_aa[vl_idx][3];
		Dnm = 0.5 * (D0[na] + D0[ma]);
		Dnm_sq = Dnm * Dnm;
		Dnm_a = a_cut_c * Dnm;
		Dnm_a_sq = Dnm_a * Dnm_a;
		Dnm_h = a_cut_h * Dnm;
		shift = (-0.25 * KAA_a / Dnm_sq) * (2. * Dnm_a_sq - 2. * Dnm * Dnm_a * (1. + a_cut_c) + Dnm_sq * (1. + a_cut_c * a_cut_c));
		// energy shift to make the adhesion potential continuous

		dx13 = xa[ms] - xa[ns];
		dx13 = pbc_separation(dx13, box_length);
		if (std::abs(dx13) < Dnm_a)
		{
			dy13 = ya[ms] - ya[ns];
			dy13 = pbc_separation(dy13, box_length);
			if (std::abs(dy13) < Dnm_a)
			{
				dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
				if (dnm < Dnm_a)
				{
					if (dnm > Dnm_h)
					{
						F = KAA_a * (Dnm_a / dnm - 1.0) / Dnm_sq;
						Eval -= 0.5 * KAA_a * Dnm_a_sq * (1.0 - dnm / Dnm_a) * (1.0 - dnm / Dnm_a) / Dnm_sq;
					}
					else if ((dnm <= Dnm_h) && (dnm > Dnm))
					{
						F = KAA_a * (1.0 - Dnm / dnm) / Dnm_sq;
						Eval += 0.5 * KAA_a * (1.0 - dnm / Dnm) * (1.0 - dnm / Dnm) + shift;
					}
					else
					{
						F = KAA * (1.0 - Dnm / dnm) / Dnm_sq;
						Eval += 0.5 * KAA * (1.0 - dnm / Dnm) * (1.0 - dnm / Dnm) + shift;
					}
					dFx = F * dx13;
					dFy = F * dy13;
					aforces[ns] += dFx;
					aforces[ns + Nhalf] += dFy;
					aforces[ms] -= dFx;
					aforces[ms + Nhalf] -= dFy;
					F_a[ns] += dFx;
					F_a[ns + Nhalf] += dFy;
					F_a[ms] -= dFx;
					F_a[ms + Nhalf] -= dFy;
				}
			}
		}
	}

	// intra-cell interaction
	double Dnm2;
	for (na = 0; na < Na; na++)
	{
		int ignore_number = 5;
		int start = idx_start[na];
		//Dnm = l0; // doesn't really matter so long as there is some exclusion principle
		Dnm = D0[na];
		Dnm_sq = Dnm * Dnm;
		for (ns = idx_start[na]; ns < idx_end[na] - ignore_number; ns++)
		{
			for (ms = ns + ignore_number+1; ms <= idx_end[na]; ms++)
			{
				int difference = ms-ns;
				int distance = Ns[na] - difference;
				if (distance<=ignore_number)
				{
					continue;
				}
				dx13 = xa[ms] - xa[ns];
				if (std::abs(dx13) < Dnm)
				{
					dy13 = ya[ms] - ya[ns];
					if (std::abs(dy13) < Dnm)
					{
						dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
						if (dnm < Dnm)
						{
							F = KAA * (1.0 - Dnm / dnm) / Dnm_sq;
							// dd = 1.0 - dnm / Dnm;
							Eval += 0.5 * KAA * (1.0 - dnm / Dnm) * (1.0 - dnm / Dnm); // cell-cell PE

							dFx = F * dx13;
							dFy = F * dy13;
							F_a[ns] += dFx;
							F_a[ns + Nhalf] += dFy;
							F_a[ms] -= dFx; // 3rd law
							F_a[ms + Nhalf] -= dFy;
						}
					}
				}
			}
		}
	}

	return Eval;
}

double Pressure_DPM(int Na, int Nhalf, double Rb, vector<int> &idx_start, vector<int> &idx_end,
					vector<double> &D0, double Kw, vector<double> &pos)
{
	double Pwall = 0.0;
	double Dnm, Dnm_sq, F;

	// adipocyte-wall force
	double r;
	for (int na = 0; na < Na; na++)
	{
		Dnm = 0.5 * D0[na];
		Dnm_sq = Dnm * Dnm;
		for (int ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			r = std::sqrt(pos[ns] * pos[ns] + pos[ns + Nhalf] * pos[ns + Nhalf]);
			if (r > Rb - Dnm)
			{
				F = Kw * (Dnm + r - Rb) / Dnm_sq / r;
				Pwall += F;
			}
		}
	}

	return Pwall / 2.0 / Rb / PI;
}
void VerletList(int Na, vector<int> &Ns, int Nhalf, int Nall,
				vector<double> &Da, double r_cut,
				vector<int> &idx_start, vector<int> &idx_end,
				vector<double> &pos_a, vector<double> &pos_save,
				vector<vector<int>> &VL_aa, vector<int> &VL_counter,
				int first_call, double box_length)
{
	// distance for making Verlet list
	double r_factor = 1.2;
	double r_cut_sq = r_cut * r_cut;
	double r_list = r_factor * r_cut;
	double r_list_sq = r_list * r_list;
	double r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double dx, dy, dr;
	int i, iy, j, jy;
	int na, ma, ns, ms, nc, mc;

	// first_call = 1; //FOR DEBUGGING

	if (first_call == 0)
	{
		double dr_max = 0.0;
		for (i = 0; i < Nhalf; i++)
		{
			iy = i + Nhalf;
			dx = pos_a[i] - pos_save[i];
			dx = pbc_separation(dx, box_length);
			dy = pos_a[iy] - pos_save[iy];
			dy = pbc_separation(dy, box_length);
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
		if (4.0 * dr_max < r_skin_sq)
		{
			return;
		}
	}

	vector<double> x_cen(Na, 0.0);
	vector<double> y_cen(Na, 0.0);
	for (na = 0; na < Na; na++)
	{
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
		{
			x_cen[na] += pos_a[ns];
			y_cen[na] += pos_a[ns + Nhalf];
		}
		x_cen[na] /= Ns[na];
		y_cen[na] /= Ns[na];
	}

	double x_cen_m, y_cen_m, x_cen_n, y_cen_n;
	double dx_cen, dy_cen;
	double dx_cen_shift, dy_cen_shift;
	double x1, y1, x3, y3;
	double lc1, lc2, dl_n, dl_m;

	VL_counter[0] = 0;

	for (na = 0; na < Na - 1; na++)
	{
		x_cen_n = x_cen[na];
		y_cen_n = y_cen[na];
		for (ma = na + 1; ma < Na; ma++)
		{
			x_cen_m = x_cen[ma];
			y_cen_m = y_cen[ma];
			dx_cen = x_cen_m - x_cen_n;
			dx_cen = pbc_separation(dx_cen, box_length);
			dy_cen = y_cen_m - y_cen_n;
			dy_cen = pbc_separation(dy_cen, box_length);
			// if (std::sqrt(dx_cen * dx_cen + dy_cen * dy_cen) > Da[na] + Da[ma] + r_list)
			// 	continue;

			for (ns = idx_start[na]; ns <= idx_end[na]; ns++)
			{
				x1 = pos_a[ns];
				y1 = pos_a[ns + Nhalf];
				for (ms = idx_start[ma]; ms <= idx_end[ma]; ms++)
				{
					dx = pos_a[ms] - x1;
					dx = pbc_separation(dx, box_length);
					dy = pos_a[ms + Nhalf] - y1;
					dy = pbc_separation(dy, box_length);
					dr = dx * dx + dy * dy;
					if (dr < r_list_sq)
					{
						VL_aa[VL_counter[0]][0] = ns;
						VL_aa[VL_counter[0]][1] = ms;
						VL_aa[VL_counter[0]][2] = na;
						VL_aa[VL_counter[0]][3] = ma;
						VL_counter[0]++;
					}
				}
			}
		}
	}

	// intra-cell interaction
	/*
	for (na = 0; na < Na; na++){
		for (ns = idx_start[na]; ns < idx_end[na] - 1; ns++){
			for (ms = ns + 2; ms <= idx_end[na]; ms++){
				if ((ns == idx_start[na]) && (ms == idx_end[na]))
					continue;
				dx = pos_a[ms] - pos_a[ns];
				dy = pos_a[ms + Nhalf] - pos_a[ns + Nhalf];
				dr = dx * dx + dy * dy;
				if (dr < r_list_sq){
					VL_aa[VL_counter[0]][0] = ns;
					VL_aa[VL_counter[0]][1] = ms;
					VL_aa[VL_counter[0]][2] = na;
					VL_aa[VL_counter[0]][3] = na;
					VL_counter[0]++;
				}
			}
		}
	}
	*/

	for (i = 0; i < Nall; i++)
		pos_save[i] = pos_a[i];
}

void restart_write(char file_name[200], double Rb, vector<double> &pos, vector<double> &Vel, double A0, vector<double> &L0,
				   double KA, double Kp, int Na, int Nall, int Nhalf, double Kb, vector<double> &D0, double box_length)
{
	// this function writes simulation configurations when needed to resume a simulation later.
	FILE *restart_file;
	restart_file = fopen(file_name, "w");
	fprintf(restart_file, "%.32e\n", Rb);

	for (int i = 0; i < Nall; i++)
	{
		fprintf(restart_file, "%.32e\n", pos[i]);
	}
	for (int i = 0; i < Nall; i++)
	{
		fprintf(restart_file, "%.32e\n", Vel[i]);
	}
	for (int na = 0; na < Na; na++)
	{
		fprintf(restart_file, "%.32e\n", A0);
	}
	for (int edg = 0; edg < Nhalf; edg++)

	{
		int nvertices = Nhalf/Na;
		int type = (int)(edg + nvertices) / nvertices;
		fprintf(restart_file, "%.32e\n", L0[type-1]);
	}
	for (int na = 0; na < Na; na++)
	{
		fprintf(restart_file, "%.32e\n", D0[na]);
	}

	fprintf(restart_file, "%.32e\n", box_length);

	fclose(restart_file);
}

double FIRE_DPM_RP_VL(int Na, vector<int> &Ns, int Nhalf, int Nall,
					  double Rb, vector<double> &D0, double A0,
					  vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
					  double Kp, double KA, double KAA, double KAA_a, double Kw,
					  vector<int> &idx_start, vector<int> &idx_end,
					  vector<int> &ift, vector<int> &jft,
					  vector<double> &pos, double Fthresh, double dt_fire, int Nt_fire,
					  vector<double> &L0, double Kb,
					  int dstep, int bending_added, double box_length, int bending_type, double l0)
{
	int n, nt = 0;
	vector<double> Acc(Nall, 0.0);
	vector<double> Vel(Nall, 0.0);
	vector<double> Kl(Na, Kp);
	double P, Vel_norm, Acc_norm, Acc_max, Acc_abs;
	// FIRE parameters
	int N_delay = 20;
	int N_pn_max = 2000;
	double f_inc = 1.1;
	double f_dec = 0.5;
	double a_start = 0.15;
	double f_a = 0.99;
	double dt_max = 10.0 * dt_fire;
	double dt_min = 0.05 * dt_fire;
	int initialdelay = 1;
	double v_rsc;
	// FIRE Initialization
	double dt = dt_fire;
	double dt_half = dt / 2.0;
	double a_fire = a_start;
	double delta_a = 1.0 - a_fire;
	int N_pp = 0, N_pn = 0;
	int dstep_mini = 10000;
	int flag_write = 0;

	// Verlet list parameters
	int first_call = 1;
	vector<vector<int>> VL_aa(20 * Nhalf, vector<int>(4, -1));
	vector<int> VL_counter_all(1, 0);
	vector<double> pos_save(Nall, 0.0);
	vector<double> aforces(Nall, 0.0);
	vector<double> lk(Nhalf, 0.0);
	vector<double> A(Na, 0.0);
	char coords[200];
	sprintf(coords, "minimize.csv");
	FILE *coords_file;

	vector<double> xaa(Nhalf, 0.0);
	vector<double> yaa(Nhalf, 0.0);
	for (int nw = 0; nw < Nhalf; nw++)
	{
		xaa[nw] = pos[nw];
		yaa[nw] = pos[nw + Nhalf];
	}

	vector<double> lxx(Nhalf, 0.0);
	vector<double> lyy(Nhalf, 0.0);
	vector<double> lk_sqq(Nhalf, 0.0);

	for (int na = 0; na < Na; na++)

	{
		// L0_a = L0[na];
		for (int nw = idx_start[na]; nw <= idx_end[na]; nw++)
		{
			int iftss = ift[nw];
			double lxii = xaa[iftss] - xaa[nw];
			double lyii = yaa[iftss] - yaa[nw];
			lxx[nw] = lxii;
			lyy[nw] = lyii;
			lk_sqq[nw] = lxii * lxii + lyii * lyii;
		}
	}

	vector<double> thetaa(Nhalf, 0.0);
	vector<double> alpha(Nhalf, 0.0);

	for (int nw = 0; nw < Nhalf; nw++)
	{
		thetaa[nw] = std::atan2(lyy[nw], lxx[nw]);
	}
	for (int na = 0; na < Na; na++)
	{
		for (int nw = idx_start[na]; nw <= idx_end[na]; nw++)
		{
			alpha[nw] = thetaa[nw] - thetaa[jft[nw]];
		}
	}

	double alphaai;
	for (int nw = 0; nw < Nhalf; nw++)
	{
		alphaai = alpha[nw];
		alpha_0[nw] = alphaai - std::round(alphaai / 2.0 / PI) * 2.0 * PI;
	}

	double r_cut = D0[0];
	for (n = 1; n < Na; n++)
		r_cut = (r_cut > D0[n]) ? r_cut : D0[n];
	r_cut *= a_cut_c;

	double Pw;

	VerletList(Na, Ns, Nhalf, Nall, Da, r_cut, idx_start, idx_end,
			   pos, pos_save, VL_aa, VL_counter_all, first_call, box_length);
	Pw = Force_DPM(Na, Ns, Nhalf, Nall, Rb, D0, A0, Da, a_cut_c,
				   Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
				   pos, VL_aa, VL_counter_all, Acc, aforces, L0, lk, A, Kb, nt, dstep, bending_added,
				   box_length, bending_type, Kl);
	/*
	Eval = Energy_DPM(Na, Ns, Nhalf, Nall, Lx, Ly, L0, D0, A0, Da,
						Kp, KA, KC, Kb, KAC, KCC, idx_start, idx_end,
						ift, jft, pos_all, Nc, Dc, Acc);
	*/
	Acc_max = std::abs(Acc[0]);
	for (n = 1; n < Nall; n++)
	{
		Acc_abs = std::abs(Acc[n]);
		Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
	}
	if (Acc_max < Fthresh)
	{
		// printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
		return Pw;
	}

	first_call = 0;
	for (nt = 1; nt < Nt_fire; nt++)
	{
		P = 0.0;
		for (n = 0; n < Nall; n++)
			P += Acc[n] * Vel[n];
		if (P > 0.0)
		{
			N_pp += 1;
			N_pn = 0;
			if (N_pp > N_delay)
			{
				dt = std::min(f_inc * dt, dt_max);
				dt_half = dt / 2.0;
				a_fire *= f_a;
				delta_a = 1.0 - a_fire;
			}
		}
		else
		{
			N_pn += 1;
			N_pp = 0;
			if (N_pn > N_pn_max)
				break;
			if ((initialdelay < 1) || (nt >= N_delay))
			{
				if (f_dec * dt > dt_min)
				{
					dt *= f_dec;
					dt_half = dt / 2.0;
				}
				a_fire = a_start;
				delta_a = 1.0 - a_fire;
				for (n = 0; n < Nall; n++)
				{
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
		for (n = 0; n < Nall; n++)
		{
			Vel_norm += Vel[n] * Vel[n];
			Acc_norm += Acc[n] * Acc[n];
		}
		v_rsc = a_fire * std::sqrt(Vel_norm / Acc_norm);
		for (n = 0; n < Nall; n++)
		{
			Vel[n] = delta_a * Vel[n] + v_rsc * Acc[n];
			pos[n] += dt * Vel[n];
		}

		VerletList(Na, Ns, Nhalf, Nall, Da, r_cut, idx_start, idx_end,
				   pos, pos_save, VL_aa, VL_counter_all, first_call, box_length);

		Pw = Force_DPM(Na, Ns, Nhalf, Nall, Rb, D0, A0, Da, a_cut_c,
					   Kp, KA, KAA, KAA_a, Kw, idx_start, idx_end, ift, jft,
					   pos, VL_aa, VL_counter_all, Acc, aforces, L0, lk, A, Kb, nt, dstep, bending_added,
					   box_length, bending_type, Kl);
		/*
		Eval = Energy_DPM(Na, Ns, Nhalf, Nall, Lx, Ly, L0, D0, A0, Da,
						Kp, KA, KC, Kb, KAC, KCC, idx_start, idx_end,
						ift, jft, pos_all, Nc, Dc, Acc);
		*/
		for (n = 0; n < Nall; n++)
			Vel[n] += dt_half * Acc[n];

		Acc_max = std::abs(Acc[0]);
		for (n = 1; n < Nall; n++)
		{
			Acc_abs = std::abs(Acc[n]);
			Acc_max = (Acc_max < Acc_abs) ? Acc_abs : Acc_max;
		}
		// if (Acc_max < Fthresh)
		// 	break;
		// if (nt % 1 == 0) printf("nt: %d  Max Force: %.5e\n", nt, Acc_max);

		if (nt % dstep_mini == 0)
		{
			if (nt == 1)
			{
				coords_file = fopen(coords, "w");

				for (int i = 0; i < Nhalf; i++)
				{
					int type = (int)(i + Ns[0]) / Ns[0];
					fprintf(coords_file, "%d, %d, %d, %f, %f, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f\n", nt, i + 1, type, pos[i], pos[i + Nhalf], 0.,
							D0[0] / 2., Nhalf, A[type - 1], A0, lk[i], L0[type-1], alpha[i], alpha_0[i], Kb, box_length);
				}
			}

			else
			{
				coords_file = fopen(coords, "a");

				for (int i = 0; i < Nhalf; i++)
				{
					int type = (int)(i + Ns[0]) / Ns[0];
					fprintf(coords_file, "%d, %d, %d, %f, %f, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f\n", nt, i + 1, type, pos[i], pos[i + Nhalf], 0.,
							D0[0] / 2., Nhalf, A[type - 1], A0, lk[i], L0[type-1], alpha[i], alpha_0[i], Kb, box_length);
				}
			}
			fclose(coords_file);
		}
	}
	// printf("Acc_max: %.4e  Iteration Number: %d\n", Acc_max, nt);
	return Pw;
}
