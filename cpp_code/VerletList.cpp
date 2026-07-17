#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <vector>

using std::vector;

void VerletList(int Na, vector<int> &Ns, int Nhalf, int Nall,
				vector<double> &Da, double r_cut,
				vector<int> &idx_start, vector<int> &idx_end,
				vector<double> &pos_a, vector<double> &pos_save,
				vector<vector <int>> &VL_aa, vector<int> &VL_counter,
				int first_call)
{
	// distance for making Verlet list
	double 		r_factor = 1.2;
	double 		r_cut_sq = r_cut * r_cut;
	double 		r_list = r_factor * r_cut;
	double 		r_list_sq = r_list * r_list;
	double 		r_skin_sq = (r_factor - 1.0) * (r_factor - 1.0) * r_cut_sq;

	double 		dx, dy, dr;
	int 		i, iy, j, jy;
	int 		na, ma, ns, ms, nc, mc;

	if (first_call == 0){
		double dr_max = 0.0;
		for (i = 0; i < Nhalf; i++){
			iy = i + Nhalf;
			dx = pos_a[i] - pos_save[i];
			dy = pos_a[iy] - pos_save[iy];
			dr = dx * dx + dy * dy;
			dr_max = (dr > dr_max) ? dr : dr_max;
		}
	    if (4.0 * dr_max < r_skin_sq){
	        return;
	    }
	}

	vector<double>	x_cen(Na, 0.0);
	vector<double> 	y_cen(Na, 0.0);
	for (na = 0; na < Na; na++){
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
		    x_cen[na] += pos_a[ns];
		    y_cen[na] += pos_a[ns + Nhalf];
		}
		x_cen[na] /= Ns[na];
		y_cen[na] /= Ns[na];
	}

	double 		x_cen_m, y_cen_m, x_cen_n, y_cen_n;
	double 		dx_cen, dy_cen;
	double 		dx_cen_shift, dy_cen_shift;
	double 		x1, y1, x3, y3;
	double 		lc1, lc2, dl_n, dl_m;

	VL_counter[0] = 0;

	for (na = 0; na < Na - 1; na++){
	    x_cen_n = x_cen[na];
	    y_cen_n = y_cen[na];
	    for (ma = na + 1; ma < Na; ma++){
	        x_cen_m = x_cen[ma];
	        y_cen_m = y_cen[ma];
	        dx_cen = x_cen_m - x_cen_n;
	        dy_cen = y_cen_m - y_cen_n;
	        if (std::sqrt(dx_cen * dx_cen + dy_cen * dy_cen) > Da[na] + Da[ma] + r_list)
	            continue;

		    for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
	            x1 = pos_a[ns]; y1 = pos_a[ns + Nhalf];
		    	for (ms = idx_start[ma]; ms <= idx_end[ma]; ms++){
	                dx = pos_a[ms] - x1;
	                dy = pos_a[ms + Nhalf] - y1;
					dr = dx * dx + dy * dy;
	                if (dr < r_list_sq){
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