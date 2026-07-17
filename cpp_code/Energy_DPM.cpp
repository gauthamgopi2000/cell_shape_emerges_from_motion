#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <string.h>
#include <vector>

#include "Energy_DPM.h"
#include "main.h"

using std::vector;

double Energy_DPM(int Na, vector<int> &Ns, int Nhalf, int Nall, double Rb,
				  vector<double> &L0, vector<double> &D0, vector<double> &A0,
				  vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				  double Kp, double KA, double KAA, double KAA_a, double Kw,
				  vector<int> &idx_start, vector<int> &idx_end,
				  vector<int> &ift, vector<int> &jft,
				  vector<double> &pos, vector<double> &force)
{
    int 		ns, ns_y, ms, na, ma;
    double  	lxi, lyi, lki, dli, dlj;
    double  	F;

    double 		Pwall 		= 0.0;

	vector<double>  x(Nhalf, 0.0);
	vector<double>  y(Nhalf, 0.0);
	vector<double>  lx(Nhalf, 0.0);
	vector<double>  ly(Nhalf, 0.0);
	vector<double>	lk_sq(Nhalf, 0.0);
	vector<double>  A(Na, 0.0);
	double  		dAa;
	vector<double>  dxA(Nhalf, 0.0);
	vector<double>  dyA(Nhalf, 0.0);
	vector<double>  lk(Nhalf, 0.0);

	int     	ifts, jfts;
	double  	d0;
	int     	countstart 	= 0;
	int     	countend 	= 0;

	for (ns = 0; ns < Nall; ns++)
        force[ns] = 0.0;

	for (na = 0; na < Na; na++){
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			x[ns] = pos[ns];
			y[ns] = pos[ns + Nhalf];
		}
	}

	// Perimeter Force
	double 		L0_a, L0_sq;
	double 		Ns_a;
	for (na = 0; na < Na; na++){
		L0_a = L0[na];
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			ifts = ift[ns];
			lxi = x[ifts] - x[ns];
			lyi = y[ifts] - y[ns];
			lx[ns] = lxi;
			ly[ns] = lyi;
			lk_sq[ns] = lxi * lxi + lyi * lyi;
			lk[ns] = std::sqrt(lk_sq[ns]);
			dli = lk[ns] / L0_a - 1.0;
			// Eval += Kp * dli * dli / double(Ns[nc]) / 2.0;
		}
	}
	for (na = 0; na < Na; na++){
		Ns_a = double(Ns[na]);
		L0_a = L0[na];
		L0_sq = L0_a * L0_a;
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			jfts = jft[ns];
	    	ns_y = ns + Nhalf;
	    	dlj = (L0_a / lk[jfts] - 1.0) / L0_sq;
	    	dli = (1.0 - L0_a / lk[ns]) / L0_sq;
	    	force[ns] += Kp / Ns_a * (dlj * lx[jfts] + dli * lx[ns]);
    		force[ns_y] += Kp / Ns_a * (dlj * ly[jfts] + dli * ly[ns]);
		}
	}
	// Area Force
	double A0_a, KA_a;
    for (na = 0; na < Na; na++){
    	A0_a = A0[na];
        KA_a = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * KA;
        // printf("Ka: %.5e\n", KA_a);
    	for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
    		ifts = ift[ns];
    		jfts = jft[ns];
    		A[na] += (x[ns] * y[ifts] - y[ns] * x[ifts]);
    		dxA[ns] = (x[ifts] - x[jfts]) / A0_a;
    		dyA[ns] = (y[jfts] - y[ifts]) / A0_a;
    	}
    	A[na] /= 2.0;
    	dAa = KA_a * (A[na] / A0_a - 1.0);
    	// Eval += 0.5 * dAc * dAc / KA;
    	for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
    		force[ns] += 0.5 * dAa * dyA[ns];
    		force[ns + Nhalf] += 0.5 * dAa * dxA[ns];
    	}
    }
    
	//Bending force
	// version 1: (\alpha - \alpha_0)^2
	/*
	vector<double>	 	theta(Nhalf, 0.0);
	vector<double>	 	d_alpha(Nhalf, 0.0);
	vector<double>	 	cos_theta(Nhalf, 0.0);
	vector<double> 		sin_theta(Nhalf, 0.0);
	double 				alphai, d_alpha_i, d_alpha_j, Kb_s;
	for (ns = 0; ns < Nhalf; ns++){
		theta[ns] = std::atan2(ly[ns], lx[ns]);
		cos_theta[ns] = lx[ns] / lk_sq[ns];
		sin_theta[ns] = ly[ns] / lk_sq[ns];
	}
	for (na = 0; na < Na; na++){
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			alphai = theta[ns] - theta[jft[ns]];
			d_alpha[ns] = alphai - std::round(alphai / 2.0 / PI) * 2.0 * PI - alpha_0[na];
			// Eval += 0.5 * Kb * d_alpha[ns] * d_alpha[ns];
		}
	}
    for (na = 0; na < Na; na++){
        for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
            Kb_s = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * Kb_all[ns];
            // printf("ns: %d  Kb: %.5e\n", ns, Kb_s);
            ifts = ift[ns];
            jfts = jft[ns];
            d_alpha_i = d_alpha[ns] - d_alpha[ifts];
            d_alpha_j = d_alpha[ns] - d_alpha[jfts];
            force[ns] -= Kb_s * (d_alpha_i * sin_theta[ns] + d_alpha_j * sin_theta[jfts]);
            force[ns + Nhalf] += Kb_s * (d_alpha_i * cos_theta[ns] + d_alpha_j * cos_theta[jfts]);
        }
	}
	*/
	
	// Cell-cell Force
    vector<double> 	xcen(Na, 0.0);
    vector<double> 	ycen(Na, 0.0);
    double 		dx_cen, dy_cen, xshift, yshift;
    double 		xcenn, ycenn, xcenm, ycenm;
    int 		i, j;
    double 		Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double 		dFx, dFy;
    double 		x1, y1, x3, y3;
    double 		dx13, dy13;
    double 		dnm, dd;
	double 		a_cut_h = 0.5 + 0.5 * a_cut_c;

    for (na = 0; na < Na; na++){
    	for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
    		xcen[na] += x[ns];
    		ycen[na] += y[ns];
    	}
    	xcen[na] /= Ns[na];
    	ycen[na] /= Ns[na];
    }
    for (na = 0; na < Na - 1; na++){
	    xcenn = xcen[na];
	    ycenn = ycen[na];
    	for (ma = na + 1; ma < Na; ma++){
    		xcenm = xcen[ma];
		    ycenm = ycen[ma];
    		dx_cen = xcenm - xcenn;
	        dy_cen = ycenm - ycenn;
	        if (std::sqrt(dx_cen * dx_cen + dy_cen * dy_cen) > (Da[na] + Da[ma]))
	            continue;

    		Dnm = (D0[na] + D0[ma]) / 2.0;
    		Dnm_sq = Dnm * Dnm;
			Dnm_a = a_cut_c * Dnm;
			Dnm_h = a_cut_h * Dnm;

		    for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
	            x1 = x[ns]; y1 = y[ns];
		    	for (ms = idx_start[ma]; ms <= idx_end[ma]; ms++){
	                dx13 = x[ms] - x1;
					if (std::abs(dx13) < Dnm_a){
	                	dy13 = y[ms] - y1;
						if (std::abs(dy13) < Dnm_a){
							dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
							if (dnm < Dnm_a){
								if (dnm > Dnm_h)
									F = KAA_a * (Dnm_a / dnm - 1.0) / Dnm_sq;
								else if ((dnm <= Dnm_h) && (dnm > Dnm))
									F = KAA_a * (1.0 - Dnm / dnm) / Dnm_sq;
								else
									F = KAA * (1.0 - Dnm / dnm) / Dnm_sq;
								dFx = F * dx13;
								dFy = F * dy13;
								force[ns]     += dFx;
								force[ns + Nhalf]   += dFy;
								force[ms]    -= dFx;
								force[ms + Nhalf]  -= dFy;
							}
						}
					}
		    	}
		    }
    	}
    }

	// intra-cell interaction
	double 		Dnm2;
	for (na = 0; na < Na; na++){
	    Dnm = L0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_start[na]; ns < idx_end[na] - 1; ns++){
	    	ns_y = ns + Nhalf;
	        for (ms = ns + 2; ms <= idx_end[na]; ms++){
	        	if ((ns == idx_start[na]) && (ms == idx_end[na]))
		    		continue;
	            dx13 = x[ms] - x[ns];
	            if (std::abs(dx13) < Dnm){
	                dy13 = y[ms] - y[ns];
	                dnm = dx13 * dx13 + dy13 * dy13;
	                if(dnm < Dnm_sq){
	                    dnm = std::sqrt(dnm);
	                    F = KAA * (1.0 - Dnm / dnm) / Dnm_sq;
						// dd = 1.0 - dnm / Dnm;
	                    // Eval += KC_half * dd * dd;  // cell-cell PE

	                    dFx = F * dx13;
                        dFy = F * dy13;
                        force[ns] 			+= dFx;
                        force[ns_y] 		+= dFy;
                        force[ms] 			-= dFx; // 3rd law
                        force[ms + Nhalf] 	-= dFy;
	                }
	            }
	        }
	    }
	}

	// cell-wall force
	double 	r;
	for (na = 0; na < Na; na++){
	    Dnm = 0.5 * D0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			r = std::sqrt(x[ns] * x[ns] + y[ns] * y[ns]);
	    	if (r > Rb - Dnm){
				F = Kw * (Dnm + r - Rb) / Dnm_sq / r;
				force[ns] -= F * x[ns];
				force[ns + Nhalf] -= F * y[ns];
				Pwall += F;
			}
	    }
	}
	
	return Pwall / 2.0 / Rb / PI;
}



double Force_DPM(int Na, vector<int> &Ns, int Nhalf, int Nall,
				double Rb, vector<double> &L0, vector<double> &D0, vector<double> &A0,
				vector<double> &Da, vector<double> &alpha_0, double a_cut_c,
				double Kp, double KA, double KAA, double KAA_a, double Kw,
				vector<int> &idx_start, vector<int> &idx_end,
				vector<int> &ift, vector<int> &jft,
				vector<double> &pos, vector<vector <int>> &VL_aa,
				vector<int> &VL_counter, vector<double> &F_a)
{
    int 		ns, ns_y, ms, na, ma;
    double  	lxi, lyi, lki, dli, dlj;
    double  	F;
	double 		Pwall 		= 0.0;

	vector<double> 	xa(Nhalf, 0.0);
	vector<double> 	ya(Nhalf, 0.0);
	vector<double> 	lx(Nhalf, 0.0);
	vector<double> 	ly(Nhalf, 0.0);
	vector<double> 	lk_sq(Nhalf, 0.0);
	vector<double> 	A(Na, 0.0);
	vector<double> 	dxA(Nhalf, 0.0);
	vector<double> 	dyA(Nhalf, 0.0);
	vector<double> 	lk(Nhalf, 0.0);
	double  	dAa;

	int     	ifts, jfts;
	double  	d0;
	int     	countstart 	= 0;
	int     	countend 	= 0;

	for (ns = 0; ns < Nhalf; ns++){
		xa[ns] = pos[ns];
		ya[ns] = pos[ns + Nhalf];
		F_a[ns] = 0.0;
		F_a[ns + Nhalf] = 0.0;
	}

	// Perimeter Force
	double 		L0_a, L0_sq;
	double 		Ns_a;
	for (na = 0; na < Na; na++){
		L0_a = L0[na];
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
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
	for (na = 0; na < Na; na++){
		Ns_a = double(Ns[na]);
		L0_a = L0[na];
		L0_sq = L0_a * L0_a;
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			jfts = jft[ns];
	    	dlj = (L0_a / lk[jfts] - 1.0) / L0_sq;
	    	dli = (1.0 - L0_a / lk[ns]) / L0_sq;
	    	F_a[ns] += Kp / Ns_a * (dlj * lx[jfts] + dli * lx[ns]);
    		F_a[ns + Nhalf] += Kp / Ns_a * (dlj * ly[jfts] + dli * ly[ns]);
		}
	}
	// Area Force
	double A0_a, KA_a;
    for (na = 0; na < Na; na++){
    	A0_a = A0[na];
        KA_a = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * KA;
        // printf("na: %d  Ka: %.5e\n", na, KA_a);
    	for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
    		ifts = ift[ns];
    		jfts = jft[ns];
    		A[na] += (xa[ns] * ya[ifts] - ya[ns] * xa[ifts]);
    		dxA[ns] = (xa[ifts] - xa[jfts]) / A0_a;
    		dyA[ns] = (ya[jfts] - ya[ifts]) / A0_a;
    	}
    	A[na] /= 2.0;
    	dAa = 0.5 * KA_a * (A[na] / A0_a - 1.0);
    	// Eval += 0.5 * dAc * dAc / KA;
    	for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
    		F_a[ns] += dAa * dyA[ns];
    		F_a[ns + Nhalf] += dAa * dxA[ns];
    	}
    }
    
	//Bending force
	// version 1: (\alpha - \alpha_0)^2
	/*
	vector<double>	 	theta(Nhalf, 0.0);
	vector<double>	 	d_alpha(Nhalf, 0.0);
	vector<double>	 	cos_theta(Nhalf, 0.0);
	vector<double> 		sin_theta(Nhalf, 0.0);
	double 				alphai, d_alpha_i, d_alpha_j, Kb_s;
	for (ns = 0; ns < Nhalf; ns++){
		theta[ns] = std::atan2(ly[ns], lx[ns]);
		cos_theta[ns] = lx[ns] / lk_sq[ns];
		sin_theta[ns] = ly[ns] / lk_sq[ns];
	}
	for (na = 0; na < Na; na++){
		for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			alphai = theta[ns] - theta[jft[ns]];
			d_alpha[ns] = alphai - std::round(alphai / 2.0 / PI) * 2.0 * PI - alpha_0[na];
			// Eval += 0.5 * Kb * d_alpha[ns] * d_alpha[ns];
		}
	}
    for (na = 0; na < Na; na++){
        for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
            Kb_s = ((double)Ns[na] / 20.0) * ((double)Ns[na] / 20.0) * Kb_all[ns];
            // printf("ns: %d  Kb: %.5e\n", ns, Kb_s);
            ifts = ift[ns];
            jfts = jft[ns];
            d_alpha_i = d_alpha[ns] - d_alpha[ifts];
            d_alpha_j = d_alpha[ns] - d_alpha[jfts];
            F_a[ns] -= Kb_s * (d_alpha_i * sin_theta[ns] + d_alpha_j * sin_theta[jfts]);
            F_a[ns + Nhalf] += Kb_s * (d_alpha_i * cos_theta[ns] + d_alpha_j * cos_theta[jfts]);
        }
	}
	*/
	
	// Cell-cell Force
    double 		Dnm, Dnm_sq, Dnm_a, Dnm_h;
    double 		dFx, dFy;
    double 		dx13, dy13;
    double 		dnm, dd;
	int 		vl_idx;
	double 		a_cut_h = 0.5 + 0.5 * a_cut_c;

	for (vl_idx = 0; vl_idx < VL_counter[0]; vl_idx++){
		ns = VL_aa[vl_idx][0];
		ms = VL_aa[vl_idx][1];
		na = VL_aa[vl_idx][2];
		ma = VL_aa[vl_idx][3];
		Dnm = 0.5 * (D0[na] + D0[ma]);
		Dnm_sq = Dnm * Dnm;
		Dnm_a = a_cut_c * Dnm;
		Dnm_h = a_cut_h * Dnm;
		
		dx13 = xa[ms] - xa[ns];
		if (std::abs(dx13) < Dnm_a){
			dy13 = ya[ms] - ya[ns];
			if (std::abs(dy13) < Dnm_a){
				dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
				if (dnm < Dnm_a){
					if (dnm > Dnm_h)
						F = KAA_a * (Dnm_a / dnm - 1.0) / Dnm_sq;
					else if ((dnm <= Dnm_h) && (dnm > Dnm))
						F = KAA_a * (1.0 - Dnm / dnm) / Dnm_sq;
					else
						F = KAA * (1.0 - Dnm / dnm) / Dnm_sq;
					dFx = F * dx13;
					dFy = F * dy13;
					F_a[ns]     	+= dFx;
					F_a[ns + Nhalf] += dFy;
					F_a[ms]    		-= dFx;
					F_a[ms + Nhalf] -= dFy;
				}
			}
		}
	}

	// intra-cell interaction
	double 		Dnm2;
	for (na = 0; na < Na; na++){
	    Dnm = L0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_start[na]; ns < idx_end[na] - 1; ns++){
	        for (ms = ns + 2; ms <= idx_end[na]; ms++){
	        	if ((ns == idx_start[na]) && (ms == idx_end[na])){
		    		continue;
		    	}
	            dx13 = xa[ms] - xa[ns];
	            if (std::abs(dx13) < Dnm){
	                dy13 = ya[ms] - ya[ns];
					if (std::abs(dy13) < Dnm){
						dnm = std::sqrt(dx13 * dx13 + dy13 * dy13);
						if(dnm < Dnm){
							F = KAA * (1.0 - Dnm / dnm) / Dnm_sq;
							// dd = 1.0 - dnm / Dnm;
							// Eval += KC_half * dd * dd;  // cell-cell PE

							dFx = F * dx13;
							dFy = F * dy13;
							F_a[ns] 			+= dFx;
							F_a[ns + Nhalf] 	+= dFy;
							F_a[ms] 			-= dFx; // 3rd law
							F_a[ms + Nhalf] 	-= dFy;
						}
					}
	            }
	        }
	    }
	}
	
	// adipocyte-wall force
	double 	r;
	for (na = 0; na < Na; na++){
	    Dnm = 0.5 * D0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (ns = idx_start[na]; ns <= idx_end[na]; ns++){
			r = std::sqrt(xa[ns] * xa[ns] + ya[ns] * ya[ns]);
	    	if (r > Rb - Dnm){
				F = Kw * (Dnm + r - Rb) / Dnm_sq / r;
				F_a[ns] -= F * xa[ns];
				F_a[ns + Nhalf] -= F * ya[ns];
				Pwall += F;
			}
	    }
	}
	
	return Pwall / 2.0 / Rb / PI;
}


double Pressure_DPM(int Na, int Nhalf, double Rb, vector<int> &idx_start, vector<int> &idx_end,
					vector<double> &D0, double Kw, vector<double> &pos)
{
	double 		Pwall 		= 0.0;
	double 		Dnm, Dnm_sq, F;

	// adipocyte-wall force
	double 	r;
	for (int na = 0; na < Na; na++){
	    Dnm = 0.5 * D0[na];
	    Dnm_sq = Dnm * Dnm;
	    for (int ns = idx_start[na]; ns <= idx_end[na]; ns++){
			r = std::sqrt(pos[ns] * pos[ns] + pos[ns + Nhalf] * pos[ns + Nhalf]);
	    	if (r > Rb - Dnm){
				F = Kw * (Dnm + r - Rb) / Dnm_sq / r;
				Pwall += F;
			}
	    }
	}
	
	return Pwall / 2.0 / Rb / PI;
}
