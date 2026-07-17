#ifndef MAIN_H_
#define MAIN_H_

const double PI = 3.14159265358979323846264338328;

void Packing(int , double , double *, double );
void restart_write(char file_name[200], double Rb, vector<double> &pos, vector<double> &Vel, double A0, vector<double> &L0i,
double KA, double Kp, int Na, int Nall, int Nhalf, double Kb, vector<double> &D0, double box_length);
double pbc_separation(double dl, double box_length);

#endif /*MAIN_H_*/
