#ifndef UTILS_H
#define UTILS_H

// Miscellaneous utility functions here -- initialization, timing, I/O

void make_sim_dir(void);

double get_a_from_snap(int snap);
double get_t_from_snap(int snap);
double get_a(double t);
double get_t(double a);
int set_snap(double a);
double get_dt(double a);

void ioSnap(int mode, int snap);

#endif
