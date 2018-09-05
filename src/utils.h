#ifndef UTILS_H
#define UTILS_H

// Miscellaneous utility functions here -- initialization, timing, I/O

void make_sim_dir(void);

double get_a_from_snap(int snap);
double get_t_from_snap(int snap);
double get_a(double t);
double get_t(double a);
double get_dt(double a);

void readIC(void);
void saveSnap(int snap);
int set_snap(double a);


#endif
