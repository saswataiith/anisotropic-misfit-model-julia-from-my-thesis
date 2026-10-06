#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include"gsl_support.h"
#define PI M_PI
#define Tolerance 1.0e-10
#define COMPERR 1.0e-8
#define Re 0
#define Im 1 


unsigned fftw_flag;

long SEED;

fftw_complex *comp, *dfdc;

fftw_complex *eta_1, *eta_2, *eta_3;

fftw_complex *dfdeta_1, *dfdeta_2, *dfdeta_3;

fftw_plan p_up, p_dn;

int **occupancy;

int num_steps, count, ncount;

int initcount, flag;

int noise_flag;

double alloycomp, noise_level; 

double sustained_noise_level_cons, sustained_noise_level_uncons_1, 
			 sustained_noise_level_uncons_2, sustained_noise_level_uncons_3;

double del_x, del_y, del_t;

double del_t1, del_t2;

double sim_time, total_time;

double err;

double sigma;

int nx, ny, nx_half, ny_half;

int n_cout, noise_steps;

double A1, A2, A41, A42, A61, A62;
double mobility, kappa;

double L1, L2, L3;
 
double kappa1, kappa2, kappa3;

double one_by_nxny;

double supersaturation;

int nucleation_state, iter;

FILE *fpout;

double epsc, eps_eta, tetra;

double mubar, nubar, Aniso;

int steps_t1, steps_t2, time_to_change;

double eigen_strain_0[3][3], eigen_strain_1[3][3], eigen_strain_2[3][3], 
			 eigen_strain_3[3][3];

double **B_00, **B_10, **B_20, **B_30, **B_11, **B_21, **B_31, **B_22, **B_32,
			 **B_33;
