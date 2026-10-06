#include"binary.h"
void
Evolve ()
{
 
void Output_Conf ();

double gsl_gaussian_sample(long *idum);

double gsl_uniform_sample(long *idum);

fftw_complex *var_elast0, *var_elast1, *var_elast2, *var_elast3;

fftw_complex *eta1sq, *eta2sq, *eta3sq;

double **random_num, **random_num1, **random_num2, **random_num3;

int i, j, loop_condition;

double del_kx, del_ky, kx, ky, kpow2, kpow4;

double rc, fp, rc_new;

double reta_1, reta_2, reta_3, eta_new_1, eta_new_2,
		eta_new_3;

double fp_reta_1, fp_reta_2, fp_reta_3;
 
double total;
 
double **tempreal;
 
double lhs, rhs;

double lhs1, lhs2, lhs3;

double rhs1, rhs2, rhs3;
	
double c, eta[3];

double error, maxerror;

double sum, mean;
 
tempreal = dmatrix (0, nx, 0, ny);
random_num = dmatrix (0, nx, 0, ny);
random_num1 = dmatrix (0, nx, 0, ny);
random_num2 = dmatrix (0, nx, 0, ny);
random_num3 = dmatrix (0, nx, 0, ny);



var_elast0 = fftw_malloc(nx * ny * sizeof(fftw_complex));
var_elast1 = fftw_malloc(nx * ny * sizeof(fftw_complex));
var_elast2 = fftw_malloc(nx * ny * sizeof(fftw_complex));
var_elast3 = fftw_malloc(nx * ny * sizeof(fftw_complex));


eta1sq = 	fftw_malloc(nx * ny * sizeof(fftw_complex));		
eta2sq = 	fftw_malloc(nx * ny * sizeof(fftw_complex));		
eta3sq = 	fftw_malloc(nx * ny * sizeof(fftw_complex));		


del_kx = 2.0 * PI / ((double) nx * del_x);
del_ky = 2.0 * PI / ((double) ny * del_y);

 
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   tempreal[i][j] = comp[j + i * ny][Re];
  }
 }
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
  
					var_elast0[j + i * ny][Re] = 0.0;
					var_elast1[j + i * ny][Re] = 0.0;
					var_elast2[j + i * ny][Re] = 0.0;
					var_elast3[j + i * ny][Re] = 0.0;
					var_elast0[j + i * ny][Im] = 0.0;
					var_elast1[j + i * ny][Im] = 0.0;
					var_elast2[j + i * ny][Im] = 0.0;
					var_elast3[j + i * ny][Im] = 0.0;

					eta1sq[j + i * ny][Re] = 0.0;
					eta2sq[j + i * ny][Re] = 0.0;
					eta3sq[j + i * ny][Re] = 0.0;
					eta1sq[j + i * ny][Im] = 0.0;
					eta2sq[j + i * ny][Im] = 0.0;
					eta3sq[j + i * ny][Im] = 0.0;
  }}

 
 loop_condition = 1;
 iter = 0;
 del_t = del_t1;
 for (count = 0; count <= num_steps; count++) {
     ncount = count + initcount;
 
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
		 comp[j + i * ny][Re] = dfdc[ j + i * ny][Re];
		 comp[j + i * ny][Im] = dfdc[ j + i * ny][Im];
		 
		 eta_1[j + i * ny][Re] = dfdeta_1[ j + i * ny][Re];
		 eta_1[j + i * ny][Im] = dfdeta_1[ j + i * ny][Im];

		 eta_2[j + i * ny][Re] = dfdeta_2[ j + i * ny][Re];
		 eta_2[j + i * ny][Im] = dfdeta_2[ j + i * ny][Im];
	
		 eta_3[j + i * ny][Re] = dfdeta_3[ j + i * ny][Re];
		 eta_3[j + i * ny][Im] = dfdeta_3[ j + i * ny][Im];
     }
	}

/************************************************************************************* 
  if ((flag == 0) && (ncount < noise_steps)){
  Nucleate ();
  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
		 
		 dfdc[j + i * ny][Re] = comp[ j + i * ny][Re];
		 dfdc[j + i * ny][Im] = comp[ j + i * ny][Im];
		 
		 dfdeta_1[j + i * ny][Re] = eta_1[ j + i * ny][Re];
		 dfdeta_1[j + i * ny][Im] = eta_1[ j + i * ny][Im];

		 dfdeta_2[j + i * ny][Re] = eta_2[ j + i * ny][Re];
		 dfdeta_2[j + i * ny][Im] = eta_2[ j + i * ny][Im];
	
		 dfdeta_3[j + i * ny][Re] = eta_3[ j + i * ny][Re];
		 dfdeta_3[j + i * ny][Im] = eta_3[ j + i * ny][Im];
     }
	 }	 
	}		
**************************************************************************************/

fftw_execute_dft (p_up, comp, comp);
fftw_execute_dft (p_up, eta_1, eta_1);
fftw_execute_dft (p_up, eta_2, eta_2);
fftw_execute_dft (p_up, eta_3, eta_3);

if (ncount > noise_steps){
for (i = 0; i < nx; ++i) {
 for (j = 0; j < ny; ++j) {
 eta1sq[j + i * ny][Re] = dfdeta_1[j + i * ny][Re] * dfdeta_1[j + i * ny][Re];
 eta2sq[j + i * ny][Re] = dfdeta_2[j + i * ny][Re] * dfdeta_2[j + i * ny][Re];
 eta3sq[j + i * ny][Re] = dfdeta_3[j + i * ny][Re] * dfdeta_3[j + i * ny][Re];
 eta1sq[j + i * ny][Im] = 0.0;
 eta2sq[j + i * ny][Im] = 0.0;
 eta3sq[j + i * ny][Im] = 0.0;
 }}

fftw_execute_dft (p_up, eta1sq, eta1sq);
fftw_execute_dft (p_up, eta2sq, eta2sq);
fftw_execute_dft (p_up, eta3sq, eta3sq);

for (i = 0; i < nx; ++i) {
 for (j = 0; j < ny; ++j) {

	var_elast0[j + i * ny][Re] = 
	B_00[i][j] * comp[j + i * ny][Re] + B_10[i][j] * eta1sq[j + i * ny][Re] 
	+ B_20[i][j] * eta2sq[j + i * ny][Re] + B_30[i][j] * eta3sq[j + i * ny][Re];

	var_elast1[j + i * ny][Re] =
	B_10[i][j] * comp[j + i * ny][Re] + B_11[i][j] * eta1sq[j + i * ny][Re]
	+ B_21[i][j] * eta2sq[j + i * ny][Re] + B_31[i][j] * eta3sq[j + i * ny][Re];

	var_elast2[j + i * ny][Re] = 
	B_20[i][j] * comp[j + i * ny][Re] + B_21[i][j] * eta1sq[j + i * ny][Re]
	+ B_22[i][j] * eta2sq[j + i * ny][Re] + B_32[i][j] * eta3sq[j + i * ny][Re];

	var_elast3[j + i * ny][Re] =
	B_30[i][j] * comp[j + i * ny][Re] + B_31[i][j] * eta1sq[j + i * ny][Re] 
	+ B_32[i][j] * eta2sq[j + i * ny][Re] + B_33[i][j] * eta3sq[j + i * ny][Re];

	var_elast0[j + i * ny][Im] = 
	B_00[i][j] * comp[j + i * ny][Im] + B_10[i][j] * eta1sq[j + i * ny][Im] 
	+ B_20[i][j] * eta2sq[j + i * ny][Im] + B_30[i][j] * eta3sq[j + i * ny][Im];

	var_elast1[j + i * ny][Im] =
	B_10[i][j] * comp[j + i * ny][Im] + B_11[i][j] * eta1sq[j + i * ny][Im]
	+ B_21[i][j] * eta2sq[j + i * ny][Im] + B_31[i][j] * eta3sq[j + i * ny][Im];

	var_elast2[j + i * ny][Im] = 
	B_20[i][j] * comp[j + i * ny][Im] + B_21[i][j] * eta1sq[j + i * ny][Im]
	+ B_22[i][j] * eta2sq[j + i * ny][Im] + B_32[i][j] * eta3sq[j + i * ny][Im];

	var_elast3[j + i * ny][Im] =
	B_30[i][j] * comp[j + i * ny][Im] + B_31[i][j] * eta1sq[j + i * ny][Im] 
	+ B_32[i][j] * eta2sq[j + i * ny][Im] + B_33[i][j] * eta3sq[j + i * ny][Im];

 }}
 
fftw_execute_dft(p_dn, var_elast1, var_elast1);
fftw_execute_dft(p_dn, var_elast2, var_elast2);
fftw_execute_dft(p_dn, var_elast3, var_elast3);


for (i = 0; i < nx; ++i) {
 for (j = 0; j < ny; ++j) {
	var_elast1[j + i * ny][Re] = var_elast1[j + i * ny][Re] * one_by_nxny;
	var_elast2[j + i * ny][Re] = var_elast2[j + i * ny][Re] * one_by_nxny;
	var_elast3[j + i * ny][Re] = var_elast3[j + i * ny][Re] * one_by_nxny;
	var_elast1[j + i * ny][Im] = var_elast1[j + i * ny][Im] * one_by_nxny;
	var_elast2[j + i * ny][Im] = var_elast2[j + i * ny][Im] * one_by_nxny;
	var_elast3[j + i * ny][Im] = var_elast3[j + i * ny][Im] * one_by_nxny;
 }}


for (i = 0; i < nx; ++i) {
 for (j = 0; j < ny; ++j) {
	
 var_elast1[j + i * ny][Re] = 
				 2.0 * dfdeta_1[j + i * ny][Re] * var_elast1[j + i * ny][Re];	

 var_elast2[j + i * ny][Re] = 
				 2.0 * dfdeta_2[j + i * ny][Re] * var_elast2[j + i * ny][Re];		

 var_elast3[j + i * ny][Re] = 
				 2.0 * dfdeta_3[j + i * ny][Re] * var_elast3[j + i * ny][Re];	

 var_elast1[j + i * ny][Im] = 0.0;
 var_elast2[j + i * ny][Im] = 0.0;
 var_elast3[j + i * ny][Im] = 0.0;
 }}


fftw_execute_dft(p_up, var_elast1, var_elast1);
fftw_execute_dft(p_up, var_elast2, var_elast2);
fftw_execute_dft(p_up, var_elast3, var_elast3);
}


	if (count > noise_steps) del_t = del_t2;
  
	if (((ncount < noise_steps) && (ncount % 100)==0)
	    || (( ncount <= 5000) && (( ncount % 1000)==0)) 
			|| (( ncount > 5000 && ncount <= 50000 ) && 
      (( ncount % 2000 ) == 0)) ||
      (( ncount > 50000 ) && (( ncount % 10000 ) == 0)) ||
      (count == num_steps) || (loop_condition == 0)) {
   
	 printf ("total_time=%lf\n", sim_time);
   printf ("SUCCESS!!!!\n");
   Output_Conf ();
  
	}
//  printf ("num_steps=%d\n", count);
  
   if (count > num_steps || loop_condition == 0)
   break;


for (i = 0; i < nx; ++i) {
 for (j = 0; j < ny; ++j) {
				
			c = dfdc[j + i * ny][Re];
			
			eta[0] = dfdeta_1[j + i * ny][Re];

			eta[1] = dfdeta_2[j + i * ny][Re];

			eta[2] = dfdeta_3[j + i * ny][Re];

dfdc[j + i * ny][Re] =  
				c - eta[0] * eta[0] - eta[1] * eta[1] - eta[2] * eta[2];
				              

dfdeta_1[j + i * ny][Re] =  
       2.0*(1.0-c)*eta[0]-6.0*eta[0]*eta[0]*eta[0]+6.0*pow((eta[0]*
eta[0]+eta[1]*eta[1]+eta[2]*eta[2]),2.0)*eta[0];
      
dfdeta_2[j + i * ny][Re] = 
       2.0*(1.0-c)*eta[1]-6.0*eta[1]*eta[1]*eta[1]+6.0*pow((eta[0]*
			 eta[0] + eta[1]*eta[1] + eta[2]*eta[2]),2.0)*eta[1];
				
dfdeta_3[j + i * ny][Re] = 
       2.0*(1.0-c)*eta[2]-6.0*eta[2]*eta[2]*eta[2]+6.0*pow((eta[0]*
			 eta[0] + eta[1]*eta[1] + eta[2]*eta[2]),2.0)*eta[2];

dfdc[j + i * ny][Im] = 0.0;
dfdeta_1[j + i * ny][Im] = 0.0;
dfdeta_2[j + i * ny][Im] = 0.0;
dfdeta_3[j + i * ny][Im] = 0.0;
}}

	 
fftw_execute_dft (p_up, dfdc, dfdc);
fftw_execute_dft (p_up, dfdeta_1, dfdeta_1);
fftw_execute_dft (p_up, dfdeta_2, dfdeta_2);
fftw_execute_dft (p_up, dfdeta_3, dfdeta_3);

for (i = 0; i < nx; i++) {
 if (i <= nx_half)
  kx = (double) i *del_kx;
  else
  kx = (double) (i - nx) * del_kx;
   
  kx = kx * kx;
   
for (j = 0; j < ny; j++) {
 if (j <= ny_half)
  ky = (double) j *del_ky;
  else
  ky = (double) (j - ny) * del_ky;
    
  ky = ky * ky;
  kpow2 = kx + ky;
  kpow4 = kpow2 * kpow2;

  lhs = 1.0 + 2.0 * mobility * kappa * kpow4 * del_t;
    
  rc = comp[j + i * ny][Re];
  fp = dfdc[j + i * ny][Re];
  rhs = rc - mobility * kpow2 * del_t * fp
			  - mobility * kpow2 * del_t * var_elast0[j + i * ny][Re];
  rc_new = rhs / lhs;
  comp[j + i * ny][Re] = rc_new;
  dfdc[j + i * ny][Re] = comp[j + i * ny][Re];

  rc = comp[j + i * ny][Im];
  fp = dfdc[j + i * ny][Im];
  rhs = rc - mobility * kpow2 * del_t * fp
			  - mobility * kpow2 * del_t * var_elast0[j + i * ny][Im];
  rc_new = rhs / lhs;
  comp[j + i * ny][Im] = rc_new;
  dfdc[j + i * ny][Im] = comp[j + i * ny][Im];
 
/* Real part */		
	reta_1 = eta_1[j + i * ny][Re];
	fp_reta_1 = dfdeta_1[j + i * ny][Re];
	lhs1 = 1.0 + 2.0 * L1 * kappa1 * kpow2 * del_t;
	rhs1 = reta_1 - L1 * del_t * fp_reta_1
				 - L1 * del_t * var_elast1[j + i * ny][Re];
	eta_new_1 = rhs1/lhs1;
	eta_1[j + i * ny][Re] = eta_new_1;
	dfdeta_1[j + i * ny][Re] = eta_1[j + i * ny][Re];

  reta_2 = eta_2[j + i * ny][Re];
	fp_reta_2 = dfdeta_2[j + i * ny][Re];
	lhs2 = 1.0 + 2.0 * L2 * kappa2 * kpow2 * del_t;
	rhs2 = reta_2 - L2 * del_t * fp_reta_2
				 - L2 * del_t * var_elast2[j + i * ny][Re];
	eta_new_2 = rhs2/lhs2;
	eta_2[j + i * ny][Re] = eta_new_2;
	dfdeta_2[j + i * ny][Re] = eta_2[j + i * ny][Re];
	
  reta_3 = eta_3[j + i * ny][Re];
	fp_reta_3 = dfdeta_3[j + i * ny][Re];
	lhs3 = 1.0 + 2.0 * L3 * kappa3 * kpow2 * del_t;
	rhs3 = reta_3 - L3 * del_t * fp_reta_3
				 - L3 * del_t * var_elast3[j + i * ny][Re];
	eta_new_3 = rhs3/lhs3;
	eta_3[j + i * ny][Re] = eta_new_3;
	dfdeta_3[j + i * ny][Re] = eta_3[j + i * ny][Re];

/* Imaginary part */   
	
	reta_1 = eta_1[j + i * ny][Im];
	fp_reta_1 = dfdeta_1[j + i * ny][Im];
	lhs1 = 1.0 + 2.0 * L1 * kappa1 * kpow2 * del_t;
	rhs1 = reta_1 - L1 * del_t * fp_reta_1
				 - L1 * del_t * var_elast1[j + i * ny][Im];
	eta_new_1 = rhs1/lhs1;
	eta_1[j + i * ny][Im] = eta_new_1;
	dfdeta_1[j + i * ny][Im] = eta_1[j + i * ny][Im];

  reta_2 = eta_2[j + i * ny][Im];
	fp_reta_2 = dfdeta_2[j + i * ny][Im];
	lhs2 = 1.0 + 2.0 * L2 * kappa2 * kpow2 * del_t;
	rhs2 = reta_2 - L2 * del_t * fp_reta_2
				 - L2 * del_t * var_elast2[j + i * ny][Im];
	eta_new_2 = rhs2/lhs2;
	eta_2[j + i * ny][Im] = eta_new_2;
	dfdeta_2[j + i * ny][Im] = eta_2[j + i * ny][Im];
	
  reta_3 = eta_3[j + i * ny][Im];
	fp_reta_3 = dfdeta_3[j + i * ny][Im];
	lhs3 = 1.0 + 2.0 * L3 * kappa3 * kpow2 * del_t;
	rhs3 = reta_3 - L3 * del_t * fp_reta_3
				 - L3 * del_t * var_elast3[j + i * ny][Im];
	eta_new_3 = rhs3/lhs3;
	eta_3[j + i * ny][Im] = eta_new_3;
	dfdeta_3[j + i * ny][Im] = eta_3[j + i * ny][Im];
	
   }
  }

/* Check for conservation of mass */
  total = dfdc[0][Re] * one_by_nxny;
//  printf ("total_A=%le total_B=%le\n", total_A, total_B);
  err = fabs (total - alloycomp);
  if (err > COMPERR) {
   printf ("ELEMENTS ARE NOT CONSERVED,SORRY!!!!\n");
   printf ("error=%lf\n", err);
   exit (0);
  }
  

  fftw_execute_dft (p_dn, dfdc, dfdc);
  fftw_execute_dft (p_dn, dfdeta_1, dfdeta_1);
  fftw_execute_dft (p_dn, dfdeta_2, dfdeta_2);
  fftw_execute_dft (p_dn, dfdeta_3, dfdeta_3);

  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    dfdc[j + i * ny][Re] = dfdc[j + i * ny][Re] * one_by_nxny;
    dfdc[j + i * ny][Im] = dfdc[j + i * ny][Im] * one_by_nxny;
	  dfdeta_1[j + i * ny][Re] *= one_by_nxny; 
	  dfdeta_2[j + i * ny][Re] *= one_by_nxny; 
	  dfdeta_3[j + i * ny][Re] *= one_by_nxny; 
	  dfdeta_1[j + i * ny][Im] *= one_by_nxny; 
	  dfdeta_2[j + i * ny][Im] *= one_by_nxny; 
	  dfdeta_3[j + i * ny][Im] *= one_by_nxny; 
    dfdc[j + i * ny][Im] = 0.0;
   }
  }

for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = 0.0;
	 random_num1[i][j] = 0.0;
	 random_num2[i][j] = 0.0;
	 random_num3[i][j] = 0.0;
	}}
/* Sustained noise starts */
if(ncount < noise_steps){
/******************************************************************************
 * Apply sustained conserved noise to the compositions                        
 * ***************************************************************************/
for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = gsl_gaussian_sample(&SEED);
  }}
sum = 0.0;
for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = random_num[i][j] * sustained_noise_level_cons;
	 sum += random_num[i][j];
	}}
  mean = sum * one_by_nxny;

for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
	 random_num[i][j] = mean - random_num[i][j];
   dfdc[j + i * ny][Re] = dfdc[j + i * ny][Re] + random_num[i][j];
  }
 }
/***************************************************************************
 * Apply sustained nonconserved noise to the order parameters             *
 ***************************************************************************/

for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
	random_num1[i][j] = gsl_gaussian_sample(&SEED);
	random_num2[i][j] = gsl_gaussian_sample(&SEED);
	random_num3[i][j] = gsl_gaussian_sample(&SEED);
  }}
for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
	random_num1[i][j] = random_num1[i][j] * sustained_noise_level_uncons;
	random_num2[i][j] = random_num2[i][j] * sustained_noise_level_uncons;
	random_num3[i][j] = random_num3[i][j] * sustained_noise_level_uncons;
  }}
  
for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
	dfdeta_1[j + i * ny][Re] = dfdeta_1[j + i * ny][Re] + random_num1[i][j];
	dfdeta_2[j + i * ny][Re] = dfdeta_2[j + i * ny][Re] + random_num2[i][j];
	dfdeta_3[j + i * ny][Re] = dfdeta_3[j + i * ny][Re] + random_num3[i][j];
	}}/* End of sustained noise */
} 


/* Check for bounds */
/*  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    if (dfdc[j + i * ny][Re] < -0.4
        || dfdc[j + i * ny][Re] > 1.4) {
     printf ("compositions out of bounds. Exiting\n");
     exit (0);
    }
   }
  }
*/
/* Check for convergence */
  maxerror = 0.0;

  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    error = fabs (tempreal[i][j] - dfdc[j + i * ny][Re]);
    if (error > maxerror)
     maxerror = error;
   }
  }
  if (maxerror <= Tolerance) {
   printf ("maxerror=%lf\tnumbersteps=%d\n", maxerror, count);
   loop_condition = 0;
  }
  sim_time = sim_time + del_t;
  printf("time=%lf\n", sim_time);
  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    tempreal[i][j] = dfdc[j + i * ny][Re];
   }
  }
 }
 free_dmatrix (random_num, 0, nx, 0, ny);
 free_dmatrix (random_num1, 0, nx, 0, ny);
 free_dmatrix (random_num2, 0, nx, 0, ny);
 free_dmatrix (random_num3, 0, nx, 0, ny);
 free_dmatrix (tempreal, 0, nx, 0, ny);
 fftw_free(eta1sq);
 fftw_free(eta2sq);
 fftw_free(eta3sq);
 fftw_free(var_elast0);
 fftw_free(var_elast1);
 fftw_free(var_elast2);
 fftw_free(var_elast3);
}
