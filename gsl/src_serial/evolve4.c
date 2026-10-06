#include"ternary.h"
void
Evolve ()
{

 double **random_numA, **random_numB;
 
 void Output_Conf ();
 
 int i, j, loop_condition, ncount;
 
 double del_kx, del_ky, kx, ky, kpow2, kpow4, kpow6, kpow8;
 
 double rc_b, rc_c, rc_a, fp_b, fp_c, rc_new_b, rc_new_c;
 
 double tmp_b, tmp_c, sumcomp;
 
 double **tempreal_b, **tempreal_c;
 
 double lhs, rhs_b, rhs_c;
 
 double error1, error2, maxerror1, maxerror2;
 
 double Pbb, Pbc, Pcb, Pcc, Qbb, Qbc, Qcb, Qcc, Tbb, Tbc, Tcb, Tcc;
 
 gsl_rng *r;

 gsl_rng_env_setup ();
 
 T = gsl_rng_taus2;
 
 r = gsl_rng_alloc (T);
 
 tempreal_b = dmatrix (0, nx, 0, ny);
 tempreal_c = dmatrix (0, nx, 0, ny);
 
 random_numA = dmatrix (0, nx, 0, ny);
 random_numB = dmatrix (0, nx, 0, ny);
 
 del_kx = 2.0 * PI / ((double) nx * del_x);
 del_ky = 2.0 * PI / ((double) ny * del_y);

 Qbb = mobil_bb*kappa_bb + mobil_bc*kappa_bc;
 Qbc = mobil_bb*kappa_bc + mobil_bc*kappa_cc;
 Qcb = mobil_bc*kappa_bb + mobil_cc*kappa_bc;
 Qcc = mobil_bc*kappa_bc + mobil_cc*kappa_cc;
 
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   tempreal_b[i][j] = comp_b[j + i * ny][Re];
   tempreal_c[i][j] = comp_c[j + i * ny][Re];
  }
 }

 fftw_execute_dft (p_up, comp_b, comp_b);
 fftw_execute_dft (p_up, comp_c, comp_c);

 loop_condition = 1;
 for (count = 0; count <= num_steps; count++) {
     ncount = count + initcount;
  if ((( ncount <= 5000 ) && (( ncount % 1000 ) == 0)) ||
      (( ncount > 5000 && ncount <= 50000 ) && 
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

  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    rc_b = dfdc_b[j + i * ny][Re];
    rc_c = dfdc_c[j + i * ny][Re];
    sumcomp = dfdc_b[j + i * ny][Re] + dfdc_c[j + i * ny][Re];
    rc_a = 1.0 - sumcomp;
    dfdc_b[j + i * ny][Re] =
     L * 2.0 * (rc_b * (rc_c * rc_c + rc_a * rc_a) -
                rc_a * (rc_b * rc_b + rc_c * rc_c)) +
     M * 2.0 * (rc_b * rc_c * rc_c * rc_a) * (rc_a - rc_b);
    dfdc_b[j + i * ny][Im] = 0.0;

    dfdc_c[j + i * ny][Re] =
     L * 2.0 * (rc_c * (rc_b * rc_b + rc_a * rc_a) -
                rc_a * (rc_b * rc_b + rc_c * rc_c)) +
     M * 2.0 * (rc_b * rc_b * rc_c * rc_a) * (rc_a - rc_c);
    dfdc_c[j + i * ny][Im] = 0.0;

   }
  }

  fftw_execute_dft (p_up, dfdc_b, dfdc_b);
  fftw_execute_dft (p_up, dfdc_c, dfdc_c);

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
    kpow6 = kpow2 * kpow4;
    kpow8 = kpow4 * kpow4;
   
    Pbb = mobil_bb * B_bb[i][j] + mobil_bc * B_bc[i][j];
    Pbc = mobil_bb * B_bc[i][j] + mobil_bc * B_cc[i][j];
    Pcb = mobil_bc * B_bb[i][j] + mobil_cc * B_bc[i][j];
    Pcc = mobil_bc * B_bc[i][j] + mobil_cc * B_cc[i][j];
 
    Tbb = 1.0 + kpow2 * del_t * Pbb + 2.0 * kpow4 * del_t * Qbb;
    Tbc = kpow2 * del_t * Pbc + 2.0 * kpow4 * del_t * Qbc;
    Tcb = kpow2 * del_t * Pcb + 2.0 * kpow4 * del_t * Qcb;
    Tcc = 1.0 + kpow2 * del_t * Pcc + 2.0 * kpow4 * del_t * Qcc;

    lhs = 1.0 + kpow2 * del_t * (Pbb + Pcc) + 2.0 * kpow4 * del_t * (Qbb +
          Qcc) + kpow4 * del_t * del_t * ( Pbb * Pcc - Pcb * Pbc) + 2.0 *
          kpow6 * del_t * del_t * (Pbb * Qcc + Pcc * Qbb - Pcb * Qbc 
          - Qcb * Pbc) + 4.0 * kpow8 * del_t * del_t * (Qbb * Qcc - Qcb *
          Qbc); 
    
    rc_b = comp_b[j + i * ny][Re];
    rc_c = comp_c[j + i * ny][Re];
    fp_b = dfdc_b[j + i * ny][Re];
    fp_c = dfdc_c[j + i * ny][Re];
    tmp_b = rc_b - kpow2 * del_t * (mobil_bb * fp_b + mobil_bc * fp_c);
    tmp_c = rc_c - kpow2 * del_t * (mobil_bc * fp_b + mobil_cc * fp_c);
    rhs_b = Tcc * tmp_b - Tbc * tmp_c;
    rhs_c = Tbb * tmp_c - Tcb * tmp_b;
    rc_new_b = rhs_b / lhs;
    rc_new_c = rhs_c / lhs;
    comp_b[j + i * ny][Re] = rc_new_b;
    comp_c[j + i * ny][Re] = rc_new_c;
    dfdc_b[j + i * ny][Re] = comp_b[j + i * ny][Re];
    dfdc_c[j + i * ny][Re] = comp_c[j + i * ny][Re];

    rc_b = comp_b[j + i * ny][Im];
    rc_c = comp_c[j + i * ny][Im];
    fp_b = dfdc_b[j + i * ny][Im];
    fp_c = dfdc_c[j + i * ny][Im];
    tmp_b = rc_b - kpow2 * del_t * (mobil_bb * fp_b + mobil_bc * fp_c);
    tmp_c = rc_c - kpow2 * del_t * (mobil_bc * fp_b + mobil_cc * fp_c);
    rhs_b = Tcc * tmp_b - Tbc * tmp_c;
    rhs_c = Tbb * tmp_c - Tcb * tmp_b;
    rc_new_b = rhs_b / lhs;
    rc_new_c = rhs_c / lhs;
    comp_b[j + i * ny][Im] = rc_new_b;
    comp_c[j + i * ny][Im] = rc_new_c;
    dfdc_b[j + i * ny][Im] = comp_b[j + i * ny][Im];
    dfdc_c[j + i * ny][Im] = comp_c[j + i * ny][Im];
   }
  }


/* Check for conservation of mass */
  total_A = dfdc_b[0][Re] * one_by_nxny;
  total_B = dfdc_c[0][Re] * one_by_nxny;
//  printf ("total_A=%le total_B=%le\n", total_A, total_B);
  err1 = fabs (total_A - alloycomp_b);
  err2 = fabs (total_B - alloycomp_c);
  if (err1 > COMPERR || err2 > COMPERR) {
   printf ("ELEMENTS ARE NOT CONSERVED,SORRY!!!!\n");
   printf ("error_a=%lf error_b=%lf\n", err1, err2);
   exit (0);
  }
  

  fftw_execute_dft (p_dn, dfdc_b, dfdc_b);
  fftw_execute_dft (p_dn, dfdc_c, dfdc_c);

  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    dfdc_b[j + i * ny][Re] = dfdc_b[j + i * ny][Re] * one_by_nxny;
    dfdc_b[j + i * ny][Im] = dfdc_b[j + i * ny][Im] * one_by_nxny;
    dfdc_c[j + i * ny][Re] = dfdc_c[j + i * ny][Re] * one_by_nxny;
    dfdc_c[j + i * ny][Im] = dfdc_c[j + i * ny][Im] * one_by_nxny;
    dfdc_b[j + i * ny][Im] = 0.0;
    dfdc_c[j + i * ny][Im] = 0.0;
   }
  }


/* Check for bounds */
  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    if (dfdc_b[j + i * ny][Re] < -0.2
        || dfdc_b[j + i * ny][Re] > 1.2
        || dfdc_c[j + i * ny][Re] < -0.2 || dfdc_c[j + i * ny][Re] > 1.2) {
     printf ("compositions out of bounds. Exiting\n");
     exit (0);
    }
   }
  }


/* Check for convergence */
  maxerror1 = 0.0;
  maxerror2 = 0.0;

  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    error1 = fabs (tempreal_b[i][j] - dfdc_b[j + i * ny][Re]);
    error2 = fabs (tempreal_c[i][j] - dfdc_c[j + i * ny][Re]);
    if (error1 > maxerror1)
     maxerror1 = error1;
    if (error2 > maxerror2)
     maxerror2 = error2;
   }
  }
  if (maxerror1 <= Tolerance && maxerror2 <= Tolerance) {
   printf ("maxerror=%lf\t%lf\tnumbersteps=%d\n", maxerror1, maxerror2,
           count);
   loop_condition = 0;
  }
  sim_time = sim_time + del_t;
  printf("time=%lf\n", sim_time);
  for (i = 0; i < nx; i++) {
   for (j = 0; j < ny; j++) {
    tempreal_b[i][j] = dfdc_b[j + i * ny][Re];
    tempreal_c[i][j] = dfdc_c[j + i * ny][Re];
   }
  }
 }
 free_dmatrix (tempreal_b, 0, nx, 0, ny);
 free_dmatrix (tempreal_c, 0, nx, 0, ny);
 free_dmatrix (random_numA, 0, nx, 0, ny);
 free_dmatrix (random_numB, 0, nx, 0, ny);
 gsl_rng_free (r);
}
