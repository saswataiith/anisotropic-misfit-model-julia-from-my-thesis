#include"binary.h"

int
main (void)
{
 void Get_Input_Spino (char *fnin, char *fnout);
 void Init_Conf_Rand ();
 void Init_Conf_File ();
 void Calculate_EigenStrain ();
 void calculate_Bn(double **B, double eigen_strain_a[3][3],double eigen_strain_b[3][3], double stress_a[3][3], double stress_b[3][3]);
 void Evolve ();

 double stress_0[3][3], stress_1[3][3], stress_2[3][3], stress_3[3][3];
 char finput[15] = "bin1ary";
 
 char fnin[15], fnout[15];
 
 FILE *fp;
 
 unsigned FLAG;
 /* int i, j; */

 if (!(fp = fopen (finput, "r"))) {
  printf ("File:%s could not be opened\n", finput);
  exit (EXIT_FAILURE);
 }
  fscanf (fp, "%s", fnin);
  fscanf (fp, "%s", fnout);
 if (!(fpout = fopen (fnout, "w"))) {
  printf ("File:%s could not be opened\n", fnout);
  exit (EXIT_FAILURE);
 }
  fclose (fp);

 Get_Input_Spino (fnin, fnout);
 
 comp = fftw_malloc (sizeof (fftw_complex) * nx * ny);
 eta_1  = fftw_malloc (sizeof (fftw_complex) * nx * ny); 
 eta_2  = fftw_malloc (sizeof (fftw_complex) * nx * ny); 
 eta_3  = fftw_malloc (sizeof (fftw_complex) * nx * ny); 
 dfdc = fftw_malloc (sizeof (fftw_complex) * nx * ny);
 dfdeta_1 = fftw_malloc (sizeof (fftw_complex) * nx * ny);
 dfdeta_2 = fftw_malloc (sizeof (fftw_complex) * nx * ny);
 dfdeta_3 = fftw_malloc (sizeof (fftw_complex) * nx * ny);
 B_00 = dmatrix(0, nx, 0, ny);
 B_10 = dmatrix(0, nx, 0, ny);
 B_20 = dmatrix(0, nx, 0, ny);
 B_30 = dmatrix(0, nx, 0, ny);
 B_11 = dmatrix(0, nx, 0, ny);
 B_21 = dmatrix(0, nx, 0, ny);
 B_31 = dmatrix(0, nx, 0, ny);
 B_22 = dmatrix(0, nx, 0, ny);
 B_32 = dmatrix(0, nx, 0, ny);
 B_33 = dmatrix(0, nx, 0, ny);


 nx_half = nx / 2;

 ny_half = ny / 2;

 one_by_nxny = 1.0 / (double) (nx * ny);


//      printf("%lf\t%lf\t%lf\n",mobil_bb,mobil_cc,mobil_bc);


 FLAG = FFTW_ESTIMATE;
 
 if(fftw_flag == 1)
 FLAG = FFTW_MEASURE;
 
 if(fftw_flag == 2)
 FLAG = FFTW_PATIENT;
 
 if(fftw_flag == 3)
 FLAG = FFTW_EXHAUSTIVE;
 
 p_up = fftw_plan_dft_2d (nx, ny, comp, comp, FFTW_FORWARD,
                          FLAG);
 p_dn = fftw_plan_dft_2d (nx, ny, comp, comp, FFTW_BACKWARD,
                          FLAG);

 if (flag == 0) {
  Init_Conf_Rand ();
 }
 else
  Init_Conf_File ();
 
 sim_time = 0.0;
 
 count = 0;

 Calculate_EigenStrain ();
 calculate_Bn(B_00, eigen_strain_0, eigen_strain_0, stress_0, stress_0); 
 calculate_Bn(B_10, eigen_strain_1, eigen_strain_0, stress_1, stress_0); 
 calculate_Bn(B_20, eigen_strain_2, eigen_strain_0, stress_2, stress_0); 
 calculate_Bn(B_30, eigen_strain_3, eigen_strain_0, stress_3, stress_0); 
 calculate_Bn(B_11, eigen_strain_1, eigen_strain_1, stress_1, stress_1); 
 calculate_Bn(B_21, eigen_strain_2, eigen_strain_1, stress_2, stress_1); 
 calculate_Bn(B_31, eigen_strain_3, eigen_strain_1, stress_3, stress_1); 
 calculate_Bn(B_22, eigen_strain_2, eigen_strain_2, stress_2, stress_2); 
 calculate_Bn(B_32, eigen_strain_3, eigen_strain_2, stress_3, stress_2); 
 calculate_Bn(B_33, eigen_strain_3, eigen_strain_3, stress_3, stress_3); 


 
 Evolve ();
 

 fclose (fpout);
 
 fftw_destroy_plan (p_up);
 fftw_destroy_plan (p_dn);
 fftw_free (comp);
 fftw_free (eta_1);
 fftw_free (eta_2);
 fftw_free (eta_3);
 fftw_free (dfdc);
 fftw_free (dfdeta_1);
 fftw_free (dfdeta_2);
 fftw_free (dfdeta_3);
 free_dmatrix(B_00, 0, nx, 0, ny);
 free_dmatrix(B_10, 0, nx, 0, ny);
 free_dmatrix(B_20, 0, nx, 0, ny);
 free_dmatrix(B_30, 0, nx, 0, ny);
 free_dmatrix(B_11, 0, nx, 0, ny);
 free_dmatrix(B_21, 0, nx, 0, ny);
 free_dmatrix(B_31, 0, nx, 0, ny);
 free_dmatrix(B_22, 0, nx, 0, ny);
 free_dmatrix(B_32, 0, nx, 0, ny);
 free_dmatrix(B_33, 0, nx, 0, ny);
 return (0);
}
