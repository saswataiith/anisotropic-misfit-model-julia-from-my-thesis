#include"binary.h"

int
main (void)
{
 void Get_Input_Spino (char *fnin, char *fnout);
 void Init_Conf_Rand ();
 void Init_Conf_File ();
 void Evolve ();

 
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
 return (0);
}
