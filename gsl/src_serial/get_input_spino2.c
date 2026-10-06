#include"binary.h"
void
Get_Input_Spino (char *fnin, char *fnout)
{
 FILE *fpin, *fpcout;
 char  param[100];

 if (!(fpcout = fopen (fnout, "w"))) {
  printf ("File:%s could not be opened \n", fnout);
  exit (1);
 }
 fprintf (fpcout, "The name of this file is : %s \n", fnout);
 fprintf (fpcout, "Input is from            : %s \n", fnin);

 if (!(fpin = fopen (fnin, "r"))) {
  printf ("File: %s could not be opened \n", fnin);
  exit (1);
 }

 fscanf (fpin, "%s%d", param,&nx);
 fscanf (fpin, "%s%d", param,&ny);
 fscanf (fpin, "%s%lf",param,&del_x);
 fscanf (fpin, "%s%lf",param,&del_y);
 fscanf (fpin, "%s%lf",param,&del_t1);
 fscanf (fpin, "%s%lf",param,&del_t2);
 fscanf (fpin, "%s%lf",param,&noise_level);
 fscanf (fpin, "%s%lf",param,&sustained_noise_level_cons);
 fscanf (fpin, "%s%lf",param,&sustained_noise_level_uncons_1);
 fscanf (fpin, "%s%lf",param,&sustained_noise_level_uncons_2);
 fscanf (fpin, "%s%lf",param,&sustained_noise_level_uncons_3);
 fscanf (fpin, "%s%ld",param,&SEED);
 fscanf (fpin, "%s%d", param,&num_steps);
 fscanf (fpin, "%s%d", param,&noise_steps);
 fscanf (fpin, "%s%lf",param,&alloycomp);
 fscanf (fpin, "%s%lf",param,&A1);
 fscanf (fpin, "%s%lf",param,&A2);
 fscanf (fpin, "%s%lf",param,&A41);
 fscanf (fpin, "%s%lf",param,&A42);
 fscanf (fpin, "%s%lf",param,&A61);
 fscanf (fpin, "%s%lf",param,&A62);
 fscanf (fpin, "%s%lf",param,&kappa);
 fscanf (fpin, "%s%lf",param,&mobility);
 fscanf (fpin, "%s%lf",param,&kappa1);
 fscanf (fpin, "%s%lf",param,&kappa2);
 fscanf (fpin, "%s%lf",param,&kappa3);
 fscanf (fpin, "%s%lf",param,&L1);
 fscanf (fpin, "%s%lf",param,&L2);
 fscanf (fpin, "%s%lf",param,&L3);
 fscanf (fpin, "%s%lf",param,&epsc);
 fscanf (fpin, "%s%lf",param,&eps_eta);
 fscanf (fpin, "%s%lf",param,&tetra);
 fscanf (fpin, "%s%lf",param,&mubar);
 fscanf (fpin, "%s%lf",param,&nubar);
 fscanf (fpin, "%s%lf",param,&Aniso);
 fscanf (fpin, "%s%d",param,&steps_t1);
 fscanf (fpin, "%s%d",param,&steps_t2);
 fscanf (fpin, "%s%d",param,&time_to_change);
 fscanf (fpin, "%s%d", param,&initcount);
 fscanf (fpin, "%s%d", param,&flag);
 fscanf (fpin, "%s%d", param,&fftw_flag);
 fclose (fpin);

 printf ("nx = %d\n", nx);
 printf ("ny = %d\n", ny);
 printf ("del_x = %lf\n", del_x);
 printf ("del_y = %lf\n", del_y);
 printf ("del_t1 = %lf\n", del_t1);
 printf ("del_t2 = %lf\n", del_t2);
 printf ("noise_level = %lf\n", noise_level);
 printf ("sustained noise level (conserved) = %lf\n", 
		sustained_noise_level_cons);
 printf ("sustained noise level (unconserved) varaiant 1 = %lf\n", 
	    sustained_noise_level_uncons_1);
 printf ("sustained noise level (unconserved) varaiant 2 = %lf\n", 
	    sustained_noise_level_uncons_2);
 printf ("sustained noise level (unconserved) varaiant 3 = %lf\n", 
	    sustained_noise_level_uncons_3);
 printf ("Seed = %ld\n", SEED);
 printf ("num_steps = %d\n", num_steps);
 printf ("noise_steps = %d\n", noise_steps);
 printf ("composition_of_alloy = %lf\n", alloycomp);
 printf ("Coefficients A1 = %lf  A2 = %lf A41 = %lf A42 = %lf A61 = %lf A62 = %lf\n",A1,A2,A41,A42,A61,A62);
 printf ("KAPPA_c=%lf\n", kappa);
 printf ("Mobility=%lf\n", mobility);
 printf ("kappa1 = %lf kappa2 = %lf kappa3 = %lf\n", kappa1, kappa2, kappa3);
 printf ("L1 = %lf L2 = %lf L3 = %lf\n",L1, L2, L3);
 printf ("Lattice expansion coefficient epsc = %lf\n", epsc);
 printf ("Lattice expansion coefficient eps_eta = %lf\n", eps_eta);
 printf ("Tetragonality = %lf\n", tetra);
 printf ("Shear modulus = %lf\n", mubar);
 printf ("Poisson's ratio = %lf\n", nubar);
 printf ("Zener Anisotropy Parameter = %lf\n", Aniso);
 printf ( "Time interval for output (before switching of noise) = %d\n", steps_t1);
 printf ( "Time interval for output (after switching of noise) = %d\n", steps_t2);
 printf ( "Time to change (after switching of noise) = %d\n", time_to_change);
 printf ("Initial count = %d\t flag = %d\n", initcount, flag);
 printf ( "FFTW_FLAG=%d\n", fftw_flag);
 
 fprintf (fpcout, "nx = %d\n", nx);
 fprintf (fpcout, "ny = %d\n", ny);
 fprintf (fpcout, "del_x = %lf\n", del_x);
 fprintf (fpcout, "del_y = %lf\n", del_y);
 fprintf (fpcout, "del_t1 = %lf\n", del_t1);
 fprintf (fpcout, "del_t2 = %lf\n", del_t2);
 fprintf (fpcout, "noise_level = %lf\n", noise_level);
 fprintf (fpcout, "sustained noise level (conserved) = %lf\n", 
		  sustained_noise_level_cons);
 fprintf (fpcout, "sustained noise level (unconserved) variant 1= %lf\n", 
	      sustained_noise_level_uncons_1);
 fprintf (fpcout, "sustained noise level (unconserved) variant 2= %lf\n", 
	      sustained_noise_level_uncons_2);
 fprintf (fpcout, "sustained noise level (unconserved) variant 3= %lf\n", 
	      sustained_noise_level_uncons_3);
 fprintf (fpcout, "Seed = %ld\n", SEED);
 fprintf (fpcout, "num_steps = %d\n", num_steps);
 fprintf (fpcout, "noise_steps = %d\n", noise_steps);
 fprintf (fpcout, "composition_of_alloy=%lf\n", alloycomp);
 fprintf (fpcout, "Coefficients A1 = %lf  A2 = %lf A41 = %lf A42 = %lf A61 = %lf A62 = %lf\n",A1,A2,A41,A42,A61,A62);
 fprintf (fpcout, "\n");
 fprintf (fpcout, "KAPPA_c=%lf\n", kappa);
 fprintf (fpcout, "Mobility=%lf\n", mobility);
 fprintf (fpcout, "kappa1 = %lf kappa2 = %lf kappa3 = %lf\n",kappa1, kappa2, kappa3);
 fprintf (fpcout, "L1 = %lf L2 = %lf L3 = %lf\n", L1, L2, L3);
 fprintf (fpcout, "Lattice expansion coefficient epsc = %lf\n", epsc);
 fprintf (fpcout, "Lattice expansion coefficient eps_eta = %lf\n", eps_eta);
 fprintf (fpcout, "Tetragonality = %lf\n", tetra);
 fprintf (fpcout, "Shear modulus = %lf\n", mubar);
 fprintf (fpcout, "Poisson's ratio = %lf\n", nubar);
 fprintf (fpcout, "Zener Anisotropy Parameter = %lf\n", Aniso);
 fprintf (fpcout, "Time interval for output (before switching of noise) = %d\n", steps_t1);
 fprintf (fpcout, "Time interval for output (after switching of noise) = %d\n", steps_t2);
 fprintf (fpcout, "Time to change (after switching of noise) = %d\n", time_to_change);
 fprintf (fpcout, "Initial count = %d\t flag = %d\n", initcount, flag);
 if (flag == 0)
  fprintf (fpcout, "Configuration is initialised by me\n");
 else
  fprintf (fpcout, "Configuration is read from file\n");
 fclose (fpcout);
}
