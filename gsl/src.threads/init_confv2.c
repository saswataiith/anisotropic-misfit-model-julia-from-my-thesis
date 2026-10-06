#include"binary.h"

void
Init_Conf_Rand ()
{
 double gsl_gaussian_sample(long *idum);
 double gsl_uniform_sample(long *idum);
 double **random_num;
 double **random_num1, **random_num2, **random_num3; 
 double sum = 0.0, mean;
 int i, j;

 random_num = dmatrix (0, nx, 0, ny);
 random_num1 = dmatrix (0, nx, 0, ny);
 random_num2 = dmatrix (0, nx, 0, ny);
 random_num3 = dmatrix (0, nx, 0, ny);

 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   comp[j + i * ny][Re] = alloycomp;
   comp[j + i * ny][Im] = 0.0;
  }
 }
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
		  
		  eta_1[j + i * ny][Re] = 0.0;
		  eta_2[j + i * ny][Re] = 0.0;
		  eta_3[j + i * ny][Re] = 0.0;
		  eta_1[j + i * ny][Im] = 0.0;
		  eta_2[j + i * ny][Im] = 0.0;
		  eta_3[j + i * ny][Im] = 0.0;
    } 
  }		  
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = gsl_uniform_sample(&SEED);
   random_num[i][j] = (2.0 * random_num[i][j]) - 1.0;
   random_num[i][j] = random_num[i][j] * noise_level * alloycomp;
   sum += random_num[i][j];
  }
 }
 mean = sum / (nx * ny);
 
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = mean - random_num[i][j];
   comp[j + i * ny][Re] = comp[j + i * ny][Re] + random_num[i][j];
  }
 }
 
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
  random_num1[i][j] = gsl_gaussian_sample(&SEED);
  random_num2[i][j] = gsl_gaussian_sample(&SEED);
  random_num3[i][j] = gsl_gaussian_sample(&SEED);
  }
 }

 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
		  eta_1[j + i * ny][Re] += random_num1[i][j] * noise_level;
		  eta_2[j + i * ny][Re] += random_num2[i][j] * noise_level;
		  eta_3[j + i * ny][Re] += random_num3[i][j] * noise_level;
  }}
          
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   dfdc[j + i * ny][Re] = comp[j + i * ny][Re];
   dfdc[j + i * ny][Im] = comp[j + i * ny][Im];
   dfdeta_1[j + i * ny][Re] = eta_1[j + i * ny][Re];
   dfdeta_2[j + i * ny][Re] = eta_2[j + i * ny][Re];
   dfdeta_3[j + i * ny][Re] = eta_3[j + i * ny][Re];
   dfdeta_1[j + i * ny][Im] = eta_1[j + i * ny][Im];
   dfdeta_2[j + i * ny][Im] = eta_2[j + i * ny][Im];
   dfdeta_3[j + i * ny][Im] = eta_3[j + i * ny][Im];
  }
 }
 free_dmatrix (random_num, 0, nx, 0, ny);
 free_dmatrix (random_num1, 0, nx, 0, ny);
 free_dmatrix (random_num2, 0, nx, 0, ny);
 free_dmatrix (random_num3, 0, nx, 0, ny);
}

void
Init_Conf_File ()
{
 double gsl_gaussian_sample(long *idum);
 double gsl_uniform_sample(long *idum);
 FILE *fpread;
 char fr[100], fr1[100], fr2[100], fr3[100];
 int i, j;
 double **random_num;
 double **random_num1, **random_num2, **random_num3; 
 double sum = 0.0, mean;

 random_num = dmatrix (0, nx, 0, ny);
 random_num1 = dmatrix (0, nx, 0, ny);
 random_num2 = dmatrix (0, nx, 0, ny);
 random_num3 = dmatrix (0, nx, 0, ny);
 sprintf (fr,  "comp.%06d", initcount);
 sprintf (fr1, "eta1.%06d", initcount);
 sprintf (fr2, "eta2.%06d", initcount);
 sprintf (fr3, "eta3.%06d", initcount);

 fpread = fopen (fr, "r");
 fread (&comp[0][0], sizeof(double), 2 * nx * ny, fpread);
 fclose (fpread);

 fpread = fopen (fr1, "r");
 fread (&eta_1[0][0], sizeof(double), 2 * nx * ny, fpread);
 fclose (fpread);
 
 fpread = fopen (fr2, "r");
 fread (&eta_2[0][0], sizeof(double), 2 * nx * ny, fpread);
 fclose (fpread);

 fpread = fopen (fr3, "r");
 fread (&eta_3[0][0], sizeof(double), 2 * nx * ny, fpread);
 fclose (fpread);
 
 
 
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = gsl_uniform_sample(&SEED);
   random_num[i][j] = (2.0 * random_num[i][j]) - 1.0;
   random_num[i][j] = random_num[i][j] * noise_level;
   sum += random_num[i][j];
  }
 }
 mean = sum / (nx * ny);


 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   random_num[i][j] = mean - random_num[i][j];
   comp[j + i * ny][Re] = comp[j + i * ny][Re] + random_num[i][j];
  }
 }

 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
  random_num1[i][j] = gsl_uniform_sample(&SEED);
  random_num2[i][j] = gsl_uniform_sample(&SEED);
  random_num3[i][j] = gsl_uniform_sample(&SEED);
  }
 }
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
	  eta_1[j + i * ny][Re] += random_num1[i][j] * noise_level;
	  eta_2[j + i * ny][Re] += random_num2[i][j] * noise_level;
	  eta_3[j + i * ny][Re] += random_num3[i][j] * noise_level;
  }}

 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
   dfdc[j + i * ny][Re] = comp[j + i * ny][Re];
   dfdc[j + i * ny][Im] = comp[j + i * ny][Im];
   dfdeta_1[j + i * ny][Re] = eta_1[j + i * ny][Re];
   dfdeta_2[j + i * ny][Re] = eta_2[j + i * ny][Re];
   dfdeta_3[j + i * ny][Re] = eta_3[j + i * ny][Re];
   dfdeta_1[j + i * ny][Im] = eta_1[j + i * ny][Im];
   dfdeta_2[j + i * ny][Im] = eta_2[j + i * ny][Im];
   dfdeta_3[j + i * ny][Im] = eta_3[j + i * ny][Im];
  }
 }
free_dmatrix (random_num, 0, nx, 0, ny);
free_dmatrix (random_num1, 0, nx, 0, ny);
free_dmatrix (random_num2, 0, nx, 0, ny);
free_dmatrix (random_num3, 0, nx, 0, ny);
}
