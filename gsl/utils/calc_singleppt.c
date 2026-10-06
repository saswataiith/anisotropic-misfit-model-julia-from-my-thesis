#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include"gsl_support.h"
#include"gsl_support.c"
#define pi M_PI
int
main (void)
{
 int nx, ny, nx_half, ny_half;
 FILE *fpread, *fpin, *fpout;
 char fr[100],fw[100];
 double del_x, del_y;
 int initial, final, steps, count;
 char new[100];
 fftw_complex *comp, *eta1, *eta2, *eta3;
 int i, j;
 double interfacial_energy, mubar, epsilon_T;
 double size, characteristic_length, radius, area_frac;
 int area_counts;
 double I_xx, I_yy, I_xy, I_yx;
 double eigen_Ixx, eigen_Iyy, theta;
 double cm_x, cm_y, weighted_sum, sum_x, sum_y, sum, shape_param1, 
				shape_param2;
 double one_by_nxny;

 fpin = fopen("SystemParameters","r");
 fscanf(fpin,"%d",&nx);
 fscanf(fpin,"%d",&ny);
 fscanf(fpin,"%d", &initial);
 fscanf(fpin,"%d", &final);
 fscanf(fpin,"%d", &steps);
 fscanf(fpin,"%lf", &del_x);
 fscanf(fpin,"%lf", &del_y);
 fscanf(fpin,"%lf", &interfacial_energy);
 fscanf(fpin,"%lf", &mubar);
 fscanf(fpin,"%lf", &epsilon_T);
 fclose(fpin);
 one_by_nxny = 1.0/(double) (nx * ny);
 nx_half = nx/2;
 ny_half = ny/2;

 comp = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta1 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta2 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta3 = fftw_malloc (nx * ny * sizeof (fftw_complex));

 fpout = fopen("Single_ppt_params", "a");
 for(count=initial;count<=final;count=count+steps){

 sprintf (fr, "comp.%06d",count);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&comp[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);

 sprintf (fr, "eta1.%06d",count);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&eta1[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);

 sprintf (fr, "eta2.%06d",count);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&eta2[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);
 
 sprintf (fr, "eta3.%06d",count);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&eta3[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);

 /* Calculate the center of mass */

 cm_x = 0.0;
 cm_y = 0.0;
 sum  = 0.0;
 for(j = 0; j < ny; ++j){
	for (i = 0; i < nx; ++i){
	 if ((comp[i + j * nx][0] >= 0.5) && (eta1[i + j * nx][0] >= 0.5)){
	 sum += comp[i + j * nx][0];
	 cm_x += comp[i + j * nx][0] * i;
   cm_y += comp[i + j * nx][0] * j;
	 }
	}}
 cm_x = cm_x/sum;
 cm_y = cm_y/sum;
 printf("%lf\t%lf\n",cm_x, cm_y);

/* Calculate the size of the particle */

 area_counts = 0;
 for(j = 0; j < ny; ++j){
	for (i = 0; i < nx; ++i){
	 if ((comp[i + j * nx][0] >= 0.5) && (eta1[i + j * nx][0] >= 0.5)){
		area_counts++;
	 }
	}}
 area_frac = area_counts * one_by_nxny * nx * ny * del_x * del_y;
 radius = sqrt(area_frac/pi);
 characteristic_length = interfacial_energy/(mubar * epsilon_T * epsilon_T);
 printf("Characteristic_length = %lf\n", characteristic_length);
 size = radius/characteristic_length;

/* Calculate moment of inertia */
 I_xx = 0.0;
 I_xy = 0.0;
 I_yy = 0.0;
 for(j = 0; j < ny; ++j){
	for (i = 0; i < nx; ++i){
	 if ((comp[i + j * nx][0] >= 0.5) && (eta1[i + j * nx][0] >= 0.5)){
					 I_xx += (j-nx_half) * (j-nx_half) * del_x * del_y;
					 I_xy += (i - ny_half) * (j - nx_half) * del_x * del_y;
					 I_yy += (i - ny_half) * (i-ny_half) * del_x * del_y;
	 }
  }}
  if ( I_xx >= I_yy){ 
  eigen_Ixx = (I_xx + I_yy)/2.0 + sqrt(((I_xx - I_yy)*(I_xx - I_yy))/4.0 + 
							 I_xy * I_xy); 
  eigen_Iyy = (I_xx + I_yy)/2.0 - sqrt(((I_xx - I_yy)*(I_xx - I_yy))/4.0 + 
							 I_xy * I_xy);
	} 
	if ( I_xx < I_yy){
  eigen_Ixx = (I_xx + I_yy)/2.0 + sqrt(((I_yy - I_xx)*(I_yy - I_xx))/4.0 + 
							 I_xy * I_xy); 
  eigen_Iyy = (I_xx + I_yy)/2.0 - sqrt(((I_yy - I_xx)*(I_yy - I_xx))/4.0 + 
							 I_xy * I_xy);
	}
  theta = ((atan(-2.0 * I_xy/(I_xx - I_yy)))/2.0) * (180.0/pi);
	shape_param1 = (eigen_Ixx - eigen_Iyy)/(eigen_Ixx + eigen_Iyy);
	shape_param2 = (eigen_Iyy - eigen_Ixx)/(eigen_Ixx + eigen_Iyy);

 fprintf(fpout,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", count, size, eigen_Ixx, eigen_Iyy, theta, shape_param1, shape_param2);
 }
 fclose(fpout);
 fftw_free(comp);
 fftw_free(eta1);
 fftw_free(eta2);
 fftw_free(eta3);
}
