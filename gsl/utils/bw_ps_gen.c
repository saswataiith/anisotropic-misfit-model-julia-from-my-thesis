#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include"gsl_support.h"
#include"gsl_support.c"
int main (void)
{

 int N;
 FILE *fpt, *fpw;
 int i, j, t1, t2, tempr, tempg, templ;
 fftw_complex *comp, *eta1;
 double comp_l, comp_r, eta1_l, eta1_r;
 int tmp_y1, tmp_y2, interpolated_y_b, interpolated_y_t;
 int temp;
 char fn1[100], fn2[100], fn3[100];
 int count;
 int initial;
 int final;
 int steps;
 char new[100];
 printf("Enter the system size (N)\n");
 scanf("%d",&N);
 printf ("Initial count (initial)\n");
 scanf ("%d", &initial);
 printf ("Final count (final) \n");
 scanf ("%d", &final);
 printf ("Steps (steps) \n");
 scanf ("%d", &steps);
 comp = fftw_malloc (sizeof (fftw_complex) * N * N);
 eta1 = fftw_malloc (sizeof (fftw_complex) * N * N);
 for (count = initial; count <= final; count = count + steps) {
  sprintf (fn1, "comp.%06d", count);
  sprintf (fn1, "eta1.%06d", count);
  sprintf (fn2, "compcont%06d.dat",count);

  fpt = fopen (fn1, "r");

  fread (&comp[0][0], sizeof(double), 2 * N * N, fpt);

  fclose (fpt);
  
	fpt = fopen (fn1, "r");

  fread (&eta1[0][0], sizeof(double), 2 * N * N, fpt);

  fclose (fpt);

	 for(j = 0; j < N; ++j){
		interpolated_y_b = 0;			 
	  for(i = 1; i < N - 1; ++i){
		comp_l = comp[j + (i-1)*N][0];
		comp_r = comp[j + i*N][0];
		eta1_l = eta1[j + (i-1)*N][0];
		eta1_r = eta1[j + i*N][0];
		if (((comp_l < 0.5 && comp_r > 0.5) && (eta1_l < 0.5 && eta1_r > 0.5))){
		tmp_y1 = i - 1;
		tmp_y2 = i;
		interpolated_y_b =(int) (tmp_y1 + (tmp_y2 - tmp_y1) * 
				 (0.5 - comp[j + tmp_y1*N][0])/
				 (comp[j + tmp_y2 * N][0] - comp[j + tmp_y1 * N][0]));
	  printf("%d\t%d\n", j, interpolated_y_b);
		}
		}}
		
	 for(j = 0; j < N; ++j){
		interpolated_y_t = 0;			 
	  for(i = 1; i < N - 1; ++i){
		comp_l = comp[j + (i-1)*N][0];
		comp_r = comp[j + i*N][0];
		eta1_l = eta1[j + (i-1)*N][0];
		eta1_r = eta1[j + i*N][0];
		if(((comp_l > 0.5 && comp_r < 0.5) && (eta1_l > 0.5 && eta1_r < 0.5))){
		tmp_y1 = i - 1;
		tmp_y2 = i;
		interpolated_y_t =(int) (tmp_y1 + (tmp_y2 - tmp_y1) * 
				 (0.5 - comp[j + tmp_y1*N][0])/
				 (comp[j + tmp_y2 * N][0] - comp[j + tmp_y1 * N][0]));
	  printf("%d\t%d\n", j, interpolated_y_t);
		}
		}}
 }
 fftw_free(comp);
 fftw_free(eta1);
 return(0);
}
