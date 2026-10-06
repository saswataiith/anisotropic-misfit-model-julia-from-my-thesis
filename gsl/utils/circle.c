#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<fftw3.h>
#define nx 1024
#define ny 1024
#define Re 0
#define Im 1

int main(void)
{
  FILE *fpt1, *fpt2, *fpt3, *fpt4;
  fftw_complex *comp;
  fftw_complex *eta1, *eta2, *eta3;
  double c_matrix;
  char fn2[100];
  int o_sites, i, j;
	double rad, del_x;
  int nx_half, ny_half;
  double comp_b, total;
  int occupancy[nx][ny];
  int term1;
  int term2;


  comp = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  eta1 = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  eta2 = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  eta3 = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      comp[j + i * ny][Re] = 0.0;
      eta1[j + i * ny][Re] = 0.0;
      eta2[j + i * ny][Re] = 0.0;
      eta3[j + i * ny][Re] = 0.0;
      comp[j + i * ny][Im] = 0.0;
      eta1[j + i * ny][Im] = 0.0;
      eta2[j + i * ny][Im] = 0.0;
      eta3[j + i * ny][Im] = 0.0;
    }
  }
  fputs("Enter the initial composition\n", stdout);
  scanf("%lf", &c_matrix);
  fputs("Enter the radius\n", stdout);
  scanf("%lf", &rad);
  fputs("Enter the grid size\n", stdout);
  scanf("%lf", &del_x);
  term1 = ceil(100 * c_matrix);
  nx_half = nx/2;
  ny_half = ny/2;
  comp_b = 1.0;
  rad = rad/del_x;

    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
      if ((double)((i-nx_half) * (i-nx_half) +  (j-ny_half) * (j-ny_half)) 
           <= (rad * rad)){
	    comp[j + i*ny][Re] = 1.0;
	    eta1[j + i*ny][Re] = 1.0;
	    eta2[j + i*ny][Re] = 1.0;
	    eta3[j + i*ny][Re] = 1.0;
	  }
      }}
  total = 0.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      total = total + comp[j + i * ny][Re];
    }
  }
  total = total / (double) (nx * ny);
  printf("%lf\n", total);
  o_sites = 1;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      if ((comp[j+i*ny][Re] < DBL_EPSILON)){
    o_sites = o_sites + 1;
    }
    }
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      if ((comp[j+i*ny][Re] < DBL_EPSILON)){
	comp[j + i * ny][Re] = c_matrix;
    }
    }
  }

  total = 0.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      total = total + comp[j + i * ny][Re];
    }
  }
  total = total / (double) (nx * ny);
  printf("%lf\n", total);
  fpt1 = fopen("comp.000000", "w");
  fpt2 = fopen("eta1.000000", "w");
  fpt3 = fopen("eta2.000000", "w");
  fpt4 = fopen("eta3.000000", "w");
  fwrite(&comp[0][0], sizeof(double), 2 * nx * ny, fpt1);
  fwrite(&eta1[0][0], sizeof(double), 2 * nx * ny, fpt2);
  fwrite(&eta2[0][0], sizeof(double), 2 * nx * ny, fpt3);
  fwrite(&eta3[0][0], sizeof(double), 2 * nx * ny, fpt4);
  fclose(fpt1);
  fclose(fpt2);
  fclose(fpt3);
  fclose(fpt4);
  fftw_free(comp);
  fftw_free(eta1);
  fftw_free(eta2);
  fftw_free(eta3);
}
