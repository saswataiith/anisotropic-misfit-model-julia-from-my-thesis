#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<fftw3.h>
#include"gsl_support.h"
#include"gsl_support.c"
#include"gsl_uniform_sample.c"
#define Re 0
#define Im 1
#define TRUE 0
#define FALSE 1

/***************************************************************************** 
* Profile generator for a random distribution of particles version 6.0
* This program generates a random distribution of precipitates inside a
* box of given dimensions from a given distribution of nuclei. The 
* compositions   of the
* nuclei (cb & cc), the radius of the nuclei, the number of nuclei with
* the same radius, and the total number of nuclei are read from the file.
* Overall alloy compositon (cb_matrix & cc_matrix) are given as input from the
* keyboard. This program ensures that the nuclei do not overlap and they are
* distributed in the box with periodic boundary conditions. 
*******************************************************************************/

int main(void)
{
  int nx, ny;
  long SEED = -1111;
  FILE *fpt1, *fpt2;
  fftw_complex *comp_b, *comp_c;
  double ca[101], cb[101], radius[101];
  int n_ppt[101], state, **occupancy;
  int tmp1, tmp2, tmp3, tmp4;
  int new1, new2;
  int iter, nparticles;
  int count, i, j, l, m, p, q, r, s;
  int ncount, flag[400], check;
  double rad[2000], compb[2000], compc[2000];
  int x[2000], y[2000];
  double dist, distance[2000][2000];
  int randnum1[400], randnum2[400];
  double alloycomp_b, alloycomp_c, total_b, total_c;
  int term1, term2, term3, term4, term5;
  double cb_matrix, cc_matrix;
  char fn1[100], fn2[100];
  int o_sites, drad[2000];
  printf("Enter the domain size nx x ny\n");
  scanf("%d%d",&nx, &ny);
  occupancy = imatrix(0,nx,0,ny);
  comp_b = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  comp_c = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      comp_b[j + i * ny][Re] = 0.0;
      comp_c[j + i * ny][Re] = 0.0;
      comp_b[j + i * ny][Im] = 0.0;
      comp_c[j + i * ny][Im] = 0.0;
    }
  }
  fputs("Enter the initial composition\n", stdout);
  scanf("%lf", &cb_matrix);
  scanf("%lf", &cc_matrix);
// new1 = ceil(100 * cb_matrix);
// new2 = ceil(100 * cc_matrix);
// sprintf(fn1, "Nucleate_compB%02d_compC%02d", new1, new2); 
// sprintf(fn2, "profile_cb%02d_cc%02d", new1, new2); 
 sprintf(fn1, "Nucleate_compB%02d_compC%02d", (int) ceil(100 * cb_matrix), 
                                              (int) ceil(100 * cc_matrix));
 sprintf(fn2, "profile_cb%02d_cc%02d",  (int) ceil(100 * cb_matrix),
                                        (int) ceil(100 * cc_matrix)); 
 fpt1 = fopen(fn1, "r");

/* Reading relevant information from the file Nucleation_cb_cc */

  for (count = 0; count <= 100; count++) {
    fread(&ca[count], sizeof(double), (size_t) 1, fpt1);
    fread(&cb[count], sizeof(double), (size_t) 1, fpt1);
    fread(&radius[count], sizeof(double), (size_t) 1, fpt1);
    fread(&n_ppt[count], sizeof(int), (size_t) 1, fpt1);
    printf("%lf\t%lf\t%lf\t%d\n", ca[count], cb[count], radius[count],
	   n_ppt[count]);
  }
  fclose(fpt1);

/* Count the total number of particles */  
  nparticles = 0;
  for (count = 0; count <= 100; count++) {
    for (ncount = 0; ncount < n_ppt[count]; ncount++) {
      rad[nparticles] = radius[count];
      compb[nparticles] = ca[count];
      compc[nparticles] = cb[count];
      nparticles++;
    }
  }

/* Remember that the radius was calculated using classical nucleation theory.
 * We increase the size of the nucleus by 15% of the calculated value */
 
  for (i = 0; i < nparticles; i++){
  rad[i] = rad[i] + 0.15 * rad[i];
  printf("%d\t%lf\t%lf\t%lf\n",  i, compb[i],
  compc[i], rad[i]);
  }
  for (i = 0; i < nparticles; i++){
  drad[i] = ceil(rad[i]);
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      occupancy[i][j] = 0;
    }
  }
 printf("%d\n", nparticles);

  l = 0;
  iter = 0;
  
/* Here is the main portion of the program. Here we generate random numbers
 * using Numeriacal Recipes routine "gsl_uniform_sample.c". 
 * We choose the coordinates of the
 * center of the nuclei using the random numbers. We also apply periodic
 * boundary conditions. We scan the box to see which positions are occupied by
 * the nuclei. We use 4x(diameter of a particle) as the width of the domain
 * in which the occupancy is assigned as 1 (i.e. occupied by the particle). 
 * We put new particles in the unoccupied sites (occupancy = 0) 
 * to ensure that there is no overlap. */

 while (l < nparticles) {
 x[l] = (int) (nx-1) * gsl_uniform_sample(&SEED);
 y[l] = (int) (ny-1) * gsl_uniform_sample(&SEED);

  tmp1 = x[l] - 2 * drad[l];
  p = x[l];
  if(tmp1 < 0){
  tmp1 = tmp1 + nx;
  p  = tmp1 + 2 * drad[l];
  }
  tmp2 = x[l] + 2 * drad[l];
  q = x[l];
  if(tmp2 > (nx-1)){
  tmp2 = tmp2 - nx;
  q  = tmp2 - 2 * drad[l];
  }
  tmp3 = y[l] - 2 * drad[l];
  r = y[l];
  if(tmp3 < 0){
  tmp3 = tmp3 + ny;
  r  = tmp3 + 2 * drad[l];
  }
  tmp4 = y[l] + 2 * drad[l];
  s = y[l];
  if(tmp4 > (ny-1)){
  tmp4 = tmp4 - ny;
  s  = tmp4 - 2 * drad[l];
  } 
    
    state = TRUE;
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
      if (  ((i-x[l]) * (i-x[l]) +  (j-y[l]) * (j-y[l])) <= 4 * drad[l] * drad[l])
           {
	  if (occupancy[i][j] != 0) {
	    state = FALSE;
	  }
	}
      if (  ((i-p) * (i-p) +  (j-r) * (j-r)) <= 4 * drad[l] * drad[l])
           {
	  if (occupancy[i][j] != 0) {
	    state = FALSE;
	  }
	}
      if (  ((i-q) * (i-q) +  (j-r) * (j-r)) <= 4 * drad[l] * drad[l])
           {
	  if (occupancy[i][j] != 0) {
	    state = FALSE;
	  }
	}
      if (  ((i-p) * (i-p) +  (j-s) * (j-s)) <= 4 * drad[l] * drad[l])
           {
	  if (occupancy[i][j] != 0) {
	    state = FALSE;
	  }
	}
      if (  ((i-q) * (i-q) +  (j-s) * (j-s)) <= 4 * drad[l] * drad[l])
           {
	  if (occupancy[i][j] != 0) {
	    state = FALSE;
	  }
	}
      }
    }
    if (state == TRUE) {
      for (i = 0; i < nx; i++) {
	for (j = 0; j < ny; j++) {
      if (  ((i-x[l]) * (i-x[l]) +  (j-y[l]) * (j-y[l])) <= 4 * drad[l] * drad[l])
           {
	    occupancy[i][j] = 1;
	  }
      if (  ((i-p) * (i-p) +  (j-r) * (j-r)) <= 4 * drad[l] * drad[l])
           {
	    occupancy[i][j] = 1;
	  }
      if (  ((i-q) * (i-q) +  (j-r) * (j-r)) <= 4 * drad[l] * drad[l])
           {
	    occupancy[i][j] = 1;
	  }
      if (  ((i-p) * (i-p) +  (j-s) * (j-s)) <= 4 * drad[l] * drad[l])
           {
	    occupancy[i][j] = 1;
	  }
      if (  ((i-q) * (i-q) +  (j-s) * (j-s)) <= 4 * drad[l] * drad[l])
           {
	    occupancy[i][j] = 1;
	  }
      if (  ((i-x[l]) * (i-x[l]) +  (j-y[l]) * (j-y[l])) <= drad[l] * drad[l])
           {
	    comp_b[j + i * ny][Re] = compb[l];
	    comp_c[j + i * ny][Re] = compc[l];
	  }
      if (  ((i-p) * (i-p) +  (j-r) * (j-r)) <= drad[l] * drad[l])
           {
	    comp_b[j + i * ny][Re] = compb[l];
	    comp_c[j + i * ny][Re] = compc[l];
	  }
      if (  ((i-q) * (i-q) +  (j-r) * (j-r)) <= drad[l] * drad[l])
           {
	    comp_b[j + i * ny][Re] = compb[l];
	    comp_c[j + i * ny][Re] = compc[l];
	  }
      if (  ((i-p) * (i-p) +  (j-s) * (j-s)) <= drad[l] * drad[l])
           {
	    comp_b[j + i * ny][Re] = compb[l];
	    comp_c[j + i * ny][Re] = compc[l];
	  }
      if (  ((i-q) * (i-q) +  (j-s) * (j-s)) <= drad[l] * drad[l])
           {
	    comp_b[j + i * ny][Re] = compb[l];
	    comp_c[j + i * ny][Re] = compc[l];
	  }
      }}
    }
    iter++;
    printf("Iter=%d\n", iter);
    if (state == TRUE) {
      l++;
      printf("%d\n", l);
    }
  } 
/* Here we calculate the total composition of the particles that have been
 * introduced in the box */  
  total_b = 0.0;
  total_c = 0.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      total_b = total_b + comp_b[j + i * ny][Re];
      total_c = total_c + comp_c[j + i * ny][Re];
    }
  }
  total_b = total_b / (double) (nx * ny);
  total_c = total_c / (double) (nx * ny);
  printf("%le\t %le\n", total_b, total_c);
  
/* We now calculate the total number of sites in which the particles are
 * absent */  
  o_sites = 1;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      if ((comp_b[j + i * ny][Re] < DBL_EPSILON) &&
       (comp_c[j + i * ny][Re] < DBL_EPSILON)){
    o_sites = o_sites + 1;
    }
    }
  }

/* Oce we know the total number of vacant sites we assign the supersaturation
 * at those points so the total composition is equal to the starting alloy
 * composition. */
  
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      if ((comp_b[j + i * ny][Re] < DBL_EPSILON) &&
       (comp_c[j + i * ny][Re] < DBL_EPSILON)){
	comp_b[j + i * ny][Re] = ((cb_matrix - total_b)*((double)(nx*ny)))/(double)
                              o_sites;
	comp_c[j + i * ny][Re] = ((cc_matrix - total_c)*((double)(nx*ny)))/(double)
                              o_sites;
    }
    }
  }

  total_b = 0.0;
  total_c = 0.0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      total_b = total_b + comp_b[j + i * ny][Re];
      total_c = total_c + comp_c[j + i * ny][Re];
    }
  }
  total_b = total_b / (double) (nx * ny);
  total_c = total_c / (double) (nx * ny);
  printf("%le\t %le\n", total_b, total_c);
  fpt2 = fopen(fn2, "w");
  
  fwrite(&comp_b[0][0], sizeof(double), 2 * nx * ny, fpt2);
  fwrite(&comp_c[0][0], sizeof(double), 2 * nx * ny, fpt2);
  fclose(fpt2);
  fftw_free(comp_b);
  fftw_free(comp_c);
  free_imatrix(occupancy,0,nx,0,ny);
  return(0);
}
