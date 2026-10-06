#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include"gsl_support.h"
#include"gsl_support.c"
int
main (void)
{
 int nx, ny;
 
 FILE *fpread, *fpout;
 char fr[100],fw[100];
 double total;
 char new[100];
 fftw_complex *comp, *eta1, *eta2, *eta3;
 double **etamap;
 int initial, final, steps, count;
 int i, j;
 
 printf("Enter the system size in nx and ny \n");
 scanf("%d%d",&nx,&ny);
 printf ("Initial count (initial)\n");
 scanf ("%d", &initial);
 printf ("Final count (final) \n");
 scanf ("%d", &final);
 printf ("Steps (steps) \n");
 scanf ("%d", &steps);
 comp = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta1 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta2 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta3 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 etamap = dmatrix(0,nx,0,ny);
 for (count = initial; count <= final; count = count + steps) {
 sprintf (fr, "comp.%06d",count);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&comp[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);
 total = 0.0;
 for (i = 0; i < nx; i++) {
  for (j = 0; j < ny; j++) {
    total += comp[j + i * ny][0];
  }
 }
 total = total/(double) (nx * ny); 
 printf("total = %le\n", total);
 sprintf (fw, "comp%06d.profile",count);
 fpout= fopen(fw,"w");
  for (j = 0; j < ny; j++){
   for (i = 0; i < nx; i++){
    fprintf(fpout,"%d\t%d\t%le\n", j, i, comp[i+j*ny][0]);
	}
	fprintf(fpout,"\n");
 }

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
 sprintf (fw, "etamap.%06d", count);
 fpout = fopen (fw, "w");
 for ( i = 0; i < nx; ++i){
  for ( j = 0; j < ny; ++j){
//  etamap[i][j] = (4.0 * eta1[j + i * ny][0] + 2.0 * eta2[j + i * ny][0] + 
//                 eta3[j + i * ny][0] + 6.0*comp[j + i * ny][0]);
  etamap[i][j] = ( eta1[j + i * ny][0] * eta1[j + i * ny][0] - 
                   eta2[j + i * ny][0] * eta2[j + i * ny][0] 
                   + 2.0 * eta2[j + i * ny][0] * eta3[j + i * ny][0]);
  fprintf(fpout, " %lf ", etamap[i][j]);
  }
  fprintf(fpout, "\n");
 }
 fclose(fpout);		 
 }
 fftw_free (comp);
 fftw_free (eta1);
 fftw_free (eta2);
 fftw_free (eta3);
 free_dmatrix(etamap,0,nx,0,ny);
 return(0);
}
