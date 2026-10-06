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
 int numsteps;
 char new[100];
 fftw_complex *comp, *eta1, *eta2, *eta3;

 printf("Enter the system size in nx and ny \n");
 scanf("%d%d",&nx,&ny);
 fputs("Enter the extension (numsteps)\n",stdout);
 scanf ("%d", &numsteps);
 comp = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta1 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta2 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 eta3 = fftw_malloc (nx * ny * sizeof (fftw_complex));
 int i, j;

 sprintf (fr, "comp.%06d",numsteps);
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
 sprintf (fw, "comp%06d.profile",numsteps);
 fpout= fopen(fw,"w");
 for (i = 910; i < nx; i++){
  for (j = 0; j < ny; j++){
   if (j == 420)
   fprintf(fpout,"%d\t%le\n ", i, comp[j+i*ny][0]);
	}
 }

sprintf (fr, "eta1.%06d",numsteps);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&eta1[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);
 sprintf (fw, "eta1%06d.profile",numsteps);
 fpout= fopen(fw,"w");
 for (i = 910; i < nx; i++){
  for (j = 0; j < ny; j++){
   if (j == 420)
   fprintf(fpout,"%d\t%le\n ", i, eta1[j+i*ny][0]);
	}
 }
// for (i = 0; i < nx; i++){
//  for (j = 0; j < ny; j++){
//   fprintf(fpout," %le ",eta1[j+i*ny][0]);
//	}
//	fprintf(fpout,"\n");
//}
  
  fclose(fpout);
 
 sprintf (fr, "eta2.%06d",numsteps);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&eta2[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);
 sprintf (fw, "eta2%06d.profile",numsteps);
 fpout= fopen(fw,"w");
 for (i = 910; i < nx; i++){
  for (j = 0; j < ny; j++){
   if (j == 420)
   fprintf(fpout,"%d\t%le\n ", i, eta2[j+i*ny][0]);
	}
 }
// for (i = 0; i < nx; i++){
//  for (j = 0; j < ny; j++){
//   fprintf(fpout," %le ",eta2[j+i*ny][0]);
//	}
//	fprintf(fpout,"\n");
// }
  
  fclose(fpout);

 sprintf (fr, "eta3.%06d",numsteps);
 printf ("%s\n", fr);
 fpread = fopen (fr, "r");
 fread (&eta3[0][0], sizeof (double), 2 * nx * ny, fpread);
 fclose (fpread);
 sprintf (fw, "eta3%06d.profile",numsteps);
 fpout= fopen(fw,"w");
 for (i = 910; i < nx; i++){
  for (j = 0; j < ny; j++){
   if (j == 420)
   fprintf(fpout,"%d\t%le\n ", i, eta3[j+i*ny][0]);
	}
 }
// for (i = 0; i < nx; i++){
//  for (j = 0; j < ny; j++){
//   fprintf(fpout," %le ",eta3[j+i*ny][0]);
//	}
//	fprintf(fpout,"\n");
// }
  fclose(fpout);

 fftw_free (comp);
 fftw_free (eta1);
 fftw_free (eta2);
 fftw_free (eta3);
 return(0);
}
