#include"binary.h"
void
Output_Conf ()
{
 FILE *fpt1, *fpt2, *fpt3, *fpt4;

// int i, j;

 char fn1[100], fn2[100], fn3[100], fn4[100];

 
 sprintf (fn1, "comp.%06d", count + initcount);
 sprintf (fn2, "eta1.%06d", count + initcount);
 sprintf (fn3, "eta2.%06d", count + initcount);
 sprintf (fn4, "eta3.%06d", count + initcount);



 fpt1 = fopen (fn1, "w");
 fwrite (&dfdc[0][0], sizeof(double), 2 * nx * ny, fpt1);
 fclose (fpt1);
 
 fpt2 = fopen(fn2, "w");
 fwrite (&dfdeta_1[0][0], sizeof(double), 2 * nx * ny, fpt2);
 fclose(fpt2);
 fpt3 = fopen(fn3, "w");
 fwrite (&dfdeta_2[0][0], sizeof(double), 2 * nx * ny, fpt3);
 fclose(fpt3);
 fpt4 = fopen(fn4, "w");
 fwrite (&dfdeta_3[0][0], sizeof(double), 2 * nx * ny, fpt4);
 fclose(fpt4);

}
