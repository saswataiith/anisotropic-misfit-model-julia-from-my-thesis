/*
 * ============================================================================
 * 
 *        Filename:  bulk_mu.c
 * 
 *     Description:  calculate dfdc_b, dfdc_c, dfdeta_i (i=1,2,3,4,5,6)
 * 
 *         Version:  1.0
 *         Created:  07/14/2006 03:20:21 PM IST
 *        Revision:  none
 *        Compiler:  gcc
 *          Author:  saswata  
 *         Company:  cmsg
 * 
 * ============================================================================
 */

#include "binary.h"

void Calculate_bulk_mu () 
{
		int i, j;	
		double c, eta[3];

		for (i = 0; i < nx; ++i) {
		 for (j = 0; j < ny; ++j) {
				
			c = dfdc[j + i * ny][Re];
			
			eta[0] = dfdeta_1[j + i * ny][Re];

			eta[1] = dfdeta_2[j + i * ny][Re];

			eta[2] = dfdeta_3[j + i * ny][Re];

dfdc[j + i * ny][Re] = 
100.0*c-0.15E2-0.12E2*eta[0]*eta[0]-0.12E2*eta[1]*eta[1]-0.12E2*eta
[2]*eta[2];

				              

dfdeta_1[j + i * ny][Re] =  
0.24E2*(0.245-c)*eta[0]-0.22E2*eta[0]*eta[0]*eta[0]+0.278E2*pow(eta
[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2],2.0)*eta[0];
 
      
dfdeta_2[j + i * ny][Re] = 
0.24E2*(0.245-c)*eta[1]-0.22E2*eta[1]*eta[1]*eta[1]+0.278E2*pow(eta
[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2],2.0)*eta[1];

				
dfdeta_3[j + i * ny][Re] = 
0.24E2*(0.245-c)*eta[2]-0.22E2*eta[2]*eta[2]*eta[2]+0.278E2*pow(eta
[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2],2.0)*eta[2];

dfdc[j + i * ny][Im] = 0.0;
dfdeta_1[j + i * ny][Im] = 0.0;
dfdeta_2[j + i * ny][Im] = 0.0;
dfdeta_3[j + i * ny][Im] = 0.0;
		 }}
}
