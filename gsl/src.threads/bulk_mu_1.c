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
		double cb, eta1, eta2, eta3;

		for (i = 0; i < nx; ++i) {
		 for (j = 0; j < ny; ++j) {
				
			cb = dfdc[j + i * ny][Re];
			
			eta1 = dfdeta_1[j + i * ny][Re];

			eta2 = dfdeta_2[j + i * ny][Re];

			eta3 = dfdeta_3[j + i * ny][Re];

dfdc[j + i * ny][Re] = 905.387 * (cb - 0.125) - (211.611/6.0) * 
				              (eta1 * eta1 + eta2 * eta2 + eta3 * eta3);

dfdeta_1[j + i * ny][Re] = (211.611/6.0) * (0.383 - cb) * (2.0 * eta1) 
				- (165.031/3.0) * ( eta2 * eta3) 
				+ (30.0/4.0) * (4.0 * eta1 * eta1 * eta1) 
				+ (5.071/4.0) * 4.0 * eta1 * (eta1 * eta1 + eta2 * eta2 + eta3 * eta3);
				                 
				
      
dfdeta_2[j + i * ny][Re] = (211.611/6.0) * (0.383 - cb) * (2.0 * eta2) 
				- (165.031/3.0) * ( eta1 * eta3) 
				+ (30.0/4.0) * (4.0 * eta2 * eta2 * eta2) 
				+ (5.071/4.0) * 4.0 * eta2 * (eta1 * eta1 + eta2 * eta2 + eta3 * eta3);
				
dfdeta_3[j + i * ny][Re] = (211.611/6.0) * (0.383 - cb) * (2.0 * eta3) 
				- (165.031/3.0) * ( eta1 * eta2) 
				+ (30.0/4.0) * (4.0 * eta3 * eta3 * eta3) 
				+ (5.071/4.0) * 4.0 * eta3 * (eta1 * eta1 + eta2 * eta2 + eta3 * eta3);
				
		
dfdc[j + i * ny][Im] = 0.0;
dfdeta_1[j + i * ny][Im] = 0.0;
dfdeta_2[j + i * ny][Im] = 0.0;
dfdeta_3[j + i * ny][Im] = 0.0;
		 }}
}
