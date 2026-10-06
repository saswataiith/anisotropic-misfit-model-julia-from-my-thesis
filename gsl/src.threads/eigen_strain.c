/*
 * =====================================================================================
 * 
 *        Filename:  eigen_strain.c
 * 
 *     Description:  
 * 
 *         Version:  1.0
 *         Created:  10/12/2006 07:15:20 PM IST
 *        Revision:  none
 *        Compiler:  gcc
 * 
 *          Author:   (), 
 *         Company:  
 * 
 * =====================================================================================
 */

#include "binary.h"

void Calculate_EigenStrain ()
{
int i, j;

for(i = 0; i < 3; ++i){
 for(j = 0; j < 3; ++j){
  if(i == j && (i+j)<=2){
  eigen_strain_0[i][j] = epsc;
	} else {
  eigen_strain_0[i][j] = 0.0;
	}
}} 
	   
eigen_strain_1[0][0] = eps_eta;
eigen_strain_1[0][1] = 0.0;
eigen_strain_1[1][0] = 0.0;
eigen_strain_1[1][1] = eps_eta * tetra;


eigen_strain_2[0][0] = eps_eta * (( 1.0 + 3.0 * tetra)/4.0);
eigen_strain_2[0][1] = eps_eta * ((sqrt(3.0) * (1.0 - tetra))/4.0);
eigen_strain_2[1][0] = eps_eta * ((sqrt(3.0) * (1.0 - tetra))/4.0) ;
eigen_strain_2[1][1] = eps_eta * ((3.0 + tetra)/4.0);


eigen_strain_3[0][0] = eps_eta * (( 1.0 + 3.0 * tetra)/4.0);
eigen_strain_3[0][1] = eps_eta * -1.0 * ((sqrt(3.0) * (1.0 - tetra))/4.0);
eigen_strain_3[1][0] = eps_eta * -1.0 * ((sqrt(3.0) * (1.0 - tetra))/4.0) ;
eigen_strain_3[1][1] = eps_eta * ((3.0 + tetra)/4.0);
}
