#ifndef PW_SP_H_
#define PW_SP_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "osu_mksa.h"

typedef struct planewave_sp{
  double lambda0;       // wavelength in vacuum 
  double n;             // refractive index of incident field
  double complex es;    // s-polarization amplitude include initial phase
  double complex ep;    // p-polarization amplitude include initial phase
  double theta;         // incident angle (0 <= theta < pi/2)
  double phi;           // parameter for incident direction (-pi < phi <= pi) 
                        // if theta=0, phi=0 then ep=E_x, es=E_y. 

  double k0;            // wave number in vacuum, same as the angular frequency in this system of units (OSU)
  double vk[3];         // wave vector, vk[0]=k_x, vk[1]=k_y, vk[2]=k_z
  
  double complex e[3];  // electric field vector, e[0]=E_x,e[1]=E_y,e[2]=E_z. 
  double complex h[3];  // magnetic field vector, h[0]=H_x,h[1]=H_y,h[2]=H_z.
}Pwsp;

void read_data_pwsp(char *filename,Pwsp *pw); // read datafile
void print_data_pwsp(Pwsp *pw);               // print data
void print_data_pwsp_mksa(Pwsp *pw);          // print data in MKSA system of units
void setup_pwsp(Pwsp *pw);                    // calculate coefficients 

void EH_pwsp(double complex *e,double complex *h,double *r,Pwsp *pw); 
// inputs
// r[0]=x, r[1]=y, r[2]=z, pointer of Pwsp object
// outputs
// e[0]=E_x, e[1]=E_y, e[2]=E_z, h[0]=H_x, h[1]=H_y, h[2]=H_z 

#endif
