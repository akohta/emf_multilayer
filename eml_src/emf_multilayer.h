/*
 * emf_multilayer.h
 *
 *  Created on: Feb 28, 2019
 *      Author: ohta
 */

#ifndef EMF_MULTILAYER_H_
#define EMF_MULTILAYER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "my_utils.h"
#include "osu_mksa.h"
#include "physical_constant.h"
#include "pw_sp.h" // incident plane wave 


#define AT_COEF 1.0e-100 // amplutude attenuation limit for avoid underflow and overflow
                         // If the amplitude falls below this value, it is considered to have been attenuated to 0. 

typedef struct multi_layer_object{
  Pwsp pw ;                 // incident plane wave
  
  int No;                   // number of layers except incident field layer. total number of layers is No+1
  int Ns;                   // the last layer id 
  
  double *zj;               // boundary position 
  double complex *nj;       // refractive index of each layer

  double complex *Zj;       // wave impedance of each layer
  double complex *kzj;      // z-component of wave vector in each layer
  double complex *cos_tj;   // kzj = k * cos_tj

  double complex **eC;      // coefficients for E-wave (p-polarization)
  double complex **mC;      // coefficients for H-wave (s-polarization)
  double Rs,Rp,Ts,Tp;       // power reflection (incident field layer) and transmission (last layer) coefficients 

  double complex **Ep,**Em; // electric field amplitude vector in each layer.
                            // Ep is the component propagating in the positive direction.
                            // Em is the component propagating in the negative direction.
  double complex **Hp,**Hm; // magnetic field amplitude vector in each layer.
                            // Hp is the component propagating in the positive direction.
                            // Hm is the component propagating in the negative direction.
}Mlpw;


void read_data_mlpw(char *fname_ipw,char *fname_ml,Mlpw *ml); // read data and allocate memory, fname_ipw:filename for incident plane wave datafile, fname_ml:filename for multilayer datafile.
void print_data_mlpw(Mlpw *ml);       // print data 
void print_data_mlpw_mksa(Mlpw *ml);  // print data in MKSA system of units 
void setup_mlpw(Mlpw *ml);            // calculate coefficients
void mfree_mlpw(Mlpw *ml);            // free allocated memory 
int layer_id_type0(double z,Mlpw *ml);// return the layer id. On the boundary, it returns the layer id of the smaller side. 
int layer_id_type1(double z,Mlpw *ml);// return the layer id. On the boundary, it returns the layer id of the larger side. 

int EH_t_mlpw(double complex *e,double complex *h,double *r,int type,Mlpw *ml); // total field
int EH_i_mlpw(double complex *e,double complex *h,double *r,int type,Mlpw *ml); // incident field
int EH_p_mlpw(double complex *e,double complex *h,double *r,int type,Mlpw *ml); // the component propagating in the positive z-axis direction 
int EH_m_mlpw(double complex *e,double complex *h,double *r,int type,Mlpw *ml); // the component propagating in the negative z-axis direction 
// outputs
// e[0]=E_x, e[1]=E_y, e[2]=E_z (electric field), h[0]=H_x, h[1]=H_y, h[2]=H_z (magnetic field), return layer id ( 0:incident field layer ~ Ns:last layer ).
// inputs
// r[0]=x, r[1]=y, r[2]=z, type=0:select layer_id_type0(), type!=0:select layer_id_type1(), pointer of Mlpw object.

void power_RT_coefficient(double *R,double *T,Mlpw *ml); // power reflection and transmission coefficients 
// outputs
// R:power reflection coefficients,   R[0] for p-polarization, R[1] for s-polarization
// T:power transmission coefficients, T[0] for p-polarization, T[1] for s-polarization 

int grad_EH_t_mlpw(double complex *e,double complex *grad_e,double complex *h,double complex *grad_h,double *r,int type,Mlpw *ml);
// outputs
// grad_e[0]=dE_x/dx, grad_e[1]=dE_x/dy, grad_e[2]=dE_x/dz,
// grad_e[3]=dE_y/dx, grad_e[4]=dE_y/dy, grad_e[5]=dE_y/dz,
// grad_e[6]=dE_z/dx, grad_e[7]=dE_z/dy, grad_e[8]=dE_z/dz.
// grad_h[0]=dH_x/dx, grad_h[1]=dH_x/dy, grad_h[2]=dH_x/dz,
// grad_h[3]=dH_y/dx, grad_h[4]=dH_y/dy, grad_h[5]=dH_y/dz,
// grad_h[6]=dH_z/dx, grad_h[7]=dH_z/dy, grad_h[8]=dH_z/dz.
// others are the same as EH_t_mlpw().

#endif 
