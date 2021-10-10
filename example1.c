#include "emf_multilayer.h"

int main(int argc,char *argv[])
{
  if(argc!=3){
    printf("Usage: %s incident_plane_wave_datafile_name multilayer_datafile_name\n",argv[0]);
    exit(0);
  }
  
  Mlpw ml;
  double complex e[3],h[3];
  double r[3],R[2],T[2];
  int type,lid;
  
  read_data_mlpw(argv[1],argv[2],&ml); // read data and memory allocation 
  print_data_mlpw(&ml);                // print data
  //print_data_mlpw_mksa(&ml);         // print data in MKSA system of units
  setup_mlpw(&ml);                     // calculation of coefficients

  printf("power reflection and transmission coefficients\n");
  power_RT_coefficient(R,T,&ml);
  printf("p-polarization Rp=% 15.14f,\tTp=% 15.14f,\tRp+Tp=%15.14f \n",R[0],T[0],R[0]+T[0]);
  printf("s-polarization Rs=% 15.14f,\tTs=% 15.14f,\tRs+Ts=%15.14f \n",R[1],T[1],R[1]+T[1]);
 
  r[0]=0.0; 
  r[1]=0.0; 
  r[2]=0.0; // r[0]=x, r[1]=y, r[2]=z
  type=0;   // select layer_id_type0() (On the boundary, it returns the layer id of the smaller side.)
  lid=EH_t_mlpw(e,h,r,type,&ml); // total field
  printf("\nelectromagnetic field at r=(%g, %g, %g), layer id=%d\n",r[0],r[1],r[2],lid);
  printf("E_x=% 15.14e %+15.14e I\n",creal(e[0]),cimag(e[0]));
  printf("E_y=% 15.14e %+15.14e I\n",creal(e[1]),cimag(e[1]));
  printf("E_z=% 15.14e %+15.14e I\n",creal(e[2]),cimag(e[2]));
  printf("H_x=% 15.14e %+15.14e I\n",creal(h[0]),cimag(h[0]));
  printf("H_y=% 15.14e %+15.14e I\n",creal(h[1]),cimag(h[1]));
  printf("H_z=% 15.14e %+15.14e I\n",creal(h[2]),cimag(h[2]));
  
  type=1;  // select layer_id_type1() (On the boundary, it returns the layer id of the larger side.)
  lid=EH_t_mlpw(e,h,r,type,&ml); // total field
  printf("electromagnetic field at r=(%g, %g, %g), layer id=%d\n",r[0],r[1],r[2],lid);
  printf("E_x=% 15.14e %+15.14e I\n",creal(e[0]),cimag(e[0]));
  printf("E_y=% 15.14e %+15.14e I\n",creal(e[1]),cimag(e[1]));
  printf("E_z=% 15.14e %+15.14e I\n",creal(e[2]),cimag(e[2]));
  printf("H_x=% 15.14e %+15.14e I\n",creal(h[0]),cimag(h[0]));
  printf("H_y=% 15.14e %+15.14e I\n",creal(h[1]),cimag(h[1]));
  printf("H_z=% 15.14e %+15.14e I\n",creal(h[2]),cimag(h[2]));
  
  mfree_mlpw(&ml);                      // free allocated memory
  return 0;
}
