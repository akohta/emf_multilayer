// Analysis example of power reflection and transmission coefficients, intensity distributions.
#include "emf_multilayer.h"

int main(int argc,char *argv[]) 
{
  if(argc!=3){
    printf("Usage: %s incident_plane_wave_datafile_name multilayer_datafile_name\n",argv[0]);
    exit(0);
  }
  
  Mlpw ml;
  FILE *fp1,*fp2;
  double complex e[3],h[3];
  double angle_min,angle_max,dt,cdr,angle,R[2],T[2],z_min,z_max,dz,r[3];
  int sn,i,j,type;

  // -- settings --
  sn=400;         // sampling number 
  angle_min=0.0;  //
  angle_max=90.0; // incident angle range, angle_min <= theta < angle_max [degree]
  z_min=-5.0;
  z_max= 5.0;     // range for intensity analysis, z_min <= z <= z_max 
  // --------------
  
  dt=(angle_max-angle_min)/(double)(sn);
  cdr=M_PI/180.0; 
  dz=(z_max-z_min)/(double)(sn-1);
  r[0]=0.0;
  r[1]=0.0; // set x y position for intensity analysis
  type=0;   // select layer_id_type0() (On the boundary, it returns the layer id of the smaller side.)
  
  read_data_mlpw(argv[1],argv[2],&ml); // read data and memory allocation
  print_data_mlpw(&ml);                // print data
  //print_data_mlpw_mksa(&ml);         // print data in MKSA system of units
  
  if((fp1=fopen("power_coefficients.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# incident_angle_in_degree Rp Rs Tp Ts");
  if((fp2=fopen("intensity_distributions.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp2,"%s\n","# incident_angle_in_degree z |E|^2 |H|^2");
  
  for(i=0;i<sn;i++){
    angle=angle_min+(double)i*dt;       // incident angle in degree 
    ml.pw.theta=angle*cdr;              // set incident angle in radian
    setup_mlpw(&ml);                    // calculate coefficients
    power_RT_coefficient(R,T,&ml);      // take power coefficients
    fprintf(fp1,"%15.14f %15.14f %15.14f %15.14f %15.14f\n",angle,R[0],R[1],T[0],T[1]);
    
    for(j=0;j<sn;j++){
      r[2]=z_min+(double)j*dz;
      EH_t_mlpw(e,h,r,type,&ml); // total field
      fprintf(fp2,"%15.14f %15.14f %15.14f %15.14f\n",angle,r[2],creal(e[0]*conj(e[0])+e[1]*conj(e[1])+e[2]*conj(e[2])),creal(h[0]*conj(h[0])+h[1]*conj(h[1])+h[2]*conj(h[2])));
    }
    fprintf(fp2,"\n");
  }  
 
  fclose(fp1);
  fclose(fp2);

  printf("Done.\n");

  mfree_mlpw(&ml);
  return 0;
}
