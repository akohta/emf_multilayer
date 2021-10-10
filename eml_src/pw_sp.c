#include "pw_sp.h"

void read_data_pwsp(char *fname,Pwsp *pw)
{
  FILE *fp;
  char buf[256]="";
  double td1,td2;
  
  if((fp=fopen(fname,"rt"))==NULL){ 
    printf("failed to read %s. Exit...\n",fname);
    exit(1);
  }  
  fgets(buf,256,fp);  fgets(buf,256,fp);
  
  fscanf(fp,"%lf",&td1);  pw->lambda0=td1;
  fscanf(fp,"%lf",&td1);  pw->n=td1;    
  fscanf(fp,"%lf",&td1);
  fscanf(fp,"%lf",&td2);  pw->ep=td1+I*td2;
  fscanf(fp,"%lf",&td1);
  fscanf(fp,"%lf",&td2);  pw->es=td1+I*td2; 
  fscanf(fp,"%lf",&td1);  pw->theta=td1;
  fscanf(fp,"%lf",&td1);  pw->phi=td1;
  
  fclose(fp);
}

void print_data_pwsp(Pwsp *pw)
{
  printf("-- incident plane wave --\n");
  printf("wavelength in vacuum                         : %-15.14g\n",pw->lambda0);
  printf("refractive index                             : %-15.14g\n",pw->n);
  printf("p-polarization amplitude                     : %-7.6e %+7.6e I\n",creal(pw->ep),cimag(pw->ep)); 
  printf("s-polarization amplitude                     : %-7.6e %+7.6e I\n",creal(pw->es),cimag(pw->es));
  printf("incident angle (theta)                 [rad] : %-15.14g\n",pw->theta);
  printf("parameter for incident direction (phi) [rad] : %-15.14g\n",pw->phi); 
  //printf("incident wave power                          : %-15.14g\n",0.5*pw->n*creal(pw->es*conj(pw->es)+pw->ep*conj(pw->ep)));
  
  printf("\n");
}

void print_data_pwsp_mksa(Pwsp *pw)
{
  printf("-- incident plane wave --\n");
  printf("wavelength in vacuum                     [m] : %-15.14g\n",OSUtoMKSA_length(pw->lambda0));
  printf("refractive index                             : %-15.14g\n",pw->n);
  printf("p-polarization amplitude               [V/m] : %-7.6e %+7.6e I\n",creal(OSUtoMKSA_ElectricField(pw->ep)),cimag(OSUtoMKSA_ElectricField(pw->ep))); 
  printf("s-polarization amplitude               [V/m] : %-7.6e %+7.6e I\n",creal(OSUtoMKSA_ElectricField(pw->es)),cimag(OSUtoMKSA_ElectricField(pw->es)));
  printf("incident angle (theta)                 [rad] : %-15.14g\n",pw->theta);
  printf("parameter for incident direction (phi) [rad] : %-15.14g\n",pw->phi); 
  printf("incident wave power                  [W/m^2] : %-15.14g\n",OSUtoMKSA_Power_per_unit_area(0.5*pw->n*creal(pw->es*conj(pw->es)+pw->ep*conj(pw->ep))));
  
  printf("\n");
}

void setup_pwsp(Pwsp *pw)
{
  double complex te[3],th[3];
  double ct,st,cp,sp,rm[9],k;
  int i,j;
  
  ct=cos(pw->theta);
  st=sin(pw->theta);
  cp=cos(pw->phi);
  sp=sin(pw->phi);

  pw->k0=2.0*M_PI/pw->lambda0;
  k=pw->k0*pw->n;
  pw->vk[0]=k*st*cp;
  pw->vk[1]=k*st*sp;
  pw->vk[2]=k*ct;
  
  for(i=0;i<3;i++){
    pw->e[i]=0.0;
    pw->h[i]=0.0;
  }
  te[0]=pw->ep;
  te[1]=pw->es;
  te[2]=0.0;
  th[0]=-te[1]*pw->n;
  th[1]= te[0]*pw->n;
  th[2]=0.0;
  
  // rotation matrix arround y-axis
  rm[0*3+0]= ct;  rm[0*3+1]=0.0;  rm[0*3+2]= st;
  rm[1*3+0]=0.0;  rm[1*3+1]=1.0;  rm[1*3+2]=0.0;
  rm[2*3+0]=-st;  rm[2*3+1]=0.0;  rm[2*3+2]= ct;

  for(j=0;j<3;j++){
    for(i=0;i<3;i++){
      pw->e[j]+=rm[j*3+i]*te[i];
      pw->h[j]+=rm[j*3+i]*th[i];
    }
  }
  
  for(i=0;i<3;i++){
    te[i]=pw->e[i];
    th[i]=pw->h[i];
    pw->e[i]=0.0;
    pw->h[i]=0.0;
  }
  // rotation matrix arround z-axis
  rm[0*3+0]= cp;  rm[0*3+1]=-sp;  rm[0*3+2]=0.0;
  rm[1*3+0]= sp;  rm[1*3+1]= cp;  rm[1*3+2]=0.0;
  rm[2*3+0]=0.0;  rm[2*3+1]=0.0;  rm[2*3+2]=1.0;
  
  for(j=0;j<3;j++){
    for(i=0;i<3;i++){
      pw->e[j]+=rm[j*3+i]*te[i];
      pw->h[j]+=rm[j*3+i]*th[i];
    }
  }
}

void EH_pwsp(double complex *e,double complex *h,double *r,Pwsp *pw)
{
  double complex ce;
  double arg;
  int i;
  
  arg=pw->vk[0]*r[0]+pw->vk[1]*r[1]+pw->vk[2]*r[2];
  ce=cos(arg)+I*sin(arg);
  
  for(i=0;i<3;i++){
    e[i]=pw->e[i]*ce;
    h[i]=pw->h[i]*ce;
  }
}
