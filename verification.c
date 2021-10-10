// verification code using boundary condition and Maxwell's equations
#include "emf_multilayer.h"
#include "time.h"

double rand_num_d(double min,double max,int resolution); 
void verification_1(Mlpw *ml); // using boundary condition 
void verification_2(Mlpw *ml); // using Maxwell's equations

int main(int argc,char *argv[])
{
  if(argc!=3 && argc!=4){
    printf("Usage: %s incident_plane_wave_datafile_name multilayer_datafile_name [incident_angle](optional)\n",argv[0]);
    exit(0);
  }
  
  Mlpw ml;
  
  read_data_mlpw(argv[1],argv[2],&ml);    // read data and memory allocation 
  if(argc==4) ml.pw.theta=atof(argv[3]); // set incident angle 
  print_data_mlpw(&ml);                   // print data
  //print_data_mlpw_mksa(&ml);            // print data in MKSA system of units
  setup_mlpw(&ml);                        // calculation of coefficients

  verification_1(&ml);
  verification_2(&ml);
  
  mfree_mlpw(&ml);                      // free allocated memory
  return 0;
}

double rand_num_d(double min,double max,int res)
{
  double dl=(max-min)/(double)res;
  static int s=0;
  if(s==0){
    srand((unsigned)time(NULL)); // random seed 
    s++;
  }
  return min+dl*(double)(rand()%(res+1));
}

void verification_1(Mlpw *ml)
{
  double complex e0[3],h0[3],e1[3],h1[3],te0,te1;
  double r[3];
  int i,lid0,lid1;
  
  printf("-- verification using boundary conditions --\n");
  
  r[0]=rand_num_d(-10.0,10.0,100);
  r[1]=rand_num_d(-10.0,10.0,100); 
  for(i=1;i<=ml->Ns;i++){
    printf("boundary %d, z=%15.14f\n",i,ml->zj[i]);
    r[2]=ml->zj[i];
 
    lid0=EH_t_mlpw(e0,h0,r,0,ml);
    printf("  electromagnetic field at r=(%g, %g, %g), layer id=%d\n",r[0],r[1],r[2],lid0);
    printf("    E_x=% 15.14e %+15.14e I\n",creal(e0[0]),cimag(e0[0]));
    printf("    E_y=% 15.14e %+15.14e I\n",creal(e0[1]),cimag(e0[1]));
    printf("    E_z=% 15.14e %+15.14e I\n",creal(e0[2]),cimag(e0[2]));
    printf("    H_x=% 15.14e %+15.14e I\n",creal(h0[0]),cimag(h0[0]));
    printf("    H_y=% 15.14e %+15.14e I\n",creal(h0[1]),cimag(h0[1]));
    printf("    H_z=% 15.14e %+15.14e I\n",creal(h0[2]),cimag(h0[2]));
    lid1=EH_t_mlpw(e1,h1,r,1,ml);
    printf("  electromagnetic field at r=(%g, %g, %g), layer id=%d\n",r[0],r[1],r[2],lid1);
    printf("    E_x=% 15.14e %+15.14e I\n",creal(e1[0]),cimag(e1[0]));
    printf("    E_y=% 15.14e %+15.14e I\n",creal(e1[1]),cimag(e1[1]));
    printf("    E_z=% 15.14e %+15.14e I\n",creal(e1[2]),cimag(e1[2]));
    printf("    H_x=% 15.14e %+15.14e I\n",creal(h1[0]),cimag(h1[0]));
    printf("    H_y=% 15.14e %+15.14e I\n",creal(h1[1]),cimag(h1[1]));
    printf("    H_z=% 15.14e %+15.14e I\n",creal(h1[2]),cimag(h1[2]));
    printf("  boundary condition : n dot D_0 = n dot D_1, n dot B_0 = n dot B_1, n=(0,0,1)\n");
    te0=ml->nj[lid0]*ml->nj[lid0]*e0[2];
    printf("    n dot D_0 = % 15.14f %+15.14f I\n",creal(te0),cimag(te0));
    te1=ml->nj[lid1]*ml->nj[lid1]*e1[2];
    printf("    n dot D_1 = % 15.14f %+15.14f I\n",creal(te1),cimag(te1));
    printf("    n dot B_0 = % 15.14f %+15.14f I\n",creal(h0[2]),cimag(h0[2]));
    printf("    n dot B_1 = % 15.14f %+15.14f I\n",creal(h1[2]),cimag(h1[2]));
    printf("  boundary condition : n times E_0 = n times E_1, n times H_0 = n times H_1\n");
    printf("    n times E_0 = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, 0)\n",creal(-e0[1]),cimag(-e0[1]),creal(e0[0]),cimag(e0[0]));
    printf("    n times E_1 = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, 0)\n",creal(-e1[1]),cimag(-e1[1]),creal(e1[0]),cimag(e1[0]));
    printf("    n times H_0 = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, 0)\n",creal(-h0[1]),cimag(-h0[1]),creal(h0[0]),cimag(h0[0]));
    printf("    n times H_1 = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, 0)\n",creal(-h1[1]),cimag(-h1[1]),creal(h1[0]),cimag(h1[0]));
    printf("\n");
  }
 
}

void verification_2(Mlpw *ml)
{
  double complex e[3],h[3],grad_e[9],grad_h[9],c1,c2,c3;
  double r[3];
  int lid;
  
  r[0]=rand_num_d(-10.0,10.0,100);
  r[1]=rand_num_d(-10.0,10.0,100);
  r[2]=rand_num_d(-10.0,10.0,100);
  
  printf("-- verification using Maxwell's equations --\n");
 
  lid=grad_EH_t_mlpw(e,grad_e,h,grad_h,r,0,ml);
  printf("electromagnetic field at r=(%g, %g, %g), layer id=%d\n",r[0],r[1],r[2],lid);
  printf("  E_x=% 15.14e %+15.14e I\n",creal(e[0]),cimag(e[0]));
  printf("  E_y=% 15.14e %+15.14e I\n",creal(e[1]),cimag(e[1]));
  printf("  E_z=% 15.14e %+15.14e I\n",creal(e[2]),cimag(e[2]));
  printf("  H_x=% 15.14e %+15.14e I\n",creal(h[0]),cimag(h[0]));
  printf("  H_y=% 15.14e %+15.14e I\n",creal(h[1]),cimag(h[1]));
  printf("  H_z=% 15.14e %+15.14e I\n",creal(h[2]),cimag(h[2]));
  
  printf("Faraday equation : rot E = -dB/dt\n");
  c1=grad_e[2*3+1]-grad_e[1*3+2]; // dE_x/dy - dE_y/dz
  c2=grad_e[0*3+2]-grad_e[2*3+0]; // dE_x/dz - dE_z/dx
  c3=grad_e[1*3+0]-grad_e[0*3+1]; // dE_y/dx - dE_x/dy
  printf("   rot E = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, % 15.14f %+15.14f I)\n",creal(c1),cimag(c1),creal(c2),cimag(c2),creal(c3),cimag(c3));
  c1=I*ml->pw.k0*h[0]; // i * omega * mu0 * H_x
  c2=I*ml->pw.k0*h[1]; // i * omega * mu0 * H_y
  c3=I*ml->pw.k0*h[2]; // i * omega * mu0 * H_z
  printf("  -dB/dt = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, % 15.14f %+15.14f I)\n",creal(c1),cimag(c1),creal(c2),cimag(c2),creal(c3),cimag(c3));
  
  printf("Ampere equation  : rot H = dD/dt\n");
  c1=grad_h[2*3+1]-grad_h[1*3+2]; // dH_x/dy - dH_y/dz
  c2=grad_h[0*3+2]-grad_h[2*3+0]; // dH_x/dz - dH_z/dx
  c3=grad_h[1*3+0]-grad_h[0*3+1]; // dH_y/dx - dH_x/dy
  printf("   rot H = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, % 15.14f %+15.14f I)\n",creal(c1),cimag(c1),creal(c2),cimag(c2),creal(c3),cimag(c3));
  c1=-I*ml->pw.k0*ml->nj[lid]*ml->nj[lid]*e[0]; // -i * omega * epsilon * E_x
  c2=-I*ml->pw.k0*ml->nj[lid]*ml->nj[lid]*e[1]; // -i * omega * epsilon * E_y
  c3=-I*ml->pw.k0*ml->nj[lid]*ml->nj[lid]*e[2]; // -i * omega * epsilon * E_z
  printf("   dD/dt = (% 15.14f %+15.14f I, % 15.14f %+15.14f I, % 15.14f %+15.14f I)\n",creal(c1),cimag(c1),creal(c2),cimag(c2),creal(c3),cimag(c3));
 
  printf("Gauss's law in this system : div D = 0, div B = 0\n");
  c1=grad_e[0*3+0]+grad_e[1*3+1]+grad_e[2*3+2]; 
  c2=grad_h[0*3+0]+grad_h[1*3+1]+grad_h[2*3+2];
  printf("   div D = % 15.14f %+15.14f I\n",creal(c1),cimag(c1));
  printf("   div B = % 15.14f %+15.14f I\n",creal(c2),cimag(c2));
}

