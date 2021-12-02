#include "emf_multilayer.h"

void read_data_mlpw(char *fname_ipw,char *fname_ml,Mlpw *ml)
{
  void read_data_ml(char *fname_ml,Mlpw *ml);
  
  // multilayer
  read_data_ml(fname_ml,ml);
  // incident plane wave 
  read_data_pwsp(fname_ipw,&(ml->pw));
}

void mfree_mlpw(Mlpw *ml)
{
  void mfree_ml(Mlpw *md);
  
  mfree_ml(ml);
}

void print_data_mlpw(Mlpw *ml)
{
  void print_data_ml(Mlpw *md);
  
  // incident plane wave 
  print_data_pwsp(&(ml->pw)); 
  // multilayer
  print_data_ml(ml);
}
 
void print_data_mlpw_mksa(Mlpw *ml)
{
  void print_data_ml_mksa(Mlpw *md);
  
  // incident plane wave 
  print_data_pwsp_mksa(&(ml->pw)); 
  // multilayer
  print_data_ml_mksa(ml);
}

void setup_mlpw(Mlpw *ml)
{
  void calc_coefficient_C (Mlpw *md);
  void calc_coefficient_EH(Mlpw *md);
  
  double complex ce;
  double k02,i_k0,kxy2,arg;
  int i;
  
  // incident plane wave
  setup_pwsp(&(ml->pw)); 
  if(ml->pw.theta<0.0 || ml->pw.theta>=0.5*M_PI){
    printf("incident angle is out of range.  theta=%g. Exit...\n",ml->pw.theta);
    exit(1);
  }
 
  k02=ml->pw.k0*ml->pw.k0;
  i_k0=1.0/ml->pw.k0;
  kxy2=ml->pw.vk[0]*ml->pw.vk[0]+ml->pw.vk[1]*ml->pw.vk[1];

  // layer 0 data
  ml->nj[0]=ml->pw.n;
  ml->kzj[0]=ml->pw.vk[2];
  ml->zj[0]=ml->zj[1];
  // E0+,H0+ (incident field amplitude)
  arg=ml->pw.vk[2]*ml->zj[0];
  ce=cos(arg)+I*sin(arg);
  for(i=0;i<3;i++){
    ml->Ep[0][i]=ml->pw.e[i]*ce;
    ml->Hp[0][i]=ml->pw.h[i]*ce;
  }
  // Z_j (wave impedance)
  for(i=0;i<=ml->No;i++) ml->Zj[i]=1.0/ml->nj[i];
  // k_zj, cos(theta_j) 
  ml->cos_tj[0]=ml->kzj[0]*i_k0*ml->Zj[0];
  for(i=1;i<=ml->No;i++){
    ml->kzj[i]=csqrt(k02*ml->nj[i]*ml->nj[i]-kxy2); // select positive (principal value of square root)
    ml->cos_tj[i]=ml->kzj[i]*i_k0*ml->Zj[i];
    if(cabs(ml->kzj[i])==0.0){
      printf("There is no solution since k_z becomes 0. Exit...\n");
      printf("layer id = %d. kx=%g, ky=%g, kz=%g %+g I\n",i,ml->pw.vk[0],ml->pw.vk[1],creal(ml->kzj[i]),cimag(ml->kzj[i]));
      exit(1);
    }
  }
  // attenuation check
  ml->Ns=ml->No;
  for(i=1;i<ml->No;i++){
    ce=cexp(I*ml->kzj[i]*(ml->zj[i+1]-ml->zj[i]));
    if(AT_COEF>cabs(ce)){
      printf("The wave attenuated to zero in the layer %d.\n",i);
      printf("Consider this layer as the last (ignore boundaries after the z_%02d boundary).\n\n",i);
      ml->Ns=i;
      break;
    }
  }
  // calc coefficient 
  calc_coefficient_C(ml);
  calc_coefficient_EH(ml);
}

int layer_id_type0(double z,Mlpw *ml)
{
  int i;
  
  if(ml->zj[1]>=z) return 0;
  for(i=2;i<=ml->Ns;i++){
    if(ml->zj[i]>=z) return i-1;
  }
  return ml->Ns;
}

int layer_id_type1(double z,Mlpw *ml)
{
  int i;
  
  if(ml->zj[1]>z) return 0;
  for(i=2;i<=ml->Ns;i++){
    if(ml->zj[i]>z) return i-1;
  }
  return ml->Ns;
}

int EH_t_mlpw(double complex *E,double complex *H,double *r,int type,Mlpw *ml)
{
  int lid;
  double tc,tz,tti;
  double complex texy,tzp,tzm,ttr;
  
  if(type==0) lid=layer_id_type0(r[2],ml);
  else        lid=layer_id_type1(r[2],ml);
  tc=ml->pw.vk[0]*r[0]+ml->pw.vk[1]*r[1];
  texy=cos(tc)+I*sin(tc);

  tz=creal(ml->kzj[lid])*(r[2]-ml->zj[lid]);
  ttr=cos(tz)+I*sin(tz);
  tti=exp(-cimag(ml->kzj[lid])*(r[2]-ml->zj[lid]));
  tzp=     ttr *tti;
  if(lid!=ml->Ns) tzm=conj(ttr)/tti;
  else tzm=0.0;

  E[0]=(ml->Ep[lid][0]*tzp+ml->Em[lid][0]*tzm)*texy;
  E[1]=(ml->Ep[lid][1]*tzp+ml->Em[lid][1]*tzm)*texy;
  E[2]=(ml->Ep[lid][2]*tzp+ml->Em[lid][2]*tzm)*texy;

  H[0]=(ml->Hp[lid][0]*tzp+ml->Hm[lid][0]*tzm)*texy;
  H[1]=(ml->Hp[lid][1]*tzp+ml->Hm[lid][1]*tzm)*texy;
  H[2]=(ml->Hp[lid][2]*tzp+ml->Hm[lid][2]*tzm)*texy;
  
  return lid;
}

int EH_i_mlpw(double complex *e,double complex *h,double *r,int type,Mlpw *ml)
{
  int lid;
  
  if(type==0) lid=layer_id_type0(r[2],ml);
  else        lid=layer_id_type1(r[2],ml);

  EH_pwsp(e,h,r,&(ml->pw)); 

  return lid;
}

int EH_p_mlpw(double complex *E,double complex *H,double *r,int type,Mlpw *ml)
{
  int lid;
  double tc,tz,tti;
  double complex texy,tzp,ttr;
  
  if(type==0) lid=layer_id_type0(r[2],ml);
  else        lid=layer_id_type1(r[2],ml);
  tc=ml->pw.vk[0]*r[0]+ml->pw.vk[1]*r[1];
  texy=cos(tc)+I*sin(tc);

  tz=creal(ml->kzj[lid])*(r[2]-ml->zj[lid]);
  ttr=cos(tz)+I*sin(tz);
  tti=exp(-cimag(ml->kzj[lid])*(r[2]-ml->zj[lid]));
  tzp=     ttr *tti;

  E[0]=ml->Ep[lid][0]*tzp*texy;
  E[1]=ml->Ep[lid][1]*tzp*texy;
  E[2]=ml->Ep[lid][2]*tzp*texy;

  H[0]=ml->Hp[lid][0]*tzp*texy;
  H[1]=ml->Hp[lid][1]*tzp*texy;
  H[2]=ml->Hp[lid][2]*tzp*texy;
  
  return lid;
}

int EH_m_mlpw(double complex *E,double complex *H,double *r,int type,Mlpw *ml)
{
  int lid;
  double tc,tz,tti;
  double complex texy,tzm,ttr;
  
  if(type==0) lid=layer_id_type0(r[2],ml);
  else        lid=layer_id_type1(r[2],ml);
  tc=ml->pw.vk[0]*r[0]+ml->pw.vk[1]*r[1];
  texy=cos(tc)+I*sin(tc);

  tz=creal(ml->kzj[lid])*(r[2]-ml->zj[lid]);
  ttr=cos(tz)+I*sin(tz);
  tti=exp(-cimag(ml->kzj[lid])*(r[2]-ml->zj[lid]));
  if(lid!=ml->Ns) tzm=conj(ttr)/tti;
  else tzm=0.0;

  E[0]=ml->Em[lid][0]*tzm*texy;
  E[1]=ml->Em[lid][1]*tzm*texy;
  E[2]=ml->Em[lid][2]*tzm*texy;

  H[0]=ml->Hm[lid][0]*tzm*texy;
  H[1]=ml->Hm[lid][1]*tzm*texy;
  H[2]=ml->Hm[lid][2]*tzm*texy;
  
  return lid;
}

void power_RT_coefficient(double *R,double *T,Mlpw *ml)
{
  R[0]=ml->Rp;
  R[1]=ml->Rs;
  T[0]=ml->Tp;
  T[1]=ml->Ts;
}

int grad_EH_t_mlpw(double complex *E,double complex *grad_E,double complex *H,double complex *grad_H,double *r,int type,Mlpw *ml)
{
  int lid,i;
  double tc,tz,tti;
  double complex texy,tzp,tzm,ttr,ce1,ce2;
  
  if(type==0) lid=layer_id_type0(r[2],ml);
  else        lid=layer_id_type1(r[2],ml);
  tc=ml->pw.vk[0]*r[0]+ml->pw.vk[1]*r[1];
  texy=cos(tc)+I*sin(tc);

  tz=creal(ml->kzj[lid])*(r[2]-ml->zj[lid]);
  ttr=cos(tz)+I*sin(tz);
  tti=exp(-cimag(ml->kzj[lid])*(r[2]-ml->zj[lid]));
  tzp=     ttr *tti;
  if(lid!=ml->Ns) tzm=conj(ttr)/tti;
  else tzm=0.0;

  E[0]=(ml->Ep[lid][0]*tzp+ml->Em[lid][0]*tzm)*texy;
  E[1]=(ml->Ep[lid][1]*tzp+ml->Em[lid][1]*tzm)*texy;
  E[2]=(ml->Ep[lid][2]*tzp+ml->Em[lid][2]*tzm)*texy;

  H[0]=(ml->Hp[lid][0]*tzp+ml->Hm[lid][0]*tzm)*texy;
  H[1]=(ml->Hp[lid][1]*tzp+ml->Hm[lid][1]*tzm)*texy;
  H[2]=(ml->Hp[lid][2]*tzp+ml->Hm[lid][2]*tzm)*texy;
  
  // d/dx, d/dy
  ce1=I*ml->pw.vk[0];
  ce2=I*ml->pw.vk[1];
  for(i=0;i<3;i++){
    grad_E[i*3+0]=E[i]*ce1;
    grad_H[i*3+0]=H[i]*ce1;
    grad_E[i*3+1]=E[i]*ce2;
    grad_H[i*3+1]=H[i]*ce2;
  }
  // d/dz
  ce1= I*ml->kzj[lid];
  ce2=-I*ml->kzj[lid];
  for(i=0;i<3;i++){
    grad_E[i*3+2]=(ce1*ml->Ep[lid][i]*tzp+ce2*ml->Em[lid][i]*tzm)*texy;
    grad_H[i*3+2]=(ce1*ml->Hp[lid][i]*tzp+ce2*ml->Hm[lid][i]*tzm)*texy;
  }
  
  return lid;
}

//////////////////////////////////////////////////////////////////////
void read_data_ml(char *fname_ml,Mlpw *md)
{
  void malloc_ml(Mlpw *md);
  
  FILE *fp;
  char buf[256]="";
  double td1,td2;
  int ti,i;
  
  if((fp=fopen(fname_ml,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",fname_ml);    exit(1);  }
  if(fgets(buf,256,fp)==NULL){
    printf("emf_multilayer.c, read_data_ml(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("emf_multilayer.c, read_data_ml(), failed to read the line. exit...\n");
    exit(1);
  }
  
  if(fscanf(fp,"%d\n",&ti)!=1){
    printf("emf_multilayer.c, read_data_ml(), failed to read the No. exit...\n");
    exit(1);
  }
  md->No=ti;
  if(fgets(buf,256,fp)==NULL){
    printf("emf_multilayer.c, read_data_ml(), failed to read the line. exit...\n");
    exit(1);
  }

  malloc_ml(md);
  for(i=1;i<=ti;i++){
    if(fscanf(fp,"%lf\n",&td1)!=1){
      printf("emf_multilayer.c, read_data_ml(), failed to read the zj. exit...\n");
      exit(1);
    }
    md->zj[i]=td1;
    if(fscanf(fp,"%lf\n",&td1)!=1){
      printf("emf_multilayer.c, read_data_ml(), failed to read the real(nj). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf\n",&td2)!=1){
      printf("emf_multilayer.c, read_data_ml(), failed to read the imag(nj). exit...\n");
      exit(1);
    }
    md->nj[i]=td1+I*td2;
  }
  
  fclose(fp);
}

void malloc_ml(Mlpw *md)
{
  int N,i;
  
  N=md->No;
  
  md->zj=(double *)m_alloc2(N+1,sizeof(double),"emf_multilayer.c, malloc_ml(), md->zj");
  md->nj=(double complex *)m_alloc2(N+1,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->nj");
  
  md->Zj    =(double complex *)m_alloc2(N+1,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->Zj");
  md->kzj   =(double complex *)m_alloc2(N+1,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->kzj");
  md->cos_tj=(double complex *)m_alloc2(N+1,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->cos_tj");
  md->eC    =(double complex **)m_alloc2(N+1,sizeof(double complex *),"emf_multilayer.c, malloc_ml(), md->eC");
  md->mC    =(double complex **)m_alloc2(N+1,sizeof(double complex *),"emf_multilayer.c, malloc_ml(), md->mC");
  md->Ep    =(double complex **)m_alloc2(N+1,sizeof(double complex *),"emf_multilayer.c, malloc_ml(), md->Ep");
  md->Em    =(double complex **)m_alloc2(N+1,sizeof(double complex *),"emf_multilayer.c, malloc_ml(), md->Em");
  md->Hp    =(double complex **)m_alloc2(N+1,sizeof(double complex *),"emf_multilayer.c, malloc_ml(), md->Hp");
  md->Hm    =(double complex **)m_alloc2(N+1,sizeof(double complex *),"emf_multilayer.c, malloc_ml(), md->Hm");
  for(i=0;i<=N;i++){
    md->eC[i]=(double complex *)m_alloc2(2,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->eC[i]");
    md->mC[i]=(double complex *)m_alloc2(2,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->mC[i]");
    md->Ep[i]=(double complex *)m_alloc2(3,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->Ep[i]");
    md->Em[i]=(double complex *)m_alloc2(3,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->Em[i]");
    md->Hp[i]=(double complex *)m_alloc2(3,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->Hp[i]");
    md->Hm[i]=(double complex *)m_alloc2(3,sizeof(double complex),"emf_multilayer.c, malloc_ml(), md->Hm[i]");
  }
}

void mfree_ml(Mlpw *md)
{
  int i;
  
  for(i=0;i<=md->No;i++){
    free(md->eC[i]);    free(md->mC[i]);
    free(md->Ep[i]);    free(md->Em[i]);
    free(md->Hp[i]);    free(md->Hm[i]);
  }
  free(md->eC);  free(md->mC);
  free(md->Ep);  free(md->Em);
  free(md->Hp);  free(md->Hm);
  free(md->zj);
  free(md->nj);
  free(md->Zj);
  free(md->kzj);
  free(md->cos_tj);
  
  md->No=0;
  md->Ns=0;
}

void print_data_ml(Mlpw *md)
{
  int i;
  
  printf("-- multilayer --\n");
  printf("number of layers : %-d\n",md->No+1);
  printf(" layer %d, refractive index n_%02d = %-7.6f %+7.6f I\n",0,0,md->pw.n,0.0);
  for(i=1;i<=md->No;i++){
    printf("z_%02d ------------------------------------------------ : %-g\n",i,md->zj[i]);
    printf(" layer %d, refractive index n_%02d = %-7.6f %+7.6f I\n",i,i,creal(md->nj[i]),cimag(md->nj[i]));
  }
  printf("\n");
}

void print_data_ml_mksa(Mlpw *md)
{
  int i;
  
  printf("-- multilayer --\n");
  printf("number of layers : %-d\n",md->No+1);
  printf(" layer %d, refractive index n_%02d = %-7.6f %+7.6f I\n",0,0,md->pw.n,0.0);
  for(i=1;i<=md->No;i++){
    printf("z_%02d ---------------------------------------------[m] : %-g\n",i,OSUtoMKSA_length(md->zj[i]));
    printf(" layer %d, refractive index n_%02d = %-7.6f %+7.6f I\n",i,i,creal(md->nj[i]),cimag(md->nj[i]));
  }
  printf("\n");
}

void calc_coefficient_C(Mlpw *md)
{
  int i,j;
  double eki0,eki1;
  double complex te[4],tm[4],ek0p,ek0m,ek1p,ek1m,i_ce,i_cm,c0e,c1e,c0m,c1m,ekr0,ekr1,**eM,**mM;
  
  eM=(double complex **)m_alloc2(md->Ns+1,sizeof(double complex *),"emf_multilayer.c, calc_coefficient_C(),eM");
  mM=(double complex **)m_alloc2(md->Ns+1,sizeof(double complex *),"emf_multilayer.c, calc_coefficient_C(),mM");
  for(i=0;i<=md->Ns;i++){
    eM[i]=(double complex *)m_alloc2(4,sizeof(double complex),"emf_multilayer.c, calc_coefficient_C(),eM[i]");
    mM[i]=(double complex *)m_alloc2(4,sizeof(double complex),"emf_multilayer.c, calc_coefficient_C(),mM[i]");
  }
  // Transfer matrix
  eM[0][0]=1.0;  eM[0][1]=0.0;
  eM[0][2]=0.0;  eM[0][3]=1.0;
  mM[0][0]=1.0;  mM[0][1]=0.0;
  mM[0][2]=0.0;  mM[0][3]=1.0;
  for(j=1;j<=md->Ns;j++){
    ekr0=cos( creal(md->kzj[j-1])*(md->zj[j]-md->zj[j-1]))+I*sin( creal(md->kzj[j-1])*(md->zj[j]-md->zj[j-1]));
    eki0=exp(-cimag(md->kzj[j-1])*(md->zj[j]-md->zj[j-1]));
    ek0p=     ekr0 *eki0;
    ek0m=conj(ekr0)/eki0;
    ekr1=1.0;
    eki1=1.0;
    ek1p=     ekr1 *eki1;
    ek1m=conj(ekr1)/eki1;
    
    // E-wave
    i_ce=1.0/(2.0*md->Zj[j]*md->cos_tj[j]);
    c0e=i_ce*(md->Zj[j]*md->cos_tj[j]+md->Zj[j-1]*md->cos_tj[j-1]);
    c1e=i_ce*(md->Zj[j]*md->cos_tj[j]-md->Zj[j-1]*md->cos_tj[j-1]);
    
    eM[j][0]=c0e*ek0p*ek1m;    eM[j][1]=c1e*ek0m*ek1m;
    eM[j][2]=c1e*ek0p*ek1p;    eM[j][3]=c0e*ek0m*ek1p;
    
    te[0]=eM[j][0]*eM[0][0]+eM[j][1]*eM[0][2];    te[1]=eM[j][0]*eM[0][1]+eM[j][1]*eM[0][3];
    te[2]=eM[j][2]*eM[0][0]+eM[j][3]*eM[0][2];    te[3]=eM[j][2]*eM[0][1]+eM[j][3]*eM[0][3];
    eM[0][0]=te[0];    eM[0][1]=te[1];
    eM[0][2]=te[2];    eM[0][3]=te[3];
    
    // H-wave
    i_cm=1.0/(2.0*md->Zj[j-1]*md->cos_tj[j]);
    c0m=i_cm*(md->Zj[j-1]*md->cos_tj[j]+md->Zj[j]*md->cos_tj[j-1]);
    c1m=i_cm*(md->Zj[j-1]*md->cos_tj[j]-md->Zj[j]*md->cos_tj[j-1]);
    
    mM[j][0]=c0m*ek0p*ek1m;    mM[j][1]=c1m*ek0m*ek1m;
    mM[j][2]=c1m*ek0p*ek1p;    mM[j][3]=c0m*ek0m*ek1p;
    
    tm[0]=mM[j][0]*mM[0][0]+mM[j][1]*mM[0][2];    tm[1]=mM[j][0]*mM[0][1]+mM[j][1]*mM[0][3];
    tm[2]=mM[j][2]*mM[0][0]+mM[j][3]*mM[0][2];    tm[3]=mM[j][2]*mM[0][1]+mM[j][3]*mM[0][3];
    mM[0][0]=tm[0];    mM[0][1]=tm[1];
    mM[0][2]=tm[2];    mM[0][3]=tm[3];
  }
  
  if(cabs(eM[0][3])==0.0 || cabs(mM[0][3])==0.0){
    printf("The amplitude coefficient diverged. Exit...\n");
    exit(1);
  }
  else {
    md->eC[     0][0]=1.0;
    md->eC[     0][1]=-eM[0][2]/eM[0][3];
    md->eC[md->Ns][0]=(eM[0][0]*eM[0][3]-eM[0][1]*eM[0][2])/eM[0][3];
    md->eC[md->Ns][1]=0.0;

    md->mC[     0][0]=1.0;
    md->mC[     0][1]=-mM[0][2]/mM[0][3];
    md->mC[md->Ns][0]=(mM[0][0]*mM[0][3]-mM[0][1]*mM[0][2])/mM[0][3];
    md->mC[md->Ns][1]=0.0;

    for(j=1;j<md->Ns;j++){
      md->eC[j][0]=eM[j][0]*md->eC[j-1][0]+eM[j][1]*md->eC[j-1][1];
      md->eC[j][1]=eM[j][2]*md->eC[j-1][0]+eM[j][3]*md->eC[j-1][1];
      md->mC[j][0]=mM[j][0]*md->mC[j-1][0]+mM[j][1]*md->mC[j-1][1];
      md->mC[j][1]=mM[j][2]*md->mC[j-1][0]+mM[j][3]*md->mC[j-1][1];
    }
    
    // power coefficients
    md->Rp=creal(md->eC[0][1]*conj(md->eC[0][1]));
    md->Rs=creal(md->mC[0][1]*conj(md->mC[0][1]));
    md->Tp=creal(md->Zj[md->Ns]*md->cos_tj[md->Ns])/creal(md->Zj[0]*md->cos_tj[0])*creal(md->eC[md->Ns][0]*conj(md->eC[md->Ns][0]));
    md->Ts=creal(md->nj[md->Ns]*md->cos_tj[md->Ns])/creal(md->nj[0]*md->cos_tj[0])*creal(md->mC[md->Ns][0]*conj(md->mC[md->Ns][0]));
  }

  for(i=0;i<=md->Ns;i++){
    free(eM[i]);    free(mM[i]);
  }
  free(eM);  free(mM);
}

void calc_coefficient_EH(Mlpw *md)
{
  double complex eHx,eHy,mEx,mEy;
  double kxy2,k;
  int j;
  
  kxy2=md->pw.vk[0]*md->pw.vk[0]+md->pw.vk[1]*md->pw.vk[1];
  k=md->pw.k0*md->pw.n;

  // E-wave H-wave decomposition 
  if(kxy2!=0.0){
    eHx= k*md->pw.vk[1]/(md->Zj[0]*kxy2)*md->Ep[0][2];
    eHy=-k*md->pw.vk[0]/(md->Zj[0]*kxy2)*md->Ep[0][2];
    mEx=-k*md->pw.vk[1]*md->Zj[0]/kxy2*md->Hp[0][2];
    mEy= k*md->pw.vk[0]*md->Zj[0]/kxy2*md->Hp[0][2];
  }
  else {
    eHx=0.0;    eHy=0.0;
    mEx=md->Ep[0][0];
    mEy=md->Ep[0][1];
  }
  // coefficients
  md->Em[0][0]=-eHy*md->eC[0][1]*md->Zj[0]*md->cos_tj[0] + mEx*md->mC[0][1];
  md->Em[0][1]= eHx*md->eC[0][1]*md->Zj[0]*md->cos_tj[0] + mEy*md->mC[0][1];
  md->Em[0][2]=-(md->pw.vk[0]*eHy-md->pw.vk[1]*eHx)*md->eC[0][1]*md->Zj[0]/(md->pw.k0*md->nj[0]);
  md->Hm[0][0]= eHx*md->eC[0][1] + mEy*md->mC[0][1]*md->cos_tj[0]/md->Zj[0];
  md->Hm[0][1]= eHy*md->eC[0][1] - mEx*md->mC[0][1]*md->cos_tj[0]/md->Zj[0];
  md->Hm[0][2]= (md->pw.vk[0]*mEy-md->pw.vk[1]*mEx)*md->mC[0][1]/(md->Zj[0]*md->pw.k0*md->nj[0]);
  for(j=1;j<=md->Ns;j++){
    md->Ep[j][0]= eHy*md->eC[j][0]*md->Zj[j]*md->cos_tj[j] + mEx*md->mC[j][0];
    md->Ep[j][1]=-eHx*md->eC[j][0]*md->Zj[j]*md->cos_tj[j] + mEy*md->mC[j][0];
    md->Ep[j][2]=-(md->pw.vk[0]*eHy-md->pw.vk[1]*eHx)*md->eC[j][0]*md->Zj[j]/(md->pw.k0*md->nj[j]);
    md->Em[j][0]=-eHy*md->eC[j][1]*md->Zj[j]*md->cos_tj[j] + mEx*md->mC[j][1];
    md->Em[j][1]= eHx*md->eC[j][1]*md->Zj[j]*md->cos_tj[j] + mEy*md->mC[j][1];
    md->Em[j][2]=-(md->pw.vk[0]*eHy-md->pw.vk[1]*eHx)*md->eC[j][1]*md->Zj[j]/(md->pw.k0*md->nj[j]);

    md->Hp[j][0]= eHx*md->eC[j][0] - mEy*md->mC[j][0]*md->cos_tj[j]/md->Zj[j];
    md->Hp[j][1]= eHy*md->eC[j][0] + mEx*md->mC[j][0]*md->cos_tj[j]/md->Zj[j];
    md->Hp[j][2]= (md->pw.vk[0]*mEy-md->pw.vk[1]*mEx)*md->mC[j][0]/(md->Zj[j]*md->pw.k0*md->nj[j]);
    md->Hm[j][0]= eHx*md->eC[j][1] + mEy*md->mC[j][1]*md->cos_tj[j]/md->Zj[j];
    md->Hm[j][1]= eHy*md->eC[j][1] - mEx*md->mC[j][1]*md->cos_tj[j]/md->Zj[j];
    md->Hm[j][2]= (md->pw.vk[0]*mEy-md->pw.vk[1]*mEx)*md->mC[j][1]/(md->Zj[j]*md->pw.k0*md->nj[j]);
  }
}
