#include "D_Cont.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define qave_EXP  0.547
#define qave_COEF 0.386
#define qdev_EXP  0.673
#define qdev_COEF 1.327
#define L_max 2000.0
#define A 5.0
int Ini_Dcont=0;
float D0;

float Compute_Dcont(float q, int L1, int L2, int *related){
  float D1, D2;
  double L12, qave, qdev, qinf, z;
  if(Ini_Dcont==0){
    Ini_Dcont=1; L12=L_max;
    qave=qave_COEF*pow(L12, -qave_EXP);
    qdev=qdev_COEF*pow(L12, -qdev_EXP);
    qinf=qave+A*qdev;
    D0=1.+A-log(qdev)+log(1.-qinf);
  }

  L12=sqrt((float)(L1*L2));
  qave=qave_COEF*pow(L12, -qave_EXP);
  qdev=qdev_COEF*pow(L12, -qdev_EXP);
  qinf=qave+A*qdev;
  z=(q-qave)/qdev;
  D2=D0-z;
  //printf("qave= %.4f qdev= %.4f qinf= %.4f D0=%.2f L=%.0f\n",
  // qave, qdev, qinf, D0, L12);
  if(q>qinf){
    D1=-log((q-qinf)/(1.-qinf));
    *related=1;
    //printf("D1=%.2f D2=%.2f\n", D1, D2);
    if((z>1)||(D1<D2))return(D1);
  }
  *related=0;
  return(D2);
}
