/* TM score = 1/L sum_i 1/[1+(di/d0)^2]
   d0=TM_coeff*(L-L_offset)^(1/3)-TM_offset
*/
#define TM_coeff 1.24
#define TM_offset 1.8
#define L_offset 15
#include "McLachlan.h"
#include "tm_score.h"
#include "math.h"
#include "stdio.h"

static int verbose=0;

float dist2(float *x1, float *x2);

float TM_score(float **xca1_store, int *ali1, int nres_1,
	       float **xca2_store, int *ali2, int nres_2,
	       int nali)
{
  int lmin=nres_1; if(nres_2<lmin)lmin=nres_2;
  if(lmin <= L_offset)return(-1);

  float d0=TM_coeff*pow(lmin-L_offset, 1./3)-TM_offset;
  if(d0 <= 0)return(-1); d0*=d0;

  // Copy coordinates of aligned residues
  int i, n=0, k=0;
  float xca1[3*lmin], xca2[3*lmin], w[lmin];
  for(i=0; i<nali; i++){
    int i1=ali1[i], i2=ali2[i];
    if((i1<0)||(i2<0))continue;
    w[n]=1; n++;
    for(int a=0; a<3; a++){
      xca1[k]=xca1_store[i1][a];
      xca2[k]=xca2_store[i2][a]; k++;
    }
  }

  double TM_max=-1;
  // Optimizing rotations
  for(int L_ini=n/8; L_ini<=n; L_ini*=2){
    for(int ini=0; ini<n; ini+=L_ini){
      int end=ini+L_ini; if(end>n)break;
      // initial superimposed fragment
      for(i=0; i<n; i++){
	if((i>=ini)&&(i<end)){w[i]=1;}else{w[i]=0;}
      }
      double TM_tmp=0;
      int discard=-1, round;

      // Optimize rotations by removing most divergent residues
      for(round=0; round<n; round++){
	rmsd_mclachlan_f(xca1, xca2, w, n);
	double TM=0, d2_max=0; int i_max=-1; k=0;
	for(i=0; i<n; i++){
	  float d2=dist2(xca1+k, xca2+k); k+=3;
	  TM+=1./(1+d2/d0);
	  if(w[i]&&(d2>d2_max)){d2_max=d2; i_max=i;}
	}
	if(i_max>=0){w[i_max]=0;}
	if(TM>TM_tmp){TM_tmp=TM;}
	else if(round){break;}
	discard++;
      }
      // Optimize rotations superimposing only if di<d0
      for(round=0; round<n; round++){
	  rmsd_mclachlan_f(xca1, xca2, w, n);
	  double TM=0; discard=0; k=0;
	  for(i=0; i<n; i++){
	    float d2=dist2(xca1+k, xca2+k); k+=3;
	    TM+=1./(1+d2/d0);
	    if(d2<d0){w[i]=1;}else{w[i]=0; discard++;}
	  }
	  if(TM>TM_tmp){TM_tmp=TM;}
	  else if(round){break;}
	} // end rounds
	if(TM_tmp>TM_max){
	  TM_max=TM_tmp;
	  if(verbose){
	    printf("TM= %.1f  d0=%.1f lmin=%d discard=%d %d rounds\n",
		   TM_max, sqrt(d0), lmin, discard, round+1);
	  }
	}
    }
  }
  if(verbose){
    printf("TM= %.1f  d0=%.1f lmin=%d\n", TM_max, sqrt(d0), lmin);
  }
  return(TM_max/lmin);
}

float dist2(float *x1, float *x2){
  float dx=x1[0]-x2[0], dy=x1[1]-x2[1], dz=x1[2]-x2[2];
  return(dx*dx+dy*dy+dz*dz);
}
