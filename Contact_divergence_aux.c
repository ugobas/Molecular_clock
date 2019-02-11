#include "Contact_divergence_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "NeedlemanWunsch.h"
int VERBOSE=1;  // Verbose: 

static int Equivalent(char aa1, char aa2);
static char *Remove_gaps(char *gapped_seq, int N_ali, int *N);

int Align_seq(int *Prot_ali, int N_ali, char *Prot_seq, char *PDB_seq, int n2)
{
  // Prot_ali: position l in the PDB aligned to i in seq
  int l=0, i, Lali=0, lmax=n2-1;
  for(i=0; i<N_ali; i++){
    Prot_ali[i]=-1;
    if(Prot_seq[i]!='-'){
      if(Equivalent(Prot_seq[i],PDB_seq[l])){
	Prot_ali[i]=l; l++; Lali++;
      }else if(l==0){
	while((l<lmax)&&(Equivalent(Prot_seq[i],PDB_seq[l]==0)))l++;
      }else if(l<n2){
	goto align;
      }
    }
  }
  if(VERBOSE)printf("%d residues perfectly aligned\n", Lali);
  return(Lali);

 align:
  Lali=0;
  // Sequence alignment parameters
  int IDE=1;  // Use identity to score alignment
  int GAP=7;  // Gap opening penalty
  int n1, nali;
  char *seq_nogap=Remove_gaps(Prot_seq, N_ali, &n1);
  char *ali1=malloc((n1+n2)*sizeof(char));
  char *ali2=malloc((n1+n2)*sizeof(char));
  int verbose=0;
  if(alignNW(seq_nogap,n1,PDB_seq,n2,verbose,IDE,GAP,ali1,ali2,&nali)==0){
    printf("ERROR, alignment failed\n"); return(-1);
  }
  int mismatch=0, nomatch[n1], lnomatch[n1];
  int l1=-1, l2=0, gap=0;;
  int *ali=malloc(n1*sizeof(int));
  for(i=0; i<nali; i++){
    if(ali1[i]!='-'){
      l1++;
      if(ali2[i]!='-'){
	ali[l1]=l2; Lali++;
	if(Equivalent(ali1[i],ali2[i])==0){
	  printf("Mismatch: %c-%c i=%d\n", ali1[i], ali2[i], i);
	  nomatch[mismatch]=i; lnomatch[mismatch]=l1;
	  mismatch++;
	}
      }else{
	ali[l1]=-1;  gap++;
      }
    }
    if(ali2[i]!='-')l2++;
  }
  l1=0;
  for(i=0; i<N_ali; i++){
    if(Prot_seq[i]=='-'){
      Prot_ali[i]=-1;
    }else{
      Prot_ali[i]=ali[l1]; l1++;
    }
  }
  if(VERBOSE){
    for(i=0; i<N_ali; i++)printf("%c", Prot_seq[i]); printf("\n");
    for(i=0; i<N_ali; i++){
      int k=Prot_ali[i];
      if(k>=0){printf("%c", PDB_seq[k]);}else{printf("-");}
    }
    printf("\n");
    printf("%d residues aligned with %d gaps and %d mismatches: ",
	   Lali, gap, mismatch);
    for(i=0; i<mismatch; i++){
      int k=nomatch[i];
      printf(" %c%d%c", ali1[k], lnomatch[i], ali2[k]);
    }
    printf("\n");
  }
  free(ali1); free(ali2); free(ali); free(seq_nogap);
  return(Lali);
}

static char *Remove_gaps(char *gapped_seq, int N_ali, int *N){
  *N=0;
  char *seq=malloc(N_ali*sizeof(char));
  for(int i=0; i<N_ali; i++){
    if(gapped_seq[i]!='-'){
      seq[*N]=gapped_seq[i]; (*N)++;
    }
  }
  return(seq);
}

static int Equivalent(char aa1, char aa2){
  if((aa1==aa2)||(aa1=='X')||(aa2=='X'))return(1);
  return(0);
}

float Seq_identity(char *seq1, char *seq2, int N_ali){
  int i, id=0, l1=0, l2=0; float si;
  for(i=0; i<N_ali; i++){
    if(seq1[i]!='-'){
      l1++; if(seq1[i]==seq2[i])id++;
    }
    if(seq2[i]!='-')l2++; 
  }
  if(l1>l2){si=(float)id/l1;}else{si=(float)id/l2;}
  return(si);
}

float Compute_overlap(short **Cont_mat_1, int N1, int *ali1,
		      short **Cont_mat_2, int N2, int *ali2,
		      int N_ali)
{
  float q=0;
  short i, *j1, *j2, i2, j1ali;
  int num_cont_1=0, num_cont_2=0, num_cont_12=0;

  // Align sequence 1 and 2
  int *ali=malloc(N1*sizeof(int));
  for(i=0; i<N1; i++)ali[i]=-1;
  for(i=0; i<N_ali; i++){
    if(ali1[i]>=0)ali[ali1[i]]=ali2[i];
  }
  for(i=0; i<N1; i++){
    j1=&Cont_mat_1[i][0];
    while(*j1>=0){num_cont_1++; j1++;}
  }
  for(i=0; i<N2; i++){
    j1=&Cont_mat_2[i][0];
    while(*j1>=0){num_cont_2++; j1++;}
  }
  for(i=0; i<N1; i++){
    i2=ali[i]; if(i2<0)continue;
    j1=&Cont_mat_1[i][0];   
    j2=&Cont_mat_2[i2][0];
    while(*j1>=0){
      j1ali=ali[*j1];
      if(j1ali >=0){
	while((*j2<j1ali)&&(*j2>=0))j2++;
	if(*j2==j1ali)num_cont_12++;
      }
      j1++;
    }
  }
  if((num_cont_1!=0)&&(num_cont_2!=0)){
    q=(float)num_cont_12/sqrt((float)(num_cont_1*num_cont_2));
  }
  return(q);

}


char Get_compression(char *pdb_name){
  char *tmp=pdb_name;
  while((*tmp!='\0')&&(*tmp!='\n')){
    if(*tmp=='.'){
      if((*(tmp+1)=='g')&&(*(tmp+2)=='z')){
	return('g');
      }else if((*(tmp+1)=='Z')&&(*(tmp+2)=='\0')){
	return('Z');
      }else if((*(tmp+1)=='z')&&(*(tmp+2)=='i')&&
	       (*(tmp+3)=='p')&&(*(tmp+4)=='\0')){
	return('z');
      }
    }
    tmp++;
  }
  return('\0');
}

int Find_name(char *name, char **names, int N, int Np)
{
  int ip;
  for(ip=Np; ip<N; ip++)if(strcmp(name, names[ip])==0)return(ip);
  for(ip=0; ip<Np; ip++)if(strcmp(name, names[ip])==0)return(ip);
  return(-1);
}

int Change_ext(char *name_out, char *name_in, char *ext){
  char *ptr1=name_in, *ptr2=name_out;
  while((*ptr1!='\0')&&(*ptr1!='.')){
    *ptr2=*ptr1; ptr1++; ptr2++;
  }
  sprintf(name_out, "%s%s", name_in, ext);
  //sprintf(name_out, "%s%s", name_out, ext);
  printf("file name: %s\n", name_out);
  return(0);
}
