/* 
   Program Clock_violations
   Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)
   ubastolla@cbm.csic.es
   Reads a multiple alignment and computes contact divergence and other
   structure comparison measures.
   Estimates or reads outgroups from a tree and computes clock violations

   INPUT: file with multiple alignment in FASTA format (at least 2)
   Protein names must be names of PDB files.
   The first line may be PDBDIR=<directory of PDB files>
   (default: current directory)

   OUTPUT: For each protein pair, structural scores printed are
   Contact_divergence, contact overlap, TM score (default no)

*/

#include "Contact_divergence_aux.h"
#include "D_Cont.h"
#include "protein.h"
#include "cont_list.h"
#include "allocate.h"
#include "tm_score.h"
#include "read_structures.h"
#include "tree.h"
#include "CV_statistics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float S0=0.06;    // for Tajima-Nei divergence
//float alpha=0.7; // Exponent for computing clock violations
int PRINT_SIM=0;  // Print matrix of similarities?
int PRINT_DIV=1;  // Print matrix of divergences?
int PRINT_CV=1;   // Print clock violations?
int PRINT_DIFF=1; // Compare sequence and structure alignments if both given?
int PRINT_AVE=1;  // Print average values as a function of d?
int PRINT_GAP=1; // Examine relationship between Triangle Inequality and gaps?
float HUGE=10;    // Very large divergence
//float TM0=0.17;   // TM score of unrelated proteins

#define EXT_DIFF ".diff" // Extension for comparison of alignments
#define EXT_DIV ".div"   // Extension for divergence output file
#define EXT_SIM ".sim"   // Extension for similarity output file
#define EXT_CV  ".cv"    // Extension for clock-violation output file

#define IJ_MIN_DEF 3
#define CONT_TYPE_DEF 'c'
#define CONT_THR_DEF 4.5

float fun_high=0.90;
float fun_low=0.70;
int diff_thr=4;

int IJ_MIN=IJ_MIN_DEF;        // Only contacts with |i-j|>=IJ_MIN
float CONT_THR=CONT_THR_DEF;
char CONT_TYPE=CONT_TYPE_DEF;  //a=alpha b=beta c=all atoms
char CONT_STRING[80];
char CODENAME[40]="main_Evol_div.c";
int CONT_DEF=1;

struct Prot_input{
  char name[80];
  char chain;
  char *seq;
};

void help(char *pname);
void Get_input(char *file_ali, char *file_ali_str, char *file_fun,
	       char *name, char *PDBDIR, char *PDBEXT, char *OUTG,
	       int *PRINT_SIM, int *PRINT_CV, int *PRINT_DIV, int *PRINT_DIFF,
	       int argc, char **argv);
int Get_alignment_old(struct Prot_input **Prot_input, int *Nali,
		      char *PDBDIR, char *PDBEXT, int *PRINT_SIM,
		      char *file_ali);
int Get_alignment(struct Prot_input **Prot_input, int *Nali, char *file_ali,
		  char *PDBDIR, char *PDBEXT);
int *Match_alignments(struct Prot_input *Prot1,
		      struct Prot_input *Prot2, int N);
int Find_prot(char *name, struct Prot_input *Prot, int *index, int N);
float **Read_function(char *file_fun, struct Prot_input *Prot, int *index, int N);

// Auxiliary
void Set_contact_type();
int Count_AA(char *seq, int N_ali);
int Seq_differences(int *id, char *seq1, char *seq2, int N_ali);
int Count_gaps(char *seq1, char *seq2, int N_ali);
float Seqid(float id, int L1, int L2);
float Min_dist(int i, int j, int **conformation, int *N_conf, float **div);
float Min_CV(int i, int j, int k, int **conformation, int *N_conf, float **div);
void Change_conformations(int **conformation, int N_seq, int *N_conf, int *ali_str);
void Get_file_name(char *name, char *file);

// Eliminate pairs with SI<SI_thr=2*S0=0.12
float SI_thr=0.10; //0.12


/*********************************************************************
                          MAIN routine
**********************************************************************/
int main(int argc, char **argv)
{
  // INPUT
  char PDB_DIR[100]="./", PDB_EXT[10]="", OUTG[10]="TN";
  char file_ali[200], file_ali_str[200]="\0", file_fun[200], name_in[80]="";
  Get_input(file_ali, file_ali_str, file_fun, name_in,
	    PDB_DIR, PDB_EXT, OUTG, &PRINT_SIM, &PRINT_CV, &PRINT_DIV, &PRINT_DIFF,
	    argc, argv);

  //strcpy(file_ali, argv[1]);
  //int N_prot=Get_alignment_old(&Prot_in, &N_ali, PDB_DIR, PDB_EXT,
  //			       &PRINT_SIM, file_ali);

  // Alignments
  int N_ali=0, N_ali_str=0, i, j;
  struct Prot_input *Prot_in, *Prot2;
  int N_prot=Get_alignment(&Prot_in, &N_ali, file_ali, PDB_DIR, PDB_EXT);
  // Read structure alignment, if any
  int Np2=Get_alignment(&Prot2, &N_ali_str, file_ali_str, PDB_DIR, PDB_EXT);
  if(Np2 && (Np2!=N_prot)){
    printf("WARNING, structure alignent contains %d proteins instead of %d\n",
	   Np2, N_prot);
    printf("Discarding structure alignment\n"); Np2=0;
  }
  int *ali_str=NULL;
  if(Np2){
    printf("Matching multiple sequence and multiple structure alignments\n");
    ali_str=Match_alignments(Prot_in, Prot2, N_prot);
  }

  // What to print ? 
  if((PRINT_CV==0)&&(PRINT_SIM==0)&&((PRINT_DIFF==0)||(ali_str==NULL)))
    PRINT_DIV=1;


  /**************   READ PROTEIN STRUCTURES  ******************/
  // Read PDB files and compute contact matrices
  Set_contact_type();
  printf("Contact type: %c Threshold: %.2f A |i-j|>%d\n",
	 CONT_TYPE, CONT_THR,IJ_MIN);
  struct protein prots[N_prot], *prot=prots;
  int N_pdb=0, index[N_prot];
  int **Prot_ali=Allocate_mat2_i(N_prot, N_ali);
  int **Prot_ali_str=NULL;
  if(ali_str)Prot_ali_str=Allocate_mat2_i(N_prot, N_ali_str);
  for(i=0; i<N_prot; i++){
    char pdb[80];
    sprintf(pdb, "%s%s", Prot_in[i].name, PDB_EXT);
    if(Read_PDB_compress(prot, pdb, &(Prot_in[i].chain), PDB_DIR)>0){
      if(Align_seq(Prot_ali[N_pdb], N_ali,
		   Prot_in[i].seq, prot->aseq, prot->len)<0)continue;
      if(Prot_ali_str && 
	 (Align_seq(Prot_ali_str[N_pdb], N_ali_str,
		    Prot2[ali_str[i]].seq, prot->aseq, prot->len)<0))ali_str[i]=-1;
      index[N_pdb]=i; N_pdb++;
      int NC=Compute_contact_list(prot, CONT_TYPE, CONT_THR, IJ_MIN);
      printf("%d contacts\n", NC); prot++;
    }
  }
  printf("%d proteins read out of %d listed in %s\n", N_pdb, N_prot, file_ali);
  if(N_pdb<2){
    printf("ERROR, fewer than 2 proteins found\n"); exit(8);
  }
  float **fun_sim=Read_function(file_fun, Prot_in, index, N_pdb);

  // Prepare output
  char name_sim[100]; FILE *file_sim=NULL;
  if(PRINT_SIM){
    Change_ext(name_sim, name_in, EXT_SIM);
    file_sim=fopen(name_sim, "w");
    fprintf(file_sim, "# Cont_overlap and TM_score obtained with multiple ");
    if(Prot_ali_str){fprintf(file_sim, "structure alignments\n");}
    else{fprintf(file_sim, "sequence alignments\n");}
    fprintf(file_sim, "#Prot1 Prot2 Seq_Id Cont_Overlap TM_Score\n");
  }
  char name_out[100]; FILE *file_out=NULL;
  if(PRINT_DIV){
    Change_ext(name_out, name_in, EXT_DIV);
    file_out=fopen(name_out, "w");
    fprintf(file_out, "# Cont_Div and TM_Div obtained with multiple ");
    if(Prot_ali_str){fprintf(file_out, "structure alignments\n");}
    else{fprintf(file_out, "sequence alignments\n");}
    fprintf(file_out, "#Prot1 Prot2 Tajima-Nei_Div Cont_Div TM_Div\n");
  }
  char name_diff[100]; FILE *file_diff=NULL;
  if(PRINT_DIFF && Prot_ali_str){
    Change_ext(name_diff, name_in, EXT_DIFF);
    file_diff=fopen(name_diff, "w");
    fprintf(file_diff, "#Prot1 Prot2 TN_Div_SqA TN_Div_StA ");
    fprintf(file_diff, " Cont_Div_SqA Cont_Div_StA ");
    fprintf(file_diff, " TM_Div_SqA TM_Div_StA\n");
  }

  // Pairwise computations only for i>j
  int L_seq[N_pdb];
  float **Seq_diff_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  float **SI_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  float **TN_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);;
  float **Cont_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  float **TM_Div_PDB=Allocate_mat2_f(N_pdb, N_pdb);
  float overlap=0, overlap_SqA=0, overlap_StA=0;
  float TM=0, TM_SqA=0, TM_StA=0, SI_StA=0, D_TM=0;
  printf("Computing pairwise similarities\n");
  int al2i=-1, al2j=-1;
  for(i=0; i<N_pdb; i++){
    int al1i=index[i]; if(ali_str)al2i=ali_str[al1i];
    L_seq[i]=Count_AA(Prot_in[al1i].seq, N_ali);
    struct protein *proti=prots+i, *protj=prots;
    for(j=0; j<i; j++){
      int homo, id, al1j=index[j]; if(ali_str)al2j=ali_str[al1j];

      // Similarities
      Seq_diff_PDB[i][j]=
	Seq_differences(&id, Prot_in[al1i].seq, Prot_in[al1j].seq, N_ali);
      SI_PDB[i][j]=Seqid(id, L_seq[i], L_seq[j]);
      //if((Prot_ali_str==0) || file_diff){
      overlap_SqA=
	Compute_overlap(proti->Cont_map, proti->len, Prot_ali[i],
			protj->Cont_map, protj->len, Prot_ali[j], N_ali);
      TM_SqA=
	TM_score(proti->xca, Prot_ali[i], proti->len,
		 protj->xca, Prot_ali[j], protj->len, N_ali);

      if((al2i>=0) && (al2j>=0)){
      
	overlap_StA=
	  Compute_overlap(proti->Cont_map, proti->len, Prot_ali_str[i],
			  protj->Cont_map, protj->len, Prot_ali_str[j],
			  N_ali_str);
	TM_StA=TM_score(proti->xca, Prot_ali_str[i], proti->len,
			protj->xca, Prot_ali_str[j], protj->len,N_ali_str);
	Seq_differences(&id, Prot2[al2i].seq, Prot2[al2j].seq, N_ali_str);
	SI_StA=Seqid(id, L_seq[i], L_seq[j]);
      }

      if(overlap_StA>overlap_SqA){overlap=overlap_StA;}
      else{overlap=overlap_SqA;}
      if(TM_StA > TM_SqA){TM=TM_StA;}
      else{TM=TM_SqA;}

      // Divergences
      if(SI_PDB[i][j]>S0){
	TN_Div_PDB[i][j]=-log((SI_PDB[i][j]-S0)/(1.-S0));
      }else{
	TN_Div_PDB[i][j]=HUGE;
      }
      Cont_Div_PDB[i][j]=
	Compute_Dcont(overlap, proti->len, protj->len, &homo);
      if(TM){D_TM=-log(TM);}else{D_TM=HUGE;}
      TM_Div_PDB[i][j]=D_TM;

      // Print
      if(file_sim){
	fprintf(file_sim, "%s\t%s\t%.3f\t%.3f\t%.3f\n",
		proti->name, protj->name, SI_PDB[i][j], overlap, TM);
      }
      if(file_out){
	fprintf(file_out, "%s\t%s\t%.3f\t%.3f\t%.3f\n",
		proti->name, protj->name,
		TN_Div_PDB[i][j], Cont_Div_PDB[i][j], TM_Div_PDB[i][j]);
      }
      if(Prot_ali_str && file_diff){
	float TN_div_StA;
	if(SI_StA>S0){TN_div_StA=-log((SI_StA-S0)/(1.-S0));}
	else{TN_div_StA=HUGE;}
	fprintf(file_diff, "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
		proti->name, protj->name,
		TN_Div_PDB[i][j], TN_div_StA,
		//SI_PDB[i][j], SI_StA,
		Compute_Dcont(overlap_SqA, proti->len, protj->len, &homo),
		Compute_Dcont(overlap_StA, proti->len, protj->len, &homo),
		//overlap_SqA, overlap_StA, 
		-log(TM_SqA), -log(TM_StA)); 
		//TM_SqA, TM_StA); 
      }
      protj++;
    }
  }
  // Symmetrize
  for(i=0; i<N_pdb; i++){
    for(j=0; j<i; j++){
      Cont_Div_PDB[j][i]=Cont_Div_PDB[i][j];
      TM_Div_PDB[j][i]=TM_Div_PDB[i][j];
      TN_Div_PDB[j][i]=TN_Div_PDB[i][j];
      SI_PDB[j][i]=SI_PDB[i][j];
    }
  }
  printf("End pairwise computations\n");

  if(file_sim){
    printf("Similarities written in file %s\n", name_sim);
    fclose(file_sim);
  }
  if(file_out){
    printf("Divergences written in file %s\n", name_out);
    fclose(file_out);
  }
  if(file_diff){
    printf("Alignment scores written in file %s\n", name_diff);
    fclose(file_diff);
  }

  // Group identical sequences
  // Structural divergence: minimum among all conformations 
  printf("Grouping identical sequences by single linkage\n");
  int **conformation, *N_conf, *rep_str, *seq_clus;
  int N_seq=
    Single_linkage(&conformation, &N_conf, &rep_str, &seq_clus,
		   Seq_diff_PDB, N_pdb, diff_thr);
  //if(ali_str)Change_conformations(conformation, N_seq, N_conf, ali_str);
  float **Seq_Id=Allocate_mat2_f(N_seq, N_seq);
  float **TN_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **TM_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **Cont_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  float **p_Div_Seq=NULL, **Poiss_Div_Seq=NULL;
  int POISS=0;
  if(POISS){
    p_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
    Poiss_Div_Seq=Allocate_mat2_f(N_seq, N_seq);
  }

  for(i=0; i<N_seq; i++){
    for(j=0; j<i; j++){
      // Look for minimum among conformations
      float TN=1000, TM=1000, CD=1000, SI;
      for(int ki=0; ki<N_conf[i]; ki++){
	int ci=conformation[i][ki];
	for(int kj=0; kj<N_conf[j]; kj++){
	  int cj=conformation[j][kj];
	  if((TN_Div_PDB[ci][cj]<TN)){
	    TN=TN_Div_PDB[ci][cj]; SI=SI_PDB[ci][cj];
	  }
	  if((TM_Div_PDB[ci][cj]<TM))TM=TM_Div_PDB[ci][cj];
	  if((Cont_Div_PDB[ci][cj]<CD))CD=Cont_Div_PDB[ci][cj];
	}
      }
      Seq_Id[i][j]=SI;
      TN_Div_Seq[i][j]=TN;
      TM_Div_Seq[i][j]=TM;
      Cont_Div_Seq[i][j]=CD;
      Seq_Id[j][i]=SI;
      TN_Div_Seq[j][i]=TN;
      TM_Div_Seq[j][i]=TM;
      Cont_Div_Seq[j][i]=CD;
      if(POISS){
	p_Div_Seq[i][j]=1-SI;
	Poiss_Div_Seq[i][j]=-log(SI);
	p_Div_Seq[j][i]=p_Div_Seq[i][j];
	Poiss_Div_Seq[j][i]=Poiss_Div_Seq[i][j];
      }
    }
  }
  Empty_matrix_f(SI_PDB, N_pdb);
  Empty_matrix_f(Seq_diff_PDB, N_pdb);
  //Empty_matrix_f(TN_Div_PDB, N_pdb);

  // Function similarity for genes
  float **fun_sim_Seq=NULL;
  int n_low=0, n_high=0;
  if(fun_sim){
    fun_sim_Seq=malloc(N_seq*sizeof(int *));
    for(i=0; i<N_seq; i++){
      fun_sim_Seq[i]=malloc(N_seq*sizeof(int));
      for(j=0; j<N_seq; j++)fun_sim_Seq[i][j]=-1;
    }
    for(i=0; i<N_pdb; i++){
      int i1=seq_clus[i];
      for(j=0; j<i; j++){
	if(fun_sim[i][j]<0)continue;
	int j1=seq_clus[j];
	float s=fun_sim[i][j];
	fun_sim_Seq[i1][j1]=s;
	fun_sim_Seq[j1][i1]=s;
	if(s<=fun_low)n_low++;
	if(s>=fun_high)n_high++;
      }
    }
    Empty_matrix_f(fun_sim, N_pdb);
    printf("%d pairs with function similarity <= %.3f\n",
	   n_low, fun_low);
    printf("%d pairs with function similarity >= %.3f\n",
	   n_high, fun_high);
  }

  // Neighbor Joining and outgroups
  /* struct treenode *nodes=Neighbor_Joining(TajNei_Div, N_seq);
  int **N_out, ***outgroup=Build_outgroups(&N_out, nodes, N_seq); */
  float **Out_Div=TN_Div_Seq;
  if(strcmp(OUTG,"CD")==0){Out_Div=Cont_Div_Seq;}
  else if(strcmp(OUTG,"TM")==0){Out_Div=TM_Div_Seq;}
  printf("Determining outgroups through neighbor joining ");
  printf("with divergence %s\n", OUTG);
  int **N_out, ***outgroup=Outgroups_NJ(&N_out, Out_Div, N_seq);

  /*********************** Clock violations ***********************/
  int Nd=3, dist; // types of distances
  if(POISS){Nd=5;}
  char *name_dist[Nd]; float **Div[Nd], **Div_PDB[Nd];
  for(dist=0; dist<Nd; dist++)name_dist[dist]=malloc(80*sizeof(char));
  name_dist[0]="Tajima-Nei"; Div[0]=TN_Div_Seq;   Div_PDB[0]=TN_Div_PDB;   
  name_dist[1]="Cont_Div";   Div[1]=Cont_Div_Seq; Div_PDB[1]=Cont_Div_PDB;
  name_dist[2]="TM_Div";     Div[2]=TM_Div_Seq;   Div_PDB[2]=TM_Div_PDB;
  if(POISS){
    name_dist[3]="p_Div"; Div[3]=p_Div_Seq;
    name_dist[4]="Poiss_Div"; Div[4]=Poiss_Div_Seq;
  }


  // Prepare output
  char head[5000], ALI[100];
  sprintf(head, "# %d sequences aligned with %d PDB files\n", N_prot, N_pdb);
  sprintf(head, "%s# %d groups of sequences with < %d intragroup substitutions\n",
	  head, N_seq, diff_thr);
  sprintf(head, "%s# Average number of structures per group: %.1f\n",
	  head, N_pdb/(float)N_seq);
  sprintf(head, "%s# Alignment: multiple sequence", head);
  if(ali_str)sprintf(head, "%s, multiple structure",head);
  sprintf(head, "%s\n# Outgroups assigned with NJ using divergence %s\n",
	  head, OUTG);
  if(ali_str){strcpy(ALI,"StAli");}else{strcpy(ALI,"SqAli");}

  Change_ext(name_out, name_in, EXT_CV);
  file_out=fopen(name_out, "w");
  fprintf(file_out, "# Output of the program %s\n", argv[0]);
  fprintf(file_out, "%s", head);
  fprintf(file_out, "# Prot1, Prot2, number_of_outgroups_C\n");
  fprintf(file_out, "# For each divergence measure d, plot:\n");
  fprintf(file_out, "# d(A,B)\n");
  fprintf(file_out, "# CV(A,B)=sum_C(d(A,C)-d(B,C))/nC*d(A,B)\n"); //^%.2f,alpha
  fprintf(file_out, "# t=|CV(A,B)|/S.E.M.(CV) (Standard Error of Mean)\n");
  fprintf(file_out, "# n_TI number of violations of Triangle Inequality\n");
  fprintf(file_out, "# number of outgroups with minority sign\n");
  fprintf(file_out, "# Note: used outgroups = num_out-num_TI\n");
  fprintf(file_out, "#Prot1 Prot2 L1-L2 num_out");
  int k=5;
  for(dist=0; dist<Nd; dist++){
    fprintf(file_out, " %s:", name_dist[dist]);
    fprintf(file_out, " %d=d %d=CV %d=t_CV %d=n_TI %d=n_sign",
	    k, k+1, k+2, k+3, k+4); k+=5;
  }
  fprintf(file_out, "\n");

  // Compute CV
  float **CV_dist[Nd], **t_dist[Nd];
  int **nout_dist[Nd], **nsign_dist[Nd], **nTIV_dist[Nd];
  if(PRINT_AVE){
    for(dist=0; dist<Nd; dist++){
      CV_dist[dist]=Allocate_mat2_f(N_seq, N_seq);
      t_dist[dist]=Allocate_mat2_f(N_seq, N_seq);
      nout_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
      nTIV_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
      nsign_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
    }
  }

  int **N_gaps=NULL, **All_gaps, **TIV_gaps, kgap_max=40, bingap=5;
  if(PRINT_GAP){
    N_gaps=Allocate_mat2_i(N_seq, N_seq);
    All_gaps=Allocate_mat2_i(Nd, kgap_max+1);
    TIV_gaps=Allocate_mat2_i(Nd, kgap_max+1);
    for(i=0; i<N_seq; i++){
      int i1=rep_str[i];
      for(j=0; j<i; j++){
	int j1=rep_str[j];
	N_gaps[i][j]=Count_gaps(Prot_in[i1].seq, Prot_in[j1].seq, N_ali);
	N_gaps[j][i]=N_gaps[i][j];
      }
    }
  }

  int SI_low=0, No_out[Nd]; long Num_out=0;
  for(i=0; i<Nd; i++)No_out[i]=0;

  // Sum over pairs of sequences
  for(i=0; i<N_seq; i++){
    for(j=0; j<i; j++){
      Num_out+=N_out[i][j];
      if(Seq_Id[i][j]<SI_thr){SI_low++; continue;}
      fprintf(file_out, "%s\t%s\t%d\t%d",
	      prots[rep_str[i]].name, prots[rep_str[j]].name,
	      L_seq[rep_str[i]]-L_seq[rep_str[j]], 
	      N_out[i][j]);

      int *outg=outgroup[i][j], used[N_seq];
      for(dist=0; dist<Nd; dist++){
	float **diver=Div[dist], d=diver[i][j];
	// outgroups, triangle inequality, diff.sign
	int nout=0, nTIV=0, nplus=0, kgap=-1, TIV=0; 
	double CV1=0, CV2=0, t=0, nindep=0; // norm=pow(d,alpha);
	for(int k1=0; k1<N_out[i][j]; k1++){
	  // Eliminate outgroups that are very far away
	  int k=outg[k1]; float CV;
	  if((Seq_Id[i][k]<SI_thr)||(Seq_Id[j][k]<SI_thr))continue;
	  if(dist<3){
	    CV=Min_CV(i, j, k, conformation, N_conf, Div_PDB[dist]);
	  }else{
	    CV=diver[i][k]-diver[j][k];
	  }
	  if((CV>d)||(CV<-d)){ // check triangle inequality TI
	    used[k1]=0; TIV=1; nTIV++; 
	  }else{
	    used[k1]=1; if(TIV)TIV=0;
	    CV1+=CV; CV2+=CV*CV; nout++; if(CV>0)nplus++;
	    float s_max=0, *Sk=Seq_Id[k];
	    for(int k2=0; k2<k1; k2++){
	      if((used[k2])&&(Sk[outg[k2]]>s_max))s_max=Sk[outg[k2]];
	    }
	    nindep+=s_max;
	  }
	  if(N_gaps){
	    kgap=N_gaps[i][j]+N_gaps[i][k]+N_gaps[j][k];
	    kgap/=bingap; if(kgap>kgap_max)kgap=kgap_max;
	    All_gaps[dist][kgap]++;
	    if(TIV)TIV_gaps[dist][kgap]++;
	  }
	} // end sum over outgroups

	if(nout){CV1/=nout;}
	else{No_out[dist]++;}
	nindep=nout-nindep;
	if(nout <= 1){t=1;}
	else{
	  CV2=(CV2-nout*CV1*CV1)/(nout-1);
	  if(CV2<=0){t=1;}
	  else{t=fabs(CV1)/sqrt(CV2/nindep);}
	}
	CV1/=d;

	int nm=nout-nplus; if(nm<nplus)nplus=nm;
	fprintf(file_out, "\t%.3f\t%.3f\t%.1f\t%d\t%d",
		d, CV1, t, nTIV, nplus); //norm
	if(PRINT_AVE){
	  CV_dist[dist][i][j]=CV1;
	  t_dist[dist][i][j]=t;
	  nout_dist[dist][i][j]=nout;
	  nsign_dist[dist][i][j]=nplus;
	  nTIV_dist[dist][i][j]=nTIV;
	}
      } // end dists
      fprintf(file_out, "\n");
    }
  } // end pairs
  fclose(file_out);

  if(PRINT_AVE){
    sprintf(head, "%s# Total number of outgroups: %ld\n", head, Num_out);
    sprintf(head, "%s# %d pairs with SI < %.2f\n", head, SI_low, SI_thr);
    sprintf(head, "%s# pairs with zero outgroups: ", head); 
    for(i=0; i<Nd; i++){
      sprintf(head, "%s %d (%s)", head, No_out[i], name_dist[i]);
    }
    sprintf(head, "%s\n", head);

    for(dist=0; dist<Nd; dist++){
      CV_statistics(name_in, head, OUTG, ALI,
		    name_dist[dist], Div[dist], CV_dist[dist],
		    t_dist[dist], nout_dist[dist],
		    nsign_dist[dist], nTIV_dist[dist], Seq_Id, SI_thr,
		    NULL, 0.00, 1.00, N_seq, dist, "AllFun");
      if(fun_sim_Seq){
	CV_statistics(name_in, head, OUTG, ALI,
		      name_dist[dist], Div[dist], CV_dist[dist],
		      t_dist[dist], nout_dist[dist],
		      nsign_dist[dist], nTIV_dist[dist], Seq_Id, SI_thr,
		      fun_sim_Seq, fun_high, 1.00, N_seq, dist, "SameFun");
	CV_statistics(name_in, head, OUTG, ALI,
		      name_dist[dist], Div[dist], CV_dist[dist],
		      t_dist[dist], nout_dist[dist],
		      nsign_dist[dist], nTIV_dist[dist], Seq_Id, SI_thr,
		      fun_sim_Seq, 0.00, fun_low, N_seq, dist, "DiffFun");
      }

      Empty_matrix_f(CV_dist[dist], N_seq);
      Empty_matrix_f(t_dist[dist], N_seq);
      Empty_matrix_i(nout_dist[dist], N_seq);
      Empty_matrix_i(nTIV_dist[dist], N_seq);
      Empty_matrix_i(nsign_dist[dist], N_seq);
    }
  }

  /*****************************************************************************/
  /* Relation between gaps and triangle inequality */
  if(PRINT_GAP){
    Change_ext(name_out, name_in, ".gaps");
    file_out=fopen(name_out, "w");
    for(int dist=0; dist<Nd; dist++){
      fprintf(file_out, "# dist=%s\n", name_dist[dist]);
      fprintf(file_out, "#ngap P(TIV) s.e. num\n");
      double norm=0;
      for(i=0; i<=kgap_max; i++)norm+=All_gaps[dist][i];
      for(i=0; i<=kgap_max; i++){
	if(All_gaps[dist][i]==0)continue;
	float p=(float)TIV_gaps[dist][i]/All_gaps[dist][i];
	fprintf(file_out, "%.1f %.3f %.3f %.3f\n", (i+0.5)*bingap,
		p, sqrt(p*(1-p)/All_gaps[dist][i]),
		All_gaps[dist][i]/norm);
      }
    }
    fclose(file_out);
    printf("Writing %s\n", name_out);
  }


  /*****************************************************************************/
  return(0);
}

void help(char *pname){
  printf("Program %s\n", pname);
  printf("Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa ");
  printf("(CSIC-UAM), Madrid, Spain\nEmail: <ubastolla@cbm.csic.es>\n");
  printf("\n");
  printf("Given a multiple sequence alignment (MSA) of proteins with known structures and optionally a multiple structure alignment, it computes and prints sequence and structural similarity measures (sequence identity SI, contact overlap q and TM-score TM) and the corresponding divergence measures (Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=%.2f, Contact_divergence CD=-log((q-q0(L))/(1-q0(L)), TM_divergence=-log(TM)) for all aligned pairs.\n", S0);
  printf("Optionally, it computes and prints for all three divergence measures the violation of the molecular clock averaged over all possible outgroups identified with the Neighbor-Joining criterion and the corresponding significance score\n");
  printf("==========================================================\n");
  printf("RUN: %s <Config file>\n", pname);
  printf("==========================================================\n");
  printf("Config file:\n");
  printf("ALI=<MSA file in FASTA format, indicating names of PDB files>\n");
  printf("Optional parameters:\n");
  printf("STR_ALI=<MSA file in FASTA format, with names of PDB files>\n");
  printf("# (multiple structure alignment, optional).\n");
  printf("FUN_SIM=<file with function similarity for pairs of PDB files, optional>\n");
  printf("OUTGROUP=<Method to assign outgroups> ");
  printf("Allowed: TN (Tajima-Nei, default) CD (Contact divergence) TM (TM-score)\n"); 
  printf("NAME= <Name of output files> (default: alignment file)\n");
  printf("PDBDIR=<directory of pdb files>  (default: current directory)\n");
  printf("PDBEXT=<extension of pdb files>  (default: none)\n");
  printf("PRINT_SIM=<0,1>   Print similarity measures? (default: 0)\n");
  printf("PRINT_DIV=<0,1>   Print divergence measures? (default: 1)\n");
  printf("PRINT_CV=<0,1>    Print clock violations? (default: 1)\n");
  printf("PRINT_DIFF=<0,1>  Compare seq and str alignments if both present? ");
  printf("  (default: 0)\n");
  printf("The protein name is the name of a PDB file, optionally followed\n");
  printf("by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A)\n\n");
  printf("==========================================================\n");
  printf("OUTPUT (for each pair of proteins):\n");
  printf("File with extension .sim (similarity):\n");
  printf("Sequence identity SI (sequence alignment)\n");
  printf("Contact overlap q (structural, str.ali if present)\n");
  printf("TM-score TM (structural, str.ali if present), Zhang & Skolnick Proteins 2004 57:702\n");
  printf("\n");
  printf("File with extension .div (divergence):\n");
  printf("Tajima-Nei divergence TN=-log((SI-S0)/(1-S0) S0=%.2f (Tajima F & Nei 1984, Mol Biol Evol 1:269).\n", S0);
  printf("Contact divergence  CD=-log((q-q0(L))/(1-q0(L)) (Pascual-García et al Proteins 2010 78:181-96)\n");
  printf("TM_divergence=-log(TM)\n");
  printf("\n");
  printf("File with extension .cv (clock violations):\n");
  printf("CV=sum_c (D(a,c)-D(b,c))/(n_c*D(a,b)) (Pascual-Garcia, Arenas & Bastolla, submitted).\n");
  printf("\n");
  exit(8);
}

int Get_alignment_old(struct Prot_input **Prot_input, int *Nali,
		      char *PDB_DIR, char *PDB_EXT, int *PRINT_SIM,
		      char *file_ali)
{
  // Open file
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", file_ali); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int dir=0, n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>'){
      n++;
    }else if(n==0){
      if(strncmp(string, "PDBDIR", 6)==0){
	sscanf(string+7,"%s", PDB_DIR);
	printf("Directory for PDB files: %s\n", PDB_DIR); dir=1;
      }else if(strncmp(string, "PDBEXT", 6)==0){
	sscanf(string+7, "%s", PDB_EXT);
      }else if(strncmp(string, "PRINT_SIM", 8)==0){
	sscanf(string+9, "%d", PRINT_SIM);
      }
    }
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in file %s\n", file_ali); exit(8);
  }
  printf("%d sequences found in %s\n", n, file_ali);

  // Allocate and read sequences
  int LMAX=10000, l=0, i;
  char chain[10], dumm[40];
  char *Seq=malloc(LMAX*sizeof(char)), *s=NULL;
  *Prot_input=malloc(n*sizeof(struct Prot_input));
  n=-1; *Nali=0;
  file_in=fopen(file_ali, "r");
  if(dir)fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(string[0]=='>'){
      n++;
      sscanf(string+1, "%s", (*Prot_input)[n].name);
      for(i=0; i<10; i++)chain[i]='\0';
      int c=sscanf(string, "%s%s\n", dumm, chain);
      if((c>1)&&(chain[0]!='\0')&&(chain[0]!='\n')){
	(*Prot_input)[n].chain=chain[0];
      }else if(string[5]=='_'){
	printf("Getting chain after _\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[6];
      }else if((string[5]>=65)&&(string[5]<=90)){ // Maiuscule
	printf("Getting chain after character 4\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[5];
      }else{
	(*Prot_input)[n].chain=' ';
      }
      printf("%s %c\n", (*Prot_input)[n].name, (*Prot_input)[n].chain);
      if((*Nali==0)&&(l)){
	*Nali=l;
	(*Prot_input)[0].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[0].seq;
	for(l=0; l<*Nali; l++){*s=Seq[l]; s++;}
      }
      if(*Nali){
	(*Prot_input)[n].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[n].seq; l=0;
      }else{
	s=Seq; l=0;
      }
    }else if(n>=0){
      char *c=string;
      while(*c!='\n'){*s=*c; l++; s++; c++;}
      if(l > LMAX){
	printf("ERROR, alignment length larger than maximum allowed %d\n", l);
	printf("Increase LMAX in code %s\n", CODENAME); exit(8);
      }
      if((*Nali)&&(l>*Nali)){
	printf("ERROR, too many column in alignment %d.",n+1);
	printf(" Expected %d, found >= %d\n", *Nali, l); exit(8); 
      }
    }
  }
  n++;
  fclose(file_in);
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);
  return(n);
}

int Get_alignment(struct Prot_input **Prot_input, int *Nali,
		  char *file_ali, char *PDB_DIR, char *PDB_EXT)
{
  // Open file
  if(file_ali[0]=='\0')return(0);
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", file_ali); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>')n++;
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in file %s\n", file_ali); exit(8);
  }
  printf("%d sequences found in %s\n", n, file_ali);

  // Allocate and read sequences
  int LMAX=10000, l=0, i;
  char chain[10], dumm[40];
  char *Seq=malloc(LMAX*sizeof(char)), *s=NULL;
  *Prot_input=malloc(n*sizeof(struct Prot_input));
  n=-1; *Nali=0;
  file_in=fopen(file_ali, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(string[0]=='>'){
      n++;
      sscanf(string+1, "%s", (*Prot_input)[n].name);
      for(i=0; i<10; i++)chain[i]='\0';
      int c=sscanf(string, "%s%s\n", dumm, chain);
      if((c>1)&&(chain[0]!='\0')&&(chain[0]!='\n')){
	(*Prot_input)[n].chain=chain[0];
      }else if(string[5]=='_'){
	printf("Getting chain after _\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[6];
      }else if((string[5]>=65)&&(string[5]<=90)){ // Maiuscule
	printf("Getting chain after character 4\n");
	(*Prot_input)[n].name[4]='\0';
	(*Prot_input)[n].chain=string[5];
      }else{
	(*Prot_input)[n].chain=' ';
      }
      printf("%s %c\n", (*Prot_input)[n].name, (*Prot_input)[n].chain);
      if((*Nali==0)&&(l)){
	*Nali=l;
	(*Prot_input)[0].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[0].seq;
	for(l=0; l<*Nali; l++){*s=Seq[l]; s++;}
      }
      if(*Nali){
	(*Prot_input)[n].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[n].seq; l=0;
      }else{
	s=Seq; l=0;
      }
    }else if(n>=0){
      char *c=string;
      while(*c!='\n'){*s=*c; l++; s++; c++;}
      if(l > LMAX){
	printf("ERROR, alignment length larger than maximum allowed %d\n", l);
	printf("Increase LMAX in code %s\n", CODENAME); exit(8);
      }
      if((*Nali)&&(l>*Nali)){
	printf("ERROR, too many column in alignment %d.",n+1);
	printf(" Expected %d, found >= %d\n", *Nali, l); exit(8); 
      }
    }
  }
  n++;
  fclose(file_in);
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);
  return(n);
}


void Get_file_name(char *name, char *file){
  char *s=file, *t=name;
  while(*s!='\0'){
    if(*s=='/'){s++; t=name;}
    else if(*s=='.'){break;}
    *t=*s; s++;
  }
}

void Set_contact_type(){

  if(CONT_TYPE=='a'){strcpy(CONT_STRING, "Alpha");}
  else if(CONT_TYPE=='b'){strcpy(CONT_STRING, "Beta");}
  else if(CONT_TYPE=='c'){strcpy(CONT_STRING, "All atoms");}
  else{
    printf("WARNING, undefined contact %c\n", CONT_TYPE);
    CONT_TYPE='c'; strcpy(CONT_STRING, "All atoms");
    printf("Using default %s\n", CONT_STRING);
  }
  // Default type of contacts?
  if((CONT_TYPE!=CONT_TYPE_DEF)||(CONT_THR!=CONT_THR_DEF)||
     (IJ_MIN!=IJ_MIN_DEF))CONT_DEF=0;
}

int Count_AA(char *seq, int N_ali){
  int L=0; char *s=seq;
  for(int i=0; i<N_ali; i++){if(*s!='-')L++; s++;}
  return(L);
}

int Seq_differences(int *id, char *seq1, char *seq2, int N_ali){
  int d=0; (*id)=0; char *s1=seq1, *s2=seq2;
  for(int i=0; i<N_ali; i++){
    if((*s1!='-')&&(*s2!='-')){
      if(*s1!=*s2){d++;}else{(*id)++;}
    }
    s1++; s2++;
  }
  return(d);
}

float Seqid(float id, int L1, int L2){
  // Normalize by the longer protein
  if(L1>L2){return(id/L1);}
  else{return(id/L2);}
}

int Count_gaps(char *seq1, char *seq2, int N_ali){
  int ngap=0, open=0; char *s1=seq1, *s2=seq2;
  for(int i=0; i<N_ali; i++){
    if((*s1!='-')&&(*s2!='-')){open=0;}
    else if(open==0){ngap++; open=1;}
    s1++; s2++;
  }
  return(ngap);
}

float Min_CV(int i, int j, int k, int **conformation, int *N_conf, float **div)
{
  float CV=10000;
  for(int k1=0; k1<N_conf[k]; k1++){
    float *dk=div[conformation[k][k1]];
    for(int i1=0; i1<N_conf[i]; i1++){
      float d_ki=dk[conformation[i][i1]];
      for(int j1=0; j1<N_conf[j]; j1++){
	float cv=d_ki-dk[conformation[j][j1]];
	if(fabs(cv)<fabs(CV)){CV=cv;}
      }
    }
  }
  return(CV);
}

void Change_conformations(int **conformation, int N_seq, int *N_conf, int *ali_str)
{
  for(int i=0; i<N_seq; i++){
    for(int k=0; k<N_conf[i]; k++){
      int l=conformation[i][k];
      conformation[i][k]=ali_str[l];
    }
  }
}

float Min_dist(int i, int j, int **conformation, int *N_conf, float **div)
{
  float d=1000; int ini=1;
  for(int i1=0; i1<N_conf[i]; i1++){
    float *d_i=div[conformation[i][i1]];
    for(int j1=0; j1<N_conf[j]; j1++){
      float d_il=d_i[conformation[j][j1]];
      if(ini){d=d_il; ini=0;}
      else if(d_il<d){d=d_il;}
    }
  }
  return(d);
}

int *Match_alignments(struct Prot_input *Prot1,
		      struct Prot_input *Prot2, int N)
{
  int *ali_str=malloc(N*sizeof(int));
  struct Prot_input *P1=Prot1;
  for(int i=0; i<N; i++){
    ali_str[i]=-1;
    for(int j=0; j<N; j++){
      if((strcmp(Prot2[j].name, P1->name)==0)&&
	 (Prot2[j].chain==P1->chain)){
	ali_str[i]=j; break;
      }
    }
    if(ali_str[i]<0){
      printf("WARNING, protein %s%c i=%d N=%d not matched\n",
	     P1->name, P1->chain, i, N);
      free(ali_str); return(NULL);
    } 
    P1++;
  }

  return(ali_str);
}

void Get_input(char *file_ali, char *file_ali_str, char *file_fun,
	       char *name_in,  char *PDB_DIR, char *PDB_EXT, char *OUTG,
	       int *PRINT_SIM, int *PRINT_CV, int *PRINT_DIV, int *PRINT_DIFF,
	       int argc, char **argv)
{
  printf("Starting %s\n", argv[0]);
  if((argc<2)||(strncmp(argv[1], "-h", 2)==0))help(argv[0]);

  FILE *file_in=fopen(argv[1], "r");
  if(file_in==NULL){
    printf("ERROR, input file %s does not exist\n\n", argv[1]);
    exit(8);
  }
  char string[100], read[80]; int ali=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "ALI=", 4)==0){
      sscanf(string+4,"%s", file_ali);
      printf("File with multiple sequence alignment: %s\n", file_ali);
      ali=1;
    }else if(strncmp(string, "STR_ALI=", 8)==0){
      sscanf(string+8,"%s", file_ali_str);
      printf("File with multiple structure alignment: %s\n", file_ali_str);
    }else if(strncmp(string, "NAME", 4)==0){
      sscanf(string+5, "%s", name_in);
    }else if(strncmp(string, "FUN_SIM", 7)==0){
      sscanf(string+8, "%s", file_fun);
    }else if(strncmp(string, "OUTGROUP", 8)==0){
      sscanf(string+9, "%s", read);
      if((strcmp(read,"TN")!=0) && 
	 (strcmp(read,"CD")!=0) && 
	 (strcmp(read,"TM")!=0)){
	printf("WARNING, outgroup method %s not implemented, using default\n",
	       read);
      }else{
	strcpy(OUTG, read);
      }
      printf("Outgroup method set to %s\n", OUTG);
    }else if(strncmp(string, "PDBDIR", 6)==0){
      sscanf(string+7,"%s", PDB_DIR);
      printf("Directory for PDB files: %s\n", PDB_DIR);
    }else if(strncmp(string, "PDBEXT", 6)==0){
      sscanf(string+7, "%s", PDB_EXT);
    }else if(strncmp(string, "PRINT_SIM", 9)==0){
      sscanf(string+10, "%d", PRINT_SIM);
    }else if(strncmp(string, "PRINT_CV", 8)==0){
      sscanf(string+9, "%d", PRINT_CV);
    }else if(strncmp(string, "PRINT_DIV", 9)==0){
      sscanf(string+10, "%d", PRINT_DIV);
    }else if(strncmp(string, "PRINT_DIFF", 10)==0){
      sscanf(string+11, "%d", PRINT_DIFF);
    }else{
      printf("WARNING, unknown command %s", string);
    }
  }
  if(ali==0){
    printf("ERROR, file with alignment not provided in %s\n", argv[1]);
    exit(8);
  }
  if(name_in[0]=='\0')Get_file_name(name_in, file_ali);

}

float **Read_function(char *file_fun, struct Prot_input *Prot, int *index, int N)
{
  if(file_fun[0]=='\0')return(NULL);
  FILE *file_in=fopen(file_fun, "r");
  if(file_in==NULL)return(NULL);

  int n=0, i, j;
  float **f_sim=malloc(N*sizeof(float *)), sim;
  for(i=0; i<N; i++){
    f_sim[i]=malloc(N*sizeof(float));
    for(j=0; j<N; j++)f_sim[i][j]=-1;
  }
  char string[1000], name1[40], name2[40];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    sscanf(string, "%s%s%f", name1, name2, &sim);
    i=Find_prot(name1, Prot, index, N);
    j=Find_prot(name2, Prot, index, N);
    if((i<0)||(j<0)){
      printf("WARNING, pair %s %s not identified in function file %s\n",
	     name1, name2, file_fun);
    }else{
      f_sim[i][j]=sim; f_sim[j][i]=sim; n++;
    }
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no pair was identified in function file %s\n",
	   file_fun);
    Empty_matrix_f(f_sim, N);
    return(NULL);
  }else{
    printf("%d protein pairs identified in %s\n", n, file_fun);
    return(f_sim);
  }
}

int Find_prot(char *name, struct Prot_input *Prot, int *index, int N){
  for(int i=0; i<N; i++){
    int i1=index[i];
    if((strncmp(name, Prot[i1].name, 4)==0)&&(name[4]==Prot[i1].chain))
      return(i);
  }
  return(-1);
}
