/* 
   Program Contact_divergence
   Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)
   ubastolla@cbm.csic.es
   Reads a multiple alignment and computes contact divergence and other
   structure comparison measures.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float S0=0.06;  // Sequence similarity of unrelated positions
#define EXT_DIV ".div"  // Extension for divergence output file
#define EXT_SIM ".sim"  // Extension for similarity output file


#define IJ_MIN_DEF 3
#define CONT_TYPE_DEF 'c'
#define CONT_THR_DEF 4.5
int IJ_MIN=IJ_MIN_DEF;        // Only contacts with |i-j|>=IJ_MIN
float CONT_THR=CONT_THR_DEF;
char CONT_TYPE=CONT_TYPE_DEF;  //a=alpha b=beta c=all atoms
char CONT_STRING[80];
char CODENAME[40]="main_Contact_divergence.c";
int CONT_DEF=1;

struct Prot_input{
  char name[80];
  char chain;
  char *seq;
};

void help(char *pname);
void Get_input_old(char *file_ali, int *PRINT_TM, int *PRINT_CONT,
		   int *PRINT_SEQ, int argc, char **argv);
void Get_input(char *file_ali, int *PRINT_TM, int *PRINT_CONT, int *PRINT_SEQ,
	       char *PDBDIR, char *PDBEXT, int argc, char **argv);
int Get_alignment(struct Prot_input **Prot_input, int *Nali,
		  char *PDBDIR, char *PDBEXT, int *PRINT_SIM, char *file_ali);
void Set_contact_type();


/*********************************************************************
                          MAIN routine
**********************************************************************/
int main(int argc, char **argv)
{
  // INPUT
  printf("Starting %s\n", argv[0]);
  if((argc<2)||(strncmp(argv[1], "-h", 2)==0))help(argv[0]);

  char file_ali[200], PDBDIR[100]="./", PDBEXT[10]="";
  int PRINT_SIM=0, PRINT_DIV=1, N_ali=0;
  struct Prot_input *Prot_in;
  strcpy(file_ali, argv[1]);
  int N_prot=Get_alignment(&Prot_in, &N_ali, PDBDIR, PDBEXT,
			   &PRINT_SIM, file_ali);

  /**************+++   READ PROTEINS  ********************/
  // Look for contact matrix and sequence files. If not found,
  // read coordinates from PDB files and print contact matrices
  Set_contact_type();
  printf("Contact type: %c Threshold: %.2f A |i-j|>%d\n",
	 CONT_TYPE, CONT_THR,IJ_MIN);
  struct protein prots[N_prot], *prot=prots;
  int N_pdb=0, index[N_prot], i, j;
  int **Prot_ali=Allocate_mat2_i(N_prot, N_ali);
  for(i=0; i<N_prot; i++){
    sprintf(Prot_in[i].name, "%s%s", Prot_in[i].name, PDBEXT);
    if(Read_PDB_compress(prot,Prot_in[i].name,&(Prot_in[i].chain),PDBDIR)>0){
      if(Align_seq(Prot_ali[N_pdb],N_ali,
		   Prot_in[i].seq,prot->aseq,prot->len)<0)continue;
      int NC=Compute_contact_list(prot, CONT_TYPE, CONT_THR, IJ_MIN);
      printf("%d contacts\n", NC);
      index[N_pdb]=i; prot++; N_pdb++;
    }
  }
  printf("%d proteins read out of %d listed in %s\n", N_pdb, N_prot, file_ali);
  if(N_pdb<2){
    printf("ERROR, zero protein pairs found\n"); exit(8);
  } 

  // Prepare output
  char name_sim[100], name_out[100];
  FILE *file_sim=NULL, *file_out=NULL;
  if(PRINT_SIM){
    Change_ext(name_sim, file_ali, EXT_SIM);
    file_sim=fopen(name_sim, "w");
    fprintf(file_sim, "#Prot1 Prot2 Seq_Id Cont_Overlap TM_Score\n");
  }
  if(PRINT_DIV){
    Change_ext(name_out, file_ali, EXT_DIV);
    file_out=fopen(name_out, "w");
    fprintf(file_out, "#Prot1 Prot2 Tajima-Nei_Div Cont_Div TM_Div\n");
  }

  
  // Pairwise computations
  for(i=0; i<N_pdb; i++){
    int i1=index[i];
    struct protein *proti=prots+i, *protj=prots;
    for(j=0; j<i; j++){
      // Similarities
      float SeqId=Seq_identity(Prot_in[i1].seq, Prot_in[index[j]].seq, N_ali);
      float overlap=Compute_overlap(proti->Cont_map, proti->len, Prot_ali[i],
				    protj->Cont_map, protj->len, Prot_ali[j],
				    N_ali);
      float tm=TM_score(proti->xca, Prot_ali[i], proti->len,
			protj->xca, Prot_ali[j], protj->len, N_ali);
      if(file_sim){
	fprintf(file_sim, "%s\t%s\t%.3f\t%.3f\t%.3f\n",
		proti->name, protj->name, SeqId, overlap, tm);
      }
      // Divergences
      int homo;
      float CD_Div=Compute_Dcont(overlap, proti->len, protj->len, &homo);
      float TM_Div; if(tm){TM_Div=-log(tm);}else{TM_Div=10;}
      float TN_Div=-log((SeqId-S0)/(1.-S0));
      if(file_out){
	fprintf(file_out, "%s\t%s\t%.3f\t%.3f\t%.3f\n",
		proti->name, protj->name, TN_Div, CD_Div, TM_Div);
      }
      protj++;
    }
  }
  if(file_sim){
    printf("Similarities written in file %s\n", name_sim);
    fclose(file_sim);
  }
  if(file_out){
    printf("Divergences written in file %s\n", name_out);
    fclose(file_out);
  }
  return(0);
}

void Get_input(char *file_ali, int *PRINT_TM, int *PRINT_CONT, int *PRINT_SEQ,
	       char *PDBDIR, char *PDBEXT, int argc, char **argv)
{
  if((argc<2)||(strncmp(argv[1], "-h", 2)==0))help(argv[0]);

  // Open file
  char filename[100];
  strcpy(filename, argv[1]);
  FILE *file_in=fopen(filename, "r");
  if(file_in==NULL){
    printf("ERROR, input file %s not found\n", filename); 
    exit(8);
  }

  // Read
  printf("Reading parameters in %s\n", filename);
  strcpy(file_ali, ""); char string[1000];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "ALI", 2)==0){
      sscanf(string+4, "%s", file_ali);
    }else if(strncmp(string, "PDBDIR", 6)==0){
      sscanf(string+7, "%s", PDBDIR);
    }else if(strncmp(string, "PDBEXT", 6)==0){
      sscanf(string+7, "%s", PDBEXT);
    }else if(strncmp(string, "CD", 2)==0){
      sscanf(string+3, "%d", PRINT_CONT);
    }else if(strncmp(string, "TM", 2)==0){
      sscanf(string+3, "%d", PRINT_TM);
    }else if(strncmp(string, "SEQ", 3)==0){
      sscanf(string+4, "%d", PRINT_SEQ);
    }else{
      printf("WARNING, unknown command %s", string);
    }
  }
}

void Get_input_old(char *file_ali, int *PRINT_TM, int *PRINT_CONT,
		   int *PRINT_SEQ, int argc, char **argv)
{
  if((argc<2)||(strncmp(argv[1], "-h", 2)==0))help(argv[0]);
  strcpy(file_ali, argv[1]);
  printf("Alignment file: %s\n", file_ali);
  int i;
  for(i=2; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){
      help(argv[0]);
    }else if(strncmp(argv[i], "-tm", 3)==0){
      *PRINT_TM=1;
    }else if(strncmp(argv[i], "-cont", 5)==0){
      *PRINT_CONT=1;
    }else if(strncmp(argv[i], "-seq", 5)==0){
      *PRINT_SEQ=1;
    }else{
      printf("WARNING, unrecognized option %s\n", argv[i]);
    }
  }
}

void help(char *pname){
  printf("Program %s\n", pname);
  printf("Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa ");
  printf("(CSIC-UAM), Madrid, Spain\nEmail: <ubastolla@cbm.csic.es>\n");
  printf("\n");
  printf("Given a multiple sequence alignment (MSA) of proteins with known structures, it computes and prints sequence and structural similarity measures (sequence identity SI, contact overlap q and TM-score TM) and the corresponding divergence measures (Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=%.2f, Contact_divergence CD=-log((q-q0(L))/(1-q0(L)), TM_divergence=-log(TM)) for all aligned pairs.\n", S0);
  printf("==========================================================\n");
  printf("RUN: %s <MSA file>\n", pname);
  printf("==========================================================\n");
  printf("INPUT: MSA in FASTA format, indicating names of PDB files\n");
  printf("Optional parameters in the MSA file (before first sequence):\n");
  printf("PDBDIR=<directory of pdb files>  (default: current directory)\n");
  printf("PDBEXT=<extension of pdb files>  (default: none)\n");
  printf("PRINT_SIM=<0,1>   Print similarity measures? (default: 0)\n");
  printf("PRINT_CV=<0,1>    Print clock violations? (default: 1)\n");
  printf("The protein name is the name of a PDB file, optionally followed\n");
  printf("by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A)\n\n");
  printf("==========================================================\n");
  printf("OUTPUT (for each pair of proteins):\n");
  printf("File with extension .sim (similarity):\n");
  printf("Sequence identity SI\n");
  printf("Contact overlap q (structural)\n");
  printf("TM-score TM (structural), Zhang & Skolnick Proteins 2004 57:702\n");
  printf("\n");
  printf("File with extension .div (divergence):\n");
  printf("Tajima-Nei divergence TN=-log((SI-S0)/(1-S0) S0=%.2f (Tajima F & Nei 1984, Mol Biol Evol 1:269.\n", S0);
  printf("Contact divergence  CD=-log((q-q0(L))/(1-q0(L)) (Pascual-García et al Proteins 2010 78:181-96)\n");
  printf("TM_divergence=-log(TM)\n");
  printf("\n");
  exit(8);
}


int Get_alignment(struct Prot_input **Prot_input, int *Nali,
		  char *PDB_DIR, char *PDB_EXT, int *PRINT_SIM, char *file_ali)
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
  fclose(file_in);
  n++;
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);
  return(n);
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
