#define MAXATOM 100
#define MAXRES 10000
#define CHARPDB 20

struct atom{
  float r[3];
  char name[4];
};
struct residue{
  int n_atom;
  struct atom atom[MAXATOM];
  char aa;
  char label[6];
};
struct sec_el{
  int ini, len;
  char type;
};

struct protein{
  int len;
  int nca;
  char name_file[20];
  char name[CHARPDB];
  char chain;
  char *ss;
  char *aseq;
  int *nseq;
  int *n_atom;
  char **pdbres;
  float **xca;
  float **vec;
  struct atom **res_atom;
  char exp_meth;
  //
  int N_cont, n_sec_el;
  char  *seq2;
  short *ss_num;
  short **Cont_map;
  short *ncont;
  float *EC;
  struct sec_el *sec_el;
};

// Amino acid names
#define AANAME1 "ACDEFGHIKLMNPQRSTVWYX"
#define AANAME3 "ALACYSASPGLUPHEGLYHISILELYSLEUMETASNPROGLNARGSERTHRVALTRPTYRXXX"
