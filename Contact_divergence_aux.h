int Align_seq(int *Prot_ali, int N_ali, char *Prot_seq, char *PDB_seq, int n2);
float Seq_identity(char *seq_i, char *seq_j, int N_ali);
float Compute_overlap(short **Cont_map1, int n1, int *ali1,
		      short **Cont_map2, int n2, int *ali2,
		      int N_ali);
char Get_compression(char *file_name);
int Find_name(char *name, char **names, int N, int Np);
int Change_ext(char *name_out, char *name_in, char *ext);
