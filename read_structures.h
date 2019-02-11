int Read_structures(struct protein **prots1, char *PDB_PATH, char *FILE_LIST,
		    char *EXT_CM, char *EXT_SEQ, int CONT_DEF,
		    char CONT_TYPE, float CONT_THR, int IJ_MIN);
int Read_PDB_compress(struct protein *prot,
		      char *Prot_name, char *chain, char *PDB_PATH);
int Read_pdb(struct protein *prot, char *filename, char *chain_to_read);
