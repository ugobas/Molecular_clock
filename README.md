# Molecular_clock
Software related with the paper: "**The molecular clock in the evolution of protein structures**" Pascual-García, Alberto, Miguel Arenas, and Ugo Bastolla. _Systematic Biology_ (2019). [doi: 10.1093/sysbio/syz022](https://doi.org/10.1093/sysbio/syz022)

* Program: Evol_div
* Author: Ugo Bastolla, (Centro de Biologia Molecular Severo Ochoa, CSIC-UAM, Spain)
* Email: <ubastolla@cbm.csic.es>
* Short description: Given a multiple sequence alignment (MSA) of proteins with known structures, for all aligned pairs it computes and prints sequence and structure similarity measures (.sim), divergence measures (.div) and violations of the molecular clock computed across all possible outgroups (.cv).
* License: Please see license file, and note that it includes the needlemanwunsch aligner developed by Dr. Andrew C. R. Martin in the Profit suite of programs, (c) SciTech Software 1993-2007.

Similarity measures:
--------------------
(1) Sequence identity SI,

(2) Contact overlap q,

(3) TM-score TM (Zhang & Skolnick Proteins 2004 57:702)

Divergence measures:
---------------------
(1) Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=0.06 (Tajima F & Nei 1984, _Mol Biol Evol_ 1:269),

(2) Contact_divergence CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al _Proteins_ 2010 78:181-96),

(3) TM_divergence=-log(TM).

Test of molecular clock
-----------------------
For all pairs a,b and all divergence measures, it computes violations of the molecular clock  across all outgroups c

CV=(1/Nc)sum_c (D(a,c)-D(b,c))/D(a,b) (Pascual-Garcia, Arenas & Bastolla, "The molecular clock in the evolution of protein structures", submitted)

t=|Mean(D(a,c)-D(b,c))|/S.E.M.(D(a,c)-D(b,c)) (S.E.M.=Standard Error of Mean across the different outgroups c).

Compile:
--------
>unzip Evol_div.zip (if you downloaded the whole repository as a zip file)

>make -f Evol_div.makefile

>cp Evol_div ~/bin/ (or whatever path directory you like)

RUN:
-----
>Evol_div <alignment file>
EXAMPLE: Evol_div Input_Evol_Div.in
(you have to modify the names of the pdb files and directory in Input_Cont_Div_50044.aln)

INPUT: 
------
MSA in FASTA format, indicating names of PDB files

Optional parameters in the MSA file (before the first sequence):
----------------------------------------------------------------
PDBDIR=<directory of pdb files>  (default: current directory)

PDBEXT=<extension of pdb files>  (default: none)

PRINT_SIM=<0,1>   Print similarity measures? (default: 0)

PRINT_CV=<0,1>    Print clock violations? (default: 1)

The protein name is the name of a PDB file, optionally followed by the chain index (Ex: >1opd.pdb A or >1opdA or >1opd_A)\n\n");

OUTPUT (one line for each pair of proteins):
-----------------------------------

For all pairs of proteins, it prints, for each pair: 

1. File with extension .sim (similarity):
 * Sequence identity SI
 * Contact overlap q
 * TM-score TM (structural), Zhang & Skolnick Proteins 2004 57:702

2. File with extension .div (divergence):

* "Tajima-Nei divergence TN=-log((SI-S0)/(1-S0) S0=0.06 (Tajima F & Nei 1984, _Mol Biol Evol_ 1:269).
* Contact divergence  CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garci­a et al _Proteins_ 2010 78:181-96)
* TM_divergence=-log(TM)

3. File with extension .cv (clock violations):

 * PDB names 
 * difference of length 
 * number of outgroups;
 * Three types of divergences: (Tajima-Nei, Contact divergence, -log(TM)).
 * CV=(1/Nc)sum_c (D(a,c)-D(b,c))/D(a,b) (Pascual-Garcia, Arenas & Bastolla, _submitted_).
 * t=|Mean(D(a,c)-D(b,c))|/S.E.M.(D(a,c)-D(b,c)) 
 * n_TI: Number of outgroups for which the triangular inequality is violated, i.e. D(a,c)>D(a,b)+D(b,c) (or same interchanging a with b)
 * n_sign: number of outgroups for which the sign of D(a,c)-D(b,c) is different from the sign of the mean.

