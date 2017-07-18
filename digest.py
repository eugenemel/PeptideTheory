#PeptideTheory: Theoretical Ion Generation
#Author: Eugene Melamud / Calico / 2017

from Protein import *
from Peptide import *

#load protein sequence from fasta file
def readFasta(filename): 
    f = open(filename)

    proteins=list()
    seq=""
    seq_id=""
    seq_desc=""

    for line in f:
        line = line.rstrip().lstrip()
        if len(line) == 0: continue

        if line[0] == ">":
            if len(seq):
                proteins.append(Protein(seq_id,seq_desc,seq))
                seq=""
            seq_id, seq_desc = line.split(" ",1)
        else:
            seq += line.replace(" ","").upper()

    #last sequence
    if len(seq):
            proteins.append(Protein(seq_id,seq_desc,seq))
    return(proteins)



proteins = readFasta("example.fasta")

for protein in proteins:

    #
    #complete tryptic digest
    #
    peptides = protein.digest(pattern="[KR]", minlength=5, maxmiss=3, maxlength=30,missprob=0)


    #
    #partial tryptic digest missing 50% of cuts
    #
    #peptides = protein.digest(pattern="[KR]", minlength=10, maxmiss=-1, maxlength=30,missprob=0.5)

    #
    # generate theoretical spectra
    #
    for peptide_sequence in peptides:
        pep = Peptide(peptide_sequence)

        #digest
        all_ions = pep.all_ions( ionseries = ['b','y','a'], 
                frg_z_list = [1,2],
                fragmentlossgains = [0, -mNH3, -mH2O ], 
                mass_limits = [80,50000]
                )

        print protein.id, peptide_sequence, "preMz: ", pep.getMZ(2)

        #print theory digest
        for ion in all_ions: 
            print "\t", ion

