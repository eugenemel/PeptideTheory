#PeptideTheory: Theoretical Ion Generation
#Author: Eugene Melamud / Calico / 2017

from Protein import * 
from Peptide import *
import optparse

#load protein sequence from fasta file
def readFastaFile(filename): 
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
            seq_id = seq_id[1:]
        else:
            seq += line.replace(" ","").upper()

     #last sequence
    if len(seq): proteins.append(Protein(seq_id,seq_desc,seq))
    return(proteins)


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest = 'inputfile')
    parser.add_option('--tmt', action="store_true", dest = 'tmt')
    parser.add_option('--iod', action="store_true", dest = 'iod', default=True)
    parser.add_option('--reverse', action="store_true", dest = 'reverse', default=False)
    parser.add_option('--noiod', action="store_false", dest = 'iod')
    parser.add_option('--oxy',   action="store_true", dest = 'oxy', default=True)
    parser.add_option('--nooxy',   action="store_false", dest = 'oxy')
    parser.add_option('--minz',   type="int", dest='minz', default = 2)
    parser.add_option('--maxz',   type="int", dest='maxz', default = 5)
    parser.add_option('--minlen', type="int", dest='minlength', default = 4)
    parser.add_option('--maxlen', type="int", dest='maxlength', default = 30)
    parser.add_option('--maxmiss',type="int", dest='maxmiss', default = 3)
    parser.add_option('--pattern', type="string", dest='pattern', default = "[KR]")
    parser.add_option('--semiTryptic', action="store_true", dest = 'semiTryptic', default=False)

    (opts, args) = parser.parse_args()

    if not opts.inputfile:
        print "please run with -i fastafile"
        exit(-1)

    proteins = readFastaFile(opts.inputfile)

    for protein in proteins:

        #complete tryptic digest
        peptides = protein.digest(pattern=opts.pattern, minlength=opts.minlength, maxmiss=3, maxlength=opts.maxlength,missprob=0,semiTryptic=opts.semiTryptic,reverse=opts.reverse)

        #name of the protein
        proteinId = protein.id
        if opts.reverse: proteinId = "DECOY_" + proteinId

        for peptide_sequence in peptides:
            for z in xrange(opts.minz,opts.maxz):
                p = Peptide(peptide_sequence,modifications=dict())

                #fixed mods
                if opts.tmt: p.addFixedModifications("TMT")
                if opts.iod:  p.addFixedModifications("Iodoacetamide")
                print("{}\t{}\t{}".format(p.modString(z),proteinId,p.getMZ(z)))
                
                if opts.oxy: 
                    #varible oxidataion
                    oxidations=[]
                    for i in range(0,len(peptide_sequence)):
                        if peptide_sequence[i] == "M": oxidations.append(i+1) #aapos

                    for pos in oxidations:
                        p.mods[pos] = "M[147]"
                        print("{}\t{}\t{}".format(p.modString(z),proteinId,p.getMZ(z)))
