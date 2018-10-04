#PeptideTheory: Theoretical Ion Generation
#Author: Eugene Melamud / Calico / 2017


import re
import random

class Protein:

    def __init__(self, id, descirption=None, sequence=None,  modifications=dict()):
        self.id = id
        self.descirption = descirption
        self.sequence = sequence

    def __str__(self):
        return self.sequence

    def reverse(self):
        return ''.join(reversed(sequence))

    @staticmethod
    def enzymePattern(enzyme_name="Trypsin"):

        D= {
            "ArgC" :            ["[R]",     1],         #c-term
            "AspN" :            ["[D]",     0],         #n-termi 
            "Chymotrypsin" :    ["[FLYWM]", 0],         #(c-term) but not before P
            "GluC" :            ["[E]",     1],         #(c-term) but not before P
            "LysC" :            ["[K]",     1],         #(c-term)
            "LysN" :            ["[K]",     1],         #(n-term)
            "Trypsin" :         ["[KR]",    1],         #(c-term)
        }

        if enzyme_name in D:
            return D[enzyme_name]
        else:
            raise "Unknown enzyme"


    def digest(self,pattern="[RK]",cutoffset=1,minlength=3,maxlength=-1,maxmiss=-1, missprob=0.0, semiTryptic=False, reverse=False):

        protein_sequence = self.sequence
        if reverse: protein_sequence = self.sequence[::-1]

        #get cut locations
        positions = []
        itr = re.finditer(pattern, protein_sequence, flags=0)
        for m in itr:
            positions.append(m.start())

        #partial digest, cuts missed with uniform frequency
        if missprob > 0:
            K = int(len(positions)*(1-missprob))
            indices = random.sample(range(len(positions)), K)
            positions = [positions[i] for i in sorted(indices)]

        #always include ends of the sequene
        positions.insert(0,0)
        positions.append(len(protein_sequence))

        #print positions

        #get paptides, filtering on length
        peptides=[]
        inside_cut=0
        for a in range(0,len(positions)):
            
            brange = len(positions);
            if maxmiss >= 0:
                brange = min(a+1+maxmiss,len(positions))

            for b in range(a+1, brange):
                if a>0: inside_cut=1
                peptide = protein_sequence[ positions[a]+inside_cut : positions[b]+cutoffset ]
                #print positions[a], positions[b], peptide #debug

                if maxlength > 0 and len(peptide) > maxlength: continue
                if len(peptide) <= minlength: continue
                peptides.append(peptide)

                #semi tryptic digest
                if semiTryptic:
                    for c in range(positions[a]+inside_cut+1,  positions[b]):
                        peptide = protein_sequence[ c : positions[b]+cutoffset ]
                        #print "\t", c, positions[b], peptide #debug
                        if maxlength > 0 and len(peptide) > maxlength: continue
                        if len(peptide) <= minlength: continue
                        peptides.append(peptide)

        return(peptides)

