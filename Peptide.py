#PeptideTheory: Theoretical Ion Generation
#Author: Eugene Melamud / Calico / 2017


from MolecularWeights import *

class Peptide:
    def __init__(self, sequence,  modifications=dict()):
        self.sequence = sequence
        self.mods = modifications       ##dictionary example { 5: "C[160]" }

    def getMZ(self, charge):

        M=0
        for aa in self.sequence:
            M += C12MASS[aa]

        if self.mods:
            for m in self.mods.values():
                if m in MODS: 
                    M += MODS[m]
                else:
                    print "Unknown modification type:", m 

        if charge >= 1:
            return((M + mH2O + mProton*charge)/charge)
        else:
            return(M + mH2O + mProton)


    def all_ions(self, ionseries=['b','y'], frg_z_list=[1,2], fragmentlossgains=[0,], mass_limits=None):

        left  =  0
        L =  len(self.sequence)
        M =  self.getMZ(0);

        if "n" in self.mods: 
            left += MODS[ self.mods["n"] ]     #n-terminal mod

        ALLFRAG=list()
        for pos in range(0,len(self.sequence)-1):
            aa = self.sequence[pos]

            #split peptide into two segmentes
            aanum  = pos+1
            aamass = C12MASS[aa]
            if aanum in self.mods: aamass += MODS[ self.mods[aanum] ]

            left = left + aamass            #mass of the peptide to the left of break point
            right = M - left                #mass of the peptide to the right of break point

            #calculate m/z for different charges, gains and loses, and type of ion
            for z in frg_z_list:
                for W in fragmentlossgains:
                    for ion in ionseries:

                        mz=0
                        ion_num = aanum

                        if ion == 'b':
                            mz = (left + z*mProton+W)/z
                        elif ion == 'y':
                            ion_num = L - aanum
                            mz =  (right+(z-1)*mProton+W)/z

                        elif ion == 'a':
                            mz = (left -mCO + z*mProton+W)/z
                        elif ion == 'x':
                            ion_num = L - aanum
                            mz =  (right +mCO -mH2 + (z-1)*mProton+W)/z

                        elif ion == 'c':
                            mz = (left +mNH3 + z*mProton+W)/z
                        elif ion == 'z':
                            ion_num = L - aanum
                            mz =  (right -mNH3 + (z-1)*mProton+W)/z

                        if mz > 0 and mz >= mass_limits[0] and mz <= mass_limits[1]:
                            #print ion, ion_num, z, int(W), mz
                            ALLFRAG.append((ion, ion_num, z, int(W), mz))
        return(ALLFRAG)

def test_peptide_class():
        #p1  = Peptide('CKDALA',modifications={1: "C[160]"})
        p1  = Peptide('CKDALA',modifications={1: "C[160]", "n": "n[230]", 2: "K[357]"})

        print "preMz: ", p1.getMZ(1)
        print "preMz: ", p1.getMZ(2)
        print "preMz: ", p1.getMZ(3)
        print "preMz: ", p1.getMZ(4)

        all_ions = p1.all_ions( ionseries = ['b','y','a'], 
                frg_z_list = [1,2],
                fragmentlossgains = [0, -mNH3, -mH2O ], 
                mass_limits = [80,3000]
                )

