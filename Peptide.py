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
            if aa in C12MASS:
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

    def hydrophobicity(self):

        SCALE = { 
           "A"  :  0.620, "R"  : -2.530,  "N"  : -0.780, "D"  : -0.900, "C"  :  0.290,
           "Q"  : -0.850, "E"  : -0.740,  "G"  :  0.480, "H"  : -0.400, "I"  :  1.380,
           "L"  :  1.060, "K"  : -1.500,  "M"  :  0.640, "F"  :  1.190, "P"  :  0.120,
           "S"  : -0.180, "T"  : -0.050,  "W"  :  0.810, "Y"  :  0.260, "V"  :  1.080
        }

        SCORE=0
        for aa in self.sequence:
            if aa in SCALE:
                SCORE += SCALE[aa]

        return SCORE
        
    def kideraFactors(self):
        # 188 physical properties of the 20 amino acids 
        #Reference Kidera, A., Konishi, Y., Oka, M., Ooi, T., & Scheraga, H. A. (1985). Statistical analysis of the physical properties of the 20 naturally occurring amino acids. Journal of Protein Chemistry, 4(1), 23-55

        SCALE = {
            "A": [ -1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.20,-0.08, 0.21,-0.48],
            "R": [  0.22, 1.27, 1.37, 1.87,-1.70, 0.46, 0.92,-0.39, 0.23, 0.93],
            "N": [  1.14,-0.07,-0.12, 0.81, 0.18, 0.37,-0.09, 1.23, 1.10,-1.73],
            "D": [  0.58,-0.22,-1.58, 0.81,-0.92, 0.15,-1.52, 0.47, 0.76, 0.70],
            "C": [  0.12,-0.89, 0.45,-1.05,-0.71, 2.41, 1.52,-0.69, 1.13, 1.10],
            "Q": [ -0.47, 0.24, 0.07, 1.10, 1.10, 0.59, 0.84,-0.71,-0.03,-2.33],
            "E": [ -1.45, 0.19,-1.61, 1.17,-1.31, 0.40, 0.04, 0.38,-0.35,-0.12],
            "G": [  1.46,-1.96,-0.23,-0.16, 0.10,-0.11, 1.32, 2.36,-1.66, 0.46],
            "H": [ -0.41, 0.52,-0.28, 0.28, 1.61, 1.01,-1.85, 0.47, 1.13, 1.63],
            "I": [ -0.73,-0.16, 1.79,-0.77,-0.54, 0.03,-0.83, 0.51, 0.66,-1.78],
            "L": [ -1.04, 0.00,-0.24,-1.10,-0.55,-2.05, 0.96,-0.76, 0.45, 0.93],
            "K": [ -0.34, 0.82,-0.23, 1.70, 1.54,-1.62, 1.15,-0.08,-0.48, 0.60],
            "M": [ -1.40, 0.18,-0.42,-0.73, 2.00, 1.52, 0.26, 0.11,-1.27, 0.27],
            "F": [ -0.21, 0.98,-0.36,-1.43, 0.22,-0.81, 0.67, 1.10, 1.71,-0.44],
            "P": [  2.06,-0.33,-1.15,-0.75, 0.88,-0.45, 0.30,-2.30, 0.74,-0.28],
            "S": [  0.81,-1.08, 0.16, 0.42,-0.21,-0.43,-1.89,-1.15,-0.97,-0.23],
            "T": [  0.26,-0.70, 1.21, 0.63,-0.10, 0.21, 0.24,-1.15,-0.56, 0.19],
            "W": [  0.30, 2.10,-0.72,-1.57,-1.16, 0.57,-0.48,-0.40,-2.30,-0.60],
            "Y": [  1.38, 1.48, 0.80,-0.56, 0.00,-0.68,-0.31, 1.03,-0.05, 0.53],
            "V": [ -0.74,-0.71, 2.04,-0.40, 0.50,-0.81,-1.07, 0.06,-0.46, 0.65],
        }

        SCORE=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for aa in self.sequence:
            if aa in SCALE:
                S = SCALE[aa]
                for kf in range(0,10):
                    SCORE[kf] += S[kf]

        for kf in range(0,10):
            SCORE[kf] /= len(self.sequence)

        return SCORE

def test_peptide_class():
        #p1  = Peptide('CKDALA',modifications={1: "C[160]"})
        #p1  = Peptide('CKDALA',modifications={1: "C[160]", "n": "n[230]", 2: "K[357]"})
        p1  = Peptide("KLKLLLLLKLK")

        print p1.sequence;

        print "HydroPhobictiy=", p1.hydrophobicity()
        print "KideraFactors=",  p1.kideraFactors()
        print "preMz: ", p1.getMZ(1)
        print "preMz: ", p1.getMZ(2)
        print "preMz: ", p1.getMZ(3)
        print "preMz: ", p1.getMZ(4)

        all_ions = p1.all_ions( ionseries = ['b','y','a'], 
                frg_z_list = [1,2],
                fragmentlossgains = [0, -mNH3, -mH2O ], 
                mass_limits = [80,3000]
                )

#test_peptide_class()
