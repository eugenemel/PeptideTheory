import sys
import csv
import re
import glob
import SpectraReader
import Peptide
import MolecularWeights

ANNOTATION=dict()

def read_specta(spectafile):
   reader = SpectraReader.SpectraReader(spectafile)

   while 1:
        #get experimental spectrum from MGF file
        scan = reader.nextMGFSpectrum()
        if not scan: break

        #normalize scanId
        scanId = scan.title
        scanId = re.sub(r'\.\d+\.\d$','', scanId)
        scanId = scanId.lower()

        #skip scans that have not been annotated
        if scanId not in ANNOTATION:
                #print "SKIPPING SCAN=", scan.title, " Has not been annotated."
                continue

        print("SCAN=", scan.title)
        for annotation in ANNOTATION[scanId]:
            peptide = annotation['peptide']
            peptide = re.sub(r'^\w\.','',peptide)
            peptide = re.sub(r'\.\w$','',peptide)

            #add modificaitons, TMT and Cys, and variable Met
            mods = { "n" : "n[230]"}
            aanum=0
            for aa in peptide:
                if aa in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":aanum += 1
                if aa == "*": mods[aanum] = "M[147]"
                if aa == "C": mods[aanum] = "C[160]"
                if aa == "K": mods[aanum] = "K[357]"

            #create peptide object and calculate theoretical mass
            pep = Peptide.Peptide(peptide, mods)
            theoryPreMz = pep.getMZ(scan.charge)
            delta = scan.preMz - theoryPreMz

            #ERROR. scan precursorMz and theoryMz are different? why?
            if (abs(delta) > 0.1 and abs(1.0073 - delta*scan.charge) > 0.1):
                print >> sys.stderr, "SKIPPING PEPTIDE.. ERROR:", delta*scan.charge, "z=", scan.charge,  scanId, annotation['peptide'], peptide, mods, "->", theoryPreMz, " vs ", scan.preMz
            else:
                #print annotation
                #print scan
                print ("\tHIT RANK=", annotation["Rank"], "PEPTIDE=",   annotation['peptide'],\
                                                         "DECOY=",     annotation['label'],\
                                                         "MODS=",     pep.modString(),\
                                                         "deltCn=",  annotation["deltCn"],\
                                                         "Xcorr=",    annotation['Xcorr'],\
                                                         "massDiff=", delta*scan.charge)

                #print out observed spectra?
                #for i in range(0,len(scan.mz)): print scan.mz[i], scan.intensity[i]

                #print theoretical spectra?
                all_ions = pep.all_ions( ionseries = ['b','y','a'], frg_z_list = [1,2], fragmentlossgains = [0, -MolecularWeights.mNH3, -MolecularWeights.mH2O ], mass_limits = [80,50000])
                #for ion in all_ions: print "\t", ion

#function to read in peptide spectrume annotation for perculator formated files with .pin ext
def read_pin_file(pinfile): 
        reader = csv.DictReader(open(pinfile, 'r'), delimiter=' ')
        d = {}
        for row in reader:
            scanId = row["id"] 
            scanId = re.sub(r'.*/','',scanId)
            scanId = re.sub(r'\.\d+\.\d$','', scanId)
            scanId = scanId.lower()
            if scanId not in d: d[ scanId ] = list()
            d[ scanId ].append(row)
        return d



if __name__ == '__main__':
    pinfile = None
    spectrafile = None

    for arg in sys.argv:
        if arg.endswith("pin"): pinfile = arg
        if arg.endswith("mgf.gz"): spectrafile = arg

    if pinfile:
        ANNOTATION=read_pin_file(pinfile)

    if spectrafile:
        read_specta(spectrafile)
