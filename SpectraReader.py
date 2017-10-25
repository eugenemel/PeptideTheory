#PeptideTheory: Theoretical Ion Generation
#Author: Eugene Melamud / Calico / 2017

import re
import gzip

class Spectra:
    def __init__(self):
        self.title=None
        self.scannum=None
        self.rt=None
        self.charge=None
        self.preMz=None
        self.mz=[]
        self.intensity=[]

    def nobs():
        return len(self.mz)

    def info(self):
        print(self.title, "SCAN=", self.scannum,\
                          "PRE=",  self.preMz,\
                          "CHARGE=",self.charge,\
                          "RT=",self.rt,\
                          "NOBS=",len(self.mz), "TIC=", sum(self.intensity))

class SpectraReader:

    def __init__(self,filename):
        self.filename = filename

        if filename.endswith(".gz"):
            self.fh = gzip.open(filename, 'rb')
        else:
            self.fh = open(filename)

    def nextMGFSpectrum(self):
        if not self.fh: return None

        s = None
        for line in self.fh:
            line = line.rstrip().lstrip().upper()
            if not line: continue

            if line.startswith("BEGIN IONS"): s = Spectra()
            if not s: continue

            if line.startswith("TITLE"):
                s.title = line.split("=",2)[1]
                scaninfo = s.title.split(".")
                if len(scaninfo) >= 2: s.scannum=int(scaninfo[1])
            elif line.startswith("PEPMASS"):
                values = re.split(r'\s+|=',line)
                if len(values) >= 1: s.preMz =float(values[1])
            elif line.startswith("RTINSECONDS"):
                s.rt =float(line.split("=",2)[1])
            elif line.startswith("CHARGE"):
                s.charge =int(line.split("=",2)[1][0])
            elif line[0].isdigit() and line[-1].isdigit():
                _mz,_int = line.split(" ")
                s.mz.append(float(_mz))
                s.intensity.append(float(_int))
            elif line.startswith("END IONS"):
                return s

def testClass():
    reader = SpectraReader("..//comet_source_2016012/yeast_aging_test/mgf.mgf")
    while 1:
        s = reader.nextMGFSpectrum()
        if not s: break
        print(s.preMz, s.charge)

#testClass()
