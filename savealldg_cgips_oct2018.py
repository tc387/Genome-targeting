# cGIPS - cluster Genome Interaction Python Santalucia
# by tc387 starting Mar 4 for Rosalind's project on genome DNAcc interactions

# This script is a Python wrapper for the NuPack software which calculates the SantaLucia interaction free energies

# adapted on May 30th 2017: save all delG strand numbers
# adapted on Oct 2018: remove any possible digits, line breaks and white spaces in the genome file.

import math
import re
import subprocess
import sys

import numpy as np
from Bio.Seq import Seq  # biopython


# This is a python script that takes a ss genome and a blob strand, and callculates Santa Lucia (nupack)
# binding free energies between the blob strand and all strands of the genome.
def main():
    if (len(sys.argv) < 8):
        print('\n')
        print('ERROR!  missing arguments. Need 7 arguments:\n -Genome1_filename \n' \
              ' -Genome2_filename \n -strand size \n -match strand length prescreening \n -Temperature \n -blob size \n -blob serial number \n')
        print('Exiting ......\n')
        quit()

    gen1file = str(sys.argv[1])  # 1st genome filename, blobs  (e.g. Escherichia_coli.gen)
    lblob = int(sys.argv[6])  # number of bp per blob
    iiblob = int(sys.argv[7])  # blob serial number
    gen2file = str(sys.argv[2])  # 2nd genome filename, strands  (e.g. rBacillus_subtilis.gen)
    lstrand = int(sys.argv[3])  # strand length - number of bp per strand
    matchslength = int(sys.argv[4])  # matching sequence length for prescreening
    Temperature = float(sys.argv[5])  # Temperature
    # outfilename_dGbs = "outdG_bgen=" + gen1file + "_lb" + str(lblob) + "_ib" + str(iiblob) + "_sgen="+gen2file + "_ls" + str(lstrand) + "_lm" + str(matchslength) + ".dat"
    outfilename_dGbs = "outdG_bgen.dat"
    # arguments to NuPack
    nu_args = "./../../pfunc -T " + str(Temperature) + " -multi -material dna"

    kcalmol_kbT_conv = 694.77 / 1.38065 / (273.15 + Temperature)  # kcal/mol to k_BT conversion factor
    # END PARAMETERS
    # ============================================================================================== #
    # read genomes
    gen1str = "".join(line.rstrip("\n") for line in open(gen1file)).lower()  # full genome string, lower case
    gen2str = "".join(line.rstrip("\n") for line in open(gen2file)).lower()  # full genome string, lower case
    # remove count digits and white spaces
    gen1str = "".join(i for i in gen1str if any(i == j for j in ['a', 'c', 'g', 't']))
    gen2str = "".join(i for i in gen2str if any(i == j for j in ['a', 'c', 'g', 't']))
    print(gen2str)
    gen1 = Genomefrac("gen1", len(gen1str), gen1str)  # convert to the Genomefrac object
    gen2 = Genomefrac("gen2", len(gen2str), gen2str)

    # cut genome
    gen1blobs = cutGenome(gen1, lblob, cutpos="begin")  # cut from begining (full blob at the begining)
    gen2strands = cutGenome(gen2, lstrand, cutpos="end")  # cut from end (full strand at the end)

    if (iiblob > len(gen1blobs)):
        raise ValueError('Error: chosen iiblob > total number of blobs in the genome 1')

    # find matching occurrences in genome2
    candLocList = getStrandOverlaps(pattern=str(Seq(str(gen1blobs[iiblob].sequence)).reverse_complement()),
                                    matchlength=matchslength, string=gen2str)
    # list of candidate strands with removed possible duplicates,
    # matchlength/2 ensures we take the middle location of the matching length to find corresponding strand
    # (lstrand - gen2strands[0].size) ensures that a correct stand index is obtained if the first strand is shorter: gen2strands[0].size
    strcandStrandList = sorted(list(
        set(np.floor((np.array(candLocList) + matchslength / 2 + (lstrand - gen2strands[0].size)) / float(lstrand)))))
    candStrandList = [int(ii) for ii in strcandStrandList]  # convert to integers

    # callculate the free energy of blob only
    nurun1b = Nuwrap(nu_args, nstrands=1, strandlist=[gen1blobs[iiblob].sequence], permutation="1")
    nurun1b.run()
    if (nurun1b.runstatus == -1):
        raise RuntimeError("Error: NuPack blob self energy crashed: \n" + nurun1b.nu_stdout)
    delGblob = float(nurun1b.delG)  # free energy of blob only

    # callculate the free energy of all (matching) strands
    iirun = 0
    nurun2 = []
    nstrands = len(gen2strands)  # number of strands in genome 2

    deltaGlist = [0] * math.ceil(nstrands)  # define deltaG array and initialise to 0 (units: k_BT)
    for iistrand in candStrandList:  # iterate over all candidate strands
        print("Finnished ", iirun / len(candStrandList))
        iirun += 1

        # run NuPack c caluclation
        # free energy of strand only
        nurun1s = Nuwrap(nu_args, nstrands=1, strandlist=[gen2strands[iistrand].sequence], permutation="1")
        nurun1s.run()
        if (nurun1s.runstatus == -1):
            raise RuntimeError("Error: NuPack strand " + str(iistrand) + " self energy crashed: \n" + nurun1s.nu_stdout)
        delGstrand = float(nurun1s.delG)
        # free energy of blob - strand complex
        nurun2.append(
            Nuwrap(nu_args, nstrands=2, strandlist=[gen1blobs[iiblob].sequence, gen2strands[iistrand].sequence],
                   permutation="1 2"))
        nurun2[iirun - 1].run()
        if (nurun2[iirun - 1].runstatus == -1):
            raise RuntimeError("Error: NuPack blob-strand " + str(iistrand) + " free energy calc crashed: \n" + nurun2[
                iirun - 1].nu_stdout)
        iideltaG = kcalmol_kbT_conv * (float(nurun2[iirun - 1].delG) - delGstrand - delGblob)
        # print("\n")
        # print(iiblob," ",iistrand)
        # print(gen1blobs[iiblob].sequence)
        # print(gen2strands[iistrand].sequence)
        # print(float(nurun2[iirun-1].delG),delGstrand, delGblob)
        # print("delGcomplex: ", iideltaG,"\n")
        # print(nurun2[iirun-1].runstatus, nurun2[iirun-1].nu_stdout)
        deltaGlist[iistrand] = iideltaG
    sumexpdelG = float(
        -np.log(np.sum(np.exp(-np.array(deltaGlist)) - 1.0) + 1.0))  # exp delG sum, i.e. partition function

    print("sumexp_deltaG: ", sumexpdelG)
    # save deltaGs to file
    outfile = open(outfilename_dGbs, 'w')
    outfile.write("#sumexp_deltaG: " + str(sumexpdelG) + "\n")
    for ind in range(nstrands):

        if (deltaGlist[ind] < 0.0):
            outfile.write(str(gen2strands[ind].genbloc) + "  " + str(np.round(deltaGlist[ind], decimals=3)) + "\n")
    outfile.close()


# ============================================================================================== #
# END MAIN


def getStrandOverlaps(pattern, matchlength, string):
    # check all possible continous overlaps of length matchlength and find their locations in string
    if type(pattern) is not str:
        raise ValueError('Error: pattern in getStrandOverlaps must be type str')
    if type(string) is not str:
        raise ValueError('Error: string in getStrandOverlaps must be string')
    if type(matchlength) is not int:
        raise ValueError('Error: matchlength in getStrandOverlaps must be integer')

    lpattern = len(pattern)  # pattern length
    lstring = len(string)  # string length
    if (matchlength > lpattern) or (matchlength > lstring):
        return None
    matchlocations = []
    for ii in range(0, lpattern - matchlength):
        iipat = '(?=' + pattern[ii:matchlength + ii] + ')'
        # [jj.start() for jj in re.finditer(iipat, string)]
        matchlocations.extend([jj.start() for jj in re.finditer(iipat, string)])
    return matchlocations


########  NuPACK Wrapper #########
class Nuwrap(object):
    def __init__(self, nu_args=None, nstrands=None, strandlist=None, permutation=None):
        self.nu_args = nu_args
        if strandlist is None:
            raise ValueError('Error: no strands supplied to NuPack wrapper Nuwrap')
        elif type(strandlist) is not list:
            raise ValueError('Error: strandlist must be a list! supplied to NuPack wrapper Nuwrap')
        else:
            self.strandlist = strandlist
        if nstrands is None:
            self.nstrands = len(strandlist)
        else:
            self.nstrands = nstrands
        if permutation is None:
            # default to ordered permutation
            self.permutation = " ".join([str(list(range(1, self.nstrands + 1))[ii]) for ii in range(self.nstrands)])
        elif type(permutation) is not str:
            raise ValueError(
                'Error: permutation argument of Nuwrap must be a string of integers, defaulting to a series')
            self.permutation = " ".join([str(list(range(1, self.nstrands + 1))[ii]) for ii in range(self.nstrands)])
        else:
            self.permutation = permutation
        self.nu_stdout = ""  # stdout captured from NuPack
        self.runstatus = 0  # flag specifying whether the subprocess is running (0), has finnished (1) or returned an error (-1)
        self.delG = 0.0  # interaction free energy

    def run(self):
        # input string to NuPack stdin
        nu_instring = str(self.nstrands) + '\n' + "".join([self.strandlist[ii] + '\n' for ii in range(self.nstrands)]) + \
                      self.permutation + '\n'

        print("This is nupack instring: \n", nu_instring, "\n")
        nu_byteinstring = nu_instring.encode()  # encode in bytes
        # run subprocess
        try:
            nu_run = subprocess.run(self.nu_args, input=nu_byteinstring, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, \
                                    check=True, shell=True)
        except subprocess.CalledProcessError as e:
            self.runstatus = -1  # return error
            self.nu_stdout = 'NuPack crashed... return code: ' + str(
                e.returncode) + '\nNuPack stdout: ' + e.stdout.decode()
        else:
            # run finished well
            self.runstatus = 1
            self.nu_stdout = nu_run.stdout.decode()

            # get free energy from NuPack stdout
            nu_stdout_lines = self.nu_stdout.split('\n')  # split by lines
            outlen = len(nu_stdout_lines)  # # of stdout lines
            self.delG = nu_stdout_lines[outlen - 3]  # get free energy


########  GENOME CLASS  ##########
class Genomefrac(object):
    def __init__(self, genname="", size=0, sequence=None, genbpos=0, genblobindex=0):
        self.genname = genname  # genome name
        self.genbloc = genbpos  # genome location - starting base location in genome
        self.genblobindex = genblobindex  # genome location - serial number of blob in genome
        self.size = size  # size in bases of the sequence
        if sequence is None:
            self.sequence = 'T' * size  # default to polyT
        else:
            self.sequence = sequence

    #       self.size=len(sequence)

    def get_sequence(self, sizef, posindex):
        if posindex + sizef <= self.size:
            return self.sequence[posindex:posindex + sizef]
        elif posindex < self.size:
            return self.sequence[posindex:self.size]
        else:
            return None


####### CUT GENOME FUNCTION ######
def cutGenome(Genome, fsize, cutpos):
    # cuts Genome into fragements of length size.
    ind = 0
    basepos = 0
    fraglist = []  # first fragmement is here
    if (cutpos == "begin"):
        fraglist.append(Genomefrac(genname=Genome.genname + '_' + str(ind), size=fsize, \
                                   sequence=Genome.get_sequence(fsize, 0),
                                   genbpos=Genome.genbloc + basepos, genblobindex=ind))
        ind += 1
        basepos += fsize
    elif (cutpos == "end"):
        offset = Genome.size % fsize
        if (offset == 0):  # cannot have a zero length strand
            offset = fsize
        fraglist.append(Genomefrac(genname=Genome.genname + '_' + str(ind), size=offset, \
                                   sequence=Genome.get_sequence(offset, 0),
                                   genbpos=Genome.genbloc + basepos, genblobindex=ind))
        ind += 1
        basepos += offset
    else:
        raise ValueError("Error with cut position in function cutGenome ")

    while basepos < Genome.size:
        fraglist.append(Genomefrac(genname=Genome.genname + '_' + str(ind), size=fsize, \
                                   sequence=Genome.get_sequence(fsize, basepos),
                                   genbpos=Genome.genbloc + basepos, genblobindex=ind))
        ind += 1
        basepos += fsize
    return fraglist


# RUN MAIN
if __name__ == "__main__":
    main()
