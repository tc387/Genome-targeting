# GSP - Genome Similarity Python tool
# calculates the similarity between two genomes based on a matching length size
# by tc387 starting Mar 22 for Rosalind's project on genome DNAcc interactions

import numpy as np
import math
import time
import sys, getopt
import os, subprocess
from Bio.Seq import Seq # biopython
import re

# This Python script looks at two genomes and callculates the simlarity. The first genome is cut into 'blobs'.
# The script searches for the number of complementary overlaps of length 'matchlength' in the second genome.
def main():

    if (len(sys.argv) < 5) :
        print ('\n')
        print ('ERROR!  missing arguments. Need 4 arguments:\n -Genome1_filename \n' + \
               ' -Genome2_filename \n -blob size \n -match strand length \n')
        print ('Exiting ......\n')
        quit()

    gen1file = str(sys.argv[1])   # 1st genome filename, blobs  (e.g. Escherichia_coli.gen)
    lblob = int(sys.argv[3])  # number of bp per blob
    gen2file = str(sys.argv[2])   # 2nd genome filename, strands  (e.g. rBacillus_subtilis.gen)
    matchslength = int(sys.argv[4])  # matching sequence length for prescreening

    if (matchslength > lblob) :
        raise ValueError("Error: matchslength cannot be larger than the blob size. ")
        
        
    
    #outfilename_dGbs = "outdG_bgen=" + gen1file + "_lb" + str(lblob) + "_ib" + str(iiblob) + "_sgen="+gen2file + "_ls" + str(lstrand) + "_lm" + str(matchslength) + ".dat"
    outfilename = "out_sims_bg" + gen1file + "_sg="+gen2file + "_lb" + str(lblob) +  "_lm" + str(matchslength) +  ".dat"

# END PARAMETERS
# ============================================================================================== #
    # read genomes
    gen1str = "".join(line.rstrip("\n") for line in open(gen1file)) #full genome string
    gen1strc = "".join(line.rstrip("\n") for line in open("r"+gen1file))  # full genome string
    gen2str = "".join(line.rstrip("\n") for line in open(gen2file)) #full genome string
    gen2strc = "".join(line.rstrip("\n") for line in open("r"+gen2file))  # full genome string
    gen1 = Genomefrac("gen1", len(gen1str),gen1str)   # convert to the Genomefrac object
    gen2 = Genomefrac("gen2", len(gen2str),gen2str)
    gen1c = Genomefrac("gen1",  len(gen1str), gen1strc)  # convert to the Genomefrac object
    gen2c = Genomefrac("gen2",  len(gen2str), gen2strc)

    # cut genome
    gen1blobs=cutGenome(gen1, lblob, cutpos="begin")   # cut from begining (full blob at the begining)
    gen1cblobs = cutGenome(gen1c, lblob, cutpos="begin")  # cut from begining (full blob at the begining)
    # gen2strands=cutGenome(gen2, lstrand, cutpos="end") # cut from end (full strand at the end)
    nblobs = math.ceil(len(gen1str) / float(lblob))

    nmatches=np.zeros(2*nblobs)
    similarity=np.zeros(2*nblobs)

    outfile = open(outfilename, 'w')

    for iiblob in range(nblobs):
        print("Finnished.. ", iiblob / float(nblobs))
        # find matching occurrences in genome2
        candLocList = getStrandOverlaps(pattern=str(Seq(str(gen1blobs[iiblob].sequence)).reverse_complement()), matchlength=matchslength, string=gen2str)
        candLocList2 = getStrandOverlaps(pattern=str(Seq(str(gen1blobs[iiblob].sequence)).reverse_complement()),
                                        matchlength=matchslength, string=gen2strc)

        #    print (str(matchslength) + " " + str(Seq(str(gen1blobs[iiblob].sequence)).reverse_complement()) + "\n")
    #    print(candLocList)
        nmatches[iiblob] = len(candLocList) + len(candLocList2)
        outfile.write(str(iiblob) + " " + str(nmatches[iiblob].astype(np.float32) / (lblob-matchslength+1)) + "\n")
        outfile.flush()
    for iiblob in range(nblobs):
        print("Finnished.. ", iiblob / float(nblobs))
        # find matching occurrences in genome2
        candLocList = getStrandOverlaps(pattern=str(Seq(str(gen1cblobs[iiblob].sequence)).reverse_complement()), matchlength=matchslength, string=gen2str)
        candLocList2 = getStrandOverlaps(pattern=str(Seq(str(gen1cblobs[iiblob].sequence)).reverse_complement()), matchlength=matchslength, string=gen2strc)

    #    print (str(matchslength) + " " + str(Seq(str(gen1blobs[iiblob].sequence)).reverse_complement()) + "\n")
    #    print(candLocList)
        nmatches[iiblob+nblobs] = len(candLocList)+ len(candLocList2)
        outfile.write(str(iiblob+nblobs) + " " + str(nmatches[iiblob+nblobs].astype(np.float32) / (lblob-matchslength+1)) + "\n")
        outfile.flush()


    similarity = nmatches.astype(np.float32) / (lblob-matchslength+1)
    outfile.write(r"#rel_overlaps " + str(sum(similarity[:])/(2*nblobs)) + "\n")


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

    lpattern=len(pattern) # pattern length
    lstring=len(string)  # string length
    if (matchlength > lpattern) or (matchlength > lstring) :
        return None
    matchlocations=[]
    for ii in range(0,lpattern-matchlength+1) :
        iipat = '(?=' + pattern[ii:matchlength+ii] + ')'
        #[jj.start() for jj in re.finditer(iipat, string)]
        matchlocations.extend([jj.start() for jj in re.finditer(iipat, string)])
    return matchlocations


########  NuPACK Wrapper #########
class Nuwrap(object):
    def __init__(self, nu_args=None, nstrands=None, strandlist=None, permutation=None):
        self.nu_args=nu_args
        if strandlist is None:
            raise ValueError('Error: no strands supplied to NuPack wrapper Nuwrap')
        elif type(strandlist) is not list:
            raise ValueError('Error: strandlist must be a list! supplied to NuPack wrapper Nuwrap')
        else :
            self.strandlist=strandlist
        if nstrands is None:
            self.nstrands = len(strandlist)
        else :
            self.nstrands=nstrands
        if permutation is None :
            # default to ordered permutation
            self.permutation = " ".join([str(list(range(1,self.nstrands+1))[ii]) for ii in range (self.nstrands)])
        elif type(permutation) is not str :
            raise ValueError('Error: permutation argument of Nuwrap must be a string of integers, defaulting to a series')
            self.permutation = " ".join([str(list(range(1,self.nstrands+1))[ii]) for ii in range (self.nstrands)])
        else :
            self.permutation = permutation
        self.nu_stdout = ""   # stdout captured from NuPack
        self.runstatus = 0  # flag specifying whether the subprocess is running (0), has finnished (1) or returned an error (-1)
        self.delG = 0.0  # interaction free energy

    def run(self):
        # input string to NuPack stdin
        nu_instring = str(self.nstrands) + '\n' + "".join([self.strandlist[ii] + '\n' for ii in range(self.nstrands)]) + \
            self.permutation + '\n'
        nu_byteinstring=nu_instring.encode()  # encode in bytes
        # run subprocess
        try:
            nu_run = subprocess.run(self.nu_args, input=nu_byteinstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, \
                                check=True, shell=True)
        except subprocess.CalledProcessError as e :
            self.runstatus = -1  # return error
            self.nu_stdout = 'NuPack crashed... return code: '+ str(e.returncode) + '\nNuPack stdout: ' + e.stdout.decode()
        else :
            # run finished well
            self.runstatus = 1
            self.nu_stdout = nu_run.stdout.decode()

            # get free energy from NuPack stdout
            nu_stdout_lines = self.nu_stdout.split('\n') # split by lines
            outlen = len(nu_stdout_lines) # # of stdout lines
            self.delG = nu_stdout_lines[outlen-3]  # get free energy

########  GENOME CLASS  ##########
class Genomefrac(object):
    def __init__(self, genname="",size=0, sequence=None, genbpos=0, genblobindex=0):
        self.genname=genname         # genome name
        self.genbloc=genbpos         # genome location - starting base location in genome
        self.genblobindex = genblobindex  # genome location - serial number of blob in genome
        self.size = size  # size in bases of the sequence
        if sequence is None :
            self.sequence='T'*size   # default to polyT
        else :
            self.sequence=sequence
     #       self.size=len(sequence)


    def get_sequence(self, sizef, posindex):
        if posindex+sizef <= self.size :
            return self.sequence[posindex:posindex + sizef]
        elif posindex < self.size :
            return  self.sequence[posindex:self.size]
        else :
            return None

####### CUT GENOME FUNCTION ######
def cutGenome(Genome, fsize, cutpos):
    # cuts Genome into fragements of length size.
    ind = 0
    basepos = 0
    fraglist = [] # first fragmement is here
    if (cutpos=="begin"):
        fraglist.append(Genomefrac(genname=Genome.genname + '_' + str(ind), size=fsize, \
                               sequence=Genome.get_sequence(fsize, 0),
                               genbpos=Genome.genbloc + basepos, genblobindex=ind))
        ind += 1
        basepos += fsize
    elif (cutpos=="end") :
        offset = Genome.size % fsize
        if (offset == 0) :  # cannot have a zero length strand
            offset = fsize
        fraglist.append(Genomefrac(genname=Genome.genname + '_' + str(ind), size=offset, \
                                   sequence=Genome.get_sequence(offset, 0),
                                   genbpos=Genome.genbloc + basepos, genblobindex=ind))
        ind += 1
        basepos += offset
    else :
        raise ValueError("Error with cut position in function cutGenome ")

    while basepos < Genome.size :

        fraglist.append(Genomefrac(genname=Genome.genname + '_' + str(ind), size=fsize, \
                                   sequence=Genome.get_sequence(fsize, basepos),
                                   genbpos=Genome.genbloc + basepos, genblobindex=ind))
        ind += 1
        basepos += fsize
    return fraglist


# RUN MAIN
if __name__ == "__main__":
    main()
