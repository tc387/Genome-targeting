# CLUSTER  JOB SUBMISION SCRIPT
#
# for cGIPS - cluster Genome Interaction Python Santalucia
# by tc387 starting Mar 4 for Rosalind's project on genome DNAcc interactions
# this script is used for massive job submission for calculation of SANTA LUCIA free energies.

import math
import subprocess
import sys

#####  PARAMETERS #####
if (len(sys.argv) < 7):
    print('\n')
    print('ERROR!  missing arguments. Need 6 arguments:\n -Genome1_filename \n' \
          ' -Genome2_filename \n -blob size \n -strand size \n -match strand length prescreening \n -Temperature \n')
    print('Exiting ......\n')
    quit()

gen1file = str(sys.argv[1])  # 1st genome filename, blobs  (e.g. Escherichia_coli.gen)
lblob = int(sys.argv[3])  # number of bp per blob
gen2file = str(sys.argv[2])  # 2nd genome filename, strands  (e.g. rBacillus_subtilis.gen)
lstrand = int(sys.argv[4])  # strand length - number of bp per strand
matchslength = int(sys.argv[5])  # matching sequence length for prescreening
Temperature = float(sys.argv[6])  # Temperature

blobsperbatch = 100

# read genome
gen1str = "".join(line.rstrip("\n") for line in open(gen1file))  # full genome string
nblobs = math.ceil(len(gen1str) / float(lblob))  ######################################################## HACK
# main directory name
gdirname = "nGen-" + gen1file[:-4] + "-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
           "_lm" + str(matchslength) + "_T" + str(round(Temperature))
scratchdirname = "/sharedscratch/tc387/cGIPSdat/" + gdirname

subprocess.run(["mkdir " + gdirname], shell=True)
subprocess.run(["mkdir " + scratchdirname + "/"], shell=True)
iblobstart = 0  # starting blob number for the batch array
jobrunargs = r"python3 ../../savealldg_cgips_oct2017.py ../../" + gen1file + r" ../../" + gen2file + " " + \
             str(lstrand) + " " + str(matchslength) + " " + str(Temperature) + " " + str(lblob) + " "

# calculate the partition function (with all strands of genome 2) over all blobs in the genome 1
for iiblob in range(nblobs):
    # job arguments
    #   arrayjobrunargs = r"../cgips.py~../" + gen1file + r"~../" + gen2file + "~" + \
    #                str(lstrand) + "~" + str(matchslength) + "~" + str(Temperature) +  "~" + str(lblob) + "~"    # make directory and submit job
    # make directory and submit job
    dirname = "blobcgips_bg-" + gen1file[:-4] + "_sg-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
              "_lm" + str(matchslength) + "_ib" + str(iiblob)
    subprocess.run(["mkdir " + gdirname + "/" + dirname], shell=True)
    #  print(jobrunargs)
    #  subprocess.run([r"cat qsubs-master.sh | sed -e 's/JOBNAME11/cGIPS_ib" + str(iiblob) + r"/' | sed -e 's~JOBRUNARGS11~" + jobrunargs + r"~' > " + dirname + r"/qsub_local.sh"], shell=True)
    #  subprocess.run(["cp " + "pfunc  cgips.py " + dirname + "/"], shell=True)
    if (iiblob % blobsperbatch == blobsperbatch - 1):  # submit blobsperbatch blobs as a single slurm job
        # subprocess.run(["sbatch --array=0-999 array-qsubs.sh " + str(iblobstart) + " " + str(arrayjobrunargs)], shell=True )
        dirname0b = "blobcgips_bg-" + gen1file[:-4] + "_sg-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(
            lstrand) + \
                    "_lm" + str(matchslength) + "_ib" + str(iblobstart)
        # batchjobrunargs = r"cd .. ; for ii in $(seq "+str(iblobstart)+" "+str(iblobstart+blobsperbatch-1) + r" ) ; do cd blobcg*_ib$ii ; " + jobrunargs + "$ii > screen.txt ; cd .. ; mv blobcg*_ib$ii " + scratchdirname + "/ ; done"
        batchjobrunargs = r"cd .. ; for ii in $(seq " + str(iblobstart) + " " + str(
            iblobstart + blobsperbatch - 1) + r" ) ; do cd blobcg*_ib$ii ; " + jobrunargs + "$ii > screen.txt ; cd .. ; done"
        subprocess.run([r"cat qsubs-master.sh | sed -e 's/JOBNAME11/cGIPS_ib" + str(iblobstart) + "-" + str(
            iiblob) + r"/' | sed -e 's@JOBRUNARGS11@" + \
                        batchjobrunargs + r"@' > " + gdirname + "/" + dirname0b + r"/qsub_local.sh"], shell=True)
        print("Batch job iblob" + str(iblobstart))
        slurmargs = ["sbatch qsub_local.sh"]
        subprocess.Popen(slurmargs, cwd=gdirname + "/" + dirname0b, shell=True)

        iblobstart = iiblob + 1
# ru narray batch job
# print("sbatch array-qsubs.sh " + str(nblobs) + " " + str(arrayjobrunargs))
# jobname = "'GenB-" + gen1file[:-4] + "-" + gen2file[:-4] + "'"
if (nblobs > iblobstart):
    dirname0b = "blobcgips_bg-" + gen1file[:-4] + "_sg-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(
        lstrand) + "_lm" + str(matchslength) + "_ib" + str(iblobstart)
    batchjobrunargs = r"cd .. ;for ii in $(seq " + str(iblobstart) + " " + str(
        nblobs - 1) + r" ) ; do cd blobcg*_ib$ii ; " + jobrunargs + "$ii > screen.txt ; cd .. ; done"

    subprocess.run([r"cat qsubs-master.sh | sed -e 's/JOBNAME11/cGIPS_ib" + str(iblobstart) + "-" + str(
        nblobs - 1) + r"/' | sed -e 's@JOBRUNARGS11@" + \
                    batchjobrunargs + r"@' > " + gdirname + "/" + dirname0b + r"/qsub_local.sh"], shell=True)
    slurmargs = ["sbatch qsub_local.sh"]
    subprocess.Popen(slurmargs, cwd=gdirname + "/" + dirname0b, shell=True)

#   subprocess.run(["sbatch --array=0-"+str(nblobs-iblobstart-1) + " array-qsubs.sh " + str(iblobstart) + " " + str(arrayjobrunargs)], shell=True )
# subprocess.run(["sbatch --array=0-"+str(nblobs-1) + " array-qsubs.sh " + str(nblobs) + " " + str(arrayjobrunargs)], shell=True )
