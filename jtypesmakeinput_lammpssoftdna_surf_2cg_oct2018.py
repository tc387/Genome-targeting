# /usr/bin/env python

import numpy as np
import math
import sys, getopt
import os, subprocess
import re

def main():
    if (len(sys.argv) < 8):
        print('\n')
        print('ERROR!  missing arguments. Need 8 arguments:\n -Genome1_filename \n' \
              ' -Genome2_filename \n -blob size \n -strand size \n -match strand length prescreening \n -stranddensity on the surface \n -nstrandtypes \n -Temperature \n ' )
        print('Exiting ......\n')
        quit()

    gen11file = str(sys.argv[1])  # 1st genome filename, blobs  (e.g. Escherichia_coli.gen)
    gen12file="r"+gen11file
    lblob = int(sys.argv[3])  # number of bp per blob
    gen2file = str(sys.argv[2])  # 2nd genome filename, strands  (e.g. rBacillus_subtilis.gen)
    lstrand = int(sys.argv[4])  # strand length - number of bp per strand
    matchslength = int(sys.argv[5])  # matching sequence length for prescreening
    stranddensity = float(sys.argv[6])
    nstrandtypes = int(sys.argv[8])
    lencutgenseg = int(sys.argv[7])
    Temperature = round(float(sys.argv[9]))  # Temperature

    # PARAMETERS
    lbox = 50000
 #   Nfree_col = 200  # 200 # number of colloids
    Nfree_polymers = 2
    blobs_per_freepolymer = 100000  # if load_dGdata, limit the number of blobs to this number, else limit by the dGdata.
    rb = 10 # * (lblob/400)**0.588  # blob radius ; good for 400bp blobs --
 #   rc = 200.  # colloid radius
 #   colmass = 100.0
    nstrands = lbox*lbox * stranddensity
    lboxxy = lbox / rb / 5
    lboxz = lbox / rb

    load_dGdata = True  # load dG data from file

    np.random.seed(12345) # random number generator seed (numpy)
    polymer_form = 'linear'
    if (gen11file in ["HSM.gen", "rHSM.gen", "rrHSM.gen"]):
        polymer_form='circular'

    # conversion of delG from nupack to H parameter gaussian./cut pairstyle in lammps simulations, maximum value is 0 --no repulsion.
    fdG_to_H_conv_cons = math.log(rb**3 / 1.66 * 55)   # conversion constants
    fdG_to_H_stpi_cons =  math.sqrt(2*math.pi) # conversion constants
  #  fdG_to_H_conv = lambda x : np.min((x  + fdG_to_H_conv_cons ) * stpi_cons, np.zeros(x.shape))

    gen11str = "".join(line.rstrip("\n") for line in open("dGdata/"+gen11file))  # full genome string
    gen12str = "".join(line.rstrip("\n") for line in open("dGdata/" + gen12file))  # full genome string
    nblobs = math.ceil(len(gen11str) / float(lblob))  ######################################################## HACK
    gen2str = "".join(line.split()[0].rstrip("\n") for line in open("dGdata/"+gen2file))  # full genome string
    nstrands = math.ceil(len(gen2str) / float(lstrand))  ######################################################## HACK


    if (blobs_per_freepolymer > nblobs):
        blobs_per_freepolymer = nblobs
    basestrandoffset = len(gen2str) % lstrand
    totnblobs = blobs_per_freepolymer * Nfree_polymers
    max_ntypes = 2 #blobs_per_freepolymer + Nfree_col  # number of 'atom' types in the simulation

    if (len(gen2str) < nstrandtypes*lstrand) :
        print("WARNNING: set number of strand types is larger than the available numberof strands in file "  + gen2file +
              "/n" + "setting nstrandtypes=len(gen2str) " + "\n")
        nstrandtypes=np.ceil(len(gen2str)/float(lstrand))

    # main directory name
    # gdirname = "nGen-" + gen1file[:-4] + "-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
    #           "_lm" + str(matchslength) + "_T" + str(Temperature)

    ibdirname1 = "dGdata/nGen-" + gen11file[:-4] + "-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
               "_lm" + str(matchslength) + "_T" + str(Temperature) + "/blobcgips_bg-" + gen11file[:-4] + "_sg-" + \
               gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
               "_lm" + str(matchslength)
    ibdirname2 = "dGdata/nGen-" + gen12file[:-4] + "-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
                 "_lm" + str(matchslength) + "_T" + str(Temperature) + "/blobcgips_bg-" + gen12file[:-4] + "_sg-" + \
                 gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
                 "_lm" + str(matchslength)

    ldirname = "surflam-2g-" + gen11file[:-4] + "-" + gen2file[:-4] + "_lb" + str(lblob) + "_ls" + str(lstrand) + \
    "_lm" + str(matchslength) + "_T" + str(Temperature) + "_sd" + str(stranddensity) + "_nst"+ str(nstrandtypes) + \
               "_lcgs" + str(lencutgenseg) + "_lbox" + str(lbox)

    datafilename="DATA.FILE"

    # outdG_bg-EC_sg-rEC_lb200_ls20_lm12.dat
    # outdG_bg-BS_sg-rEC_lb200_ls20_lm12.dat
    # outdG_bg-EC_O157_sakai_sg-rEC_K12W3110_lb200_ls20_lm12.dat
    # outdG_bg-EC_CFT073_sg-rEC_K12W3110_lb200_ls20_lm12.dat
    # outdG_bg-EC_K12MG1655_sg-rEC_K12W3110_lb200_ls20_lm12.dat
    # outdG_bg-SE_EST_sg-rEC_K12W3110_lb200_ls20_lm12.dat
    # outdG_bg-PA_PAO1_sg-rEC_K12W3110_lb200_ls20_lm12.dat
    genfile_dG = "dGdata/outdG_bg-BS_sg-rEC_lb200_ls20_lm12_T60.dat"

    load_bnei_data = False  # load neigbour blob interaction data data from file
    neifile_dG = "dGdata/out_bneidG_BS_lb200_T60.0.dat"

    if load_dGdata :
        dG_matrix1 = np.zeros((nblobs, nstrands))
        dG_matrix2 = np.zeros((nblobs, nstrands))
        dGt_matrix1 = np.zeros((nblobs, int(np.ceil(len(gen2str)/float(lstrand)))))
        dGt_matrix2 = np.zeros((nblobs, int(np.ceil(len(gen2str) / float(lstrand)))))
       # print ("Loading dG data...")
        for ib in range(nblobs):

            ibfilename = ibdirname1 + "_ib" + str(ib) + "/outdG_bgen.dat"
            with open(ibfilename,"r") as ibf:
                firstline=True
                icount=0
                for line in ibf:
                    if firstline :
                        firstline=False
                    else :
                        sline = line.split()
                       # print (sline)


                        istrand = int(np.floor((int(sline[0]) - (basestrandoffset if basestrandoffset != 0 else lstrand)) / float(lstrand) ) + 1)  # find the
                        dG_matrix1[ib,istrand] = float(sline[1])
                        dGt_matrix1[ib, icount] = float(sline[1])
                        icount+=1
                if (icount < nstrandtypes) :
                    raise  ValueError("ERROR: nstrandtypes is larger than the dGdata data length... ")
                    exit(555)

        for ib in range(nblobs):

            ibfilename = ibdirname2 + "_ib" + str(ib) + "/outdG_bgen.dat"
            with open(ibfilename, "r") as ibf:
                firstline = True
                icount = 0
                for line in ibf:
                    if firstline:
                        firstline = False
                    else:
                        sline = line.split()
                        # print (sline)


                        istrand = int(np.floor((int(sline[0]) - (
                        basestrandoffset if basestrandoffset != 0 else lstrand)) / float(
                            lstrand)) + 1)  # find the
                        dG_matrix2[ib, istrand] = float(sline[1])
                        dGt_matrix2[ib, icount] = float(sline[1])
                        icount += 1
                if (icount < nstrandtypes):
                    raise ValueError("ERROR: nstrandtypes is larger than the dGdata data length... ")
                    exit(555)


       # print ("Finished  loading dG data.")

        # convert from NuPack dG to lammps interaction parameter H
     #   freepolymer_energy_seq = (np.sum(np.exp(-dG_matrix[:,:])-1.0e0, axis=1) + math.log(4*math.pi*math.pow(rc,2) * rb / 1.66)) * math.sqrt(2*math.pi)

         # chose random strands and build colloids
        colstrands = []
        # get genome occurance
        gen2str_counts=[]
       # with open("dGdata/" + gen2file[:-11] + "+counts_oct2018.gen", "r") as gen2read :
        with open("dGdata/" + gen2file, "r") as gen2read :
            lcount = 0
            for line in gen2read  :
                if (lcount < nstrandtypes):
                    gen2str_counts.append(math.exp(float(line.split()[1].rstrip("\n"))))
                lcount+=1
        # normalise
        sumgen2str_counts=sum(gen2str_counts)
        ngen2str_counts=[]
        for ii in gen2str_counts:
            ngen2str_counts.append(float(ii)/sumgen2str_counts)
       #  for icol in range(Nfree_col) :
       #      # chose a number of random strands
       #     # icbasevector = np.floor(nstrandtypes*lstrand*np.random.uniform(size=nstrandspercol)) # choose a random base
       #
       #      rndn=np.random.uniform(size=nstrandspercol)
       #      icbasevector=np.zeros(nstrandspercol)
       #      # chose base according to probabilities
       #      for ist in range(nstrandspercol):
       #          com=0
       #          for jst in range(len(ngen2str_counts))  :
       #              com+=ngen2str_counts[jst]
       #              if ( rndn[ist <= com]):
       #                  icbasevector.append
       #
       #
       #
       #
       #      icbasevector=np.floor(nstrandtypes * lstrand*icbasevector)
       #      # end chose base according to probabilities
       #      icsvector = np.floor( (icbasevector - (basestrandoffset if basestrandoffset != 0 else lstrand)) / float(lstrand) ) + 1 # random strand vector
       #      if (any(icsvector > nstrands)) :
       #          raise ValueError("ERROR: chosen strand number larger than the total number of strands, exiting..")
       #      colstrands.append(icsvector)
       #  colstrands = np.array(colstrands, dtype=int)
       #
       #  print(str(nstrandtypes) +" " +str(colstrands) +"\n")
       #     # print (icsvector)
       #     # exit()
       #
       # # print ("Converting delG matrix and obtaining col - blob interactions ...")
       #  # get interaction between all colloids and blobs
       #  cbdG_matrix = np.zeros((blobs_per_freepolymer, Nfree_col))

        # blob_index, col_index = np.mgrid[0:nblobs, 0:Nfree_col]
        # ind=1
        #
        # for icol in range(Nfree_col):
        #    # print(icol/Nfree_col)
        # #    print(float(icol) / Nfree_col)
        #     strand_indices = colstrands[icol]
        #     for iblob in range(blobs_per_freepolymer):
        #
        #     # get effective interaction, calculate the partition function of each blob to the nstrandpercol strands
        #         #Qib = - math.log(np.sum(np.exp([-dG_matrix[iblob,int(xx)] for xx in colstrands[icol]]) - 1.0) + 1.0)
        #         Qib = -np.log((np.exp(-dG_matrix[iblob,strand_indices])-1.0).sum()+1.0)
        #         #print(Qib)
        #         cbdG_matrix[iblob,icol] = Qib
        #       #  print (Qib)



    # get cbdG_matrix strands according to strand probabilities
        cbdG_matrix = np.zeros((blobs_per_freepolymer, 2))

        for iblob in range(blobs_per_freepolymer):

        #    print(str(dGt_matrix[iblob,0:nstrandtypes])+ " " + str(ngen2str_counts))
        #    print(" ")
        #    print(str(np.exp(-dGt_matrix[iblob,0:nstrandtypes])*ngen2str_counts))
        #    exit()



            Qib1 = -np.log(((np.exp(-dGt_matrix1[iblob,0:nstrandtypes])-1.0)*ngen2str_counts*stranddensity*rb*rb).sum()+1.0)
            #print(Qib)
            cbdG_matrix[iblob,0] = Qib1
          #  print (Qib)
            Qib2 = -np.log(((np.exp(-dGt_matrix2[iblob, 0:nstrandtypes]) - 1.0) * ngen2str_counts * stranddensity * rb * rb).sum() + 1.0)
        # print(Qib)
            cbdG_matrix[iblob, 1] = Qib2


        #    print("BBBBB " + str((np.exp(-dGt_matrix[iblob,0:nstrandtypes])-1.0)*ngen2str_counts*nstrands) + " " + str((np.exp(-dGt_matrix[iblob,0:nstrandtypes])-1.0)*ngen2str_counts) + " " + str((np.exp(-dGt_matrix[iblob,0:nstrandtypes])-1.0)))
        #    print ("CCCC " + str(nstrands))
      #  cbH_matrix = fdG_to_H_conv(cbdG_matrix)
        cbH_matrix = - np.minimum((cbdG_matrix + fdG_to_H_conv_cons) * fdG_to_H_stpi_cons, np.zeros(cbdG_matrix.shape))
       # print ("Finished matrix conversion to CG units.")
       # print(cbH_matrix[0,:])
       # print(cbH_matrix.shape)

    else :
        raise ValueError("ERROR: Not imlemented, exiting... ")
        exit(555)
        eneoffset = -5
        enevar = 2
        freepolymer_energy_seq = np.around(np.random.normal(size=blobs_per_freepolymer, scale=enevar) + eneoffset,
                                           decimals=5)

    #freepolymer_sequence = np.array(range(0, blobs_per_freepolymer)) + 1
    freepolymer_sequence = [1]*blobs_per_freepolymer

    if load_bnei_data:

        bneidG = []
        with open(neifile_dG, "r") as ndgf:
            for line in ndgf:
                sline = line.split()
                bneidG.append(float(sline[1]))
        if (len(bneidG)+1 != len(blobdG)) :
            print(str(len(bneidG)) + "  " + str(len(blobdG)))
            raise ValueError("ERROR: length of blob neighbour array must match the blob colloid delG array:" + \
                             str(len(bneidG)) + "  " + str(len(blobdG)) + "  exiting...  ")
        # convert from NuPack dG to lammps interaction parameter H
        bnei_energy_seq = (np.array(bneidG) + math.log(math.pow(rb, 3) / 1.66) + 2.0) * math.sqrt(0.8/0.5)


    # FREE COLLOIDS
    #rcol=np.around(rc/rb, decimals=2)     # 300 bp blobs -> 10nm Rg blobs, 1000bp blobs -> 21nm blobs ####  colloid size should be scaled accordingly




    # ############################################################################
    # #####   PRINT   #####   MAKE DATA.FILE    ##################################
    # with open(datafilename, "w") as outdatafile:
    #
    #     if polymer_form=='circular':
    #       #  print (print_lammps_head(totnblobs+Nfree_col, (blobs_per_freepolymer)*Nfree_polymers,
    #       #                  max_ntypes, 2,lbox))
    #         outdatafile.write(print_lammps_head(totnblobs+Nfree_col, (blobs_per_freepolymer)*Nfree_polymers,
    #                       max_ntypes, 2,lbox))
    #     else:
    #      #   print (print_lammps_head(totnblobs+Nfree_col, (blobs_per_freepolymer-1)*Nfree_polymers,
    #      #                   max_ntypes, 2,lbox))
    #         outdatafile.write(print_lammps_head(totnblobs+Nfree_col, (blobs_per_freepolymer-1)*Nfree_polymers,
    #                         max_ntypes, 2,lbox))
    #
    #     tsome_polymers_and_cols(Nfree_polymers, blobs_per_freepolymer, freepolymer_sequence, Nfree_col, rcol, \
    #                             lbox, form=polymer_form, outdatafile=outdatafile)



    # write to the coefficients PARM.FILE
    with open("coefficients.dat", "w") as cparmfile :
        for ipoly in range(Nfree_polymers):
            for iblob in range(blobs_per_freepolymer) :
                sline = str(np.around(cbH_matrix[iblob,ipoly],decimals=2))
                cparmfile.write(sline + "\n")
    with open("dGdata.dat", "w") as dgparmfile :
        for ipoly in range(Nfree_polymers):
            for iblob in range(blobs_per_freepolymer) :
                jcol=0
                sline = str(np.around(cbdG_matrix[iblob,ipoly],decimals=2))
                dgparmfile.write(sline + "\n")

    subprocess.run([r"cp dGdata.dat dGdata_" +  str(ldirname) + ".dat"], shell=True)
    

    # Update ndraw_tcl vmd script
   # subprocess.run([r"cat draw_tsoftDNA_vmd.tcl | sed -e 's/EXTRAVMDTCLLINE11/mol modselect 1 0 type > " \
   #                 + str(blobs_per_freepolymer) + r"/' > draw_softDNA_vmd.tcl"], shell=True)


    # make lammps input file:
    rcommand = r"python3 get_input_2cg.py " + str(blobs_per_freepolymer) + " " + str(lencutgenseg) + " " + str(lboxxy) + " " + str(lboxz)  + " coefficients.dat > data.genome"
    print("making lammps input with command: " + rcommand)
    subprocess.run([rcommand],shell=True)

    ##   make a directory and copy everything there
    subprocess.run([r"mkdir " + str(ldirname)], shell=True)
    subprocess.run([r"cp data.genome dGdata.dat coefficients.dat in.genome jtypesmake* " + str(ldirname)], shell=True)

    ######## END ##################

class Counter(object):
    def __init__(self):
        self.count = 1

    def increment(self, incr=0):
        self.count += incr
        return self.count - incr


class LAMMPSCounters(object):
    def __init__(self):
        self.atoms = Counter()
        self.bonds = Counter()
        self.angles = Counter()
        self.dihedrals = Counter()
        self.molecules = Counter()
        self.polytype = Counter()
        self.coltype = Counter()

    def increment_atoms(self, incr): return self.atoms.increment(incr)

    def increment_bonds(self, incr): return self.bonds.increment(incr)

    def increment_angles(self, incr): return self.angles.increment(incr)

    def increment_dihedrals(self, incr): return self.dihedrals.increment(incr)

    def increment_molecules(self, incr): return self.molecules.increment(incr)

    def increment_polytype(self, incr): return self.polytype.increment(incr)

    def increment_coltype(self, incr): return self.coltype.increment(incr)


def vec_random_ndim(n):
    """n-dimensional uniform random unit vector"""
    v = np.random.normal(size=n)
    v /= np.linalg.norm(v)
    return v


class Colloid(object):
    def __init__(self, counters, radius=10.0):
        self.atom_id = counters.increment_atoms(incr=1)
        self.molecule_id = counters.increment_molecules(incr=1)
        self.counters = counters
        self.radius = radius
        self.coordinates = np.zeros([1, 3])
        self.polymers = []

    def add_polymers(self, num_polymers, chain_length, sequence):
        for polymer in range(num_polymers):
            self.add_polymer(chain_length, sequence)

    def add_polymer(self, chain_length, sequence):
        polymer = Polymer(self.counters, chain_length, colloid_radius=self.radius, sequence=sequence, molecule_id=self.molecule_id)
        polymer.displace(self.coordinates.flatten())
        self.polymers.append(polymer)

    def print_atoms(self):
        colloid_repr = ""
        atom_type = 2
        #atom_type = self.atom_id   ### HACK
        # for natom, coordinate in enumerate(self.coordinates):
        #     print coordinate
        #     print self.coordinates
        #     colloid_repr += "{:5n}{:5n}{:5n}{:10.5f}{:10.5f}{:10.5f}\n".format(self.atom_id + natom, self.molecule_id,
        #                                                                        atom_type, * coordinate)
        colloid_repr += "{:7n}{:7n}{:7n}{:11.3f}{:11.3f}{:11.3f}\n".format(self.atom_id, self.molecule_id, atom_type,
                                                                           *self.coordinates)
        for polymer in self.polymers:
            colloid_repr += polymer.print_atoms()
        return colloid_repr

    def print_bonds(self):
        colloid_repr = ""
        for polymer in self.polymers:
            colloid_repr += polymer.print_bonds(grafted=True)
        return colloid_repr

    def displace(self, displacement):
        self.coordinates += displacement
        for polymer in self.polymers:
            polymer.displace(displacement)


class Polymer(object):
    def __init__(self, counters, chain_length, colloid_radius=0, sequence=None, molecule_id=None, form='linear'):
        self.counters = counters
        # self.chain_length = chain_length
        # self.colloid_radius = colloid_radius
        self.coordinates = np.zeros([chain_length, 3])
        self.form = form
        self.bonds = []
        self.atom_id = counters.increment_atoms(incr=chain_length)
        if self.form == 'circular' :
            self.bond_id = counters.increment_bonds(incr=chain_length)
        else :
            self.bond_id = counters.increment_bonds(incr=chain_length - 1)
        if molecule_id is None:
            molecule_id = counters.increment_molecules(incr=1)
        self.molecule_id = molecule_id
        self.build_polymer(chain_length, colloid_radius,form)
        if sequence is None:
           # print '# Error: blob type sequence length not provided, defautling to type1..'
            sequence = chain_length * [1]
        self.sequence = sequence
        if len(sequence) < chain_length :
           # print '# Error: blob type sequence length is too short: len(sequence) < chain_length, defaulting to type 1...'
            for i in range(chain_length - len(sequence)) :
                sequence.append(1)


    def build_polymer(self, chain_length, colloid_radius, form):

        if (form=='linear') :
            unit = np.array([1.0,0,0])
            uvector =  unit / chain_length*math.pow(chain_length,0.588)

            self.coordinates[:, :] = (np.arange(chain_length) * np.vstack((chain_length) * [uvector]).T).T

            CM = np.sum(self.coordinates[:, :], axis=0) / chain_length  # centre of mass)
            self.coordinates[:, :] = self.coordinates[:, :] - CM
           # self.coordinates += colloid_radius * uvector

        elif (form =='circular') :
            unit = np.array([1.0, 0, 0])
            uvector = unit / chain_length * math.pow(chain_length, 0.588)

            self.coordinates[:, :] = 2*(np.arange(chain_length) * np.vstack((chain_length) * [uvector]).T).T


            CM = np.sum(self.coordinates[:, :], axis=0) / chain_length  # centre of mass)
            self.coordinates[:, :] = self.coordinates[:, :] - CM
            # self.coordinates += colloid_radius * uvector

            ## make a closed circle
            self.coordinates[round(chain_length / 2):, 0] = -self.coordinates[round(chain_length / 2):, 0]


            #   print str(unit)
         #   print str(vec_random_ndim(3))
        elif (form=='randomwalk'):
            self.coordinates[0, :] = np.zeros(3)
            for ii in range(1,chain_length) :
                self.coordinates[ii,:] = self.coordinates[ii-1,:] + vec_random_ndim(3)
            scale_length = 1/1 * math.pow(chain_length, 0.088)  # difference between ideal and self-awoiding walk radius of gyration
            CM = np.sum(self.coordinates[:, :], axis=0) / chain_length # centre of mass
            self.coordinates[:, :] = scale_length*(self.coordinates[:, :] - CM )
    # move the polymer and put the centre of mass in [0,0,0]

        #for ii in range(chain_length) : #randomise blob positions a bit & shift polymer to the centre of box
        #    self.coordinates[ii,:] = self.coordinates[ii,:]+2*vec_random_ndim(3) - median_blob_pos

        #



        for (u, v) in zip(range(chain_length-1), range(1, chain_length)):
            self.bonds.append((u + self.atom_id, v + self.atom_id))

        # add a last bond if we have a circular polymer
        if self.form == 'circular':
            self.bonds.append((chain_length-1 + self.atom_id, self.atom_id))



    def print_atoms(self):
        polymer_repr = ""
        for natom, (coordinate, atom_type) in enumerate(zip(self.coordinates, self.sequence)):
            polymer_repr += "{:7n}{:7n}{:7n}{:11.3f}{:11.3f}{:11.3f}\n".format(self.atom_id + natom, self.molecule_id,
                                                                               atom_type, *coordinate)
        return polymer_repr

    def print_bonds(self, grafted=False):
        polymer_repr = ""
        if grafted:
            bond_type = 1
        else:
            bond_type = 2
        for nbond, (u, v) in enumerate(self.bonds):
            polymer_repr += "{:7n}{:7n}{:7n}{:7n}\n".format(self.bond_id + nbond, bond_type, u, v)
            bond_type = 2
        return polymer_repr

    def displace(self, displacement):
        self.coordinates += displacement


def dna_coated_colloid(num_colloids, polymers_per_colloid, blobs_per_polymer, polymer_sequence):
    counters = LAMMPSCounters()

    colloids = []
    for i in range(num_colloids):
        colloid = Colloid(counters, 3.0)
        colloid.add_polymers(polymers_per_colloid, blobs_per_polymer, polymer_sequence)
        colloid.displace(np.random.rand(3)*blobs_per_polymer*6)
        colloids.append(colloid)
    print ("Atoms\n")
    for colloid in colloids:
        print (colloid.print_atoms().rstrip())
    print ("\nBonds\n")
    for colloid in colloids:
        print (colloid.print_bonds().rstrip())


def tsome_polymers_and_cols(num_polymers,  blobs_per_polymer, polymer_sequence,  num_cols, rcol, lbox, form, outdatafile):
    counters = LAMMPSCounters()

    polymers = []
    colloids = []
    for i in range(num_polymers):
        polymer = Polymer(counters=counters, chain_length=blobs_per_polymer, sequence=polymer_sequence, form=form)

        polymers.append(polymer)


    if (polymers[0].coordinates[0,0] < -lbox) : # check x-axis of the first blob 
        raise ValueError("ERROR: BOX size too small to fit the polymer, increase lbox.. \
        exiting... \n")
        
                

    for i in range(num_cols):
        iColloid = Colloid(counters=counters, radius=rcol)

        overlap = True
        while (overlap) :
            newpos = (np.random.uniform(size=3) * 2.0 - 1.0) * lbox
            overlap = check_overlap(colloids, newpos, rcol, lbox) # check for overlap

        iColloid.coordinates = newpos
        colloids.append(iColloid)

#    print ("Atoms\n")
    outdatafile.write("\n Atoms \n \n")
    for polymer in polymers:
     #   polymer.displace(polymer.coordinates[])
 #       print (polymer.print_atoms().rstrip())
         outdatafile.writelines(polymer.print_atoms().rstrip()+"\n")

    for colloid in colloids :
    #    print (colloid.print_atoms().rstrip())
        outdatafile.writelines(colloid.print_atoms().rstrip()+"\n")

   # print ("\nBonds\n")
    outdatafile.write("\n \n Bonds\n \n")
    for polymer in polymers:
    #    print (polymer.print_bonds().rstrip())
         outdatafile.write(polymer.print_bonds().rstrip() + "\n")
def print_lammps_head(natoms,nbonds,natom_types,nbond_types,lbox) :
    header = ""
    header += "Lammps data file generated by JDF&TC python script for softDNA\n \n"
    header += "{:<7n}{:<7}{:<7n}{:<7}\n\n".format(natoms, ' atoms\n', nbonds, ' bonds')
    header += "{:<7n}{:<10}{:<7n}{:<10}\n".format(natom_types, ' atom types\n', nbond_types, ' bond types\n')
    header += "{:<8.3f}{:1}{:<8.3f}{:<8}\n".format(-lbox, ' ',lbox, ' xlo xhi')
    header += "{:<8.3f}{:1}{:<8.3f}{:<8}\n".format(-lbox, ' ',lbox, ' ylo yhi')
    header += "{:<8.3f}{:1}{:<8.3f}{:<8}\n\n".format(-lbox, ' ',lbox, ' zlo zhi')
    return header


def check_overlap(colloids, newpos, rcol, lbox) :

    ncol = len(colloids)

    # check xy plane (where polymers lie)
    if (np.abs(newpos[2]) < rcol + 2) :
        return True

    rcol2=rcol*rcol
    for icol in range (ncol) : # check all colloids
        distij = np.array(newpos) - np.array(colloids[icol].coordinates)
        distij = distij - 2*lbox * (np.rint(distij/(2*lbox))) # periodic boundary
        distij2 = sum (distij * distij)
        if (distij2 < 4*rcol2) :
            return True

    return False

if __name__ == "__main__":
    main()
