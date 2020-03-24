# Genome-targeting
The folder genomes/ contains all the genomes used in this work
The folder USER-SOFT-BLOB contains the LAMMPS implementation of the soft blob potential

genome_overlap.py computes the overlap score for some set of probes and genomes. To use it first run

    python3 setup.py build_ext -i

to compile the cython extensions, and then run genome_overlap.py with the appropriate command-line arguments.
