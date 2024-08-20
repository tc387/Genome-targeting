# Genome-targeting
The folder genomes/ contains all the genomes used in this work
The folder USER-SOFT-BLOB contains the LAMMPS implementation of the soft blob potential
genome_overlap.py computes the overlap score for some set of probes and genomes.

### Compiling
To use compute the overlap score, first run

    python3 setup.py build_ext -i

to compile the cython extensions, and then run genome_overlap.py with the appropriate command-line arguments.

### All probes or some probes?
Computing the score function for every length-10nt probe (as we did) is feasible for most bacterial genomes.

Computing the score function for every length-20nt probe (as we did NOT do) is both impractical (there are 1,099,511,627,776 of them) and unnecessary;
the best scoring probe is most likely to be found among the most commonly-occuring subsequences of the genome
(option --number-of-probes).

For example, of the top 5000 most commonly-occuring length-10 subsequences of _E. coli_ strain bl21-de3 genome, the top-scoring sequence is the 14th most common.

However, of the top 5000 most commonly-occuring length-20 subsequences in the same genome, the top-scoring sequence is **the** most common.

The fraction of all possible probes you need to test (probably) decreases with increasing probe length and increasing genome length (don't quote me on that). 

### Example usage
Running:

    python genome_overlap.py 10 8 genomes/rEC_bl21-de3.gen -c -n 20

should yield:

    JOB INFO:
        number of threads: 8
        probe length: 10
        number of probes: 20 most common
        genomes: ['genomes/rEC_bl21-de3.gen']
        score function: power_four_score (i.e., the one used in Curk et al, PNAS (2020))
        reverse complement?: yes
        output file: stdout
        
    probe              rEC_bl21-de3
    GCCGGATGCG        20.8797256088
    CGCATCCGGC        20.8797256088
    CTGGCGCTGG        21.2309825603
    CCAGCGCCAG        21.2309825603
    GCATCCGGCA        20.9125335234
    TGCCGGATGC        20.9125335234
    GCTGGCGCTG        21.2032711522
    CAGCGCCAGC        21.2032711522
    GCCGCATCCG        20.8121753470
    CGGATGCGGC        20.8121753470
    GGATGCGGCG        20.8510136192
    CGCCGCATCC        20.8510136192
    GGCGCTGGCG        21.2615020644
    CGCCAGCGCC        21.2615020644
    GATGCGGCGT        20.8617403144
    ACGCCGCATC        20.8617403144
    GCCTGATGCG        20.8592037187
    CGCATCAGGC        20.8592037187
    CCGGATGCGG        20.7754196406
    CCGCATCCGG        20.7754196406
