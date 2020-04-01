import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="calculate the overlap score between a genome and all length l probes")
    parser.add_argument('length', type=int, help='probe length')
    parser.add_argument('threads', type=int, help='number of openmp threads to use', default=1)
    parser.add_argument('-n', type=int, default=None, help=('only consider the n most common subsequences'
                                                            ' in the genomes'))
    parser.add_argument('-r', '--rate', type=float, default=1, help='rate parameter in the exponential weight')
    parser.add_argument('--compare', action='store_true', help='compare the scores of all pairs of genomes')
    parser.add_argument('--file', default=None, help='output file')
    parser.add_argument('genomes', nargs="+", help=('one or more genbank records (.gb)'
                                                    ' or one-line files containing genomes (.gen)'))
    return parser


def main():
    from genometargeting import exponential_score_factory, Genome, Computer

    parser = get_parser()
    args = parser.parse_args()

    gen_files = args.genomes
    length = args.length
    threads = args.threads
    n = args.n
    rate = args.rate
    compare_flag = args.compare
    file = args.file
    score_function = exponential_score_factory(rate, length)

    genomes = *map(Genome.read, gen_files),
    compute = Computer(genomes, length, score_function)
    compute(threads, n, compare_flag, file)


if __name__ == '__main__':
    main()
