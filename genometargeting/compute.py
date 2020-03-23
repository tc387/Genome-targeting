import sys
from collections import abc
from itertools import combinations
from itertools import product

import numpy as np
from genometargeting._compute import compute_full, compute_only

from genometargeting import ints_to_string, reverse_complement


class Computer:
    def __init__(self, genomes, length, score_function):
        if not isinstance(genomes, abc.Iterable):
            genomes = [genomes]
        self.genomes = genomes
        self.length = length
        self.table = self.make_table(length, score_function)
        ngen = len(genomes)
        self.compute_fmt_header = f'{{:{str(length)}}}' + ngen * ' {:>20}' + '\n'
        self.compute_fmt = f'{{:{str(length)}}}' + ngen * ' {:20.10f}' + '\n'
        self.compare_fmt_header = f'{{:{str(length)}}}' + (ngen * (ngen - 1) // 2) * ' {:>20}' + '\n'
        self.compare_fmt = f'{{:{str(length)}}}' + (ngen * (ngen - 1) // 2) * ' {:20.10f}' + '\n'

    @staticmethod
    def make_table(length, score_function):
        """
        for all possible comparisons between length "length" strands, compute and store
        the associated score
        """
        table = np.zeros(2 ** length, dtype=np.int64)
        for h, bits in enumerate(product([0, 1], repeat=length)):
            table[h] = score_function(bits)
        return table

    def __call__(self, threads, n=None, compare=False, file=None):
        if n is not None:
            computers = self.get_computers_top(threads, n)
        else:
            computers = self.get_computers_full(threads)
        if compare:
            return self.compare(computers, file)
        else:
            return self.compute(computers, file)

    def get_computers_full(self, threads):
        length = self.length
        transcribed = [g.transcribe(length) for g in self.genomes]
        return [compute_full(length, *gs, self.table, threads) for gs in transcribed]

    def get_computers_top(self, threads, n):
        length = self.length
        transcribed = [g.transcribe(length) for g in self.genomes]
        top = *map(np.copy, np.hstack([g.most_common(length, n) for g in self.genomes])),
        return [compute_only(length, *gs, *top, self.table, threads) for gs in transcribed]

    def compute(self, computers, file=None):
        with open(file, 'w') if file is not None else sys.stdout as fh:
            fh.write(self.compute_fmt_header.format('probe', *self.genomes))
            for results in zip(*computers):
                pu, pl, _ = results[0]
                probe = ints_to_string(pu, pl, self.length)
                scores = *(np.log(r[2]) for r in results),
                fh.write(self.compute_fmt.format(probe, *scores))
                fh.write(self.compute_fmt.format(reverse_complement(probe), *scores))

    def compare(self, computers, file):
        genomes = *combinations(self.genomes, 2),
        with open(file, 'w') if file is not None else sys.stdout as fh:
            fh.write(self.compare_fmt_header.format('probe', *tuple(n[0] for n in genomes)))
            fh.write(self.compare_fmt_header.format('', *tuple(n[1] for n in genomes)))
            for results in zip(*computers):
                pu, pl, _ = results[0]
                probe = ints_to_string(pu, pl, self.length)
                scores = *(np.log(r[2]) for r in results),
                scores = *((s0 - s1) for (s0, s1) in combinations(scores, 2)),
                fh.write(self.compare_fmt.format(probe, *scores))
                fh.write(self.compare_fmt.format(reverse_complement(probe), *scores))
