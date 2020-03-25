import numpy as np
from joblib import Memory

from .utils import transcribe, most_common

location = './cachedir'
memory = Memory(location, verbose=0)


class Genome:
    def __init__(self, genome_str, name):
        self.string = genome_str
        self.name = name

    @classmethod
    def read(cls, gen_file, filetype=None):
        readers = {'gen': cls.read_gen,
                   'gb': cls.read_genbank}
        if filetype is None:
            filetype = gen_file.strip().split('.')[-1]
        try:
            return readers[filetype](gen_file)
        except IndexError:
            raise ValueError("only genbank records (.gb) and single-line files (.gen) are supported")

    @classmethod
    def read_gen(cls, gen_file, name=None):
        genome_str = open(gen_file, 'r').readline().strip().upper()
        if name is None:
            name = gen_file.split('/')[-1].split('.')[0]
        return Genome(genome_str, name)

    @classmethod
    def read_genbank(cls, genbank_file, name=None):
        with open(genbank_file, 'r') as f:
            line = ""
            genome_str = ""
            while line.strip().upper() != "ORIGIN":
                line = f.readline()
            line = f.readline()
            while line.strip() != "//":
                genome_str += "".join(line.strip().split()[1:]).upper()
                line = f.readline()
        if name is None:
            name = genbank_file.split('/')[-1].split('.')[0]
        return Genome(genome_str, name)

    def __len__(self):
        return len(self.string)

    def __hash__(self):
        return hash(self.string)

    def __getitem__(self, item):
        return self.string[item]

    def __str__(self):
        return f'{self.name}'

    def __eq__(self, other):
        return self.string == other.string

    def __format__(self, format_spec):
        return format(self.name, format_spec)

    def transcribe(self, length):
        _transcribe = memory.cache(transcribe)
        gu, gl = _transcribe(self, length).T
        return np.copy(gu), np.copy(gl)

    def most_common(self, length, n):
        _most_common = memory.cache(most_common)
        return _most_common(self, length, n)


def main():
    gen_file = '../genomes/rBS_QB928.gen'
    reader = Genome.read_gen(gen_file)
    reader.transcribe(10)


if __name__ == '__main__':
    main()
