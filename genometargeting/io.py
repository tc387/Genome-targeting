import numpy as np
from joblib import Memory

from genometargeting import transcribe, most_common

location = './cachedir'
memory = Memory(location, verbose=0)


class Genome:
    def __init__(self, genome_str, name):
        self.string = genome_str
        self.name = name

    @classmethod
    def read_gen(cls, gen_file, name=None):
        genome_str = open(gen_file, 'r').readline().strip().upper()
        if name is None:
            name = gen_file.split('/')[-1].split('.')[0]
        return Genome(genome_str, name)

    def __len__(self):
        return len(self.string)

    def __hash__(self):
        return hash(self.string)

    def __getitem__(self, item):
        return self.string[item]

    def __str__(self):
        return f'{self.name}'

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
