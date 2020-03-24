from collections import deque, Counter
from functools import reduce
from itertools import chain
from itertools import product
from operator import add

import numpy as np

BASE_MAP = {'A': (1, 0), 'T': (0, 1), 'C': (0, 0), 'G': (1, 1)}
REV_BASE_MAP = {(1, 0): 'A', (0, 1): 'T', (0, 0): 'C', (1, 1): 'G'}
RC_MAP = dict(zip("ATCG", "TAGC"))
VERBOSE = False


def string_to_ints(string):
    """
    converts a base dna string into upper + lower bit integers
    """
    upper, lower = *zip(*map(BASE_MAP.get, string)),
    return int(reduce(add, map(str, upper)), base=2), int(reduce(add, map(str, lower)), base=2)


def ints_to_string(upper, lower, length):
    return "".join(ints_to_string_helper(upper, lower, length))


def ints_to_string_helper(upper, lower, length):
    if length == 0:
        return []
    return ints_to_string_helper(upper >> 1, lower >> 1, length - 1) + [REV_BASE_MAP[(upper & 1, lower & 1)]]


def make_table(length, score_function):
    """
    for all possible comparisons between length "length" strands, compute and store
    the associated score
    """
    table = np.zeros(2 ** length, dtype=np.int64)
    for h, bits in enumerate(product([0, 1], repeat=length)):
        table[h] = score_function(bits)
    return table


def reverse_complement(dna_string):
    """
    return the reverse complement of a string representing a DNA sequence
    """
    return "".join(map(RC_MAP.get, dna_string))[::-1]


def frame_iteration(sequence, length):
    """
    given an iterable sequence, efficiently generate "length" long
    overlapping subsequences
    """
    head = deque(sequence[:length])
    tail = sequence[length:]
    total = len(sequence) - length + 1
    d = total // 100
    yield string_to_ints(head)
    for i, element in enumerate(tail):
        head.popleft()
        head.append(element)
        yield string_to_ints(head)
        if VERBOSE and i % d == 0:
            print(f'{i}, {i // d:.2f}% complete')


def transcribe(genome, length):
    gen = frame_iteration(genome, length)
    genrc = frame_iteration(reverse_complement(genome), length)
    all_subs = np.array(list(chain(gen, genrc)), dtype=np.int32)
    all_subs.flags.writeable = False
    return all_subs


def most_common(genome, length, n):
    seqs = Counter(zip(*genome.transcribe(length))).most_common(n)
    return np.array([a for a, _ in seqs], dtype=np.int32).T
    # return Counter(map(tuple, genome.transcribe(length))).most_common(n)
