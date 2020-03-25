"""
Computes the exponential of the overlap score of probe p on genome G.

Each comparison between a length l probe p and a similar length
subsequence g of G can be represented as a binary integer m,
one for a match and zero for a miss, e.g.

p =   ATCCGTCGAA
g =   AGGTCGCGAG
m = 0b1000001110

0 <= m < 2**l .

The array "scores" is a length l array containing every possible score,
"compute" is agnostic of the scoring function used, assuming the
overlap score of p on G is a sum of scores of p on the {g}.

To compute m efficiently, p is expressed as two integers pu and pl,
the upper and lower bits of the binary representation of each base
in the probe, where

C = 00; T = 01; A = 10; G = 11

e.g. for p = ATCCGTCGAA, g = AGGTCGCGAG

pu = 0b1000100111 = 551
pl = 0b0100110100 = 308

gu = 0b1110010111 = 919
gl = 0b0111010101 = 469

To get the score between p and a g in G from these integers, m is computed
by computing the logical equality of pu and gu, implemented as bitwise XOR
followed by negation,

!(pu ^ gu) =  0b1111111111 -
             (0b1000100111 ^
              0b1110010111)
           =  0b1001001111

The same is done for pl and gl,

!(pl ^ gl) =  0b1111111111 -
             (0b0100110100 ^
              0b0111010101)
           =  0b1100011110

m is then the logical conjunction of the two intermediate results,

m = 0b1001001111 &
    0b1100011110
  = 0b1000001110

The overlap score between p and g is then the m^th element of "scores",
and the total overlap score between p and G the sum of the scores over
the {g}.

"""


cimport cython
from libcpp cimport bool
from cython.parallel import prange


cdef int reverseCompBits(int n, int length) nogil:
    cdef int rev = 0
    while length > 0:
        rev <<= 1
        if n & 1 == 0:
            rev |= 1
        n >>= 1
        length -= 1
    return rev


cdef bint lessThanRevComp(int pu, int pl, int length) nogil:
    cdef int ru = reverseCompBits(pu, length)
    cdef int rl = reverseCompBits(pl, length)
    return (pu << length) + pl <= (ru << length) + rl


@cython.boundscheck(False)
@cython.wraparound(False)
def compute_full(int length, int [::1,] gu , int [::1,] gl, long [::1,] scores, int omp_threads=1):
    cdef Py_ssize_t subseq = gu.size
    cdef Py_ssize_t ss, pu, pl
    cdef long score
    cdef int bmax = 2**length-1
    pu = pl = 0
    score = 0

    for pu in range(bmax+1):
        for pl in range(bmax+1):
            if lessThanRevComp(pu, pl, length):
                score = 0
                for ss in prange(subseq, nogil=True, num_threads=omp_threads, schedule="static"):
                    score += scores[(bmax - (pu ^ gu [ss])) & (bmax - (pl ^ gl[ss]))]
                yield pu, pl, score


@cython.boundscheck(False)
@cython.wraparound(False)
def compute_only(int length, int [::1,] gu, int [::1,] gl, int [::1,] ppu, int [::1,] ppl, long [::1,] scores, int omp_threads=1):
    cdef Py_ssize_t subseq = gu.size
    cdef Py_ssize_t psubseq = ppu.size
    cdef Py_ssize_t ssi, ss, pu, pl
    cdef long score
    cdef int bmax = 2**length-1
    pu = pl = 0
    score = 0

    for ssi in range(psubseq):
        pu = ppu[ssi]
        pl = ppl[ssi]
        if lessThanRevComp(pu, pl, length):
            score = 0
            for ss in prange(subseq, nogil=True, num_threads=omp_threads, schedule="static"):
                score += scores[(bmax - (pu ^ gu[ss])) & (bmax - (pl ^ gl[ss]))]
            yield pu, pl, score
