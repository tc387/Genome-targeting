from operator import mul

FAC = [sum([4 ** (m - j + 1) * j for j in range(1, m + 1)]) for m in range(1, 30)]


def get_counts(bits):
    """
    determines the numbers of length m runs of consecutive ones
    in a list of bits

    the result is returned as a dictionary keyed by m <= len(bits)
    """
    scores = {a: 0 for a in range(1, len(bits) + 1)}
    counter = 0
    for b in bits:
        if b:
            counter += 1
        elif counter > 0:
            scores[counter] += 1
            counter = 0
    if b:
        scores[counter] += 1
    return scores


def power_four_score(bits):
    """
    given a list of bits representing the logical conjunction of
    two objects of equal length, return the overlap score defined in
    <doi to paper>
    """
    counts = get_counts(bits)
    return sum(map(mul, counts.values(), FAC))
