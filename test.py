import math


test_profile = {
        'A': {0: 0.5, 1: 0.25, 2: 0.95},
        'T': {0: 0.1, 1: 0.25, 2: 0.03},
        'C': {0: 0.1, 1: 0.25, 2: 0.01},
        'G': {0: 0.3, 1: 0.25, 2: 0.01}
    }

def entropy(profile: dict[str, dict[int, float]], k) -> int:
    e_total = 0
    for pos in range(k):
        for base in ('A', 'T', 'C', 'G'):
            prob = profile[base][pos]
            e_total += -prob * math.log2(prob)

    return e_total

print(entropy(test_profile,3))