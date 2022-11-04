import random
nucs = ['A', 'C', 'G', 'T']
def dna(n):
    return ''.join(random.choices(nucs, k=n))

def spectrum(seq, k, n):

    spectrum = []
    for i in range(n - (k-1)):
        spectrum.append(seq[i:i+k])

    return spectrum

def negative_errors(spectrum, error_count):

    for i in range(error_count):
        spectrum.remove(random.choice(spectrum))

    return spectrum


def positive_errors(spectrum, error_count, k):

    for i in range(error_count):
        while random_k_mer := random.choices(nucs, k=k):
            if not spectrum.count(random_k_mer):
                break
        spectrum.append(''.join(random_k_mer))

    return spectrum