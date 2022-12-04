import random
nucs = ['A', 'C', 'G', 'T']

def save_to(filename ,to_save):
    f1 = open(filename, 'w')
    for stuff in to_save:
        f1.writelines(stuff + '\n')
    f1.close()

def dna(n):
    return ''.join(random.choices(nucs, k=n))

def spectrum(seq, k, n):

    spectrum = []
    for i in range(n - (k-1)):
        spectrum.append(seq[i:i+k])

    return spectrum

def negative_errors(spectrum, error_count):
    removed_oligos = []  # storing oligonucleotides removed from spectrum

    for i in range(error_count):
        removed_oligo = random.choice(spectrum)
        removed_oligos.append(removed_oligo)
        spectrum.remove(removed_oligo)

    print(f'Negative: {removed_oligos}')
    save_to('negative_errors.txt', removed_oligos)

    return spectrum


def positive_errors(spectrum, error_count, k):
    added_oligos = []  # storing oligonucleotides added to spectrum

    for i in range(error_count):
        while random_k_mer := random.choices(nucs, k=k):  # actually a valid operator that allows for assignment of variables within expressions
            if not spectrum.count(random_k_mer):
                break
        added_oligos.append(''.join(random_k_mer))
        spectrum.append(''.join(random_k_mer))

    print(f'Positive: {added_oligos}')
    save_to('positive_errors.txt', added_oligos)

    return spectrum