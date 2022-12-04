import generate

#parametry
n = 20
k = 3
ne = 2
pe = 2
#

seq = generate.dna(n)
spectrum = generate.spectrum(seq, k, n)
print(''.join(seq))
# sortowanie w celu pomieszania sekwnecji
start_oligo = spectrum[0]
fake_list = []  # save_to() operates on lists so we fool it with a fake one to save one element
fake_list.append(start_oligo)
generate.save_to('start.txt', fake_list)
spectrum.sort()
print(f'Oligonukleotyd startowy: {start_oligo}')
print(spectrum)
# dodanie błędów
# błędy negatywne
generate.negative_errors(spectrum, ne)
generate.positive_errors(spectrum, pe, k)
spectrum.sort()

print(set(spectrum))



