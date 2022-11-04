import generate

#parametry
n = 20
k = 4
ne = 0
pe = 2
#

seq = generate.dna(n)
spectrum = generate.spectrum(seq, k, n)
print(''.join(seq))
# sortowanie w celu rozpierdolenia sekwnecji
spectrum.sort()
print(spectrum)
# dodanie błędów

generate.negative_errors(spectrum, ne)
generate.positive_errors(spectrum, pe, k)
spectrum.sort()

print(set(spectrum))



