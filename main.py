import gen

#parametry
n = 100 #długość sekwencji
k = 10   #długość podciągu
ne = 0 #błędy negatywne
pe = 0  #błędy pozytywne
pop = 10 #populacja
#
oligo_count = n - (k - 1)

#seq = gen.dna(n) # to jest losowe
#print(seq)
seq = "ACCCGATGTGTCTAATAAGCTGTACAGTGTCCATTGCTTCGGACTTCCGGTTCGGCATTGGACAGTCGGAATCTTTATAGACTCATTACAGGCGAGGCCA" # a to nie
spectrum = gen.spectrum(seq, k, n)
original_seq_spectrum = spectrum
# print(''.join(seq))
# dodanie błędów
# błędy negatywne
gen.negative_errors(spectrum, ne)
# błędy pozytywne
gen.positive_errors(spectrum, pe, k)

start_oligo = original_seq_spectrum[0]
fake_list = []  # save_to() opgens on lists so we fool it with a fake one to save one element
fake_list.append(start_oligo)
gen.save_to('start.txt', fake_list)
print(f'Oligonukleotyd startowy: {start_oligo}')
# print(spectrum)
# sortowanie w celu pomieszania sekwnecji
spectrum.sort()

oligo_counts_dict = gen.oligo_counter(spectrum)

print(f"To nasz super słownik policzony{oligo_counts_dict}")


gen.save_to('spectrum.txt', set(spectrum))
# print(set(spectrum))

#graph itd
graph = gen.create_graph(spectrum)
paths = []
path = []
paths.append(gen.create_paths(graph, start_oligo, path, oligo_count, oligo_counts_dict))
# paths.append(gen.create_path_using_weights(graph, start_oligo, path, oligo_count, oligo_counts_dict))
print(f"To nasz super słownik policzony2{oligo_counts_dict}")

first_path = paths[0]
for x in range(1, pop):
     paths.append(gen.mutate_path(graph, first_path, oligo_count, oligo_counts_dict))

print(f'Przewidziane sciezki\n')

for solution in paths:
    print(f'{solution}')
    solution_seq = gen.build_sequence(graph, solution)
    print(f'{solution_seq}')
    coverage_sum = gen.evaluate(graph, solution)
    similarity = gen.compare(seq, solution_seq)
    print(f'Suma pokrycia: {coverage_sum}\t'
          f'Średnie pokrycie: {round(coverage_sum/oligo_count, 2)}\t'
          f'Podobieństwo do oryginalnej sekwencji: {similarity}\n')





