import gen

#parametry
n = 30
k = 5
ne = 2
pe = 2
#

#seq = gen.dna(n) # to jest losowe
seq = "AAATTAAATTGACCCAAAATCTCCGTGAGT" # a to nie
spectrum = gen.spectrum(seq, k, n)
# print(''.join(seq))
# dodanie błędów
# błędy negatywne
gen.negative_errors(spectrum, ne)
# błędy pozytywne
gen.positive_errors(spectrum, pe, k)

start_oligo = spectrum[0]
fake_list = []  # save_to() opgens on lists so we fool it with a fake one to save one element
fake_list.append(start_oligo)
gen.save_to('start.txt', fake_list)
# print(f'Oligonukleotyd startowy: {start_oligo}')
# print(spectrum)
# sortowanie w celu pomieszania sekwnecji
spectrum.sort()

oligo_counts_dict = gen.oligo_counter(spectrum)

# print(f"To nasz super słownik policzony{oligo_counts_dict}")


gen.save_to('spectrum.txt', set(spectrum))
# print(set(spectrum))

#graph itd
graph = gen.create_graph(spectrum)
oligo_count = len(spectrum)
path = gen.create_paths(graph, start_oligo, oligo_count)


print(f'Przewidziana sciezki (cud jak sie nie wypierniczy)\n{path}')





