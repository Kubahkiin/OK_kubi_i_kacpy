import gen

#parametry
n = 28
k = 5
ne = 2
pe = 2
pop = 5
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
paths = []
path = []
paths.append(gen.create_paths(graph, start_oligo, path, oligo_count))

first_path = paths[0]

for x in range(1, pop):
    immature_path = first_path[:x]
    resume_vertex = gen.add_random_edge(graph, immature_path[x-1], immature_path) # wylosuj wybór krawędzi dla wierzchołka num i dodaj do ścieżki
    paths.append(gen.create_paths(graph, resume_vertex, immature_path, oligo_count))


print(f'Przewidziane sciezki\n')

for solution in paths:
    print(f'{solution}\n')
    print(f'Suma pokrycia: {gen.evaluate(graph, solution)}')





