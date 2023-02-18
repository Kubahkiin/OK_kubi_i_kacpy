import gen

#parametry
n = 100 #długość sekwencji
k = 4   #długość podciągu
ne = 0 #błędy negatywne
pe = 0  #błędy pozytywne
pop = 10 #populacja
gen_iter = 1
mutation_chance = 0.05
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

# print(f"To nasz super słownik policzony{oligo_counts_dict}")


gen.save_to('spectrum.txt', set(spectrum))
# print(set(spectrum))

#graph itd
graph = gen.create_graph(spectrum)
paths = []
path = []
paths.append(gen.create_paths(graph, start_oligo, path, oligo_count, oligo_counts_dict))
# paths.append(gen.create_path_using_weights(graph, start_oligo, path, oligo_count, oligo_counts_dict))

first_path = paths[0]
for x in range(1, pop):
    random_start = gen.random_starting_node(graph)
    path = []
    paths.append(gen.create_paths(graph, random_start, path, oligo_count, oligo_counts_dict))

print(f'Przewidziane sciezki\n')

solutions = []
for solution in paths:
    solution_seq = gen.build_sequence(graph, solution)
    coverage_sum = gen.evaluate(graph, solution)
    mean_coverage = round(coverage_sum/len(solution), 2)

    s = gen.Solution(coverage_sum, mean_coverage, solution, solution_seq)
    similarity = gen.compare(seq, s.sequence)
    solutions.append(s)

#gen.display_solutions(solutions, seq)
#print("Ścieżki posortowane względem pokrycia")


for generation in range(gen_iter):
    print(f'######################################\n'
          f'########## GENERACJA {generation} #############\n'
          '######################################')

    gen.ranking(solutions)
    gen.display_solutions(solutions, seq)


# print("Test krzyżowania")
# print(paths[0])
#
# solution_seq = gen.build_sequence(graph, paths[0])
# print(f'{solution_seq}')
# coverage_sum = gen.evaluate(graph, paths[0])
# similarity = gen.compare(seq, solution_seq)
# print(f'Suma pokrycia: {coverage_sum}\t'
#      f'Średnie pokrycie: {round(coverage_sum/oligo_count, 2)}\t'
#      f'Podobieństwo do oryginalnej sekwencji: {similarity}\n')
#
# print(paths[1])
#
# solution_seq = gen.build_sequence(graph, paths[1])
# print(f'{solution_seq}')
# coverage_sum = gen.evaluate(graph, paths[1])
# similarity = gen.compare(seq, solution_seq)
# print(f'Suma pokrycia: {coverage_sum}\t'
#      f'Średnie pokrycie: {round(coverage_sum/oligo_count, 2)}\t'
#      f'Podobieństwo do oryginalnej sekwencji: {similarity}\n')
#
# print('Znalezione miejsce krzyżowania (LCS): ')
# lcs = gen.find_lcs(paths[0], paths[1])
# print(lcs)
# print("Ścieżki powstałe w wyniku krzyżowania: ")
# new_path = gen.crossover(paths[0], paths[1])
# print(new_path)
#
# solution_seq = gen.build_sequence(graph, new_path)
# print(f'{solution_seq}')
# coverage_sum = gen.evaluate(graph, new_path)
# similarity = gen.compare(seq, solution_seq)
# print(f'Suma pokrycia: {coverage_sum}\t'
#      f'Średnie pokrycie: {round(coverage_sum/oligo_count, 2)}\t'
#      f'Podobieństwo do oryginalnej sekwencji: {similarity}\n')
#






