import gen

#parametry
n = 300 #długość sekwencji
k = 7   #długość podciągu
ne = 150 #błędy negatywne
pe = 200  #błędy pozytywne
pop = 50 #populacja
gen_iter = 50 #liczba iteracji metaheurystyki
mutation_chance = 5 #szansa na losowe mutacje
top = 5 # ilość wybranych najlepszych rozwiązań do krzyżowania
lucky_chance = 5 # szansa na przejście do krzyżowania dla pozostałych
max_population = 100 # maksymalna liczba populacji na generację
#
oligo_count = n - (k - 1)

seq = gen.load_from("dna_300.txt")
numerek = 1
for steps in range(1):
    print(f'INSTANCJA NR {numerek}\n'
          f'{seq}')
    numerek += 1
    spectrum = gen.spectrum(seq, k, n)
    original_seq_spectrum = spectrum
    # dodanie błędów
    # błędy negatywne
    gen.negative_errors(spectrum, ne)
    # błędy pozytywne
    gen.positive_errors(spectrum, pe, k)

    start_oligo = original_seq_spectrum[0]
    print(f'Oligonukleotyd startowy: {start_oligo}')

    # sortowanie w celu pomieszania sekwnecji
    spectrum.sort()

    oligo_counts_dict = gen.oligo_counter(spectrum)

    # print(f"To nasz super słownik policzony{oligo_counts_dict}")

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


    original_solution = gen.path_to_solution(graph, first_path)
    best_solutions = [original_solution]
    #best_solutions[7]
    for generation in range(gen_iter):
        print(f'######################################\n'
              f'########### GENERACJA {generation+1} ##############\n'
              '######################################')
        # Eliminacja zbyt krótkich ścieżek
        killed = gen.kill(solutions, n)
        print(f'Eliminacja za krótkich sekwencji: {killed}')
        # Sortowanie ścieżek po funkcji celu
        gen.ranking(solutions)
        #gen.display_solutions(solutions, seq)
        # Wybór ścieżek do crossover
        pick = gen.pick_for_crossover(solutions, top, lucky_chance)
        solutions_to_cross = pick[0]
        lucky = pick[1]
        print(f'Wybrano {top} najlepszych + {lucky} miało szczęście i też zostaną skrzyżowane')
        #gen.display_solutions(solutions_to_cross, seq)
        # Crossover, jeśli nie można, to wymuszona mutacja
        baby_count = gen.make_babies(graph, solutions_to_cross, n, solutions, oligo_count, oligo_counts_dict)
        print(f'W wyniku krzyżowania powstało {baby_count} nowych rozwiązań')
        # Mutacja wszystkich ścieżek z prawdopodobieństwem jakimś
        mutant_count = gen.random_mutations(solutions, mutation_chance, graph, oligo_count, oligo_counts_dict, n)
        print(f'W wyniku losowych mutacji z prawdopodobieństwem {mutation_chance}% powstało {mutant_count} nowych rozwiązań')
        # Ponowny ranking
        gen.ranking(solutions)
        # Eliminacja ostatnich pozycji w rankingu
        solutions = gen.natural_selection(solutions, max_population)
        # Wybór obecnej najlepszej ścieżki
        print("Najlepsze rozwiązanie w tej generacji: ")
        gen.display_solutions([solutions[0]], seq)
        best_solutions.append(solutions[0])


    gen.display_solutions(best_solutions, seq)