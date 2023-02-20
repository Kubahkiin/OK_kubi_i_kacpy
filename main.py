import gen

# parametry
n = 300  # długość sekwencji
k = 7  # długość podciągu
ne = 0  # błędy negatywne
pe = 0  # błędy pozytywne
pop = 50  # populacja
generations = 10  # liczba iteracji metaheurystyki
mutation_chance = 5  # szansa na losowe mutacje
top = 5  # ilość wybranych najlepszych rozwiązań do krzyżowania
lucky_chance = 5  # szansa na przejście do krzyżowania dla pozostałych
max_population = 100  # maksymalna liczba populacji na generację
instances = 10  # liczna instancji
#
oligo_count = n - (k - 1)

grouped_results = []

seq = gen.load_from("dna_300.txt")
numerek = 1
for steps in range(instances):
    print(f'INSTANCJA NR {numerek}\n')
    numerek += 1
    spectrum = gen.spectrum(seq, k, n)
    original_seq_spectrum = spectrum
    # dodanie błędów
    # błędy negatywne
    gen.negative_errors(spectrum, ne)
    # błędy pozytywne
    gen.positive_errors(spectrum, pe, k)

    start_oligo = original_seq_spectrum[0]
    # print(f'Oligonukleotyd startowy: {start_oligo}')

    # sortowanie w celu pomieszania sekwnecji
    spectrum.sort()

    oligo_counts_dict = gen.oligo_counter(spectrum)

    # print(f"To nasz super słownik policzony{oligo_counts_dict}")

    # graph itd
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

    # print(f'Przewidziane sciezki\n')

    solutions = []
    for solution in paths:
        solution_seq = gen.build_sequence(graph, solution)
        coverage_sum = gen.evaluate(graph, solution)
        mean_coverage = round(coverage_sum / len(solution), 2)

        s = gen.Solution(coverage_sum, mean_coverage, solution, solution_seq)
        similarity = gen.compare(seq, s.sequence)
        solutions.append(s)

    # gen.display_solutions(solutions, seq)
    # print("Ścieżki posortowane względem pokrycia")

    original_solution = gen.path_to_solution(graph, first_path)
    best_solutions = [
        (original_solution.mean_coverage, original_solution.coverage, gen.compare(seq, original_solution.sequence))]
    # best_solutions[7]
    for generation in range(generations):
        print(f'GENERACJA {generation + 1}\n')
        # Eliminacja zbyt krótkich ścieżek
        solutions = gen.kill_copies(solutions)
        killed = gen.kill(solutions, n)
        # print(f'Eliminacja za krótkich sekwencji: {killed}')
        # Sortowanie ścieżek po funkcji celu
        gen.ranking(solutions)
        # gen.display_solutions(solutions, seq)
        # Wybór ścieżek do crossover
        pick = gen.pick_for_crossover(solutions, top, lucky_chance)
        solutions_to_cross = pick[0]
        lucky = pick[1]
        # print(f'Wybrano {top} najlepszych + {lucky} miało szczęście i też zostaną skrzyżowane')
        # gen.display_solutions(solutions_to_cross, seq)
        # Crossover, jeśli nie można, to wymuszona mutacja
        baby_count = gen.make_babies(graph, solutions_to_cross, n, solutions, oligo_count, oligo_counts_dict)
        # print(f'W wyniku krzyżowania powstało {baby_count} nowych rozwiązań')
        # Mutacja wszystkich ścieżek z prawdopodobieństwem jakimś
        mutant_count = gen.random_mutations(solutions, mutation_chance, graph, oligo_count, oligo_counts_dict, n)
        # print(f'W wyniku losowych mutacji z prawdopodobieństwem {mutation_chance}% powstało {mutant_count} nowych rozwiązań')
        solutions = gen.kill_copies(solutions)
        # Ponowny ranking
        gen.ranking(solutions)
        # Eliminacja ostatnich pozycji w rankingu
        solutions = gen.natural_selection(solutions, max_population)
        # Wybór obecnej najlepszej ścieżki
        # print("Najlepsze rozwiązanie w tej generacji: ")
        # gen.display_solutions([solutions[0]], seq)
        best_solution = solutions[0]
        best_solutions.append(
            (best_solution.mean_coverage, best_solution.coverage, gen.compare(seq, best_solution.sequence)))

    # gen.display_solutions_light(best_solutions)
    grouped_results.append(best_solutions)

mean_from_instances = []
for generation in range(generations):
    pokrycie = 0
    avg_pokrycie = 0
    podobienstwo = 0
    for instance in range(instances):
        avg_pokrycie += grouped_results[instance][generation][0]
        pokrycie += grouped_results[instance][generation][1]
        podobienstwo += grouped_results[instance][generation][2]
    mean_from_instances.append((avg_pokrycie / instances, pokrycie / instances, podobienstwo / instances))
##############################################
# RAPORT
# print(grouped_results)
print(f'Parametry:\tk = {k}\tpe = {pe}\tne = {ne}\tmutation_chance = {mutation_chance}')
for mean_result in mean_from_instances:
    print(f'{mean_result[0]},{mean_result[1]},{mean_result[2]}')
