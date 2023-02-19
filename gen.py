import random
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path
import copy

nucs = ['A', 'C', 'G', 'T']

out_dir = Path("DATA")
out_dir.mkdir(parents=True, exist_ok=True)


def save_to(filename, to_save):
    f1 = open(out_dir.joinpath(filename), 'w')
    for stuff in to_save:
        f1.writelines(stuff + '\n')
    f1.close()

def load_from(filename):
    f1 = open(out_dir.joinpath(filename), 'r')
    line = f1.read()
    f1.close()
    return line


def dna(n):
    return ''.join(random.choices(nucs, k=n))


def spectrum(seq, k, n):
    spectrum = []
    for i in range(n - (k - 1)):
        spectrum.append(seq[i:i + k])

    return spectrum

def oligo_counter(spectrum):
    counted_dict = {oligo:spectrum.count(oligo) for oligo in spectrum}

    return counted_dict


def negative_errors(spectrum, error_count):
    removed_oligos = []  # storing oligonucleotides removed from spectrum
    spectrum = spectrum[1:] # pierwszy jest zawsze znany, więc nie możemy go usunąć nawet błędem

    for i in range(error_count):
        removed_oligo = random.choice(spectrum)
        removed_oligos.append(removed_oligo)
        spectrum.remove(removed_oligo)

    return spectrum


def positive_errors(spectrum, error_count, k):
    added_oligos = []  # storing oligonucleotides added to spectrum

    for i in range(error_count):
        """ := is actually a valid operator that allows for assignment of variables within expressions"""
        while random_k_mer := random.choices(nucs,
                                             k=k):
            if not spectrum.count(random_k_mer):
                break
        added_oligos.append(''.join(random_k_mer))
        spectrum.append(''.join(random_k_mer))

    #print(f'Positive: {added_oligos}')

    return spectrum


def get_coverage(str_1, str_2):
    oligo_1 = list(str_1)
    oligo_2 = list(str_2)

    k = len(oligo_1)
    shift = 0

    for i in range(1, k):
        if oligo_1[i] == oligo_2[shift]:
            shift += 1
        else:
            shift = 0

    coverage = shift
    return coverage


def create_graph(spectrum):
    G = nx.MultiDiGraph()
    G.add_nodes_from(spectrum)

    for node in G.nodes():
        for other_node in G.nodes():
            if node == other_node:  # if this is the same then dont
                continue
            coverage = get_coverage(node, other_node)
            if coverage > 0:
                G.add_edge(node, other_node, weight = coverage)
            #print(f'Porównanie {node} z {other_node} = pokrywa sie na poziomie {coverage}')

    nx.write_gexf(G, 'graph.xml')

    nx.draw(G,with_labels=True)
    #plt.show()

    return G

def create_paths(G, start, path, max_length, oligo_counts_dict_original):
    oligo_counts_dict = copy.deepcopy(oligo_counts_dict_original)
    forbidden = []
    curr = start
    path.append(curr)
    oligo_counts_dict[curr] -= 1
    # display list of successors of the starting node
    # jakas petla while nie jest długość taka jak n - (k - 1)
    #print(path)
    #print(len(path))
    for step in range(max_length):
        #print(path)
        #print(len(path))
        chosen_one = 0
        chosen_weight = 0
        for successor in G.successors(curr):
            # print(f"Waga krawędzi prowadzącej z {curr} do następnika {successor}: {G.edges[curr, successor, 0]['weight']}")
            # funkcja ktora wybiera nastepny wierzcholek do sciechy
            successor_weight = G.edges[curr, successor, 0]['weight']

            if oligo_counts_dict[successor] == 0:
                continue

            if successor in forbidden:
                continue

            if chosen_one == 0 or successor_weight > chosen_weight:
                chosen_one = successor
                chosen_weight = successor_weight

        curr = chosen_one

        if curr == 0:
            return path
            # forbidden_one = path.pop()
            # forbidden.append(forbidden_one)
            # oligo_counts_dict[forbidden_one] += 1
            # curr = path[-1]
            # continue

        path.append(curr)
        oligo_counts_dict[curr] -= 1
    return path

def random_starting_node(G):
    nodes = list(G.nodes())
    return random.choice(nodes)


def add_random_edge(G, vertex, oligo_count_immature):
    successors = []

    for successor in G.successors(vertex):
        if oligo_count_immature[successor] == 0:
            continue
        successors.append(successor)

    random_vertex = random.choice(successors)
    return random_vertex

def mutate_path(G, original_path, max_length, oligo_counts_dict_original):
    oligo_count_immature = copy.deepcopy(oligo_counts_dict_original)
    mutation_position = random.randint(1, len(original_path) - 1)
    immature_path = original_path[:mutation_position]
    for oligo in immature_path:
        oligo_count_immature[oligo] -= 1

    new_max_length = max_length - len(immature_path)
    resume_vertex = add_random_edge(G, immature_path[mutation_position - 1], oligo_count_immature)  # wylosuj następnika dla wierzchołka w którym ma dojść do mutacji
    return create_paths(G, resume_vertex, immature_path, new_max_length, oligo_count_immature)

def build_sequence(G, path):
    sequence = path[0]
    for i in range(1, len(path)):
        coverage = G.edges[path[i-1], path[i], 0]['weight'] # ile odciąć
        sequence += path[i][coverage:]
    return sequence

def evaluate(G, path):
    coverage_sum = 0
    for i in range(len(path)-1):
        A = path[i]
        B = path[i + 1]
        coverage_sum += G.edges[A, B, 0]['weight']
    return coverage_sum

def compare(seq1, seq2):
    score = 0
    loop_length = min(len(seq1), len(seq2))
    for pos in range(0, loop_length):
        if seq1[pos] == seq2[pos]:
            score += 1

    return score

class Solution:
    def __init__(self, coverage, mean_coverage, oligo_list, sequence):
        self.coverage = coverage
        self.mean_coverage = mean_coverage
        self.oligo_list = oligo_list
        self.sequence = sequence

def display_solutions(solutions, original_seq):
    for s in solutions:
        similarity = compare(original_seq, s.sequence)
        print(f'{s.oligo_list}')
        print(f'{s.sequence}')
        print(f'Suma pokrycia: {s.coverage}\t'
              f'Średnie pokrycie: {s.mean_coverage}\t'
              f'Podobieństwo do oryginalnej sekwencji: {similarity}\n')

def display_solutions_light(solutions, original_seq):
    for s in solutions:
        similarity = compare(original_seq, s.sequence)
        print(f'Suma pokrycia: {s.coverage}\t'
              f'Średnie pokrycie: {s.mean_coverage}\t'
              f'Podobieństwo do oryginalnej sekwencji: {similarity}\n')

def ourKey(s):
  return (s.mean_coverage, s.coverage)

def ranking(solutions):
    solutions.sort(key=ourKey, reverse=True)

def kill(solutions, n):
    to_kill = []
    killed = 0
    for solution in solutions:
        if len(solution.sequence) < n or len(solution.sequence) > 1.5*n:
            to_kill.append(solution)
            killed += 1

    for dead in to_kill:
        solutions.remove(dead)

    return len(to_kill)

def pick_for_crossover(solutions, top, lucky_chance):
    chosen = []
    chosen.extend(solutions[:top])

    for lucky_solution in solutions[top:]:
        if random.randint(0, 100) < lucky_chance:
            chosen.append(lucky_solution)
            lucky_chance -= 1

    return (chosen, len(chosen))

def make_babies(G, solutions_to_cross, n, solutions, oligo_count, oligo_counts_dict_original):
    babies = []
    for mom in solutions_to_cross:
        for dad in solutions_to_cross:
            if mom == dad:
                continue
            baby_path = crossover(mom.oligo_list, dad.oligo_list, G, oligo_count, oligo_counts_dict_original)
            baby = path_to_solution(G, baby_path)
            babies.append(baby)

            baby_path = crossover(dad.oligo_list, mom.oligo_list, G, oligo_count, oligo_counts_dict_original)
            baby = path_to_solution(G, baby_path)
            babies.append(baby)
    kill(babies, n)
    solutions.extend(babies)

    return len(babies)

def path_to_solution(G, path):
    solution_seq = build_sequence(G, path)
    coverage_sum = evaluate(G, path)
    mean_coverage = round(coverage_sum / len(path), 2)

    return Solution(coverage_sum, mean_coverage, path, solution_seq)

def random_mutations(solutions, mutation_chance, G, oligo_count, oligo_counts_dict_original, n):
    mutants = []
    for solution in solutions:
        if random.randint(0, 100) < mutation_chance:
            mutant_path = mutate_path(G, solution.oligo_list, oligo_count, oligo_counts_dict_original)
            mutant = path_to_solution(G, mutant_path)
            mutants.append(mutant)
    kill(mutants, n)
    solutions.extend(mutants)

    return len(mutants)

def natural_selection(solutions, max_population):
    new_solutions = solutions
    if(len(solutions) > max_population):
        new_solutions = solutions[:max_population]

    return new_solutions


##############################
#### Wygenerowane funkcje ####
##############################

def find_lcs(seq1, seq2):
    """
    Znajduje najdłuższy wspólny sufiks (LCS) dwóch sekwencji oligonukleotydów.
    seq1 i seq2 to listy oligonukleotydów.
    Funkcja zwraca LCS w postaci listy oligonukleotydów.
    """
    m = len(seq1)
    n = len(seq2)
    lcs_length = 0
    end = 0

    # inicjalizacja tablicy LCS
    lcs = [[0 for j in range(n+1)] for i in range(m+1)]

    # wypełnienie tablicy LCS
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                lcs[i][j] = lcs[i-1][j-1] + 1
                if lcs[i][j] > lcs_length:
                    lcs_length = lcs[i][j]
                    end = i
            else:
                lcs[i][j] = 0

    # wyciągnięcie LCS
    start = end - lcs_length
    lcs_seq = seq1[start:end]

    return lcs_seq

def crossover(seq1, seq2, G, oligo_count, oligo_counts_dict_original):
    """
    Realizuje krzyżowanie (crossover) dwóch sekwencji oligonukleotydów.
    seq1 i seq2 to listy oligonukleotydów.
    Funkcja zwraca nową sekwencję oligonukleotydów powstałą przez krzyżowanie.
    """
    lcs = find_lcs(seq1, seq2)

    if not lcs:
        #print("Nie można zrobić crossover, wymuszona mutacja")
        return mutate_path(G, seq1, oligo_count, oligo_counts_dict_original)

    idx = seq1.index(lcs[0])
    new_seq = seq1[:idx] + seq2[seq2.index(lcs[0]):]

    return new_seq