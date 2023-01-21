import random
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path

nucs = ['A', 'C', 'G', 'T']

out_dir = Path("DATA")
out_dir.mkdir(parents=True, exist_ok=True)


def save_to(filename, to_save):
    f1 = open(out_dir.joinpath(filename), 'w')
    for stuff in to_save:
        f1.writelines(stuff + '\n')
    f1.close()


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

    for i in range(error_count):
        removed_oligo = random.choice(spectrum)
        removed_oligos.append(removed_oligo)
        spectrum.remove(removed_oligo)

    # print(f'Negative: {removed_oligos}')
    save_to('negative_errors.txt', removed_oligos)

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

    # print(f'Positive: {added_oligos}')
    save_to('positive_errors.txt', added_oligos)

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

def create_paths(G, start, path, desired_path_length):
    curr = start
    path.append(curr)
    # display list of successors of the starting node
    # jakas petla while nie jest długość taka jak n - (k - 1)
    while len(path) != desired_path_length:
        # print(f'Lista następników dla wierzchołka {curr}:')
        # Wyświetlanie wag krawędzi prowadzących do następników wierzchołka curr
        chosen_one = 0
        chosen_weight = 0
        for successor in G.successors(curr):
            # print(f"Waga krawędzi prowadzącej z {curr} do następnika {successor}: {G.edges[curr, successor, 0]['weight']}")
            # funkcja ktora wybiera nastepny wierzcholek do sciechy
            successor_weight = G.edges[curr, successor, 0]['weight']

            if chosen_one == 0 or successor_weight > chosen_weight:
                chosen_one = successor
                chosen_weight = successor_weight
        curr = chosen_one
        path.append(curr)
    return path

def add_random_edge(G, vertex, path):
    successors = []

    for successor in G.successors(vertex):
        successors.append(successor)

    random_vertex = random.choice(successors)
    return random_vertex

def evaluate(G, path):
    coverage_sum = 0
    for i in range(len(path)-1):
        A = path[i]
        B = path[i + 1]
        coverage_sum += G.edges[A, B, 0]['weight']
    return coverage_sum

#TODO zaimplementować pamięć występowania danych oligo
# (ile razy można wrócić do wierzchołka w grafie) korzystając z naszego super słownika
