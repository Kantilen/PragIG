#!/usr/bin/python
import sys
import networkx as nx
import os
import argparse as args

import input_parser
from ig_info import Intermediate_Genome as IG
import calculate_probability
from model import Genome


__author__ = 'klamkiewicz'

def get_dcj_distance_from_BP(genome, candidate):
    breakpoint_graph = IG(genome, candidate)
    breakpoint_graph.create_circular_graph()
    breakpoint_graph = breakpoint_graph.circular_breakpoint
    no_cycles = nx.number_connected_components(breakpoint_graph)
    return (candidate.length() - no_cycles, breakpoint_graph)

parser = args.ArgumentParser(description="Calculates the probability of the real ancestor being the ancestor.")
parser.add_argument('A', metavar='ANCESTOR', type=str, help="Path to the genome file containing the ancestors.")
parser.add_argument('L', metavar='LEAVES', type=str, help="Path to genome file containing the leaves.")
parser.add_argument('T', metavar='TREE', type=str, help="Path to the file containg the tree.")
parser.add_argument('AL', metavar='ALPHA', type=float, help="Alpha value that was used for the analysis.")
parser.add_argument('O', metavar='OUTPUT', type=str, help="Path to the folder where the results are saved.")

arguments = parser.parse_args()

input = input_parser.Input(arguments.A, arguments.T, False)
leaves = input.read_genomes(arguments.L)
alpha = arguments.AL

all_genomes = {}
probabilities = {}

for leaf, content in leaves.items():
    all_genomes[leaf] = Genome(leaf, content)

for ancestor, content in input.genomes.items():
    all_genomes[ancestor] = Genome(ancestor, content)

pairwise_genomes = input.find_pairwise_leaves(input.tree[0])

while pairwise_genomes:
    pair = pairwise_genomes.pop()
    names = pair[0]
    lca = input.tree[0].common_ancestor(names)

    ancestor = all_genomes[lca.name]
    children = [all_genomes[x] for x in names]

    probability = 0

    for child in children:
        tree_distance = input.tree[0].distance(child.name, ancestor.name)
        #inter_info = ig_info.Intermediate_Genome(child, ancestor)
        distance, breakpoint_graph = get_dcj_distance_from_BP(child, ancestor)
        expected_distance = child.distance_to_genome(ancestor)
        lower_bound = distance * alpha
        upper_bound = (2 - alpha) * expected_distance

        sorting_scen = calculate_probability.optimal_scenarios(breakpoint_graph)
        all_scen = calculate_probability.all_scenarios(child.adj_length(), distance)
        probability = calculate_probability.calculate_prob_ancestor(probability, sorting_scen, all_scen,
                                                                    alpha,
                                                                    distance, tree_distance,
                                                                    expected_distance,
                                                                    lower_bound, upper_bound)



    probabilities[ancestor.name] = probability
    input.tree[0].collapse(names[0])
    input.tree[0].collapse(names[1])
    all_genomes.pop(names[0])
    all_genomes.pop(names[1])
    pairwise_genomes = input.find_pairwise_leaves(input.tree[0])

folder = "%s/ancestor_probs_A%d" % (arguments.O, alpha)

if not os.path.exists(arguments.O):
    os.mkdir(arguments.O)

output = open(folder,'w')

for ancestor, prob in probabilities.items():
    output.write('%s;%d\n' % (ancestor, prob))