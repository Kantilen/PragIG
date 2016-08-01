#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys
from collections import Counter

from input_parser import Input
from genome_sampler import Genome_Sampler
from model import Genome
from ig_info import Intermediate_Genome as IG
import calculate_probability
#################################


__author__ = 'klamkiewicz'

# Commandline arguments
parser = args.ArgumentParser(description="Enter two genomes in the input format of Unimog")
parser.add_argument('G', metavar='GENOMES', type=str, help="Path to the file that contains the information of extant genomes")
parser.add_argument('T', metavar='TREE', type=str, help="Path to the file that contains the NEWICK tree")
parser.add_argument('-r', '--repetition', default=100, type=int, help="Number of sampled genomes for each ancestor, default: 100")


arguments = parser.parse_args()

# Genome content and tree file are read
input = Input(arguments.G, arguments.T)
# Get all pairwise 'siblings' in the tree

pairwise_genomes = input.find_pairwise_leaves(input.tree[0])
max_length = int(input.tree[1])

potential_ancestors = {}

#TODO: The 1000 is the number of genes in the genomes. This is hard-coded atm!
calculate_probability.preprocess_transitions(max_length, 1000)

# main iteration

while pairwise_genomes:
    pair = pairwise_genomes.pop()
    names = pair[0]
    distances = pair[1]
    first_content = input.genomes[names[0]]
    second_content = input.genomes[names[1]]

    # Create instance of Intermediate_Genome with the two current sibling-genomes.
    first_genome = Genome(names[0], first_content)
    second_genome = Genome(names[1], second_content)

    inter_info = IG(first_genome, second_genome)

    # check if content of genomes is identical.
    is_valid = inter_info.validate_input()

    # if not, quit the process
    if not is_valid[0]:
        print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
        print >> sys.stderr, " ".join(["%s" % x for x in is_valid[1]])
        sys.exit(1)

    # Create the circular breakpoint graph of the two genomes
    inter_info.create_circular_graph()
    highest_prob = -1000000000000000000000000000000000000000
    extant_adjacencies = set(first_genome.adjacency_set).union(set(second_genome.adjacency_set))

    ancestor = None
    print distances
    if distances[0] == 0.0:
        ancestor = first_genome
        print "HELLO1"
    elif distances[1] == 0.0:
        ancestor = second_genome
        print "HELLO2"
    else:
    # Sample genomes from the breakpoint graph
        for i in range(arguments.repetition):
        #print i
            pot_ancestor = Genome_Sampler(inter_info.circular_breakpoint).sampled_genomes

        #sampled_genomes = Genome_Sampler(inter_info.circular_breakpoint, 100).sampled_genomes


    #probabilities = {}

    #for pot_ancestor in sampled_genomes:
            ancestral_adjacencies = list(extant_adjacencies.union(pot_ancestor.adjacency_set))
            first_binary = first_genome.create_binary_vector(ancestral_adjacencies)
            second_binary = second_genome.create_binary_vector(ancestral_adjacencies)
            ancestor_binary = pot_ancestor.create_binary_vector(ancestral_adjacencies)

            prob = calculate_probability.calculate_probability(first_binary, second_binary, ancestor_binary, distances)

            if prob > highest_prob:
                highest_prob = prob
                ancestor = pot_ancestor
        #probabilities.update({prob:pot_ancestor})

    #ancestor = probabilities[max(probabilities.keys())]

    ancestor.name = "%s%s" % (names[0], names[1])

    clade = input.tree[0].common_ancestor(names[0],names[1])
    clade.name = ancestor.name
    input.tree[0].collapse(names[0])
    input.tree[0].collapse(names[1])

    input.genomes.update({ancestor.name:ancestor.content})
    pairwise_genomes = input.find_pairwise_leaves(input.tree[0])
    print pairwise_genomes
    print ancestor.name
    print ancestor.content