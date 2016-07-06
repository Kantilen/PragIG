#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import re
import sys

from input_parser import Input
from adjacency_creation import Adjacency_Set
import build_cBP
#from find_intermediate_adjacencies import Inter_Adjacencies
import find_intermediate_adjacencies
import binary_vector
#################################


__author__ = 'klamkiewicz'

def validate_input(first_genome, second_genome):
    '''
    At the moment no duplications or indel events are allowed. Therefore the set of the genome content
    has to be identical. This is checked here.
    :param first_genome: Sequence of the first genome
    :param second_genome:  Sequence of the second genome
    :return: Boolean variable
    '''

    # remove signs and chromosomes
    first_genome = [re.sub('-','',x) for x in first_genome]
    second_genome = [re.sub('-', '', x) for x in second_genome]

    # symmetrical difference. Just take elements that are unique in one
    # of the sets
    first_genome = set(first_genome)
    second_genome = set(second_genome)
    difference = first_genome ^ second_genome

    return (len(difference) == 0, difference)

# Commandline arguments
parser = args.ArgumentParser(description="Enter two genomes in the input format of Unimog")
parser.add_argument('G', metavar='GENOMES', type=str, help="Path to the file that contains the information of extant genomes")
parser.add_argument('T', metavar='TREE', type=str, help="Path to the file that contains the NEWICK tree")

arguments = parser.parse_args()

input = Input(arguments.G, arguments.T)
pairwise_genomes = []

for clade in input.tree.find_clades():
    if clade.count_terminals() == 2:
        leaves = clade.find_clades()
        leaves = [x.name for x in leaves if x.name]
        pairwise_genomes.append(leaves)

for pair in pairwise_genomes:
    first = input.genomes[pair[0]]
    second = input.genomes[pair[1]]

    same_content = validate_input(first, second)
    if not same_content[0]:
        print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
        print >> sys.stderr, " ".join(["%s" % x for x in same_content[1]])
        sys.exit(1)

    # Create adjacency sets of the two genomes
    adjacency_setA = Adjacency_Set(first)
    adjacency_setB = Adjacency_Set(second)

    # Create the circular breakpoint graph of the two genomes
    #circular_breakpoint = Inter_Adjacencies(adjacency_setA.adjacencies, adjacency_setB.adjacencies)
    circular_breakpoint = build_cBP.connect_adjacencies(adjacency_setA.adjacencies, adjacency_setB.adjacencies)
    intermediate_adj = find_intermediate_adjacencies.find_all_adjacencies(circular_breakpoint)

    print len(intermediate_adj)
    #bin_vector_A = binary_vector.create_vector_for_genome(adjacency_setA.adjacencies, intermediate_adj)
    #bin_vector_B = binary_vector.create_vector_for_genome(adjacency_setB.adjacencies, intermediate_adj)