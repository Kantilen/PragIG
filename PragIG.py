#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import re
import sys
import networkx as nx


import adjacency_creation
import build_cBP
import find_intermediate_adjacencies
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
    first_genome = set(re.sub('[-)|]', '', first_genome))
    second_genome = set(re.sub('[-)|]', '', second_genome))

    # symmetrical difference. Just take elements that are unique in one
    # of the sets
    difference = first_genome ^ second_genome

    return (len(difference) == 0, difference)

# Commandline arguments
# TODO: This will be done via file reading soon...
parser = args.ArgumentParser(description="Enter two genomes in the input format of Unimog")
parser.add_argument('genomeA', metavar='genomeA', type=str, help="Sequence of the first genome")
parser.add_argument('genomeB', metavar='genomeB', type=str, help="Sequence of the second genome")

arguments = parser.parse_args()

same_content = validate_input(arguments.genomeA, arguments.genomeB)
if not same_content[0]:
    print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
    print >> sys.stderr, " ".join(["%s" % x for x in same_content[1]])
    sys.exit(1)

# Create adjacency sets of the two genomes
adjacency_setA = adjacency_creation.create_adjacency_set(arguments.genomeA)
adjacency_setB = adjacency_creation.create_adjacency_set(arguments.genomeB)

# Create the circular breakpoint graph of the two genomes
circular_breakpoint = build_cBP.connect_adjacencies(adjacency_setA, adjacency_setB)
intermediate_adj = find_intermediate_adjacencies.find_all_adjacencies(circular_breakpoint)
print intermediate_adj, len(intermediate_adj)