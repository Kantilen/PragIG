#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import re
import sys

from input_parser import Input
from Intermediate_Genome import Intermediate_Genome

#################################


__author__ = 'klamkiewicz'

def validate_input(first_genome, second_genome):
    '''
    At the moment no duplications or indel events are allowed. Therefore the set of the genome content
    has to be identical. This is checked here.
    :param first_genome: Sequence of the first_content genome
    :param second_genome:  Sequence of the second_content genome
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

# Genome content and tree file are read
input = Input(arguments.G, arguments.T)
# Get all pairwise 'siblings' in the tree
pairwise_genomes = input.find_pairwise_leaves(input.tree)

potential_ancestors = {}

# main iteration
#TODO: Do a while-loop, pop the first element and update the list with the next (resolved) level in the tree.
for pair in pairwise_genomes:
    first_content = input.genomes[pair[0]]
    second_content = input.genomes[pair[1]]

    # content validation; no singletons allowed
    same_content = validate_input(first_content, second_content)
    if not same_content[0]:
        print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
        print >> sys.stderr, " ".join(["%s" % x for x in same_content[1]])
        sys.exit(1)

    # Create adjacency sets of the two genomes
    inter_info = Intermediate_Genome(first_content, second_content)
    inter_info.create_adjacency_sets()
    # Create the circular breakpoint graph of the two genomes
    inter_info.connect_adjacencies()
    # Find all intermediate genomes
    inter_info.get_all_inter_adj()


    potential_ancestors.update({(pair[0],pair[1]) : (inter_info.circular_breakpoint, inter_info.inter_adj)})

for key,value in potential_ancestors.items():
    print key, len(value[1])