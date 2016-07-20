#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys

from input_parser import Input
from genome_sampler import Genome_Sampler
from model import Genome
from ig_info import Intermediate_Genome as IG
#################################


__author__ = 'klamkiewicz'

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
    print pair
    first_content = input.genomes[pair[0]]
    second_content = input.genomes[pair[1]]


    # Create instance of Intermediate_Genome with the two current sibling-genomes.
    first_genome = Genome(pair[0], first_content)
    second_genome = Genome(pair[1], second_content)

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
    # Find all intermediate genomes
    #inter_info.get_all_inter_adj()
    # Create binary vector
    #binary_of_first = first_genome.create_binary_vector(inter_info.inter_adj)
    #binary_of_second = second_genome.create_binary_vector(inter_info.inter_adj)
    #binaries = {pair[0]:binary_of_first, pair[1]:binary_of_second}
    # Update everything

    sampled_genomes = Genome_Sampler(inter_info.circular_breakpoint, 10).sampled_genomes

    #sampled_genomes = Genome_Sampler(inter_info.circular_breakpoint, 10, len(first_genome.adjacency_set)).sampled_genomes

    #Genome.genome_from_adjacencies("s1", sampled_genomes[0])

    potential_ancestors.update({(pair[0], pair[1]): (inter_info.circular_breakpoint, sampled_genomes)})

for key,value in potential_ancestors.items():
    print key, [len(x.adjacency_set) for x in value[1]]