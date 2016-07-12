#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys

from input_parser import Input
from intermediate_genome import Intermediate_Genome
from genome_sampler import Genome_Sampler
from model import Genome
from intermediate_genome_2 import Intermediate_Genome as IG
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
    first_content = input.genomes[pair[0]]
    second_content = input.genomes[pair[1]]


    # Create instance of Intermediate_Genome with the two current sibling-genomes.
    inter_info = Intermediate_Genome(pair[0], pair[1], first_content, second_content)

    # check if content of genomes is identical.
    is_valid = inter_info.validate_input()

    # if not, quit the process
    if not is_valid[0]:
        print >> sys.stderr, "Content of the genomes is not equal. Check your input!"
        print >> sys.stderr, " ".join(["%s" % x for x in is_valid[1]])
        sys.exit(1)

    test_genome = Genome(pair[0],first_content)
    other_genome = Genome(pair[1],second_content)

    test_inter = IG(test_genome, other_genome)

    #print len(test_genome.create_adjacency_set()), test_genome.create_adjacency_set()[0].first_ex
    #print len(other_genome.create_adjacency_set()), other_genome.create_adjacency_set()[0].first_ex

    # Adjacency sets are created
    inter_info.create_adjacency_sets()
    #print (inter_info.adjacencies[pair[0]][0])
    #print(inter_info.adjacencies[pair[1]][0])

    # Create the circular breakpoint graph of the two genomes
    test_inter.create_circular_graph()
    # Find all intermediate genomes
    test_inter.get_all_inter_adj()
    # Create binary vector
    test_inter.create_binary_vector()
    # Update everything
    potential_ancestors.update({(pair[0],pair[1]) : (test_inter.circular_breakpoint, test_inter.inter_adj, test_inter.binaries)})

    sample_genomes = Genome_Sampler(potential_ancestors[(pair[0],pair[1])], 1000)

for key,value in potential_ancestors.items():
    print key, value[0].edges(), len(value[1]), value[2]