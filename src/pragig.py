#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys
import networkx as nx


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
parser.add_argument('-o', '--output_file', type=str, help="If defined, output is saved in the given file")


arguments = parser.parse_args()

# Genome content and tree file are read
input = Input(arguments.G, arguments.T)
# Get all pairwise 'siblings' in the tree

pairwise_genomes = input.find_pairwise_leaves(input.tree[0])
all_leaves = input.find_all_leaves(input.tree[0])
max_length = int(input.tree[1])

all_genomes = {}

for leaf in all_leaves:
    all_genomes[leaf.name] = Genome(leaf.name, input.genomes[leaf.name])

gene_number = all_genomes.values()[0].adj_length()
calculate_probability.preprocess_transitions(2*max_length, gene_number)

potential_ancestors = {}


# main iteration
while pairwise_genomes:
    pair = pairwise_genomes.pop()
    names = pair[0]
    lca = input.tree[0].common_ancestor(names)

    print >> sys.stderr, names

    distances = {}
    for genome_name in all_genomes.keys():
        distances[genome_name] = input.tree[0].distance(lca, genome_name)

    first_genome = all_genomes[names[0]]
    second_genome = all_genomes[names[1]]


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

    if distances[names[0]] == 0.0:
       ancestor = first_genome
    elif distances[names[1]] == 0.0:
        ancestor = second_genome
    else:
        ancestor = []
        sampler = Genome_Sampler(inter_info.circular_breakpoint)

        #all_IGs = []
        highest_prob = None

        ##############################################
        ##############################################
        # TUESDAY 13 SEPTEMBER. STARTING ALL OVER :) #
        ##############################################
        ##############################################
        #i = 0
        #while i != arguments.repetition:
        for i in range(arguments.repetition):
            candidate = sampler.enumerate_vertices()
            expected_distances = {}
            probability = 0

            for identifier, genome in all_genomes.items():
                breakpoint_graph = IG(genome, candidate)
                breakpoint_graph.create_circular_graph()
                breakpoint_graph = breakpoint_graph.circular_breakpoint

                no_cycles = nx.number_connected_components(breakpoint_graph)
                distance = candidate.length() - no_cycles

                if distance < distances[identifier]:
                    break

                sorting_scen = calculate_probability.optimal_scenarios(breakpoint_graph)
                all_scen = calculate_probability.all_scenarios(genome.adj_length(), distances[identifier])
                probability += (sorting_scen - all_scen)

            # Apparently this is only called if the for-loop did not break
            # This seems to be very fancy!
            else:
                #i += 1
                if probability > highest_prob:
                    highest_prob = probability
                    ancestor = candidate


    ancestor.name = "%s%s" % (names[0], names[1])

    clade = input.tree[0].common_ancestor(names[0],names[1])
    clade.name = ancestor.name
    input.tree[0].collapse(names[0])
    input.tree[0].collapse(names[1])
    all_genomes.pop(names[0])
    all_genomes.pop(names[1])
    input.genomes.update({ancestor.name:ancestor.content})
    all_genomes[ancestor.name] = ancestor
    pairwise_genomes = input.find_pairwise_leaves(input.tree[0])

###############################################################################################



def write_output():
    if arguments.output_file:
        output = open(arguments.output_file, 'w')

    for result_name, result_content in input.genomes.items():
        result = Genome(result_name,result_content)
        chromosome = result.chr_number()

        if not arguments.output_file:
            print '>',result_name
        else:
            output.write('>%s\n' % (result_name))
        for i in range(chromosome):
            if not arguments.output_file:
                print '# chr%d' % (i+1)
            else:
                output.write('#chr%d\n' % (i + 1))
            for index, gene in enumerate(result_content):
                if not arguments.output_file:
                    print gene,
                else:
                    output.write("%s " % (gene))
                if gene == ')' or gene == '$':
                    result_content = result_content[index+1:]
                    if not arguments.output_file:
                        print
                    else:
                        output.write('\n')
                    break

write_output()