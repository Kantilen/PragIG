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
from model import Adjacency
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

    distances = {}
    for genome_name in all_genomes.keys():
        distances[genome_name] = input.tree[0].distance(lca, genome_name)


    #distances = pair[1]
    #first_content = input.genomes[names[0]]
    #second_content = input.genomes[names[1]]

    # Create instance of Intermediate_Genome with the two current sibling-genomes.
    #first_genome = Genome(names[0], first_content)
    #second_genome = Genome(names[1], second_content)

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

        for component in nx.connected_component_subgraphs(inter_info.circular_breakpoint):
            if len(component.nodes()) == 2:
                ancestor.append(Adjacency(component.nodes()[0], component.nodes()[1]))
                continue

            adj_set_first_genome = []
            adj_set_second_genome = []
            for edge in component.edges():
                colors = component.get_edge_data(*edge).values()
                adj_set_first_genome.extend([Adjacency(edge[0],edge[1]) for color in colors if color['color'] == 'A' ])
                adj_set_second_genome.extend([Adjacency(edge[0], edge[1]) for color in colors if color['color'] == 'B'])
            extant_adjacencies = set(adj_set_first_genome).union(set(adj_set_second_genome))

            highest_candidate = None
            highest_prob = None

            if not len(component.nodes() <= 14):

                for i in range(arguments.repetition):

                    pot_ancestor = Genome_Sampler(component).intermediate_cycle
                    all_adjacencies = list(extant_adjacencies.union(pot_ancestor))

                    binaries = {}
                    for genome in all_genomes.items():
                        binaries[genome[0]] = genome[1].create_binary_vector(all_adjacencies, inter_info.circular_breakpoint)

                #print all_adjacencies
                #print pot_ancestor

                    ancestor_binary = [1 if adj in pot_ancestor else 0 for adj in all_adjacencies]
                #print ancestor_binary
                #print "\n"
                    prob = calculate_probability.calculate_probability(binaries, ancestor_binary, distances)

                    if prob > highest_prob:
                        highest_prob = prob
                        highest_candidate = pot_ancestor

                ancestor.extend(highest_candidate)
            else:
                trolo = Genome_Sampler.get_all(component)
                print trolo
                sys.exit(0)
        ancestor = Genome.genome_from_adjacencies("", ancestor)
            #print ancestor

            #highest_prob = None
            #extant_adjacencies = set(first_genome.adjacency_set).union(set(second_genome.adjacency_set))

            #ancestor = None


            #else:
            # Sample genomes from the breakpoint graph
            #    for i in range(arguments.repetition):
            #        pot_ancestor = Genome_Sampler(inter_info.circular_breakpoint).sampled_genomes

            #for pot_ancestor in sampled_genomes:
            #ancestral_adjacencies = list(extant_adjacencies.union(pot_ancestor.adjacency_set))

            #binaries = {}
            #for genome in all_genomes.items():
            #    binaries[genome[0]] = genome[1].create_binary_vector(ancestral_adjacencies, inter_info.circular_breakpoint)

            #ancestor_binary = pot_ancestor.create_binary_vector(ancestral_adjacencies, inter_info.circular_breakpoint)

            #prob = calculate_probability.calculate_probability(binaries, ancestor_binary, distances)

            #if prob > highest_prob:
            #    highest_prob = prob
            #    ancestor = pot_ancestor





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
