#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys
import networkx as nx
import math
import bigfloat

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

        print distances
        print names
        print all_genomes
        print input.tree

        for single_component in nx.connected_component_subgraphs(inter_info.circular_breakpoint):
            if len(single_component.nodes()) == 2:
                continue
            else:
                size = len(single_component.nodes())
                single_distance = (size/2) - 1

                #print [(x,single_component.get_edge_data(*x).values()[0]['color']) for x in single_component.edges()]

                all_A_adjacencies = [x for x in single_component.edges() if
                                     single_component.get_edge_data(*x).values()[0]['color'] == 'A']

                all_B_adjacencies = [x for x in single_component.edges() if
                                     single_component.get_edge_data(*x).values()[0]['color'] == 'B']

                print
                print all_A_adjacencies
                print all_B_adjacencies

                enumerated_vertices = {}
                first_adj = single_component.edges()[0]
                long_path = [x for x in nx.all_simple_paths(single_component, first_adj[0], first_adj[1])][-1]
                #assign value for each vertex from 0 to n
                for index, vertex in enumerate(long_path):
                    enumerated_vertices[index] = vertex

                for i in range(arguments.repetition):
                    inter_comp = sampler.create_adjacency_from_cycle(enumerated_vertices.values())

                    

        break

   #     for i in range(arguments.repetition):

   #         sampled_genome = sampler.enumerate_vertices()
   #         expected_DCJ_first = sampled_genome.distance_to_genome(first_genome)
   #         expected_DCJ_second = sampled_genome.distance_to_genome(second_genome)
            #print
            #print i
            #print expected_DCJ_first, distances[names[0]]
            #print expected_DCJ_second, distances[names[1]]
            #print expected_DCJ_first + expected_DCJ_second, distances[names[0]]+distances[names[1]]

            #if expected_DCJ_first < distances[names[0]] or expected_DCJ_second < distances[names[1]]:
            #    print "DISCARDING :("
            #    #continue

            #else:
            #    print "CALCULATING! :)"


   #         all_first = ((first_genome.adj_length() * (first_genome.adj_length()+1)) / 2) ** int(expected_DCJ_first)

            #first_IG = IG(first_genome, sampled_genome)
            #first_IG.create_circular_graph()

            #second_IG = IG(second_genome, sampled_genome)
            #second_IG.create_circular_graph()

            #first_lengths = {}
            #for index,component in enumerate(nx.connected_component_subgraphs(first_IG.circular_breakpoint)):
            #    if len(component.nodes()) == 2:
            #        continue
            #    first_lengths[index] = (len(component.nodes()) / 2) - 1

            #upper = 0
            #lower = 1
            #prod = 1
            #for comp, dist in first_lengths.items():
            #    upper += dist
            #    lower *= math.factorial(dist)
            #    prod *= (dist + 1) ** (dist - 1)

            #first_scenarios = (math.factorial(upper) / lower) * prod
            #prob_first = math.log10(first_scenarios) - math.log10(all_first)

            #all_second = ((second_genome.adj_length() * (second_genome.adj_length() + 1)) / 2) ** int(expected_DCJ_second)


            #second_distances = {}
            #for index, component in enumerate(nx.connected_component_subgraphs(second_IG.circular_breakpoint)):
             #   if len(component.nodes()) == 2:
             #       continue
             #   second_distances[index] = (len(component.nodes()) / 2) - 1

            #upper = 0
            #lower = 1
            #prod = 1
            #for comp, dist in second_distances.items():
             #   upper += dist
             #   lower *= math.factorial(dist)
             #   prod *= (dist + 1) ** (dist - 1)

            #second_scenarios = (math.factorial(upper) / lower) * prod
            #prob_second = math.log10(second_scenarios) - math.log10(all_second)

            #prob = prob_first + prob_second
            #if prob > highest_prob:
             #   highest_prob = prob
             #   ancestor = sampled_genome

        #print ancestor.content
        #sys.exit(0)

        #for component in nx.connected_component_subgraphs(inter_info.circular_breakpoint):
        #    if len(component.nodes()) == 2:
        #        ancestor.append(Adjacency(component.nodes()[0], component.nodes()[1]))
        #        continue

        #    adj_set_first_genome = []
        #    adj_set_second_genome = []
        #    for edge in component.edges():
        #        colors = component.get_edge_data(*edge).values()
        #        adj_set_first_genome.extend([Adjacency(edge[0],edge[1]) for color in colors if color['color'] == 'A' ])
        #        adj_set_second_genome.extend([Adjacency(edge[0], edge[1]) for color in colors if color['color'] == 'B'])
        #    extant_adjacencies = set(adj_set_first_genome).union(set(adj_set_second_genome))

        #    highest_candidate = None
        #    highest_prob = None

        #    cycle = []

        #    enumerated_vertices = {}

         #   first_adj = component.edges()[0]
         #   long_path = [x for x in nx.all_simple_paths(component, first_adj[0], first_adj[1])][-1]
            # assign value for each vertex from 0 to n
         #   for index, vertex in enumerate(long_path):
         #       enumerated_vertices[index] = vertex

            #if not len(component.nodes()) <= 14:
                #print "SAMPLE"
         #   all_IGs = []
         #   for i in range(arguments.repetition):
         #       all_IGs.append(Genome_Sampler.create_adjacency_from_cycle(enumerated_vertices.values()))
            #else:
            #    print "ALL"
            #    print enumerated_vertices.values()
            #    all_IGs = Genome_Sampler.get_all(enumerated_vertices.values(),"INIT")
                #print "ENDRESULT:", all_IGs

         #   for pot_ancestor in all_IGs:
                #print pot_ancestor
                    #cycle.append(Genome_Sampler.create_adjacency_from_cycle(enumerated_vertices.values()))
                    #pot_ancestor = [x for y in cycle for x in y]
                    #pot_ancestor = Genome_Sampler(component).intermediate_cycle
          #      all_adjacencies = list(extant_adjacencies.union(pot_ancestor))

           #     binaries = {}
           #     for genome in all_genomes.items():
           #         binaries[genome[0]] = genome[1].create_binary_vector(all_adjacencies, inter_info.circular_breakpoint)


           #     ancestor_binary = [1 if adj in pot_ancestor else 0 for adj in all_adjacencies]

           #     prob = calculate_probability.calculate_probability(binaries, ancestor_binary, distances)
           #     if prob > highest_prob:
           #         highest_prob = prob
           #         highest_candidate = pot_ancestor

           # ancestor.extend(highest_candidate)

#

        #ancestor = Genome.genome_from_adjacencies("", ancestor)

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

#write_output()