#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys
import math

import networkx as nx

from input_parser import Input
from genome_sampler import Genome_Sampler
from model import Genome
from ig_info import Intermediate_Genome as IG
import calculate_probability
from Bio import Phylo
import copy
#################################

__author__ = 'klamkiewicz'

def find_ancestral_weights(tree, extant_genomes):
    ''' TEST FUNCTION '''
    #print >> sys.stderr, extant_genomes
    #print >> sys.stderr, tree
    anc_weights = {}
    all_internal_nodes = tree.get_nonterminals()

    for ancestor in all_internal_nodes:

        #print >> sys.stderr, ancestor.name
        label = ancestor.name
        t = copy.deepcopy(tree)
        #print tree
        #Phylo.draw_ascii(tree)

        root = tree.get_nonterminals(order="level")[0]

        if ancestor != root:
            current = tree.find_clades(ancestor, order="level").next()
            children = list(current.find_clades(order="postorder"))
            children = [x.name for x in children if x != ancestor]
            for child in children:
                t.collapse(child)
            t.root_with_outgroup({'name' : current.name})

        weights = {}
        #for node in t.find_clades(order="postorder"):
        for idx, node in enumerate(t.find_clades(order="postorder")):
            #print idx, node, ancestor
            if node.is_terminal() and node.name != ancestor.name:
                weights[node.name] = {adj : 1 for adj in extant_genomes[node.name].adjacency_set}
                continue

            d = 0
            all_adj = set()
            children = list(node.find_clades())
            children = [x for x in children if x != node]



            for child in children:
                all_adj.update(weights[child.name].iterkeys())
                d += child.branch_length
            #print d, len(all_adj)
            if d == 0:
                d = 0.1
            node_adj_weights = {}
            for adj in all_adj:
                children_w = [weights[child.name][adj] * (d - child.branch_length) if adj in weights[child.name] else 0
                              for child in children]
                node_adj_weights[adj] = sum(children_w) / (d * (len(children) -1 ))

            #if len(node_adj_weights) == 0:
                #print "HELLO"

            weights[node.name] = node_adj_weights
        root = t.find_clades().next()
     #   print root
        anc_weights[label] = weights[root.name]
    #sys.exit(0)
    #print [(x,len(anc_weights[x])) for x in anc_weights.keys()]
    return anc_weights

# Commandline arguments
parser = args.ArgumentParser(description="Enter two genomes in the input format of Unimog")
parser.add_argument('G', metavar='GENOMES', type=str, help="Path to the file that contains the information of extant genomes")
parser.add_argument('T', metavar='TREE', type=str, help="Path to the file that contains the NEWICK tree")
parser.add_argument('-r', '--repetition', default=100, type=int, help="Number of sampled genomes for each ancestor, default: 100")
parser.add_argument('-o', '--output_file', type=str, help="If defined, output is saved in the given file")
parser.add_argument('-a', '--alpha', default=1.0, type=float, help="Tolerance for different tree_distances in the calculation. Closer to 1 equals 0 tolerance")
# 17.11.2016; a new parameter (hype). Just want to see what happens with higher tree-scale if I use different epsilons
parser.add_argument('-e', '--epsilon', default=0.05, type=float, help="Epsilon Parameter. Assigns weights to intermediate adjacencies that are not observed in extant genomes; default:0.05")


arguments = parser.parse_args()

if not 0 <= arguments.alpha <= 1:
    arguments.a = 1
    print >> sys.stderr, "Invalid alpha parameter; Alpha set to 1.\nUse an alpha between 0 and 1."

if not 0 <= arguments.epsilon <= 0.1:
    arguments.e = 0.05
    print >> sys.stderr, "Invalid epsilon parameter; Epsilon set to 0.05.\nUse an epsilon between 0 and 0.1."

# Genome content and tree file are read
input = Input(arguments.G, arguments.T, True)
# Get all pairwise 'siblings' in the tree

pairwise_genomes = input.find_pairwise_leaves(input.tree[0])
all_leaves = input.find_all_leaves(input.tree[0])
max_length = int(input.tree[1])


all_genomes = {}

for leaf in all_leaves:
    all_genomes[leaf.name] = Genome(leaf.name, input.genomes[leaf.name])

anc_weights = find_ancestral_weights(input.tree[0], all_genomes)

gene_number = all_genomes.values()[0].adj_length()
calculate_probability.preprocess_transitions(2*max_length, gene_number)

potential_ancestors = {}


# main iteration
while pairwise_genomes:
    pair = pairwise_genomes.pop()
    names = pair[0]
    lca = input.tree[0].common_ancestor(names)

    print >> sys.stderr, names

    tree_distances = {}
    for genome_name in all_genomes.keys():
        tree_distances[genome_name] = input.tree[0].distance(lca, genome_name)

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

    if tree_distances[names[0]] == 0.0:
       ancestor = first_genome
    elif tree_distances[names[1]] == 0.0:
        ancestor = second_genome
    else:
        ancestor = []
        sampler = Genome_Sampler(inter_info.circular_breakpoint, anc_weights[lca.name], arguments.epsilon)

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

                #19.11 - Small change. Since outgroup information is included in the
                #weighting scheme, it is not needed here anymore.
                if not identifier in names:
                    continue

                breakpoint_graph = IG(genome, candidate)
                breakpoint_graph.create_circular_graph()
                breakpoint_graph = breakpoint_graph.circular_breakpoint

                no_cycles = nx.number_connected_components(breakpoint_graph)
                distance = candidate.length() - no_cycles
                # sometimes occurs if the distances are very short
                # however, a distance of 0 doesn't make much sense
                # in terms of... well.. everything. Number of optimal sorting scenarios of length 0?
                # Number of ALL sorting scenarios of length 0? Both 1?
                if distance == 0:
                    break

                # 13.11 - MOST STUPID thing ever. Switched the tree distance and genomic distance
                # Fix now, write, check results tomorrow.
                # 14.11 - Another fix; the probability is P(G_S) for d_k <= t_k
                # and drops if d_k > t_k. Not the other way around.
                # 17.11 - NOW I got it...

                expected_distance = genome.distance_to_genome(candidate)
                if not isinstance(expected_distance, int):
                    break
                lower_bound = distance*arguments.alpha
                upper_bound = (2-arguments.alpha)*expected_distance

                if tree_distances[identifier] <= lower_bound or tree_distances[identifier] >= upper_bound:
                    break

                sorting_scen = calculate_probability.optimal_scenarios(breakpoint_graph)
                all_scen = calculate_probability.all_scenarios(genome.adj_length(), distance)
                prob_sampled_genomes = sorting_scen - all_scen

                #try:
                if arguments.alpha == 1.0 or distance <= tree_distances[identifier] <= expected_distance:
                    probability += prob_sampled_genomes
                else:
                    if tree_distances[identifier] <= distance:
                        probability += prob_sampled_genomes + math.log10(tree_distances[identifier] - lower_bound) - \
                                        math.log10(distance - lower_bound)
                    else:
                        probability += prob_sampled_genomes + math.log10(upper_bound - tree_distances[identifier]) - \
                                       math.log10(expected_distance*(1-arguments.alpha))
                #except ValueError:
                #    print >> sys.stderr, identifier, tree_distances[identifier], distance, arguments.alpha, distance*arguments.alpha
                #    sys.exit(0)



            # Apparently this is only called if the for-loop did not break
            # This seems to be very fancy!
            else:
                if probability > highest_prob:
                    highest_prob = probability
                    ancestor = candidate


    try:
        ancestor.name = "%s%s" % (names[0], names[1])
    except AttributeError:
        print >> sys.stderr, "No ancestor found. Sorry!"
        continue

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
