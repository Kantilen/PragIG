#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import argparse as args
import os
import sys
import copy
from collections import defaultdict

import networkx as nx

from input_parser import Input
from genome_sampler import Genome_Sampler
from model import Genome
from ig_info import Intermediate_Genome as IG
import calculate_probability
#################################

#################################
# FUNCTIONS                     #
#################################

def write_output():
    '''
    This function writes the output of PragIG to the specified folder.
    If the folder does not exist, it is created.
    :return:
    '''
    print >> sys.stderr, "Writing output files"

    folder = "%s/PragIG_A%s_E%s_R%d/" % (arguments.O, str(arguments.alpha), str(arguments.epsilon), arguments.repetition)

    if not os.path.exists(folder):
        os.mkdir(folder)

    ancestor_file = folder + "pragig_ancestral.txt"
    probability_file = folder + "pragig_probabilities.txt"

    ancestor = open(ancestor_file, 'w')
    probability = open(probability_file, 'w')

    for result_name, result_content in input.genomes.items():
        result = Genome(result_name, result_content)
        chromosome = result.chr_number()

        ancestor.write('>%s\n' % (result_name))
        for i in range(chromosome):
            ancestor.write('#chr%d\n' % (i + 1))
            for index, gene in enumerate(result_content):
                ancestor.write('%s ' % (gene))
                if gene == ')' or gene == '$':
                    result_content = result_content[index+1:]
                    ancestor.write('\n')
                    break

        if all_probabilities.has_key(result_name):
            probability.write('>%s\n' % (result_name))
            probability.write(' '.join(map(lambda n: '%.8f'%n ,all_probabilities[result_name])))
            probability.write('\n')


def find_ancestral_weights(tree, extant_genomes):
    '''
    Finds weights for adjacencies present in extant_genomes for each internal node in tree.
    Weighting scheme is adapted from Feijao and Araujo, 2016 (RINGO).
    :param tree: Newick tree that holds information of the relation between genomes
    :param extant_genomes: leaf genomes with conserved adjacencies
    :return: dictionary with all ancestral adjacency weights
    '''
    anc_weights = {}
    all_internal_nodes = tree.get_nonterminals()

    for ancestor in all_internal_nodes:

        label = ancestor.name
        # deepcopy is needed, otherwise the original tree gets modified
        t = copy.deepcopy(tree)

        # find the root of original tree
        root = tree.get_nonterminals(order="level")[0]

        # rerooting step in order to include outgroup information during the weighting process
        if ancestor != root:
            current = tree.find_clades(ancestor, order="level").next()
            children = list(current.find_clades(order="postorder"))
            children = [x.name for x in children if x != ancestor]
            for child in children:
                t.collapse(child)
            t.root_with_outgroup({'name' : current.name})

        weights = {}
        for idx, node in enumerate(t.find_clades(order="postorder")):
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
            if d == 0:
                d = 0.1
            node_adj_weights = {}
            for adj in all_adj:
                children_w = [weights[child.name][adj] * (d - child.branch_length) if adj in weights[child.name] else 0
                              for child in children]
                node_adj_weights[adj] = sum(children_w) / (d * (len(children) -1 ))


            weights[node.name] = node_adj_weights
        root = t.find_clades().next()
        anc_weights[label] = weights[root.name]
    return anc_weights

def get_dcj_distance_from_BP(genome, candidate):
    '''
    Creates the circular breakpint graph between genome and candidate and calculates
    the DCJ distance between the two genomes:
        d = n - c
    :param genome: First genome (object)
    :param candidate: Second genome (object, sampled ancestor)
    :return: DCJ distance between genome and candidate
    '''
    breakpoint_graph = IG(genome, candidate)
    breakpoint_graph.create_circular_graph()
    breakpoint_graph = breakpoint_graph.circular_breakpoint
    no_cycles = nx.number_connected_components(breakpoint_graph)
    return (candidate.length() - no_cycles, breakpoint_graph)

#################################

# Commandline arguments
parser = args.ArgumentParser(description="Reconstructs ancestor genomes from a given extant genome file (GRIMM) and a tree (NEWICK)")
parser.add_argument('G', metavar='GENOMES', type=str, help="Path to the file that contains the information of extant genomes. GRIMM format supported.")
parser.add_argument('T', metavar='TREE', type=str, help="Path to the file that contains the NEWICK tree")
parser.add_argument('O', metavar='OUTPUT', type=str, help="Path to folder, where results are saved.")
parser.add_argument('-r', '--repetition', default=100, type=int, help="Number of sampled genomes for each ancestor, default: 100")
parser.add_argument('-a', '--alpha', default=1.0, type=float, help="Strictness of the filter that discards sampled genomes. Closer to 0 means no genomes are discarded. Default = 1.0")
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

all_genomes = {}
for leaf in all_leaves:
    all_genomes[leaf.name] = Genome(leaf.name, input.genomes[leaf.name])

# find weights for internal nodes
anc_weights = find_ancestral_weights(input.tree[0], all_genomes)

all_probabilities = defaultdict(list)

# main iteration, as long as siblings are found
while pairwise_genomes:
    pair = pairwise_genomes.pop()
    names = pair[0]
    lca = input.tree[0].common_ancestor(names)

    print >> sys.stderr, "Resolving %s%s..." % (names[0], names[1])

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
    probs = []
    if tree_distances[names[0]] == 0.0:
       ancestor = first_genome

    elif tree_distances[names[1]] == 0.0:
        ancestor = second_genome

    else:
        ancestor = []
        sampler = Genome_Sampler(inter_info.circular_breakpoint, anc_weights[lca.name], arguments.epsilon)

        highest_prob = None

        # sample different genomes as potential ancestors
        for i in range(arguments.repetition):
            candidate = sampler.enumerate_vertices()
            probability = 0

            for identifier, genome in all_genomes.items():
                #19.11 - Small change. Since outgroup information is included in the
                #weighting scheme, it is not needed here anymore.
                if not identifier in names:
                    continue

                distance, breakpoint_graph = get_dcj_distance_from_BP(genome, candidate)
                # sometimes occurs if the distances are very short
                # however, a distance of 0 doesn't make much sense
                # in terms of... well.. everything. Number of optimal sorting scenarios of length 0?
                # Number of ALL sorting scenarios of length 0? Both 1?
                if distance == 0:
                    break

                expected_distance = genome.expected_distance_to_genome(candidate)
                # whenenver the calculation raises an ValueError
                # it is caught here and the next iteration begins
                if not isinstance(expected_distance, int):
                    break

                # lower and upper bound for the probability function
                lower_bound = distance*arguments.alpha
                upper_bound = (2-arguments.alpha)*expected_distance

                if tree_distances[identifier] <= lower_bound or tree_distances[identifier] >= upper_bound:
                    break

                # number of optimal sorting scenarions between two genomes
                sorting_scen = calculate_probability.optimal_scenarios(breakpoint_graph)
                # number of all sorting scenarios of given distance
                all_scen = calculate_probability.all_scenarios(genome.adj_length(), distance)
                # new probability for the ancestor
                probability = calculate_probability.calculate_prob_ancestor(probability, sorting_scen, all_scen, arguments.alpha,
                                                                            distance, tree_distances[identifier], expected_distance,
                                                                            lower_bound, upper_bound)
            # Apparently this is only called if the for-loop did not break
            # This seems to be very fancy!
            else:
                probs.append(probability)
                if probability > highest_prob:
                    highest_prob = probability
                    ancestor = candidate


    try:
        ancestor.name = "%s%s" % (names[0], names[1])
    except AttributeError: # if all sampled genomes are discarded, ancestor = [];
        print >> sys.stderr, "No ancestor found. Sorry!"
        continue

    # Tree and Genomes get updated for the next iteration
    all_probabilities[ancestor.name] = probs

    # rename the internal node
    clade = input.tree[0].common_ancestor(names[0],names[1])
    clade.name = ancestor.name

    # cut the current leaves from the tree and remove them from the genome list
    input.tree[0].collapse(names[0])
    input.tree[0].collapse(names[1])
    all_genomes.pop(names[0])
    all_genomes.pop(names[1])

    # update the list of input genomes to find new pairwise genomes
    input.genomes.update({ancestor.name:ancestor.content})
    all_genomes[ancestor.name] = ancestor
    pairwise_genomes = input.find_pairwise_leaves(input.tree[0])
###############################################################################################

write_output()
