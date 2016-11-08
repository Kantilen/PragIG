#!/usr/bin/python

#################################
# Import section                #
#################################
import argparse as args
import sys

from model import Genome
from input_parser import Input
#################################

__author__ = 'klamkiewicz'

parser = args.ArgumentParser(description="Calculates the TP, FP and FN for one input file. Input has to be in the GRIMM format.")
parser.add_argument('I', metavar='INPUT', type=str, help="Path to the txt-file produced by pragig.py")
parser.add_argument('R', metavar='REFERENCE', type=str, help="Path to the reference file with simulated ancestral genomes")
parser.add_argument('T', metavar='TREE', type=str, help="Path to the file that contains the complete tree, with ancestral genome annotations")

arguments = parser.parse_args()

reference = Input(arguments.R, arguments.T, False)

calculated_ancestors = reference.read_genomes(arguments.I)
leaves = reference.find_pairwise_leaves(reference.tree[0])

no_internal_nodes = len(reference.tree[0].get_nonterminals())
no_leaves = len(reference.tree[0].get_terminals())
no_nodes = no_leaves + no_internal_nodes
not_reconstructed = no_nodes - len(calculated_ancestors)

#print not_reconstructed, no_nodes, no_leaves, no_internal_nodes, len(calculated_ancestors)

print "Node;Dist;TP;FP;NP;TP%;FP%;NP%"

for i in range(not_reconstructed):
    print "NA;NA;0;0;0;0;0;0"

while leaves:
    pair = leaves.pop()
    current_leaves = pair[0]
    distance = pair[1][0] + pair[1][1]
    old_key = reference.tree[0].common_ancestor(current_leaves).name
    content = reference.genomes[old_key]
    new_key = "%s%s" % (current_leaves[0], current_leaves[1])
    reference.genomes.pop(old_key)
    reference.tree[0].common_ancestor(current_leaves).name = new_key
    reference.genomes[new_key] = content

    genome_reference = Genome("ref", content)
    try:
        genome_calculated = Genome("calc", calculated_ancestors[new_key])
    except KeyError:
        #print "%s;%d;0;0;%d;0;0;1" % (new_key, distance, genome_reference.adj_length())
        continue

    true_positives = sum(1 for adj in genome_calculated.adjacency_set if adj in genome_reference.adjacency_set)
    false_positives = sum(1 for adj in genome_calculated.adjacency_set if not adj in genome_reference.adjacency_set)
    false_negatives = sum(1 for adj in genome_reference.adjacency_set if not adj in genome_calculated.adjacency_set)

    length_reference = genome_reference.adj_length()

    print "%s;%d;%d;%d;%d;%f;%f;%f" % \
    (new_key, distance, true_positives, false_positives, false_negatives, float(true_positives)/length_reference, float(false_positives)/length_reference, float(false_negatives)/length_reference)

    reference.tree[0].collapse(current_leaves[0])
    reference.tree[0].collapse(current_leaves[1])

    leaves = reference.find_pairwise_leaves(reference.tree[0])