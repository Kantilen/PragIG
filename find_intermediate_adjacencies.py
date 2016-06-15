#!/usr/bin/python

#################################
# Import section                #
#################################
import networkx as nx
#################################

__author__ = 'klamkiewicz'

#################################
# Global definitions            #
#################################
intermediate_adjacencies = []
#################################

def perform_DCJ(graph):

    all_edges = graph.edges(data=True)

    just_a_edges = [x for x in all_edges if x[2]['color'] == 'A']
    just_b_edges = [x for x in all_edges if x[2]['color'] == 'B']
    print just_a_edges
    print
    print
    print just_b_edges
    # TODO: Fix telomeres
    # TODO: Make set out of list.
    for first, second, color in just_a_edges:
        intermediate_adjacencies.append('%s%s' % (first, second))
    for first, second, color in just_b_edges:
        intermediate_adjacencies.append('%s%s' % (first, second))

    print intermediate_adjacencies