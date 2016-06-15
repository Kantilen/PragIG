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

    for first, second, color in all_edges:
        if first.startswith('Telo') and not second.startswith('Telo'):
            intermediate_adjacencies.append('%s' % (second))
            continue
        if second.startswith('Telo') and not first.startswith('Telo'):
            intermediate_adjacencies.append('%s' % (first))
            continue
        if first.startswith('Telo') and second.startswith('Telo'):
            continue
        intermediate_adjacencies.append('%s%s' % (first, second))

    just_a_edges = [x for x in all_edges if x[2]['color'] == 'A']
    print just_a_edges
    print
    print
    print set(intermediate_adjacencies)