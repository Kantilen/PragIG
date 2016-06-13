#!/usr/bin/python

import re
import networkx as nx

__author__ = 'klamkiewicz'

adjacencies_a = {}
adjacencies_b = {}

def connect_adjacencies(adjA, adjB):

    wrapper = [(adjA, adjacencies_a), (adjB, adjacencies_b)]

    for genome,adjacencies in wrapper:
        for adjacency in genome:
            single_extremity = re.search('([\d*\w*]+[t|h])([\d*\w*]+[t|h])', adjacency)
            if not single_extremity:
                adjacencies[adjacency] = None
            else:
                adjacencies[single_extremity.group(1)] = single_extremity.group(2)
                adjacencies[single_extremity.group(2)] = single_extremity.group(1)

    create_circular_graph()

def create_circular_graph():
    G = nx.MultiGraph()
    G.add_nodes_from(adjacencies_a.keys())
    visited = []
    for key, value in adjacencies_a.items():
        if key in visited or value in visited:
            continue
    #TODO: Different telomeres! Currently there is only on telomere node.
        if not value:
            if not adjacencies_b[key]:
                G.add_edge(key, 'Telomere')
                G.add_edge('Telomere', key)
            else:
                G.add_edge(key, 'Telomere')
        else:
            visited.append(value)
            G.add_edge(key,value)
            if adjacencies_b[key]:
                G.add_edge(key, adjacencies_b[key])
            else:
                G.add_edge(key, 'Telomere')
            G.add_edge(value, adjacencies_b[value])
        visited.append(key)

    print G.nodes()
    print G.edges()
    print [x for x in nx.connected_components(G)]