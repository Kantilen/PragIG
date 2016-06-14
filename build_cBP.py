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
        visited.append(key)
        if value:
            visited.append(value)
            G.add_edge(key,value, color='A')

            if adjacencies_b[key] == value or adjacencies_b[value] == key:
                G.add_edge(key, adjacencies_b[key], color='B')
                continue

            if adjacencies_b[key]:
                G.add_edge(key, adjacencies_b[key], color='B')
            if adjacencies_b[value]:
                G.add_edge(value, adjacencies_b[value], color='B')


    components = [x for x in nx.connected_component_subgraphs(G)]

    print [x for x in nx.connected_components(G)]

    index = 0

    for comp in components:

        degree_of_comp = comp.degree()
        #print degree_of_comp
        telomeres = [vertex for vertex,degree in degree_of_comp.items() if degree == 1]

        if len(telomeres) == 2:
            G.add_node("Telo%d" % (index))
            G.add_node("Telo%d" % (index + 1))
            if G.edge[telomeres[0]][telomeres[1]][0]['color'] == 'A':
                G.add_edge("Telo%d" % (index), "Telo%d" % (index + 1), color='A')
                G.add_edge("Telo%d" % (index), telomeres[0], color='B')
                G.add_edge("Telo%d" % (index + 1), telomeres[1], color='B')
            else:
                G.add_edge("Telo%d" % (index), "Telo%d" % (index + 1), color='B')
                G.add_edge("Telo%d" % (index), telomeres[0], color='A')
                G.add_edge("Telo%d" % (index + 1), telomeres[1], color='A')

        elif len(telomeres) == 1:
            print "Do one telomere"

    print [x for x in nx.connected_components(G)]
    print G.edges()