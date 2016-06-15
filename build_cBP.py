#!/usr/bin/python

import re
import networkx as nx

__author__ = 'klamkiewicz'

adjacencies_a = {}
adjacencies_b = {}
colors = {'A' : 'B',
          'B' : 'A'}

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

    return create_circular_graph()

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
    index = 0

    for comp in components:

        if len(comp) == 1:
            G.add_node("Telo%d" % (index))
            G.add_edge("Telo%d" % (index), comp.nodes()[0], color='A')
            G.add_edge("Telo%d" % (index), comp.nodes()[0], color='B')
            index += 1
            continue

        degree_of_comp = comp.degree()
        telomeres = [vertex for vertex,degree in degree_of_comp.items() if degree == 1]

        if len(telomeres) == 2:
            data_telo = G.edges(telomeres, True)
            if len(data_telo) == 1:
                first_color = second_color = data_telo[0][2]['color']
                first_telo = data_telo[0][0]
                second_telo = data_telo[0][1]
            else:
                first_color = data_telo[0][2]['color']
                second_color = data_telo[1][2]['color']
                first_telo = data_telo[0][0]
                second_telo = data_telo[1][0]

            if first_color == second_color:
                G.add_node("Telo%d" % (index))
                G.add_node("Telo%d" % (index + 1))
                G.add_edge("Telo%d" % (index), "Telo%d" % (index + 1), color=first_color)
                G.add_edge(first_telo, "Telo%d" % (index), color=colors[first_color])
                G.add_edge(second_telo, "Telo%d" % (index + 1), color=colors[first_color])
                index += 2
                continue
            else:
                G.add_node("Telo%d" % (index))
                for telo,neighbor,color in data_telo:
                    G.add_edge(telo, "Telo%d" % (index), color=colors[color['color']])
                index += 1
                continue
    return G