#!/usr/bin/python

#################################
# Import section                #
#################################
import re
import networkx as nx
#################################

__author__ = 'klamkiewicz'

#################################
# Global definitions            #
#################################
adjacencies_a = {}
adjacencies_b = {}
colors = {'A' : 'B',
          'B' : 'A'}
#################################

def connect_adjacencies(adjA, adjB):
    '''
    From a given adjacency set the single extremities are created and stored into dicts.
    Each extremity points to its adjacency (or None, if it is a telomere)
    :param adjA: Adjacency Set of the first genome
    :param adjB: Adjacency Set of the second genome
    :return: Returns the value of create_circular_graph()
    '''
    wrapper = [(adjA, adjacencies_a), (adjB, adjacencies_b)]

    for genome,adjacencies in wrapper:
        for adjacency in genome:
            # TODO: This regex might get problems, if gene marker are called 't' or 'h'
            single_extremity = re.search('([\d*\w*]+[t|h])([\d*\w*]+[t|h])', adjacency) # Regex to find single extremities
            if not single_extremity:
                adjacencies[adjacency] = None # telomere
            else:
                adjacencies[single_extremity.group(1)] = single_extremity.group(2)
                adjacencies[single_extremity.group(2)] = single_extremity.group(1)

    return create_circular_graph() # call function to create the cBP

def create_circular_graph():
    '''
    Uses the networkx package to create the circular breakpoint graph of two given genomes
    :return: the circular breakpoint graph
    '''
    G = nx.MultiGraph()

    # Since no Indel or Duplications are allowed, the keys of adjacencies_a/b are identical
    G.add_nodes_from(adjacencies_a.keys()) #each adjacency will be a vertex
    visited = []
    # iterate through all extremities
    for key, value in adjacencies_a.items():
        if key in visited or value in visited:
            continue
        visited.append(key)
        if value: # if it is a real adjacency
            visited.append(value)
            G.add_edge(key,value, color='A') # A-edge is added

            if adjacencies_b[key] == value or adjacencies_b[value] == key: # 2-cycle
                G.add_edge(key, adjacencies_b[key], color='B') # B-edge is added as well
                continue

            if adjacencies_b[key]: # If the extremity is no telomere in the other genome
                G.add_edge(key, adjacencies_b[key], color='B') # connect it with B-edge
            if adjacencies_b[value]: # If neighbor of current extremity is no telomere in other genome
                G.add_edge(value, adjacencies_b[value], color='B') # connect it with B-edge

    # This part is to create the CIRCULAR breakpoint graph
    components = [x for x in nx.connected_component_subgraphs(G)] # get all components in a list
    index = 0

    for comp in components:

        if len(comp) == 1: # telomere in both genomes
            G.add_node("Telo%d" % (index))
            G.add_edge("Telo%d" % (index), comp.nodes()[0], color='A')
            G.add_edge("Telo%d" % (index), comp.nodes()[0], color='B')
            index += 1
            continue

        degree_of_comp = comp.degree()
        telomeres = [vertex for vertex,degree in degree_of_comp.items() if degree == 1]

        if telomeres: # if False, there are no telomeres --> already cycle
            data_telo = G.edges(telomeres, True) # get colors of the telomeres

            # This is a tricky networkx case. If both telomeres have the same edge, it returns only
            # a one-element list. Therefore I have this small hack here...
            if len(data_telo) == 1:
                first_color = second_color = data_telo[0][2]['color']
                first_telo = data_telo[0][0]
                second_telo = data_telo[0][1]

            # This is the "normal case"
            else:
                first_color = data_telo[0][2]['color']
                second_color = data_telo[1][2]['color']
                first_telo = data_telo[0][0]
                second_telo = data_telo[1][0]

            # Same color of edges means I have to add 2 telomeres, since the length of the component is even
            if first_color == second_color:
                G.add_node("Telo%d" % (index))
                G.add_node("Telo%d" % (index + 1))
                G.add_edge("Telo%d" % (index), "Telo%d" % (index + 1), color=first_color)
                G.add_edge(first_telo, "Telo%d" % (index), color=colors[first_color])
                G.add_edge(second_telo, "Telo%d" % (index + 1), color=colors[first_color])
                index += 2
                continue

            # otherwise it is odd --> add one telomere
            else:
                G.add_node("Telo%d" % (index))
                for telo,neighbor,color in data_telo: # both extremities are added to the telomere
                    G.add_edge(telo, "Telo%d" % (index), color=colors[color['color']]) # take the opposite color
                index += 1
                continue
    return G