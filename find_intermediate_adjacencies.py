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
intermediate_adjacencies = set()
#################################

def perform_DCJ(graph):

    components = [x for x in nx.connected_component_subgraphs(graph)]

    for component in components:
        all_edges = component.edges(data=True)

        for first, second, color in all_edges:
            if first.startswith('Telo') and not second.startswith('Telo'):
                intermediate_adjacencies.add('%s' % (second))
                continue
            if second.startswith('Telo') and not first.startswith('Telo'):
                intermediate_adjacencies.add('%s' % (first))
                continue
            if first.startswith('Telo') and second.startswith('Telo'):
                continue
            if not '%s%s' % (second, first) in intermediate_adjacencies:
                intermediate_adjacencies.add('%s%s' % (first, second))

        if len(all_edges) == 2:
            continue

        just_a_edges = [x for x in all_edges if x[2]['color'] == 'A']

        while just_a_edges:
            p,q,color_current = just_a_edges.pop()
            for r,s,color in just_a_edges:

                new_graph = component.copy()
                new_graph.remove_edge(p,q)
                new_graph.remove_edge(r,s)
                new_graph.add_edge(p,r,color='A')
                new_graph.add_edge(q,s,color='A')
                if len([x for x in nx.connected_components(new_graph)]) > 1:
                    perform_DCJ(new_graph)

                new_graph = component.copy()
                new_graph.remove_edge(p, q)
                new_graph.remove_edge(r, s)
                new_graph.add_edge(p, s, color='A')
                new_graph.add_edge(q, r, color='A')
                if len([x for x in nx.connected_components(new_graph)]) > 1:
                    perform_DCJ(new_graph)
    return intermediate_adjacencies

def get_all_inter_adj(graph):
    for component in nx.connected_component_subgraphs(graph):
        if len(component.nodes()) == 2:
            continue

        enumerated_vertices = {}

        for index, vertex in enumerate(component.nodes()):
            enumerated_vertices[vertex] = index+1

        odd_adjacencies = [(x,y) for x in enumerated_vertices.keys() for y in enumerated_vertices.keys()
                      if (abs(enumerated_vertices[x] - enumerated_vertices[y]) % 2 == 1)]


        for first,second in odd_adjacencies:
            if first.startswith('Telo'):
                intermediate_adjacencies.add(second)
                continue
            if second.startswith('Telo'):
                intermediate_adjacencies.add(first)
                continue
            if first.startswith('Telo') and second.startswith('Telo'):
                continue
            if not '%s%s' % (second,first) in intermediate_adjacencies:
                intermediate_adjacencies.add('%s%s' % (first,second))
    #return intermediate_adjacencies

def find_all_adjacencies(graph):
    get_all_inter_adj(graph)
    return intermediate_adjacencies