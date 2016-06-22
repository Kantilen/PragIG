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
    '''
    Obsolete function! But it does the right thing, therefore I am able to implement good stuff....
    :param graph:
    :return:
    '''
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
    '''
    This function enumerates all vertices and finds the intermediate adjacencies.
    :param graph: circular breakpoint graph
    :return: set of all intermediate adjacenices
    '''

    # each component can be solved seperately
    for component in nx.connected_component_subgraphs(graph):
        #if len(component.nodes()) == 2: # cycles of length 2 are fixed!
        #    continue

        enumerated_vertices = {}

        # assign value for each vertex from 1 to (n+1)
        for index, vertex in enumerate(component.nodes()):
            enumerated_vertices[vertex] = index+1

        # This is tricky. two vertices whose difference are odd form an intermediate adjacency.
        # I have to ask Pedro for the theoretical background here.
        odd_adjacencies = [(x,y) for x in enumerated_vertices.keys() for y in enumerated_vertices.keys()
                      if (abs(enumerated_vertices[x] - enumerated_vertices[y]) % 2 == 1)]
        #print odd_adjacencies
        # This for loops kicks the artifical telomeres.
        for first,second in odd_adjacencies:
            if first.startswith('Telo') and second.startswith('Telo'):
                continue
            if first.startswith('Telo'):
                intermediate_adjacencies.add(second)
                continue
            if second.startswith('Telo'):
                intermediate_adjacencies.add(first)
                continue
            # some nasty set issue. Since we are dealing with strings, there is a difference
            # between 1h2t and 2t1h. We might change this at some point.
            if not '%s%s' % (second,first) in intermediate_adjacencies:
                intermediate_adjacencies.add('%s%s' % (first,second))

def find_all_adjacencies(graph):
    '''
    Wrapping function that is called by the mainflow script
    :param graph: circular breakpoint graph
    :return: set of all intermediate adjacencies.
    '''
    get_all_inter_adj(graph)
    return intermediate_adjacencies