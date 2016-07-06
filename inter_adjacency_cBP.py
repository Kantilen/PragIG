#!/usr/bin/python

#################################
# Import section                #
#################################
import networkx as nx
import re
#################################

__author__ = 'klamkiewicz'

class Inter_Adjacencies():
    '''
    Given a (networkx) graph, this class creates the set of ALL intermediate adjacencies
    that arise between to genomes represented by the circular breakpoint graph.
    DEPENDENCY: networkx
    '''
    def __init__(self, graph):
        '''
        Initialization. Also calls the method to create the intermediate-adjacency set
        :param graph: circular breakpoint graph
        '''
        self.graph = graph
        self.intermediate_adjacencies = set()
        self.get_all_inter_adj(self.graph)

    def get_all_inter_adj(self, graph):
        '''
        This function enumerates all vertices and finds the intermediate adjacencies.
        :param graph: circular breakpoint graph
        :return: set of all intermediate adjacencies
        '''
        # each component can be solved seperately
        for component in nx.connected_component_subgraphs(graph):

            enumerated_vertices = {}

            first_adj = component.edges()[0]
            long_path = [x for x in nx.all_simple_paths(component, first_adj[0], first_adj[1])][-1]
            # assign value for each vertex from 1 to (n+1)
            for index, vertex in enumerate(long_path):
                enumerated_vertices[vertex] = index+1

            # This is tricky. two vertices whose difference are odd form an intermediate adjacency.
            # I have to ask Pedro for the theoretical background here.
            odd_adjacencies = [(x,y) for x in enumerated_vertices.keys() for y in enumerated_vertices.keys()
                          if (abs(enumerated_vertices[x] - enumerated_vertices[y]) % 2 == 1)]

            # This for loops kicks the artifical telomeres.
            for first,second in odd_adjacencies:
                if first.startswith('Telo') and second.startswith('Telo'):
                    continue
                if first.startswith('Telo'):
                    self.intermediate_adjacencies.add(second)
                    continue
                if second.startswith('Telo'):
                    self.intermediate_adjacencies.add(first)
                    continue
                # some nasty set issue. Since we are dealing with strings, there is a difference
                # between 1h2t and 2t1h. We might change this at some point.
                if not '%s%s' % (second,first) in self.intermediate_adjacencies:
                    self.intermediate_adjacencies.add('%s%s' % (first,second))

class Circular_Breakpoint():
    '''
    This class creates the circular breakpoint of to given genomes.
    Note that you have to input the adjacency sets and not the genomes themselves.
    DEPENDENCY: networkx
    '''
    def __init__(self, adjA, adjB):
        '''
        Initialization. Also calls the method to create the circular breakpoint graph
        :param adjA: list of adjacencies of the first genome
        :param adjB: list of adjacencies of the second genome
        '''
        self.adjA = adjA
        self.adjB = adjB

        self.adjacencies_a = {}
        self.adjacencies_b = {}

        # This dictionary is used to create "fictional" edges between two added telomeres
        self.colors = {'A': 'B',
                       'B': 'A'}

        self.graph = self.connect_adjacencies(self.adjA, self.adjB)

    def connect_adjacencies(self, adjA, adjB):
        '''
        From a given adjacency set the single extremities are created and stored into dicts.
        Each extremity points to its adjacency (or None, if it is a telomere)
        :param adjA: Adjacency Set of the first_content genome
        :param adjB: Adjacency Set of the second_content genome
        :return: Returns the value of create_circular_graph()
        '''
        wrapper = [(adjA, self.adjacencies_a), (adjB, self.adjacencies_b)]

        for genome, adjacencies in wrapper:
            for adjacency in genome:
                # TODO: This regex might get problems, if gene marker are called 't' or 'h'
                single_extremity = re.search('([\d*\w*]+[t|h])([\d*\w*]+[t|h])',
                                             adjacency)  # Regex to find single extremities
                if not single_extremity:
                    adjacencies[adjacency] = None  # telomere
                else:
                    adjacencies[single_extremity.group(1)] = single_extremity.group(2)
                    adjacencies[single_extremity.group(2)] = single_extremity.group(1)

        return self.create_circular_graph()  # call function to create the cBP

    def create_circular_graph(self):
        '''
        Uses the networkx package to create the circular breakpoint graph of two given genomes
        :return: the circular breakpoint graph
        '''
        G = nx.MultiGraph()

        # Since no Indel or Duplications are allowed, the keys of adjacencies_a/b are identical
        G.add_nodes_from(self.adjacencies_a.keys())  # each adjacency will be a vertex
        visited = []
        # iterate through all extremities
        for key, value in self.adjacencies_a.items():
            if key in visited or value in visited:
                continue
            visited.append(key)
            if value:  # if it is a real adjacency
                visited.append(value)
                G.add_edge(key, value, color='A')  # A-edge is added

                if self.adjacencies_b[key] == value or self.adjacencies_b[value] == key:  # 2-cycle
                    G.add_edge(key, self.adjacencies_b[key], color='B')  # B-edge is added as well
                    continue

                if self.adjacencies_b[key]:  # If the extremity is no telomere in the other genome
                    if not G.get_edge_data(key, self.adjacencies_b[key]):
                        G.add_edge(key, self.adjacencies_b[key], color='B')  # connect it with B-edge
                if self.adjacencies_b[value]:  # If neighbor of current extremity is no telomere in other genome
                    if not G.get_edge_data(value, self.adjacencies_b[value]):
                        G.add_edge(value, self.adjacencies_b[value], color='B')  # connect it with B-edge

        # This part is to create the CIRCULAR breakpoint graph
        components = [x for x in nx.connected_component_subgraphs(G)]  # get all components in a list
        index = 0

        for comp in components:

            if len(comp) == 1:  # telomere in both genomes
                G.add_node("Telo%d" % (index))
                G.add_edge("Telo%d" % (index), comp.nodes()[0], color='A')
                G.add_edge("Telo%d" % (index), comp.nodes()[0], color='B')
                index += 1
                continue

            degree_of_comp = comp.degree()
            telomeres = [vertex for vertex, degree in degree_of_comp.items() if degree == 1]

            if telomeres:  # if False, there are no telomeres --> already cycle
                data_telo = G.edges(telomeres, True)  # get colors of the telomeres

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
                    G.add_edge(first_telo, "Telo%d" % (index), color=self.colors[first_color])
                    G.add_edge(second_telo, "Telo%d" % (index + 1), color=self.colors[first_color])
                    index += 2
                    continue

                # otherwise it is odd --> add one telomere
                else:
                    G.add_node("Telo%d" % (index))
                    for telo, neighbor, color in data_telo:  # both extremities are added to the telomere
                        G.add_edge(telo, "Telo%d" % (index), color=self.colors[color['color']])  # take the opposite color
                    index += 1
                    continue
        return G