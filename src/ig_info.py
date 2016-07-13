    #!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import re
import networkx as nx
import numpy as np
from model import Adjacency
#################################

class Intermediate_Genome():
    '''
    This class contains all information of two given (sibling) genomes in the NEWICK tree.
    It stores the content of the genomes, the adjacency sets, the circular breakpoint graph between the genomes
    and all possible intermediate adjacencies.
    '''
    def __init__(self, genomeA, genomeB):
        '''
        Initialization. Nothing real exciting here.
        :param genomeA: gene content of the first genome (list of identifier)
        :param genomeB: gene content of the second genome (list of identifier)
        '''

        self.genomeA = genomeA.content
        self.genomeB = genomeB.content

        self.first = genomeA.name
        self.second = genomeB.name

        self.genomes = {self.first:self.genomeA, self.second:self.genomeB}

        self.adjacencies = {self.first:genomeA.adjacency_set, self.second:genomeB.adjacency_set}

        self.circular_breakpoint = None
        self.inter_adj = set()

        self.binaries = {}

    def validate_input(self):
        '''
        At the moment no duplications or indel events are allowed. Therefore the set of the genome content
        has to be identical. This is checked here. The two contents of the global variables genomeA and genomeB
        are evaluated.
        :return: Boolean variable
        '''
        # remove signs and chromosomes
        first_genome = [re.sub('-', '', x) for x in self.genomeA]
        second_genome = [re.sub('-', '', x) for x in self.genomeB]

        # symmetrical difference. Just take elements that are unique in one
        # of the sets
        first_genome = set(first_genome)
        second_genome = set(second_genome)
        difference = first_genome ^ second_genome

        return (len(difference) == 0, difference)

    def create_circular_graph(self):
        '''
        Uses the networkx package to create the circular breakpoint graph of two given genomes.
        First a dictionary for each adjacency set is created. Every extremity points to its adjacent extremity.
        Afterwards the cBP is built with the help of these dictionaries.
        The final graph is stored in the global variable circular_breakpoint.
        '''
        adjacencies_a = {} # adjacency dict
        adjacencies_b = {}

        adj_set_A = self.adjacencies[self.first]
        adj_set_B = self.adjacencies[self.second]

        # wrapper for the iteration
        wrapper = [(adj_set_A, adjacencies_a), (adj_set_B, adjacencies_b)]

        for genome, adjacencies in wrapper:
            for adjacency in genome:
                if adjacency.is_telomere():
                    adjacencies[adjacency.first_ex] = None
                else:
                    adjacencies[adjacency.first_ex] = adjacency.second_ex
                    adjacencies[adjacency.second_ex] = adjacency.first_ex

        colors = {'A': 'B',
                  'B': 'A'}
        G = nx.MultiGraph()
        # Since no Indel or Duplications are allowed, the keys of adjacencies_a/b are identical
        G.add_nodes_from(adjacencies_a.keys())  # each adjacency will be a vertex
        visited = []
        # iterate through all extremities
        for key, value in adjacencies_a.items():
            if key in visited or value in visited:
                continue
            visited.append(key)
            if value:  # if it is a real adjacency
                visited.append(value)
                G.add_edge(key, value, color='A')  # A-edge is added

                if adjacencies_b[key] == value or adjacencies_b[value] == key:  # 2-cycle
                    G.add_edge(key, adjacencies_b[key], color='B')  # B-edge is added as well
                    continue

                if adjacencies_b[key]:  # If the extremity is no telomere in the other genome
                    if not G.get_edge_data(key, adjacencies_b[key]):
                        G.add_edge(key, adjacencies_b[key], color='B')  # connect it with B-edge
                if adjacencies_b[value]:  # If neighbor of current extremity is no telomere in other genome
                    if not G.get_edge_data(value, adjacencies_b[value]):
                        G.add_edge(value, adjacencies_b[value], color='B')  # connect it with B-edge

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
                    G.add_edge(first_telo, "Telo%d" % (index), color=colors[first_color])
                    G.add_edge(second_telo, "Telo%d" % (index + 1), color=colors[first_color])
                    index += 2
                    continue

                # otherwise it is odd --> add one telomere
                else:
                    G.add_node("Telo%d" % (index))
                    for telo, neighbor, color in data_telo:  # both extremities are added to the telomere
                        G.add_edge(telo, "Telo%d" % (index),
                                   color=colors[color['color']])  # take the opposite color
                    index += 1
                    continue
        self.circular_breakpoint = G

    def get_all_inter_adj(self):
        '''
        This function enumerates all vertices and finds the intermediate adjacencies.
        :param graph: circular breakpoint graph
        :return: set of all intermediate adjacencies
        '''
        # each component can be solved seperately
        for component in nx.connected_component_subgraphs(self.circular_breakpoint):

            enumerated_vertices = {}

            first_adj = component.edges()[0]
            long_path = [x for x in nx.all_simple_paths(component, first_adj[0], first_adj[1])][-1]
            # assign value for each vertex from 0 to n
            for index, vertex in enumerate(long_path):
                enumerated_vertices[index] = vertex
            # This is tricky. two vertices whose difference are odd form an intermediate adjacency.
            # I have to ask Pedro for the theoretical background here.
            for i in range(len(enumerated_vertices.keys())):
                for j in range(i+1,len(enumerated_vertices.keys()),2):
                    first = enumerated_vertices[i]
                    second = enumerated_vertices[j]

                    if first.startswith('Telo') and second.startswith('Telo'):
                        continue
                    if first.startswith('Telo'):
                        self.inter_adj.add(Adjacency(second,None))
                        continue
                    if second.startswith('Telo'):
                        self.inter_adj.add(Adjacency(first,None))
                        continue

                    new_adj = Adjacency(first,second)
                    if not new_adj.is_in_list(self.inter_adj):
                        self.inter_adj.add(new_adj)
        self.inter_adj = list(self.inter_adj)

    def create_binary_vector(self):
        '''
        This function creates an numpy array that represents the binary vector of
        adjacencies in the two genomes that have to compared.
        '''
        for indent,adj_set in self.adjacencies.items():
            binary = np.zeros([1,len(self.inter_adj)], dtype=int)


