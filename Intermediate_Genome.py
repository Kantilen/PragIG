#!/usr/bin/python

__author__ = 'klamkiewicz'

import re
import networkx as nx

class Intermediate_Genome():
    def __init__(self, genomeA, genomeB):
        self.genomeA = genomeA
        self.genomeB = genomeB

        self.genomes = {'A':self.genomeA, 'B':self.genomeB}

        self.adj_set_A = set()
        self.adj_set_B = set()

        self.adjacencies = {'A':self.adj_set_A, 'B':self.adj_set_B}

        self.circular_breakpoint = None
        self.inter_adj = set()

    def create_adjacency_sets(self):
        '''
        This function reads the genome content and creates the adjacency set of it.
        Note that the genome should be represented in the UniMog input notation.
        :param genome_content: UniMog representation of the genome. Type=List
        :return: Adjacency set of the genome. Type=List (of strings)
        '''

        for ident,genome_content in self.genomes.items():
            adjacencies = []  # returned value
            current_chromosome = []  # this is needed if there are more chromosomes

            # if not sign is given, circular chromosomes are default
            if not (genome_content[-1] == ')' or genome_content[-1] == '|'):
                genome_content.append(')')

            # main iteration
            for index, gene in enumerate(genome_content):
                if not (re.match('[\d*\w*]+', gene) or gene.startswith('-')):  # chromosome sign detected

                    # the first_content and last extremity have to be considered separately
                    first = current_chromosome[0]
                    last = current_chromosome[-1]
                    current_chromosome.remove(first)
                    current_chromosome.remove(last)

                    # connecting adjacencies
                    current_chromosome = iter(current_chromosome)
                    adjacencies.extend([ex + next(current_chromosome) for ex in current_chromosome])

                    if gene == '|':  # linear chromosome needs telomeres
                        adjacencies.insert(0, first)
                        adjacencies.append(last)
                    if gene == ')':  # circular chromosomes needs another adjacency
                        adjacencies.append('%s%s' % (last, first))

                    current_chromosome = []  # clear the chromosome
                    continue

                # orientation of the gene is considered here
                if gene.startswith('-'):  # negative sign is ignored, but extremities are switched
                    current_chromosome.append('%sh' % gene[1:])
                    current_chromosome.append('%st' % gene[1:])
                else:
                    current_chromosome.append('%st' % gene)
                    current_chromosome.append('%sh' % gene)

            if ident == 'A':
                self.adj_set_A = adjacencies
            else:
                self.adj_set_B = adjacencies


    def connect_adjacencies(self):
        '''
        From a given adjacency set the single extremities are created and stored into dicts.
        Each extremity points to its adjacency (or None, if it is a telomere)
        :param adjA: Adjacency Set of the first_content genome
        :param adjB: Adjacency Set of the second_content genome
        :return: Returns the value of create_circular_graph()
        '''
        adjacencies_a = {}
        adjacencies_b = {}
        wrapper = [(self.adj_set_A, adjacencies_a), (self.adj_set_B, adjacencies_b)]

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

        return self.create_circular_graph(adjacencies_a, adjacencies_b)  # call function to create the cBP

    def create_circular_graph(self, adjacencies_a, adjacencies_b):
        '''
        Uses the networkx package to create the circular breakpoint graph of two given genomes
        :return: the circular breakpoint graph
        '''
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
            # assign value for each vertex from 1 to (n+1)
            for index, vertex in enumerate(long_path):
                #enumerated_vertices[vertex] = index + 1
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
                        self.inter_adj.add(second)
                        continue
                    if second.startswith('Telo'):
                        self.inter_adj.add(first)
                        continue
                    # some nasty set issue. Since we are dealing with strings, there is a difference
                    # between 1h2t and 2t1h. We might change this at some point.
                    if not '%s%s' % (second, first) in self.inter_adj:
                        self.inter_adj.add('%s%s' % (first, second))