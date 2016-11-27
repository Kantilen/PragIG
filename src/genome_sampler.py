#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import random
import sys

import networkx as nx
import numpy as np

import model
#################################

class Genome_Sampler():
    '''
    This class generates sampled genomes from the given data.
    It returns a list of potential ancestral genomes.
    '''
    def __init__(self, graph, weights, epsilon):
        self.graph = graph
        self.weights = weights
        self.epsilon = epsilon

    def weighted_choice(self, choices):
        '''
        Own weighted random choice method. Obsolete, since I use the numpy version now.
        :param choices: dict with {choice : weight}
        :return: choice based on the weights
        '''
        rnd = random.uniform(0,1)
        tmp = 0

        for choice, weight in choices.items():
            if tmp + weight >= rnd:
                return choice
            tmp += weight


    def shift(self, cycle):
        '''
        Puts the first element of a list to the end of this list
        :param cycle: list of extremities
        :return: modified list
        '''
        return cycle[1:] + cycle[:1]

    def create_adjacency_from_cycle(self, cycle):
        '''
        Samples one intermediate adjacencies from cycle. Then, a recursive call is made
        to proceed with the two new sub-cycles.
        :param cycle: list of extremities in one component
        :return: list of intermediate adjacencies
        '''
        if not cycle:
            return []
        try:
            assert(len(cycle) % 2 == 0)
        except AssertionError:
            sys.exit(1)


        if len(cycle) == 2:
            return [model.Adjacency(cycle[0],cycle[1])]

        first_ex = cycle[0]

        # If the first extremity is a Telomere, the procedure won't work, since telomeres are not
        # part of the ancestral weights.
        while first_ex.startswith("Telo"):
            cycle = list(np.roll(cycle,-1))
            first_ex = cycle[0]

        # Find possible intermediate adjacencies
        inter_ex = [x for x in cycle if (cycle.index(x) + cycle.index(first_ex)) % 2 != 0]
        # Find all adjacencies that contain the fixed extremity in the ancestral weight dict
        ancestral_adj = [x for x in self.weights.keys() if x.contains_extremity(first_ex)]
        ancestral_ex = set()
        for adj in ancestral_adj:
            ancestral_ex.update(adj.get_extremities())
        ancestral_ex.remove(first_ex)

        # Filter the observed adjacencies with possible intermediate adjacencies.
        possible_ex = [x for x in ancestral_ex if x in inter_ex]

        # If there are more intermediate adjacencies possible then observed,
        # these adjacencies are added as well with weight epsilon.
        if len(inter_ex) != len(possible_ex):
            not_observed = [ex for ex in inter_ex if ex not in possible_ex]
            possible_ex.extend(not_observed)

        weights = []
        for ex in possible_ex:
            try:
                if ex.startswith("Telo"):
                    weights.append(self.weights[model.Adjacency(first_ex,None)])

                else:
                    weights.append(self.weights[model.Adjacency(first_ex,ex)])
            except KeyError:
                weights.append(self.epsilon)

        prob_sum = sum(weights)
        weights = [x/prob_sum for x in weights]

        # sampled a random second extremity based on the weights
        second_ex = np.random.choice(possible_ex, 1, p=weights)[0]

        adj = cycle.index(second_ex)
        # recursive call to sample adjacencies from the sub-cycles.
        return self.create_adjacency_from_cycle(cycle[1:adj]) + [model.Adjacency(first_ex, second_ex)] + \
               self.create_adjacency_from_cycle(cycle[adj + 1:])

    def enumerate_vertices(self):

        #for i in range(self.iteration):
            cycle = []
            for component in nx.connected_component_subgraphs(self.graph):

               enumerated_vertices = {}

               first_adj = component.edges()[0]
               long_path = [x for x in nx.all_simple_paths(component, first_adj[0], first_adj[1])][-1]
               # assign value for each vertex from 0 to n
               for index, vertex in enumerate(long_path):
                   enumerated_vertices[index] = vertex

               cycle.append(self.create_adjacency_from_cycle(enumerated_vertices.values()))
            cycle = [x for y in cycle for x in y]
            #return cycle
            return model.Genome.genome_from_adjacencies('',cycle)