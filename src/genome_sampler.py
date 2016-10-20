#!/usr/bin/python
import itertools

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
from collections import defaultdict
import random
import sys
import networkx as nx
import model
import numpy as np
#################################

class Genome_Sampler():
    '''
    This class generates sampled genomes from the given data.
    It returns a list of potential ancestral genomes.
    '''
    def __init__(self, graph, weights):
        #self.iteration = iter

        self.graph = graph
        self.weights = weights
        #self.conflicting_adjacencies = defaultdict(set)
        #self.preprocess_conflicts()
        #self.sampled_genomes = []
        #self.create_genomes()
        #self.intermediate_cycle = self.enumerate_vertices()
        #print self.create_adjacency_from_cycle(self.graph)

    def weighted_choice(self, choices):
        rnd = random.uniform(0,1)
        tmp = 0
        #print choices, rnd

        for choice, weight in choices.items():
            if tmp + weight >= rnd:
                return choice
            tmp += weight


    #@staticmethod
    def create_adjacency_from_cycle(self, cycle):
        if not cycle:
            return []
        try:
            assert(len(cycle) % 2 == 0)
        except AssertionError:
            print cycle
            return

        if len(cycle) == 2:
            return [model.Adjacency(cycle[0],cycle[1])]

        first_ex = cycle[0]
        i = 0
        while first_ex.startswith("Telo"):
            i += 1
            first_ex = cycle[i]

        inter_ex = [x for x in cycle if (cycle.index(x) + cycle.index(first_ex)) % 2 != 0]



        ancestral_adj = [x for x in self.weights.keys() if x.contains_extremity(first_ex)]
        ancestral_ex = set()
        for adj in ancestral_adj:
            ancestral_ex.update(adj.get_extremities())
        ancestral_ex.remove(first_ex)

        possible_ex = [x for x in ancestral_ex if x in cycle]

        telomeres = [telo for telo in cycle if telo.startswith("Telo") and telo in inter_ex]

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
                weights.append(0.05)

        #if telomeres:
        #    print first_ex, cycle
        #    print telomeres
        #    print ancestral_adj
        #    for x in ancestral_adj:
        #        print x, self.weights[x]
        #    print possible_ex, weights
        #    sys.exit(0)

        prob_sum = sum(weights)
        weights = [x/prob_sum for x in weights]

        second_ex = np.random.choice(possible_ex, 1, p=weights)[0]
        adj = cycle.index(second_ex)

        return self.create_adjacency_from_cycle(cycle[1:adj]) + [model.Adjacency(first_ex, second_ex)] + \
               self.create_adjacency_from_cycle(cycle[adj + 1:])


        #sys.exit(0)

        #adj = random.randint(0, len(cycle)/2 -1) * 2 + 1
        #return self.create_adjacency_from_cycle(cycle[1:adj]) + [model.Adjacency(cycle[0],cycle[adj])] + \
        #       self.create_adjacency_from_cycle(cycle[adj+1:])

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

    @staticmethod
    def get_all(component, text):

        if not component:
            return [None]

        if len(component) == 2:
            return [model.Adjacency(component[0],component[1])]

        assert(len(component) % 2 == 0)

        ig = []
        for k in range(1,len(component),2):
            adj = model.Adjacency(component[0], component[k])
            for left in Genome_Sampler.get_all(component[1:k], "LEFT%d"%k):
                all_found = []
                if not left:
                    left = []
                if type(left) == type([]) and left:
                    all_found = [[adj] + [item] for item in left]
                else:
                    all_found = [adj,left]

                all_found = [x for x in all_found if x]


                for right in Genome_Sampler.get_all(component[k+1:], "RIGHT%d"%k):
                    print all_found
                    if not right:
                        right = []
                    if type(right) == type([]) and right:
                        if len(all_found) == 1:
                            all_found = [[item] + all_found for item in right]
                        else:
                            all_found = [[item] + [element] for item in right for element in all_found]
                    else:
                        print right
                        all_found.extend([right])

                    all_found = [x for x in all_found if x]
                    ig.append(all_found)
        return ig


    def preprocess_conflicts(self):
        for adj in self.graph:

            if adj.first_ex:
                self.conflicting_adjacencies[adj.first_ex].add(adj)
            if adj.second_ex:
                self.conflicting_adjacencies[adj.second_ex].add(adj)

    def create_genomes(self):
        for i in range(self.iteration):
            print "Step %d" % (i)
            sampled_genome = []
            conflicts = dict.fromkeys((self.conflicting_adjacencies.keys()))

            while conflicts:
                new_adj = random.sample(conflicts.keys(), 2)
                sampled_adj = self.conflicting_adjacencies[new_adj[0]].intersection(self.conflicting_adjacencies[new_adj[1]])
                if not sampled_adj:
                    continue

                conflicts.pop(new_adj[0])
                conflicts.pop(new_adj[1])
                sampled_genome.append(list(sampled_adj)[0])
            self.sampled_genomes.append(sampled_genome)