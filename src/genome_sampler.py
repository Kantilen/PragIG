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
#################################

class Genome_Sampler():
    '''
    This class generates sampled genomes from the given data.
    It returns a list of potential ancestral genomes.
    '''
    def __init__(self, graph):
        #self.iteration = iter

        self.graph = graph
        #self.conflicting_adjacencies = defaultdict(set)
        #self.preprocess_conflicts()
        #self.sampled_genomes = []
        #self.create_genomes()
        #self.intermediate_cycle = self.enumerate_vertices()
        #print self.create_adjacency_from_cycle(self.graph)

    #@staticmethod
    def create_adjacency_from_cycle(self, cycle):
        if not cycle:
            return []
        assert(len(cycle) % 2 == 0)

        adj = random.randint(0, len(cycle)/2 -1) * 2 + 1
        return self.create_adjacency_from_cycle(cycle[1:adj]) + [model.Adjacency(cycle[0],cycle[adj])] + \
               self.create_adjacency_from_cycle(cycle[adj+1:])

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