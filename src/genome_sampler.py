#!/usr/bin/python

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
    def __init__(self, graph, iter):
        self.iteration = iter

        self.graph = graph
        #self.conflicting_adjacencies = defaultdict(set)
        #self.preprocess_conflicts()
        self.sampled_genomes = []
        #self.create_genomes()
        self.enumerate_vertices()
        #print self.create_adjacency_from_cycle(self.graph)

    def create_adjacency_from_cycle(self, cycle):
        if not cycle:
            return []
        if len(cycle) == 2:
            return [model.Adjacency(cycle[0],cycle[1])]
        assert(len(cycle) % 2 == 0)

        adj = random.randint(1, len(cycle)/2 -1) * 2 + 1
        return self.create_adjacency_from_cycle(cycle[1:adj]) + [model.Adjacency(cycle[0],cycle[adj])] + self.create_adjacency_from_cycle(cycle[adj+1:])

    def enumerate_vertices(self):

        for i in range(self.iteration):
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
            self.sampled_genomes.append(model.Genome.genome_from_adjacencies(i,cycle))



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