#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
from collections import defaultdict
import random
import sys
#################################

class Genome_Sampler():
    '''
    This class generates sampled genomes from the given data.
    It returns a list of potential ancestral genomes.
    '''
    def __init__(self, data, iter, size):
        self.iteration = iter
        self.size = size

        self.adjacencies = data
        self.conflicting_adjacencies = defaultdict(set)
        self.preprocess_conflicts()
        self.sampled_genomes = []
        self.create_genomes()

    def preprocess_conflicts(self):
        for adj in self.adjacencies:

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