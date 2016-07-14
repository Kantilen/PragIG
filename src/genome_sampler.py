#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
from collections import defaultdict
import random
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
            sampled_genome = []
            conflicts = set()
            for k in range(self.size-1):
                new_adj = random.choice(self.adjacencies)
                while new_adj in conflicts:
                    new_adj = random.choice(self.adjacencies)
                conflicts = conflicts.union(self.conflicting_adjacencies[new_adj.first_ex].union(self.conflicting_adjacencies[new_adj.second_ex]))
                sampled_genome.append(new_adj)
            self.sampled_genomes.append(sampled_genome)
