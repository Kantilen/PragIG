#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import numpy as np
from model import Genome
from collections import defaultdict
import random
#################################

class Genome_Sampler():
    '''
    This class generates sampled genomes from the given data.
    It returns a list of potential ancestral genomes.
    '''
    def __init__(self, data, iter, size):
        self.all_adjacencies = data
        self.iteration = iter
        self.size = size

        self.conflicting_adjacencies = defaultdict(list)
        self.preprocess_conflicts()
        self.sampled_genomes = []
        self.create_genomes()

    def preprocess_conflicts(self):
        for adj in self.all_adjacencies:
            self.conflicting_adjacencies[adj.first_ex].append(adj)
            self.conflicting_adjacencies[adj.second_ex].append(adj)

    def create_genomes(self):
        for i in range(self.iteration):
            sampled_genome = []
            adj_copy = self.all_adjacencies[:]
            for k in range(self.size):
                new_adj = random.choice(adj_copy)
                conflicts = list(set(self.conflicting_adjacencies[new_adj.first_ex] + self.conflicting_adjacencies[new_adj.second_ex]))
                try:
                    adj_copy.remove(conflicts)
                except ValueError:
                    None
                sampled_genome.append(new_adj)
