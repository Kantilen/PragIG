#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import numpy as np
from model import Genome
#################################

class Genome_Sampler():
    '''
    This class generates sampled genomes from the given data.
    It returns a list of potential ancestral genomes.
    '''
    def __init__(self, data, iter):
        self.breakpoint_graph = data[0]
        self.all_adjacencies = data[1]
        self.iteration = iter

        self.sampled_genomes = []

        self.create_genomes(self.all_adjacencies, self.iteration)

    def create_genomes(self, adjacency_set, iteration):
        print #type(adjacency_set)

