#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import numpy as np
from model import Genome
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

        self.sampled_genomes = []

        #self.create_genomes()
