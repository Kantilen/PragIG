#!/usr/bin/python

__author__ = 'klamkiewicz'

#################################
# Import section                #
#################################
import numpy as np
#################################

class Genome_Sampler():

    def __init__(self, data, iter):
        self.breakpoint_graph = data[0]
        self.all_adjacencies = data[1]
        self.iteration = iter
