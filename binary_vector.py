#!/usr/bin/python

import numpy as np

__author__ = 'klamkiewicz'

def create_vector_for_genome(genome, adjacencies):
    vector = np.empty([1,len(adjacencies)], dtype=int)
    for index,adj in enumerate(adjacencies):
        if adj in genome:
            np.put(vector, [index], [1])
        else:
            np.put(vector, [index], [0])
    return vector