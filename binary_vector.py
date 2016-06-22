#!/usr/bin/python

import numpy as np

__author__ = 'klamkiewicz'

def create_vector_for_genome(genome, adjacencies):
    vector = np.empty([1,len(adjacencies)], dtype=int)
    for adj in adjacencies:
        if adj in genome:
            print 'Hello', adj

    return vector