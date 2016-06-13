#!/usr/bin/python

import re

__author__ = 'klamkiewicz'

adjacencies_a = {}
adjacencies_b = {}

def connect_adjacencies(adjA, adjB):

    wrapper = [(adjA, adjacencies_a), (adjB, adjacencies_b)]

    for genome,adjacencies in wrapper:
        for adjacency in genome:
            single_extremity = re.search('([\d*\w*]+[t|h])([\d*\w*]+[t|h])', adjacency)
            if not single_extremity:
                adjacencies[adjacency] = None
            else:
                adjacencies[single_extremity.group(1)] = single_extremity.group(2)
                adjacencies[single_extremity.group(2)] = single_extremity.group(1)